#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages/";
use lib "$Bin/../../../GenBo/lib/GenBoDB/writeDB/";
use lib "$Bin/packages";
use Data::Dumper;
use Getopt::Long;
use Carp;
use IO::Prompt;
use Term::Menus;
use colored;
use Cwd 'abs_path';
use Net::SFTP;
use Archive::Tar;
use File::Util;
use Parallel::ForkManager;

use GBuffer;

my $project_name;
my $project_wgs;
my $patient_names;
my $force;
my $no_exec;
my $fork = 1;
my $help;

GetOptions(
	'project|diabetome_project=s'	=> \$project_name,
	'genome_project=s'				=> \$project_wgs,
	'patients=s'					=> \$patient_names,
	'force'							=> \$force,
	'no_exec'						=> \$no_exec,
	'fork=i'						=> \$fork,
	'help'							=> \$help,
) or (warn("\nError in command line arguments\n") && usage());
usage() if ($help);
die('Enter a diabetome project name') unless ($project_name);
die('Enter the name of the corresponding genome project') unless ($project_wgs);


my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name );
my @captures = map {$_->name} @{$project->getCaptures()};
confess("Not all patients are capture 'DIABETomeV1_hg38': ".Dumper \@captures) unless (scalar @captures == 1 and $captures[0] eq 'DIABETomeV1_hg38');
warn $project_name.': '.$project->description;
my $patients = $project->get_only_list_patients($patient_names);
die("No patient in project ".$project_name."\n") unless ($patients);
@$patients = sort {$a->name cmp $b->name} @$patients;

my $buffer_genome = new GBuffer;
my $project_genome = $buffer_genome->newProject( -name => $project_wgs );
warn $project_wgs.': '.$project_genome->description;

my $samtools = $buffer->software("samtools");
my $bcftools = $buffer->software("bcftools");
my $gatk4 = $buffer->software("gatk4");
my $ref = $project->genomeFasta();

# Remove SV methods 
my @SVmethods = map{values @{$_->callingSVMethods}} @$patients;
if (grep {/canvas|manta|wisecondor/} @SVmethods and not $project->isGenome) {
	warn("Removing SV methods for genomes");
	system("del_calling_method.sh -project=$project_name -methods=canvas,manta,wisecondor") unless ($no_exec);
}

my $pm   = new Parallel::ForkManager($fork);
my $erreur_fork = 0;
my @patients_in_error;
$pm->run_on_finish(sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $block) = @_;
    if ($exit_signal || $exit_code != 0) {
        $erreur_fork = 1;
        push @patients_in_error, $ident;
		warn "ERROR: Patient '$ident' failed (PID: $pid, Exit: $exit_code)\n";
	} elsif (not $no_exec) {
		warn "OK: Patient '$ident' completed successfully\n";
	}
});

$buffer->disconnect;
foreach my $pat (@$patients) {
	my $pid = $pm->start($pat->name) and next;
	warn "\n" . $pat->name;
	my $buffer_genome = new GBuffer;
	my $project_genome = $buffer_genome->newProject( -name => $project_wgs );
	my $pat_genome = $project_genome->getPatient($pat->name);
	my $bed = $pat->getCaptureFile;
	$bed =~ s/\.gz$//;
	confess ("'$bed' does not exist. Unzip the bed file to use the ungzipped bed to extract regions from vcf") unless (-e $bed);
	warn $bed;
	
	# CRAM
	my $cram_dest = $pat->getCramFileName;
	warn ("Diabetome cram '$cram_dest' already exists") if (-e $cram_dest);
	my $cram_genome = $pat_genome->getBamFile;
	confess ("No genome cram found: $cram_genome") unless (-e $cram_genome);
	my $idxstats = $cram_dest =~ s/\.cram$/\.idxstats/r;
	my $cmd = "$samtools view $cram_genome -L $bed --cram -o $cram_dest && $samtools index $cram_dest && $samtools idxstats $cram_dest > $idxstats" unless (-e $cram_dest and -e $cram_dest.'.crai' and not $force);
	warn $cmd if ($cmd);
	system($cmd) unless ($no_exec);
	confess("Error: $cram_dest") unless (-e $cram_dest and -e $cram_dest.'.crai' or $no_exec);
	
	# Variants/VCF
	# Haplotypecaller4
	my $vcf_dest = $pat->getVariationsFileName('haplotypecaller4');
	warn ("Diabetome vcf '$vcf_dest' already exists") if (-e $vcf_dest);
	my $vcf_genome = $pat_genome->getVariationsFile('haplotypecaller4');
	confess ("No genome vcf found: $vcf_genome") unless (-e $vcf_genome);
	my $cmd = "$bcftools view $vcf_genome -T $bed -O z -o $vcf_dest && $bcftools index -tf $vcf_dest" unless (-e $vcf_dest and -e $vcf_dest.'.tbi' and not $force);
	warn $cmd if ($cmd);
	system($cmd) unless ($no_exec);
	confess("Error: $vcf_dest") unless (-e $vcf_dest or $no_exec);
	
	my $gvcf_dest = $pat->gvcfFileName('haplotypecaller4');
	warn ("Diabetome gvcf '$gvcf_dest' already exists") if (-e $vcf_dest);
	my $gvcf_genome = $pat_genome->getGvcfFile('haplotypecaller4');
	confess ("No genome gvcf found: $gvcf_genome") unless (-e $gvcf_genome);
	my $cmd = "$gatk4 SelectVariants -V $gvcf_genome -L $bed -O $gvcf_dest -R $ref && $bcftools index -tf $gvcf_dest" unless (-e $gvcf_dest and -e $gvcf_dest.'.tbi' and not $force);
	warn $cmd if ($cmd);
	system($cmd) unless ($no_exec);
	confess("Error: $gvcf_dest") unless (-e $gvcf_dest or $no_exec);
	
	# Octopus
	my $vcf_dest = $pat->getVariationsFileName('octopus');
	warn ("Diabetome vcf '$vcf_dest' already exists") if (-e $vcf_dest);
	my $vcf_genome = $pat_genome->getVariationsFileName('octopus');
	confess ("No genome vcf found: $vcf_genome") unless (-e $vcf_genome);
	my $cmd = "$bcftools view $vcf_genome -T $bed -O z -o $vcf_dest && $bcftools index -tf $vcf_dest" unless (-e $vcf_dest and -e $vcf_dest.'.tbi' and not $force);
	warn $cmd if ($cmd);
	system($cmd) unless ($no_exec);
	confess("Error: $vcf_dest") unless (-e $vcf_dest or $no_exec);
	
#	# Dragen-calling
#	my $vcf_dest = $pat->getVariationsFileName('dragen-calling');
#	my $vcf_genome = $pat_genome->getVariationsFile('dragen-calling');
#	my $cmd = "$bcftools view $vcf_genome -T $bed -O z -o $vcf_dest && $bcftools index -tf $vcf_dest" unless (-e $vcf_dest and -e $vcf_dest.'.tbi' and not $force);
#	warn $cmd if ($cmd);
#	system($cmd) unless ($no_exec);
#	confess("Error: $vcf_dest") unless (-e $vcf_dest or $no_exec);
	
	$pm->finish;
}
$pm->wait_all_children();
if ($erreur_fork) {
    warn "\n========================================\n";
    warn "ERROR: Processing failed for " . scalar(@patients_in_error) . " patient(s):\n";
    warn "  - " . join("\n  - ", @patients_in_error) . "\n";
    warn "========================================\n";
}

unless ($no_exec) {
	print("\n----------DONE----------\n\n");
	my $cmd_pipeline = "$Bin/../../bds_pipeline.pl -project $project_name -steps coverage,binary_depth";
	$cmd_pipeline .= " -patients $patient_names" if ($patient_names);
	$cmd_pipeline .= " -force 1" if ($force);
	my $cmd_dude = "$Bin/../../bds_calling.pl -project $project_name -patient all -steps dude" ;
	$cmd_dude .= " -force 1" if ($force);
	my $cmd_cache = "$Bin/../../bds_cache.pl -project $project_name" ;
	$cmd_cache .= " -force 1" if ($force);
	print("Now, run coverage, binary_depth, dude and cache on the project:\n");
	print($cmd_pipeline."\n");
	print($cmd_dude."\n");
	print($cmd_cache."\n");
}
print "\n";


sub usage {
	print "
$0
Extracts regions corresponding to the capture of <project> from the genome files of <genome_project>
-----------------	
Mandatory arguments
	-diabetome_project <s>			project name
	-genome_project <s>			project name of the corresponding genome project where the files are extracted
	
Optional arguments
	-patients <s>			patient names separated with a comma
	-force				overwrite files if existing
	-fork <i>			number of forks to use in parallele
	-no_exec			do not execute the commands
	-help				display this help message and exit

Don't forget to run coverage, lmdb_depth, dude and cache after.
";
	exit(1);
}
