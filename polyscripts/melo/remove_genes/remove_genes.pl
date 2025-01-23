#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages/";
use lib "$Bin/../../../GenBo/lib/GenBoDB/writeDB/";
use lib "$Bin/packages";use Data::Dumper;
use Getopt::Long;
use Carp;
use IO::Prompt;
use colored;
use Cwd 'abs_path';
use Net::SFTP;
use Archive::Tar;
use Text::CSV;
use File::Find;

use GBuffer;
my $buffer = new GBuffer;


my $project_name;	# 'NGS2024_8349'
my $version;
my $patient_names;	# 110018,110002,110019
my $bed;
my $extend = 300;
my $no_exec;
my $help;

GetOptions(
	'project=s'		=> \$project_name,		# Projet
	'version=s'		=> \$version,
	'patients=s'	=> \$patient_names,		# Patient(s)
	'bed=s'			=> \$bed,
	'extend=i'		=> \$extend,
	'no_exec'		=> \$no_exec,
	'help'			=> \$help,		
);

usage() unless ($project_name and $bed);
usage() if ($help);

my $project = $buffer->newProject(-name=>$project_name, -version=>$version);
$version = $project->version unless ($version);
warn $version;
my $patients = $project->get_only_list_patients($patient_names);
warn("No patient in project ".$project->name."\n") unless $patients;

# Extend bed
unless($extend==0){
	my $ext_bed = $bed =~ s/\.bed$/_extended\.bed/r;
	open(EXT_BED, '>', $ext_bed);
	open(BED, '<', $bed) or die "Can't read file '$bed' [$!]\n";
	while (my $line = <BED>) {
	    chomp $line;
	    my ($chr, $start, $end, $name) = split("\t", $line);
	    $chr =~ s/^chr// if ($version eq 'HG19_CNG');
	    print EXT_BED join("\t", ($chr, $start-$extend, $end+$extend, $name)) ."\n";
	}
	close(EXT_BED);
	close(BED);
	$bed = $ext_bed;
}
warn $bed;

my $dir = $project->getProjectPath;
my $samtools = $buffer->software("samtools");
my $bcftools = $buffer->software("bcftools");
my $cmd_filter;
my $cmd;
my @all_bam;
my @all_vcf;

foreach my $pat (@$patients) {
	my $pname = $pat->name;
	print "$pname\n";
	
#	# Get the bam and vcf files from the methods in the database
#	my @bam_files = @{$pat->getBamFiles};
#	my @vcf_files = @{$pat->getVariationsFiles};
#	my $gvcf = $pat->getGvcfFile;
#	push (@vcf_files,$gvcf) if ($gvcf);
#	print join("\n",@bam_files)."\n" if (scalar @bam_files);
#	print join("\n",@vcf_files)."\n" if (scalar @vcf_files);

	# Get all the bam and vcf files (independantly of the alignment/calling methods)
	# bam
	my @bam_files;
	my $find_bam = sub {
		push (@bam_files, $File::Find::name) if (/$pname\.bam$/);
	};
	find($find_bam, $dir);
	push (@all_bam, @bam_files);
	print scalar(@bam_files)." bam file(s):\n";
	print join("\n",@bam_files)."\n" if (scalar @bam_files);
	
	# vcf
	my @vcf_files;
	my $find_vcf = sub {
		push (@vcf_files, $File::Find::name) if (/$pname.*\.vcf(\.gz)?$/) ;
	};
	find($find_vcf, $dir);
	push (@all_vcf, @vcf_files);
	print scalar(@vcf_files)." vcf file(s):\n";
	print join("\n",@vcf_files)."\n" if (scalar @vcf_files);
	print "\n";
}


# Commandes pour filtrer et indexer les bam/vcf
foreach my $bam (@all_bam) {
#	next if ($bam =~ /old/);
	if (-e $bam.'.back') {
#	if (-e $bam.'.back' and -e $bam.'.back2') {
		warn "'$bam.back' already exists, '$bam' not filtered";
		next;
	};
	my $filtered_bam = $bam =~ s/\.bam$/\.filtered\.bam/r;
	$cmd_filter .= "$samtools view $bam -L $bed -U $filtered_bam -o /dev/null -@ 20 ; $samtools index -@ 20 $filtered_bam ; touch $filtered_bam.ok"."\n";
	$cmd->{$bam} = "mv $bam $bam.back ; mv $bam.bai $bam.bai.back ; mv $filtered_bam $bam ; mv $filtered_bam.bai $bam.bai ; rm $filtered_bam.ok" unless (-e "$bam.back");
#	$cmd->{$bam} = "mv $bam $bam.back2 ; mv $bam.bai $bam.bai.back2 ; mv $filtered_bam $bam ; mv $filtered_bam.bai $bam.bai" if (-e "$bam.back");
}
foreach my $vcf (@all_vcf) {
	if (-e $vcf.'.back') {
#	if (-e $vcf.'.back' and -e $vcf.'.back2') {
		warn "'$vcf.back' already exists, '$vcf' not filtered";
		next;
	}
	my $filtered_vcf = $vcf =~ s/\.vcf\.gz$/\.filtered\.vcf\.gz/r;
	$cmd_filter .= "$bcftools view $vcf -T ^$bed -o $filtered_vcf -l 6 --threads 20 ; bcftools index -t $filtered_vcf ; touch $filtered_vcf.ok"."\n";
	$cmd->{$filtered_vcf} = "mv $vcf $vcf.back ; mv $vcf.tbi $vcf.tbi.back ; mv $filtered_vcf $vcf ; mv $filtered_vcf.tbi $vcf.tbi ; rm $filtered_vcf.ok" unless (-e "$vcf.back");
#	$cmd->{$filtered_vcf} = "mv $vcf $vcf.back2 ; mv $filtered_vcf $vcf" if (-e "$vcf.back");
}
print "\n";
print $cmd_filter."\n";
print join("\n", values (%$cmd))."\n";

# Run les commandes
chomp $cmd_filter;
system("echo \"$cmd_filter\" | run_cluster.pl -cpu=20") unless ($no_exec);

# Vérifie que les fichiers de sortie existent et les renomme
my $errors;
print("Renaming and indexing files...\n") unless ($no_exec);
foreach my $file (keys %$cmd) {
	$errors->{$file} = "Error while filtering file '".$file=~ s/\.filtered\./\./r."'" unless (-e $file.'.ok' or $no_exec);
	system($cmd->{$file}) if (-e $file.'.ok' and not $no_exec);
}
die (join("\n",values(%$errors))) if ($errors);
print("######\nDone !\n######\n") unless ($errors or $no_exec);



sub usage {
	print "
remove_genes.pl
-----------------	
Mandatory arguments
	-project <s>		project name containing the complementary bam
	-patients <s>		samples names separated with a comma
	-bed <s>		bed file with genes to remove
	
Optional arguments
	-extend <i>		size of the bed extension, default 300
	-version <s>		project version
	-no_exec		do not run the commands
	-help				

";
	exit(1);
}