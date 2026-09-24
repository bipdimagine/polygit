#!/usr/bin/env perl
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
use colored;
use Parallel::ForkManager;
use autodie qw(system open);
use GBuffer;

my $project_name;
my $project_wgs;
my $patient_names;
my $force;
my $no_exec;
my $threads = 20;
$threads = 32 if (`hostname` =~ /^master/);
my $help;

GetOptions(
	'project|diabetome_project=s'	=> \$project_name,
	'genome_project=s'				=> \$project_wgs,
	'patients=s'					=> \$patient_names,
	'force'							=> \$force,
	'no_exec'						=> \$no_exec,
	'fork|cpu=i'					=> \$threads,
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

my @cluster_jobs;

foreach my $pat (@$patients) {
	warn "Preparation des commandes pour : " . $pat->name . "\n";

	my $pat_genome	 = $project_genome->getPatient($pat->name);

	my $bed = $pat->getCaptureFile;
	$bed =~ s/\.gz$//;
	confess("'$bed' does not exist.") unless -e $bed;

	my @patient_cmds; # Stocke les commandes sequentially pour CE patient

	# -------------------------------------------------------------------------
	# 1. CRAM
	# -------------------------------------------------------------------------
	my $cram_dest   = $pat->getCramFileName;
	my $cram_genome = $pat_genome->getBamFile;
	confess("No genome cram found: $cram_genome") unless -e $cram_genome;

	my $idxstats = $cram_dest =~ s/\.cram$/\.idxstats/r;
	my $cram_ok  = (-e $cram_dest && -e "$cram_dest.crai" && -e $idxstats);

	if ($force || !$cram_ok) {
		push @patient_cmds, "$samtools view -@ $threads $cram_genome -L $bed --cram -o $cram_dest && $samtools index  -@ $threads $cram_dest && $samtools idxstats  -@ $threads $cram_dest > $idxstats";
	}

	# -------------------------------------------------------------------------
	# 2. Variants (VCF)
	# -------------------------------------------------------------------------
	foreach my $caller ('haplotypecaller4', 'octopus') {
		my $vcf_dest   = $pat->getVariationsFileName($caller);
		my $vcf_genome = $pat_genome->getVariationsFileName($caller);
		confess("No genome vcf found for $caller: $vcf_genome") unless -e $vcf_genome;

		my $vcf_ok = (-e $vcf_dest && -e "$vcf_dest.tbi");
		if ($force || !$vcf_ok) {
			push @patient_cmds, "$bcftools view --threads $threads $vcf_genome -T $bed -O z -o $vcf_dest && $bcftools index  --threads $threads -tf $vcf_dest";
		}
	}

	# -------------------------------------------------------------------------
	# 3. GVCF (HaplotypeCaller4)
	# -------------------------------------------------------------------------
	my $gvcf_dest   = $pat->gvcfFileName('haplotypecaller4');
	my $gvcf_genome = $pat_genome->getGvcfFile('haplotypecaller4');
	confess("No genome gvcf found: $gvcf_genome") unless -e $gvcf_genome;

	my $gvcf_ok = (-e $gvcf_dest && -e "$gvcf_dest.tbi");
	if ($force || !$gvcf_ok) {
		push @patient_cmds, "$gatk4 SelectVariants -V $gvcf_genome -L $bed -O $gvcf_dest -R $ref && $bcftools index -tf $gvcf_dest";
	}

	# -------------------------------------------------------------------------
	# Aggrégation pour le cluster
	# -------------------------------------------------------------------------
	if (@patient_cmds) {
#		# Tout enchaîner sur 1 seule ligne par patient
#		push @cluster_jobs, join(' && ', @patient_cmds);
		push @cluster_jobs, @patient_cmds;
  	} else {
		warn "Tous les fichiers existent deja pour " . $pat->name . " (Skip)\n";
	}
}

# -----------------------------------------------------------------------------
# Soumission globale au cluster
# -----------------------------------------------------------------------------
if (@cluster_jobs) {
	warn "\nEnvoi de " . scalar(@cluster_jobs) . " job(s) au cluster...\n";
	
	if ($no_exec) {
		warn "[DRY-RUN] Commandes qui auraient été envoyées au cluster :\n";
		print "$_\n" for @cluster_jobs;
	} else {
		# Ouverture d'un pipe d'écriture directement vers run_cluster.pl
		open(my $cluster_fh, '|-', "run_cluster.pl -cpu=$threads");
		
		foreach my $job (@cluster_jobs) {
			print $cluster_fh "$job\n";
		}
		
		# autodie intercepte les erreurs lors de la fermeture si run_cluster.pl échoue
		close($cluster_fh); 
	}
} else {
	warn "Aucune commande à exécuter.\n";
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
	print("Now, run coverage, binary_depth, (dude) and cache on the project:\n");
	print($cmd_pipeline."\n");
#	print($cmd_dude."\n");
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
