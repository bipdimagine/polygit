#! /usr/bin/env perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use GBuffer; 
use Getopt::Long;
use List::Util qw(sum);
use autodie qw(system open);
use Carp;


my $project_name;
my $patient_name;
my $force;
my $fork = 0;
GetOptions(
	'project=s'		=> \$project_name,
	"patient=s"		=> \$patient_name,
	"fork=s"		=> \$fork,
	"force=s"		=> \$force,
) || confess("Error in command line arguments\n");


my $buffer 				= GBuffer->new();
my $project 			= $buffer->newProject( -name => $project_name );
my $patient 			= $project->getPatient($patient_name);
my $profile 			= $patient->getSampleProfile;
my $capture_name		= $patient->getCapture->name;

my $ubam_file;
my $preset;
if ($profile =~ /^singlecell/ or $capture_name eq 'transcriptome_10X') {
	# SC Kinnex ou SC Isoseq
	$ubam_file = $project->getAlignmentPipelineDir('isoseq_groupdedup') . $patient->name.'.dedup.bam';
	$preset = 'ISOSEQ';
}
elsif ($project->isRnaSeq) {
	# bulk Kinnex ou Isoseq
	$ubam_file = $project->getAlignmentPipelineDir('isoseq_cluster2') . $patient->name.'.clustered.bam';
	$preset = 'ISOSEQ';
}
else {
	# not RNAseq
	$ubam_file = $patient->uBam;
	$preset = 'HIFI';
}
warn $ubam_file;
confess ("Error: Input bam '$ubam_file' not found") unless (-f $ubam_file);


my $dir_out = $project->getAlignmentPipelineDir($patient->name);
my $bam_out = $dir_out."/".$patient->name.".bam";
if (-e  $bam_out){
	if ($force){
		unlink $bam_out;
		unlink $bam_out."bai";
	}
	else {
		warn "Already done: $bam_out";
		exit(0);
	}
}
my $bam_prod = $patient->getBamFileName("pbmm2");
if (-e $bam_prod) {
	if ($force){
		unlink $bam_prod;
		unlink $bam_prod."bai";
	}
	else {
		warn "Already done: $bam_prod";
		exit(0);
	}
}
my $cram_prod = $patient->getCramFileName("pbmm2");
if (-e $cram_prod) {
	if ($force){
		unlink $cram_prod;
		unlink $cram_prod."crai";
		unlink $cram_prod =~ s/cram$/idxstats/r;
	}
	else {
		warn "Already done: $cram_prod";
		exit(0);
	}
}

my $ref_index = eval { $project->getGenomeIndex("pbmm2").'/all.mmi' };
my $ref_fasta = $project->genomeFasta;
$ref_index = $ref_fasta if ( ! -f $ref_index or $@ );
confess("Error: reference not found: '$ref_index'") unless (-e $ref_index);
undef $@;

my $cmd = qq{cd $dir_out && exec_singularity.sh pacbio_analysis.sif pbmm2 align $ref_index $ubam_file $bam_out --sort --strip --preset $preset -j $fork};
print $cmd."\n";
system("$cmd");
die("Error running pbmm2: no output '$bam_out'") unless (-e $bam_out);

# todo garder les bam pour les RNAseq, isoseq ne supporte pas les cram
if (grep {/isoseq/} @{$patient->getCallingMethods}) {
	my $dir_align = $project->getAlignmentDir('pbmm2');
	system ("rsync -a $bam_out $bam_out.bai $dir_align");
	exit(0);
}

my $cmd2 = "cd $dir_out && samtools view -C -T $ref_fasta -o $cram_prod $bam_out -@ $fork && samtools index $cram_prod -@ $fork";
system("$cmd2") unless -e $cram_prod;
die("Error converting bam to cram: no output '$cram_prod'") unless -e $cram_prod;

my $idxstats = $cram_prod =~ s/cram$/idxstats/r;
system("samtools idxstats $cram_prod  -@ $fork > $idxstats");
die("Error running idxstats: no output '$idxstats'") unless -e $idxstats;


exit(0);




