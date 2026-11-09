#! /usr/bin/env perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use GBuffer; 
use Getopt::Long;
use Carp;
use autodie qw(system);


my $project_name;
my $patient_name;
my $force;
my $fork = 0;
GetOptions(
	'project=s'	=> \$project_name,
	'patient=s'	=> \$patient_name,
	'fork=i'	=> \$fork,
	'force'		=> \$force,
) || confess("Error in command line arguments\n");

die("Error: -project is mandatory") unless ($project_name);
die("Error: -patient is mandatory") unless ($patient_name);

my $buffer			= GBuffer->new();
my $project			= $buffer->newProject( -name => $project_name );
my $patient			= $project->getPatient($patient_name);
my $profile 		= $patient->getSampleProfile;
my $capture_name	= $patient->getCapture->name;
my $kit 			= $patient->getChemistry;

my $bam_in;
if ($profile =~ /^singlecell/ or $capture_name eq 'transcriptome_10X') {
	$bam_in = $project->getAlignmentPipelineDir('isoseq_tag') . $patient->name.'.flt.bam';
}
else {
	$bam_in = $project->getAlignmentPipelineDir('lima_isoseq') . $patient->name.'.fl.5p--3p.bam'
}
warn $bam_in;
confess ("Error: Input bam '$bam_in' not found") unless (-f $bam_in);

my $dir_out			= $project->getAlignmentPipelineDir('isoseq_refine');
my $bam_out			= $dir_out . $patient->name.'.flnc.bam';

if (-e  $bam_out){
	if ($force){
		unlink $bam_out;
		unlink $bam_out.'pbi';
	}
	else {
		warn "Already done: $bam_out";
		exit(0);
	}
}

my $primers;
if ($profile =~ /^singlecell/ or $capture_name eq 'transcriptome_10X') {
	if ($kit =~ /^SC3P/) {
		$primers = $buffer->config_path('root','public_data')."pacbio/10x_3p_primers.fa";
	}
	elsif ($kit =~ /^SC5P/) {
		$primers = $buffer->config_path('root','public_data')."pacbio/10x_5p_primers.fa";
	}
	else {
		confess("Error determinig the primers fasta for chemistry/kit '$kit'");
	}
}
else {
	$primers = $buffer->config_path('root','public_data')."pacbio/isoseq_primers.fa";
}
confess("Error: '$primers' doesn't exists.") unless (-f $primers);


my $cmd = qq{cd $dir_out && exec_singularity.sh lima-isoseq.sif isoseq refine $bam_in $primers $bam_out -j $fork};
print $cmd."\n";
system($cmd);
die("Error running isoseq refine: no output '$bam_out'") unless -e $bam_out;

exit(0);
