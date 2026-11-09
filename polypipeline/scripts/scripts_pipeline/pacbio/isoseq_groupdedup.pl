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
my $dir				= $project->getAlignmentPipelineDir('isoseq_correct');
my $bam_in			= $dir . $patient->name.'.flnc-bccorr.bam';
confess ("Error: Input bam '$bam_in' not found") unless (-f $bam_in);

my $bam_sorted		= $dir . $patient->name.'.flnc-bccorr.sorted.bam';

my $dir_out			= $project->getAlignmentPipelineDir('isoseq_groupdedup');
my $bam_out			= $dir_out . $patient->name.'.dedup.bam';


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

my $cmd = qq{cd $dir_out && samtools sort -@ 32 -t CB $bam_in -o $bam_sorted && exec_singularity.sh lima-isoseq.sif isoseq groupdedup --keep-non-real-cells $bam_sorted $bam_out -j $fork };
print $cmd."\n";
system($cmd);
die("Error running isoseq groupdedup: no output '$bam_out'") unless -e $bam_out;

exit(0);
