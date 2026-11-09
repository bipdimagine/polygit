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

#my $dir				= $project->getAlignmentPipelineDir('pbmm2');
my $bam_in			= $patient->getAlignmentFile('pbmm2');
confess ("Error: Input bam '$bam_in' not found") unless (-f $bam_in);
confess ("Error: need bam file, got cram: $bam_in") if ($bam_in =~ /cram$/);


my $dir_out			= $project->getCallingPipelineDir('isoseq_collapse');
my $gff_out			= $dir_out . $patient->name.'.collapsed.gff';


if (-e  $gff_out){
	if ($force){
		unlink $gff_out;
		unlink $gff_out.'pbi';
	}
	else {
		warn "Already done: $gff_out";
		exit(0);
	}
}

my $cmd;
if ($profile =~ /^singlecell/ or $capture_name eq 'transcriptome_10X') {
	# SC
	# isoseq collapse <mapped.bam> <collapsed.gff>
	$cmd = qq{cd $dir_out && run_singularity.sh lima-isoseq.sif isoseq collapse $bam_in $gff_out -j $fork };
}
else {
	# Bulk
	my $bam_flnc = $project->getAlignmentPipelineDir('isoseq_refine') . $patient->name.'.flnc.bam';
	# isoseq collapse --do-not-collapse-extra-5exons <mapped.bam> <flnc.bam> <collapsed.gff>
	$cmd = qq{cd $dir_out && exec_singularity.sh lima-isoseq.sif isoseq collapse --do-not-collapse-extra-5exons $bam_in $bam_flnc $gff_out -j $fork };
}
print $cmd."\n";
system($cmd);
die("Error running isoseq collapse: no output '$gff_out'") unless -e $gff_out;

exit(0);
