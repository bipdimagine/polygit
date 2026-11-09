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
my $bam_in			= $project->getAlignmentPipelineDir('lima_isoseq') . $patient->name.'.fl.5p--3p.bam';
confess ("Error: Input bam '$bam_in' not found") unless (-f $bam_in);

my $dir_out			= $project->getAlignmentPipelineDir('isoseq_tag');
my $bam_out			= $dir_out . $patient->name.'.flt.bam';

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

my $kit				= $patient->getChemistry;
my $design			= '16B-12U-13X-T' if ($kit eq 'SC5P-R2-v3');
$design				= 'T-12U-16B' if ($kit eq 'SC3Pv3' or $kit eq 'SC3Pv4');
confess ("Error: unknown desgin for kit '$kit'") unless ($design);
# T indicates the transcript position and is mandatory. It is the anchor to determine if tags are located on the 3? or 5? side.
# U as in UMI must be preceeded by the length of the UMI.
# B as in cell barcode must be preceed by the length of the cell barcode.
# G is a special PacBio 5? UMI tag with a GGG suffix.
# X must be must be preceeded by the length of an extra sequence to clip such as a TSO sequence.


my $cmd = qq{cd $dir_out && exec_singularity.sh lima-isoseq.sif isoseq tag --design $design $bam_in $bam_out -j $fork};
print $cmd."\n";
system($cmd);
die("Error running isoseq tag: no output '$bam_out'") unless -e $bam_out;

exit(0);




