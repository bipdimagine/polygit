#! /usr/bin/env perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
#use Set::IntSpan;
use GBuffer; 
use Getopt::Long;
use List::Util qw(sum);



 my $project_name;
 my $patient_name;
 #my $low_calling;
 my $method;
 my $force;
 
my $fork = 5;
GetOptions(
	'project=s'   => \$project_name,
	"patient=s" => \$patient_name,
	"fork=s" =>\$fork,
	"force=s" =>\$force,
);


my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatient($patient_name);
my $dir_ubam = $patient->getSequencesDirectory();
my $ubam_file = $patient->uBam;
my $dir_out= $project->getAlignmentPipelineDir($patient->name);
my $bam_out = $dir_out."/".$patient->name.".bam";
if (-e  $bam_out){
	if ($force){
		unlink $bam_out;
		unlink $bam_out."bai";
	}
	else {
		warn "already done";
		exit(0);
	}
}
my $bam_prod = $patient->getBamFileName("pbmm2");

my $ref_index = $project->getGenomeIndex("pbmm2").'/all.mmi';
unless (-e $bam_out){
my $cmd = qq{cd $dir_out && run_singularity.sh pacbio_analysis.sif pbmm2 align $ref_index $ubam_file $bam_out --sort --strip --preset HIFI -j $fork};

print $cmd."\n";
system("$cmd") unless -e $bam_out;
die() unless -e $bam_out;
}

exit(0);




