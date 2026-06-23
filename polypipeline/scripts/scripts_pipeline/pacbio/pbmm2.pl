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
 
 
my $fork = 5;
GetOptions(
	'project=s'   => \$project_name,
	"patient=s" => \$patient_name,
	"fork=s" =>\$fork,
);


my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatient($patient_name);
my $dir_ubam = $patient->getSequencesDirectory();
my $ubam_file = $patient->uBam;
my $dir_out= $project->getAlignmentPipelineDir($patient->name);
my $bam_out = $dir_out."/".$patient->name.".bam";
exit(0) if -e $bam_out;

my $bam_prod = $patient->getBamFileName("pbmm2");

my $ref_index = $project->getGenomeIndex("pbmm2").'/all.mmi';

my $cmd = qq{cd $dir_out && run_singularity.sh pacbio_analysis.sif pbmm2 align $ref_index $ubam_file $bam_out --sort --strip --preset HIFI}
;

print $cmd."\n";
system("$cmd") unless -e $bam_out;
die() unless -e $bam_out;
exit(0);
my $ref               = $project->genomeFasta();
my $cram_out = $patient->getCramFileName("pbmm2");
my $cmd2 = "cd $dir_out && samtools -C -T $ref -o $cram_out $bam_out -@ $fork && samtools index  $cram_out -@ $fork";
system("$cmd2");
exit(0);



