#! /usr/bin/env perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
#use Set::IntSpan;
use GBuffer; 
use Getopt::Long;



 my $bam_out;
 my $cram_out;
 #my $low_calling;
 my $method;
 my $project_name;
 my $patient_name;
my $fork = 64;
my $force;
GetOptions(
	'bam=s'   => \$bam_out,
	"cram=s" => \$cram_out,
	"project=s" => \$project_name,
	"patient=s" => \$patient_name,
	"force=s" => \$force,
);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatient($patient_name);
my $dir_out= $project->getAlignmentPipelineDir($patient->name);
my $bam_out = $dir_out."/".$patient->name.".bam";
my $cram_out = $patient->getCramFileName("pbmm2");
if (-e $cram_out){
	if ($force){
		unlink $cram_out;
		unlink $cram_out.".cai";
	}
	else {
		warn "already";
		exit(0);
	}
}
die($bam_out) unless -e $bam_out;
my $ref               = $project->genomeFasta();

$fork =30 if $fork > 40;
my $cmd2 = "cd $dir_out && samtools view -C -T $ref -o $cram_out $bam_out -@ $fork && samtools index  $cram_out -@ $fork";
system("$cmd2") unless -e $cram_out;
die() unless -e $cram_out;
 my $filename = $cram_out;
 $filename =~ s/cram/idxstats/;
system("samtools idxstats $cram_out  -\@ $fork  > $filename");


exit(0);