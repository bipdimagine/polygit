#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
#use Set::IntSpan;
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use Carp;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use String::ProgressBar;
use List::Util qw(sum);



 my $project_name;
 my $fork;
 my $callable_intspan_file;
 my $patient_name;
 #my $low_calling;
 my $method;
 my $version;
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"patient=s" => \$patient_name,
	"version=s" => \$version,
);


 my $methods = {
			"manta"				=> { "method" => sub{calling_SV::manta(@_)}, priority =>1}, 
};

my $buffer = GBuffer->new();

my $project = $buffer->newProject( -name => $project_name,-version=>$version);
my $canvas  = $project->getSoftware('canvas');
my $samtools  = $project->getSoftware('samtools');
my $tabix  = $project->getSoftware('tabix');
my $ref =  $project->genomeFasta();
my $patient = $project->getPatient($patient_name);
 my $bam = $patient->getBamFile() ;
  my $dirout= $project->getCallingPipelineDir("canvas-".$patient->name);
 if ($bam =~ /\.cram/){
 	#samtools view -b -T ref.fa -o output_bam.bam input_cram.cram
 	my $cram = $bam;
 	$bam = $dirout."/".$patient->name.".bam";
 	warn "$samtools view -b -T $ref -o $bam $cram --thread $fork && $samtools index $bam --threads $fork";
 	system("$samtools view -b -T $ref -o $bam $cram -\@ $fork && $samtools index $bam -\@ $fork");
 	die() unless -e $bam.".bai";
 }
 
my $dd .="$dirout/".$patient->name().".".time;
system("mkdir -p $dirout");
warn $dd;
if (-e"$dd" ){
	system("rm -r $dd/*");
}
my $dir_canvas=  $project->dirGenome()."/canvas/";
my $genome = $project->genomeFasta();
#dotnet /software/distrib/CANVAS/Canvas-1.40.0.1613+master_x64/Canvas.dll Germline-WGS --bam=/data-xfs/sequencing/ngs/NGS2019_2359/HG19c/align/bwa/DIN_DRI.bam  
#-r /data-xfs/public-data/HG19/canvas/genome.fa -g /data-xfs/public-data/HG19/canvas --filter-bed=/data-xfs/public-data/HG19/canvas/filter13.bed -n toto  -o toto --population-b-allele-vcf=/data-xfs/public-data/HG19/canvas/dbsnp.vcf --ploidy-vcf=./ploidy.vcf
my $ploidy = "$dir_canvas/ploidy_sample_female.vcf";
$ploidy = "$dir_canvas/ploidy_sample_male.vcf" if $patient->isMale();
 my $cmd = " $canvas Germline-WGS --bam=$bam -r $genome -g $dir_canvas --filter-bed=$dir_canvas/filter13.bed -n $patient_name -o $dd --ploidy-vcf=$ploidy --population-b-allele-vcf=$dir_canvas//dbsnp.vcf";
 warn $cmd;
 system($cmd);
 my $final_temp_vcf = "$dirout/$patient_name/CNV.vcf.gz";
 die($final_temp_vcf) unless -e $final_temp_vcf;
 
  my $final_dir =  $project->getVariationsDir("canvas") ;
 my $cmd2 = "mv  $final_temp_vcf $final_dir/$patient_name.vcf.gz && $tabix -p vcf $final_dir/$patient_name.vcf.gz";
 warn $cmd2;
 system("$cmd2");

