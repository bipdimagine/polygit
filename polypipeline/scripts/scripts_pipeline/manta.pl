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
 my $version;
 my $method;
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"patient=s" => \$patient_name,
	"version=s" => \$version,
);

$fork =10 unless $fork;
my $buffer = GBuffer->new();

my $project = $buffer->newProject( -name => $project_name,-version=>$version );
my $manta  = $project->getSoftware('manta');
my $tabix  = $project->getSoftware('tabix');
my $ref =  $project->genomeFasta();
 #$ref = qq{/data-isilon/public-data/genome/HG19_without_MT/fasta/all.fa};
my $patient = $project->getPatient($patient_name);
 my $bam = $patient->getBamFile() ;
 my $dirout= $project->getCallingPipelineDir("manta");

my $dd .="$dirout/".$patient->name().".".time;
system("mkdir -p $dirout");
warn $dd;

if (-e"$dd" ){
	system("rm  $dd/*");
}
my $cmd = "$manta --bam $bam --referenceFasta $ref --runDir $dd  ";
unless ($patient->isGenome){
	$cmd = "$manta --config /software/distrib/MANTA/manta-1.6.0.centos6_x86_64/bin/configManta.py.ini1 --bam $bam --referenceFasta $ref --runDir $dd --exome ";# --callRegions $bed2 ";
}
 system($cmd);
 
 my $cmd2 = "$dd/runWorkflow.py -m local -j $fork";
 system($cmd2);
  warn $dd;
my $final_dir= $project->getVariationsDir("manta") ;
 my $out_vcf = $final_dir."/".$patient->name.".vcf.gz";
 my $cmd3;
 if ($project->isGenome()){
 	 $cmd3 = qq{mv $dd//results/variants/diploidSV.vcf.gz $out_vcf};
 }
 else {
 	my $bed_file = $dirout."/".$patient->name.".bed"; 
 	$project->writeCaptureBedFile(200,$bed_file);
 	$bed_file = $buffer->gzip_tabix($bed_file,"bed");
 	my $bcftools = $buffer->software("bcftools");
 	 $cmd3 = qq{ $bcftools view $dd//results/variants/diploidSV.vcf.gz -R $bed_file -O z > $out_vcf};
 	
 }

 warn $cmd3;
 system($cmd3);
 $buffer->gzip_tabix($out_vcf,"vcf");
 exit(0);
 