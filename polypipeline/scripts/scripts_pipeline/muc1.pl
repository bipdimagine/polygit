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


my $buffer = GBuffer->new();
warn $patient_name;
my $project = $buffer->newProject( -name => $project_name,-version=>$version );
my $patient = $project->getPatient($patient_name);
my $bam = $patient->getBamFile();
my $fileout =  $project->getVariationsDir("muc1")."/".$patient_name.".bed";
my $fileout2 =  $project->getVariationsDir("muc1")."/".$patient_name.".aln";
my $fileout4 =  $project->getVariationsDir("muc1")."/".$patient_name.".log";
my $fileout3 =  $project->getVariationsDir("muc1")."/".$patient_name.".error";
my $bam_dir = $patient->getProject->getAlignmentDir("bwa");
 my $dirout= $project->getCallingPipelineDir("muc1");

 #exit(0) if -e $fileout;
 my $dir_working = "$dirout/$patient_name"  ;
  if (-e $dir_working){
 	system("rm $dirout/*.log*");
 }
 else {
 system ("mkdir $dir_working");
 }
my $singularity = $buffer->software("singularity");
my $image = "/software/distrib/ADVNTR/SINGLARITY/advntr_2718.sif";
my $db = "/data-isilon/public-data/repository/HG19/vntr/";
#my $cmd = qq{$singularity run --pwd /DATA/adVNTR/ -B $dir_working:/TMP -B $db:/DB/ -B $dirout:/OUT -B $bam_dir:/BAM $image  advntr genotype -a /BAM/$patient_name.bam -fs -vid 25561 -t 4  -of bed -o /OUT/$patient_name.bed --min_read_length 100 -m /DB/hg19_genic_VNTRs.db --working_directory  /TMP/};
#my $cmd = qq{$singularity run --pwd /DATA/adVNTR/ -B $dir_working:/TMP -B $db:/DB/ -B $dirout:/OUT -B $bam_dir:/BAM $image  advntr genotype -a /BAM/$patient_name.bam -fs -vid 25561 -t 4  -o /OUT/$patient_name.bed --min_read_length 100 -m /DB/hg19_genic_VNTRs.db --working_directory  /TMP/};
my $cmd = qq{$singularity run --pwd /DATA/adVNTR/ -B $dir_working:/TMP -B $db:/DB/ -B $dirout:/OUT -B $bam_dir:/BAM $image  advntr genotype -a /BAM/$patient_name.bam -fs -vid 25561 -t 4  -o /OUT/$patient_name.bed  -m /DB/hg19_genic_VNTRs.db --working_directory  /TMP/};

warn $cmd;

system($cmd." && touch $dirout/$patient_name.ok");



system("cp $dirout/$patient_name.bed $fileout");
#die() unless -e $dir_working."/log_".$patient_name.".bam.log.aln";
my $file_aln = $dir_working."/log_".$patient_name.".bam.log.aln";
my $file_log = $dir_working."/log_".$patient_name.".bam.log";
system("cp $file_aln $fileout2") if -e $file_aln;
system("cp $file_log $fileout4") if -e $file_log;
#system("touch $fileout3") unless -e "$dirout/$patient_name.ok";
die() unless -e "$dirout/$patient_name.ok";
