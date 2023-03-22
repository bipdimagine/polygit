#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages/";
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

use file_util;


 my $project_name;
 my $fork;
 my $callable_intspan_file;
 my $patient_name2;
 #my $low_calling;
 my $version;
 my $method;
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"patient=s" => \$patient_name2,
	"version=s" => \$version,
);


my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name,-version=>$version );

if ($patient_name2 =~ / /){
	my $query = $buffer->getQuery();
	my $res   = $query->getPatients( $project->id );
	my ($a,$b) = split(" ",$patient_name2);
	my (@find) = grep {$_->{name} =~/$a/ && $_->{name} =~ /$b/} @$res;
	die() if scalar(@find) ne 1;
	$patient_name2 = $find[0]->{name};
	warn $find[0]->{name};
}

#exit(0) if -e "/data-isilon/sequencing/muc1/renome/".$patient_name2.".vcf";
my $patients = $project->get_only_list_patients($patient_name2);
die() unless scalar(@$patients);
warn "coucou ".scalar(@$patients);
 my $dir_pipeline = "/tmp/pipeline/$project_name.".time."/";#

 
 system("mkdir -p $dir_pipeline") unless -e $dir_pipeline;
foreach my $patient (@$patients) {
my $patient_name = $patient->name;
my $bam = $patient->getBamFile();
my $bam_dir = $patient->getProject->getAlignmentDir("bwa");
  my $tmp_dir = $project->getCallingPipelineDir("muc1-vntyper.".$patient->name);
  $dir_pipeline = $tmp_dir;
my $fileoutx=  $project->getVariationsDir("muc1_vntyper6")."/".$patient_name.".vcf";
 my $dir_pipeline2 = $dir_pipeline ."/$patient_name/";
 system("rm -r $dir_pipeline2") if -e $dir_pipeline2;

 #exit(0) if -e $fileout;

 my $dir_working = "$dir_pipeline/$patient_name"  ;
  if (-e $dir_working){
 	system("rm $dir_pipeline/*.log*");
 }
 else {
 system ("mkdir $dir_working");
 }
my $singularity = $buffer->software("singularity");
#my $image = "/software/distrib/ADVNTR/SINGLARITY/vntyper_v7_20221209.sif";
my $image = "/software/distrib/ADVNTR/SINGLARITY/vntyper.sif";
my $db = "/data-isilon/public-data/repository/HG19/vntr/";
system ("mkdir $dir_pipeline/temp") unless -e "$dir_pipeline/temp";

#my $scommand ="python3 /SOFT/VNtyper/VNtyper.py -t 5 -k 20 -ref_VNTR /tmp/pipeline/MUC1-VNTR.fa  -f /SOFT/VNtyper/Files/ -ref /tmp/pipeline/MUC1-VNTR.fa -p Scripts/ -w $dir_pipeline -r1 $fastq1 -r2 $fastq2 -o $patient_name -m Files/vntr_data/hg19_genic_VNTRs.db --ignore_advntr  -p /SOFT/VNtyper/";
my $scommand ="python3 /SOFT/VNtyper/VNtyper.py   -t 5    --bam  -ref chr1.fa -ref_VNTR /SOFT/VNtyper/Files/MUC1-VNTR.fa -p /SOFT/VNtyper/Scripts/ -w $dir_pipeline -a $bam   -o $patient_name  -m /vntr_data/hg19_genic_VNTRs.db --ignore_advntr  -p /SOFT/VNtyper/";
my $cmd = qq{$singularity run --pwd /DATA/adVNTR/ -B /data-isilon:/data-isilon -B /tmp/:/tmp/ -H $tmp_dir  $image $scommand};
warn $cmd;
system($cmd." && touch $dir_pipeline/$patient_name.ok");

unless (-e "$dir_pipeline/$patient_name.ok"){
 #	system("rm -r $dir_pipeline");
 	die();
 }
warn "ok ************";
warn $dir_pipeline2;
my $dir_prod = $project->getVariationsDir("vntyper")."/muc1/";
my $tsv = $dir_prod."/".$patient->name."_Final_result.tsv";
if (-e $tsv ){
	system("$Bin/../rm_vcf.pl $tsv");
}
system("mkdir -p $dir_prod;chmod a+rwx $dir_prod") unless -e $dir_prod;
system("cp $dir_pipeline2/".$patient->name."* $dir_prod");
}
exit(0);
