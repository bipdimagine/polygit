#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
#use Set::IntSpan;
use GenBoNoSqlLmdb;
use Carp;
use strict;
use Set::IntSpan::Fast::XS ;
use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Carp;
 use JSON::XS;
 use List::MoreUtils qw(natatime uniq);
 use Tabix;
 
 my $buffer = new GBuffer;
my $project_name= "NGS2017_1534";
my $fork;
my $fastq1;
my $fastq2;
my $patient_name;
my $bamout;
my $version;
GetOptions(
	'project=s' => \$project_name,
	'patient=s' => \$patient_name,
	'file1=s' => \$fastq1,
	'file2=s' => \$fastq2,
	'bamout=s'=> \$bamout,
	'fork=s' => \$fork,
	'version' => \$version,
	
);
die("hey man,  no fork ") unless $fork;
 my @f1 = split(",",$fastq1);
 my @f2 = split(",",$fastq2);
 my $project = $buffer->newProject( -name 			=> $project_name );
 if ($version){
 	$project->genome_version("$version");
 	$project->version("$version");
 }
 my $patient = $project->getPatient($patient_name);
 my $bwa2 = $buffer->software("bwa2");
 my $samtools = $buffer->software("samtools");
 my $ref =  $project->getGenomeFasta;
 my $ref_bwa =  $project->dirGenome().$project->buffer->index("bwa2");

 my $dir_out = $project->getAlignmentPipelineDir("bwa2")."/$patient_name";
 
 system("mkdir -p $dir_out");
 
  my $fr1 = $dir_out."/$patient_name.R1.fastq.gz";
  my $fr2 = $dir_out."/$patient_name.R2.fastq.gz";
  unlink $fr1 if -e $fr1;
  unlink $fr2 if -e $fr2;
  my $readgroup = qq{\@RG\\tID:$patient_name\\tSM:$patient_name\\tPL:ILLUMINA\\tDS:$patient_name\\tPG:bwa-mem2};
  
unless (-e $fr1){
 if (scalar(@f1) > 1){
 	system("cat ".join(" ",@f1)." > $fr1");
 	warn "cat ".join(" ",@f1)." > $fr1";
 }
else {
	system("ln -s ".$f1[0]."  $fr1");
} 
}
unless (-e $fr2){
die()  if (scalar(@f2) ne scalar(@f1));
 if (scalar(@f2) > 1){
 	system("cat ".join(" ",@f2)." > $fr2");
 }
else {
	system("ln -s ".$f2[0]."  $fr2");
} 
}
my $cmd = qq{$bwa2 mem -R '$readgroup' -t $fork $ref_bwa $fr1 $fr2 | $samtools view -Sb - -@}.qq{ $fork > $bamout };
warn $cmd;
  system($cmd);
   unlink $fr1;
  unlink $fr2;
 exit(0) if -e $bamout;
 die();
