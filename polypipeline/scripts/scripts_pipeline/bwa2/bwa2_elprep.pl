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
my $dirin;
my $patient_name;
GetOptions(
	'project=s' => \$project_name,
	'patient=s' => \$patient_name,
	'file1=s' => \$fastq1,
	'file2=s' => \$fastq2,
#	'bam=s'=> \$bam,
	'dir=s' => \$dirin,
	'fork=s' => \$fork,
);
die("hey man,  no fork ") unless $fork;
 my @f1 = split(",",$fastq1);
 my @f2 = split(",",$fastq2);

 
 
 my $project = $buffer->newProject( -name 			=> $project_name );
 $project->genome_version("HG38_CNG");
 $project->version("HG38_CNG");
 my $patient = $project->getPatient($patient_name);
 my $bwa2 = $buffer->software("bwa2");
 my $elprep5 = $buffer->software("elprep5");
  my $bamsormadup = $buffer->software("bamsormadup");
    my $samtools = $buffer->software("samtools");
 my $ref =  $project->getGenomeFasta;
 my $ref_bwa =  $project->dirGenome().$project->buffer->index("bwa2");
 my $ref_elprep =  $project->dirGenome().$project->buffer->index("elprep");

 my $dir_out = $project->getCallingPipelineDir("genome")."/$patient_name";
 system("mkdir -p $dir_out");
 
  my $fr1 = $dir_out."/$patient_name.R1.fastq.gz";
  my $fr2 = $dir_out."/$patient_name.R2.fastq.gz";
 # unlink $fr1 if -e $fr1;
 # unlink $fr2 if -e $fr2;
  my $readgroup = qq{\@RG\\tID:$patient_name\\tSM:$patient_name\\tPL:ILLUMINA\\tDS:$patient_name\\tPG:bwa-mem2};
  my $fileout = $project->getAlignmentPipelineDir("bwa2") . "/" .$patient_name . ".align.bam";
  

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
 my $dir_tmp = $dir_out."/tmp/";
 system("mkdir -p $dir_tmp");
  my $gvcf_file =$dir_out."/$patient_name.g.vcf.gz";
  my $cram_file = $dir_out."/$patient_name.cram";
  my $cmd = qq{$bwa2 mem -R '$readgroup' -t $fork $ref_bwa $fr1 $fr2 };
   $cmd = qq{ | $bamsormadup threads=$fork inputformat=sam  > $fileout && $samtools index $fileout -@ $fork};
 
# $cmd .= qq{| $elprep5 sfm /dev/stdin /dev/stdout --tmp-path $dir_tmp --mark-duplicates --haplotypecaller $gvcf_file  --sorting-order coordinate --log-path $dir_out --intermediate-files-output-prefix $patient_name}.qq{_intermediate --intermediate-files-output-type bam --replace-read-group \"ID:group1 LB:lib1 PL:illumina PU:$project_name SM:$patient_name\" --reference $ref_elprep};
 
  $cmd .= qq{| $elprep5 sfm $fileout /dev/stdout --tmp-path $dir_tmp  --haplotypecaller $gvcf_file   --log-path $dir_out --intermediate-files-output-prefix $patient_name}.qq{_intermediate --intermediate-files-output-type bam --reference $ref_elprep >/dev/null};
 
# $cmd .= qq{ | samtools view -C -T $ref -o $cram_file -@}.$fork;
 warn $cmd;
 system($cmd);
 exit(0) if -e $fileout;
 die();
