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
my $bamin;
my $bamout;
my $version;
GetOptions(
	'project=s' => \$project_name,
	'version=s' => \$version,
);

 my $project = $buffer->newProject( -name 			=> $project_name,-version=>$version );
 
 foreach my $chr (@{$project->getChromosomes}){
 	print $chr->fasta_name."\t"."1\t".$chr->length."\n";
 }
 exit(0);
 my $patient = $project->getPatient($patient_name);
 my $bwa2 = $buffer->software("bwa2");
 my $elprep5 = $buffer->software("elprep5");
 my $ref =  $project->getGenomeFasta;
 my $ref_elprep =  $project->dirGenome().$project->buffer->index("elprep");

 my $dir_out = $project->getCallingPipelineDir("genome")."";
 my $gvcf_file = $project->getCallingPipelineDir("genome")."/$patient_name.g.vcf.gz";
 #--haplotypecaller $gvcf_file
 my $dir_tmp = "/tmp/";
 my $cmd = qq{$elprep5 sfm $bamin $bamout --tmp-path $dir_tmp --mark-duplicates  --sorting-order coordinate --log-path $dir_out  --intermediate-files-output-prefix $patient_name}.qq{_intermediate --intermediate-files-output-type bam  --reference $ref_elprep};
  
 
 system($cmd);
# die() if -s $bamout;
 exit(0) if -e $bamout;
 die();