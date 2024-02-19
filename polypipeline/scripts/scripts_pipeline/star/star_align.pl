#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/GenBoDB";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages";
use lib "$Bin/../../../dragen/scripts";
use dragen_util; 
use Getopt::Long;
use Data::Dumper;
use GBuffer;
use file_util;
my $projectName;
my $version;
my $patient_name;
my $fork = 20;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patient_name,
	'version=s' =>\$version,
	'fork=s' =>\$fork,
);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName , -version =>$version);
my $dir = $project->getJunctionsDir("star");


#exit(0) if -e "$dir/$patient_name.SJ.tab";
my $patient = $project->getPatient($patient_name);
my $sjprod = "$dir/".$patient->name.".SJ.tab";
my $bam_prod = $patient->getBamFileName("star");

#exit(0) if -e $sjprod && -e $bam_prod;
#my $patient = $project->getPatient($patient_name);

my $dir_pipeline = $project->getAlignmentPipelineDir("star-".$patient->name);
my ($fastq1,$fastq2) = dragen_util::get_fastq_file($patient,$dir_pipeline);
my $ref_root =  $project->dirGenome()."/star/";

my $star = "/software/distrib/star/STAR-2.7.5a/bin/Linux_x86_64_static/STAR";
my $cmd = "$star --runThreadN $fork  --genomeDir $ref_root --readFilesIn $fastq1 $fastq2  --readFilesCommand zcat --outFileNamePrefix $dir_pipeline/$patient_name."." --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate ";
warn $cmd;
system($cmd);
unlink $fastq1;
unlink $fastq2;
my $fileout = "$dir_pipeline/$patient_name".".Aligned.sortedByCoord.out.bam";
if (-e $fileout) {
	system("samtools index   -@ $fork "."$fileout");
	warn $bam_prod;
}
else {
	die();
}
#111
my $sjfile= "$dir_pipeline/$patient_name".".SJ.out.tab";
die($sjfile) unless -e $sjfile;


system("mv $sjfile $sjprod");
my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");
system ("$bgzip -f $sjprod &&  $tabix -f -p bed $sjprod".".gz");
$sjprod = $sjprod.".gz";
warn $sjprod;
warn $sjfile;
die("$sjprod") unless -e "$sjprod";
die($bam_prod) unless -e $fileout.".bai";
exit(0);
