#!/usr/bin/env perl

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
use autodie qw(system);

my $projectName;
my $version;
my $patient_name;
my $fork = 20;
GetOptions(
	'project=s'		=> \$projectName,
	'patients=s'	=> \$patient_name,
	'version=s'		=>\$version,
	'fork=s'		=>\$fork,
);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName , -version =>$version);
my $patient = $project->getPatient($patient_name);

my $dir_pipeline = $project->getAlignmentPipelineDir("star-".$patient->name);
my ($fastq1,$fastq2) = dragen_util::get_fastq_file($patient,$dir_pipeline);
my $ref_root =  $project->dirGenome()."/star/";
my $star = "exec_singularity.sh star.sif STAR";

my $cmd = "ulimit -n 8192 && $star --runThreadN $fork --genomeDir $ref_root --readFilesIn $fastq1 $fastq2  --readFilesCommand zcat  --twopassMode Basic --outFileNamePrefix $dir_pipeline$patient_name."
	." --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --bamRemoveDuplicatesType UniqueIdentical "; # --outBAMsortingThreadN 16
system($cmd);

unlink $fastq1;
unlink $fastq2;

my $fileout = "$dir_pipeline/$patient_name.Aligned.sortedByCoord.out.bam";
die($fileout) unless (-e $fileout);

my $sjfile = "$dir_pipeline/$patient_name".".SJ.out.tab";
die($sjfile) unless -e $sjfile;

my $dir = $project->getJunctionsDir("star");
my $sjprod = "$dir/".$patient->name.".SJ.tab";
system("mv $sjfile $sjprod");
my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");
system ("$bgzip -f $sjprod &&  $tabix -f -p bed $sjprod".".gz");
$sjprod = $sjprod.".gz";
die("$sjprod") unless -e "$sjprod";

die($fileout) unless -e $fileout;
my $bam_prod = $patient->getBamFileName("star");
system("mv $fileout $bam_prod");
system("samtools index -@ $fork $bam_prod");
die($bam_prod) unless -e $bam_prod;

exit(0);
