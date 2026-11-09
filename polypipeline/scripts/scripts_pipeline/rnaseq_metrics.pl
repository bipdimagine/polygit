#!/usr/bin/env perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../GenBo/lib/GenBoDB";
use lib "$Bin/../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages";
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

my $filein = $patient->getAlignmentFile();
my $refFlat = $project->refFlat_file();
warn $refFlat;
unless ( -e $refFlat ) {
	die("can't find $refFlat\nSee https://open.bioqueue.org/home/knowledge/showKnowledge/sig/ucsc-gtftogenepred#");
	# /software/distrib/ucsc_util/gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons gencode.v46.annotation.gtf /dev/stdout | awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > refFlat.txt
}
my $rRNA_file = $project->rRNA_file();
my $opt = "";
$opt = "RIBOSOMAL_INTERVALS=$rRNA_file" if (-e $rRNA_file);

my $dir_out = $project->getCountingDir("featureCounts") . "/metrics/";
system("mkdir $dir_out --mode 777 $dir_out") unless (-e $dir_out);
my $fileout = $dir_out . "$patient_name.metrics";

my $picard = $project->buffer->software("picard");

my $cmd = "$picard CollectRnaSeqMetrics I=$filein O=$fileout REF_FLAT=$refFlat STRAND=FIRST_READ_TRANSCRIPTION_STRAND $opt ";
system($cmd);

die($fileout) unless (-e $fileout);

exit(0);


