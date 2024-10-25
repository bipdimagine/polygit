#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/";
use lib "$Bin/../../../GenBo/lib/GenBoDB";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../../../GenBo/lib/kyoto/";
use lib "$Bin/../../../GenBo/lib/GenBoDB/writeDB";
use lib "$Bin/../../packages";
use lib "$Bin";
use dragen_util; 
use Logfile::Rotate;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
use GBuffer;
use GenBoProject;
use colored; 
use Config::Std;
use Text::Table;
use file_util;
use File::Temp qw/ tempfile tempdir /;; 

use Term::Menus;
use Proc::Simple;
use Storable;
use JSON::XS;
use Net::SSH::Perl; 

 

my $bin_cecile=qq{$Bin/scripts/scripts_db_polypipeline};
my $bin_script_pipeline = qq{$Bin/scripts/scripts_pipeline};


my $projectName;
my $filename;

my $name;
my $patients_name;
my $steps_name;
my $force;
my $type ="";
my $fastq_ext;
my $exclude_patients;
my $max_cpu ;
my $bds;
 my @running_steps;
 my $predef_steps;
 my $nocluster = 0;
my $low_calling;
my $predef_type;
my $define_steps;
my $step;
my $umi;
my $lpad;


my $limit;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'step=s'=> \$step,
	'type=s' => \$type,
	'umi=s' => \$umi,
	'padding=s' =>\$lpad,
	#'low_calling=s' => \$low_calling,
);


my $user = system("whoami");
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );

my $tm = "/staging/tmp/";

if ($project->isGenome){
	warn "coucou";
}
#system ("mkdir -p $dir_dragen/".$project->name );

my $patient = $project->getPatient($patients_name);
my $dir_pipeline = $patient->getDragenDir("pipeline");
my $prefix = $patient->name;
my $bam_prod = $patient->gvcfFileName("dragen-calling");
#warn $bam_prod if -e $bam_prod;
#exit(0) if -e $bam_prod;

my $bam_pipeline = $dir_pipeline."/".$prefix.".bam";

run_pipeline($patient) unless -e $bam_pipeline;

die() unless  -e $bam_pipeline;


exit(0);



sub run_pipeline {
my ($fastq1,$fastq2) = dragen_util::get_fastq_file($patient,$dir_pipeline);
warn "\n\n end copy";
my $ref_dragen = $project->getGenomeIndex("dragen");


my $tmp = "/staging/tmp";

my $capture_file  = dragen_util::get_capture_file($patient,$dir_pipeline."/".$patient->name.".".time.".bed");


my $runid = $patient->getRun()->id;

my $gvcf_pipeline = $dir_pipeline."/".$prefix.".hard-filtered.gvcf.gz";


#--vc-enable-vcf-output true
my $cmd = qq{dragen -f -r $ref_dragen --intermediate-results-dir $tmp --output-directory $dir_pipeline --output-file-prefix $prefix -1 $fastq1 -2 $fastq2 --RGID $runid  --RGSM $prefix  --vc-emit-ref-confidence GVCF --enable-variant-caller true   --enable-map-align-output true  --vc-target-bed $capture_file  --enable-cnv true --cnv-enable-self-normalization true --cnv-target-bed $capture_file};
#warn $cmd;
#die();
$lpad = 150 unless $lpad;
my $padding = " --vc-target-bed-padding $lpad ";
$cmd .= qq{ $padding};


if ($project->isGenome){
	$cmd = qq{dragen -f -r $ref_dragen --intermediate-results-dir $tmp --output-directory $dir_pipeline --output-file-prefix $prefix -1 $fastq1 -2 $fastq2 --RGID $runid  --RGSM $prefix  --vc-emit-ref-confidence GVCF --enable-variant-caller true   --enable-map-align-output true   --enable-cnv true --cnv-enable-self-normalization true };
	
}
if ($umi){
	$cmd .= qq{ --umi-enable true   --umi-library-type random-simplex  --umi-min-supporting-reads 1 --vc-enable-umi-germline true};
}
else {
	$cmd .= qq{  --enable-duplicate-marking true };
	 }
my $exit = system(qq{$Bin/../run_dragen.pl -cmd=\"$cmd\"}) ;#unless -e $f1;

#system("ssh pnitschk\@10.200.27.109 ". $cmd." >$dir_pipeline/dragen.stdout 2>$dir_pipeline/dragen.stderr");
#system("ssh pnitschk\@10.200.27.109 rm $fastq1 $fastq2");
#my ($out,$err,$exit) = $ssh->cmd("$cmd") ;#unless -e $bam_pipeline;
die if $exit != 0;
}




1;



