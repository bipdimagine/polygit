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

my $spipeline;

my $limit;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'step=s'=> \$step,
	'type=s' => \$type,
	'umi=s' => \$umi,
	'command=s'=>\$spipeline,
	#'low_calling=s' => \$low_calling,
);

my $pipeline;
foreach my $l (split(",",$spipeline)){
	$pipeline->{$l} ++;
}
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

run_pipeline($pipeline);# unless -e $bam_pipeline;

#die() unless  -e $bam_pipeline;


exit(0);



sub run_pipeline {
my ($pipeline) = @_;
my $param_align;
my $ref_dragen = $project->getGenomeIndex("dragen");
	my $param_umi = "";
if (exists $pipeline->{align}){
my ($fastq1,$fastq2) = dragen_util::get_fastq_file($patient,$dir_pipeline);
	my $runid = $patient->getRun()->id;
	$param_align = " -1 $fastq1 -2 $fastq2 --RGID $runid  --RGSM $prefix --enable-map-align-output true ";

	if (exists $pipeline->{umi}){
		$param_align .= qq{ --umi-enable true   --umi-library-type random-simplex  --umi-min-supporting-reads 1 --vc-enable-umi-germline true};
	}
	else {
		$param_align .= qq{ --enable-duplicate-marking true };
	 }
}
else {
	my $bam = $patient->getBamFile();
	$param_align = qq{ --bam-input $bam --enable-map-align false };
}
my $param_gvcf = "";
my $tmp = "/staging/tmp";
	my $capture_file  = $patient->getCapture->gzFileName();
my $cmd_dragen = qq{dragen -f -r $ref_dragen --output-directory $dir_pipeline --intermediate-results-dir $tmp --output-file-prefix $prefix };
if (exists $pipeline->{gvcf}){
	
	$param_gvcf = qq{--vc-emit-ref-confidence GVCF --enable-variant-caller true } ;
	unless ($project->isGenome) {
	
		$param_gvcf .= qq{ --vc-target-bed $capture_file --vc-target-bed-padding 150 };
	}
}

my $param_cnv = "";
if (exists $pipeline->{cnv}){
	$param_cnv = qq{ --enable-cnv true --cnv-enable-self-normalization true };
	unless ($project->isGenome) {
	 $param_cnv .= " --cnv-target-bed $capture_file ";
	}
}

my $param_sv = "";
if (exists $pipeline->{sv}){
	$param_sv = qq{   --enable-sv true };
	unless ($project->isGenome) {
	 $param_sv .= " --sv-exome true ";
	}
}


$cmd_dragen .= $param_umi." ".$param_align." ".$param_gvcf." ".$param_cnv." ".$param_sv;

my $exit = system(qq{$Bin/../run_dragen.pl -cmd=\"$cmd_dragen\"}) ;#unless -e $f1;

#system("ssh pnitschk\@10.200.27.109 ". $cmd." >$dir_pipeline/dragen.stdout 2>$dir_pipeline/dragen.stderr");
#system("ssh pnitschk\@10.200.27.109 rm $fastq1 $fastq2");
#my ($out,$err,$exit) = $ssh->cmd("$cmd") ;#unless -e $bam_pipeline;
die if $exit != 0;
}




1;



