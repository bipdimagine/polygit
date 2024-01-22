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

 
my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
 my $url = qq{$username\@10.200.27.109};
my $project_name;
my $patient_name;
GetOptions(
	'project=s' => \$project_name,
	'patients=s' => \$patient_name,
);

my $buffer = GBuffer->new();
my $project = $buffer->newProjectCache( -name 			=> $project_name );
my $patient = $project->getPatient($patient_name);
my $run = $patient->getRun();
my $hps =  $run->getAllPatientsInfos();
my %contr_projects;

map {$contr_projects{$_->{project}} ++} @$hps;
my $dir_pipeline_bam = $patient->getDragenDir("pipeline");
 my $dir_pipeline = $patient->getDragenDir("pipeline");
my $ref_dragen = $project->getGenomeIndex("dragen");
my $f1= $dir_pipeline."/".$patient_name.".sv.vcf.gz";
my $bam = $patient->getBamFileName();
 my $bamin = "$dir_pipeline_bam/$patient_name.bam";
unless (-e $bam){
	die($bamin) unless -e $bamin;
}
else {
system("rsync -rav $bam*  $dir_pipeline_bam/") unless -e $f1;
}
my $capture_file  = $patient->getCapture->gzFileName();
my $cmd = qq{dragen -f -r $ref_dragen --output-directory $dir_pipeline --bam-input $bamin --output-file-prefix $patient_name --enable-map-align false --vc-emit-ref-confidence GVCF --enable-variant-caller true  --vc-target-bed $capture_file --vc-target-bed-padding 150 --enable-cnv true --cnv-enable-self-normalization true };
if ($project->isGenome){
	my $cmd = qq{dragen -f -r $ref_dragen --output-directory $dir_pipeline --bam-input $bamin --output-file-prefix $patient_name --enable-map-align false --vc-emit-ref-confidence GVCF --enable-variant-caller true --enable-cnv true --cnv-enable-self-normalization true };
}

my $exit =0;
my $gvcf_pipeline = "$dir_pipeline/".$patient_name.".hard-filtered.gvcf.gz";
$exit = system(qq{$Bin/../run_dragen.pl -cmd=\"$cmd\"}) unless -e $gvcf_pipeline;

warn qq{$Bin/../run_dragen.pl -cmd=\"$cmd"};

die($cmd) unless $exit == 0;
move_gvcf($gvcf_pipeline,$patient);

exit(0);

sub move_target {
	my ($t1,$t2,$patient) = @_;
	my $dir = $patient->project->getTargetCountDir();
	system("rsync -rav  $url:$t1 $dir/");
	system("rsync -rav  $url:$t2 $dir/");
}


sub move_cnv {
	my ($t1,$patient) = @_;
	my $dir = $patient->project->getVariationsDir("dragen-cnv");
	system("rsync -rav  $url:$t1 $dir/");
	system("rsync -rav  $url:$t1.tbi $dir/");
}


sub move_gvcf {
	my ($gvcf,$patient) = @_;
	my $prod = $patient->gvcfFileName("dragen-calling");
	system("rsync -rav  $url:$gvcf $prod");
	system("rsync -rav  $url:$gvcf.tbi $prod.tbi");
}
