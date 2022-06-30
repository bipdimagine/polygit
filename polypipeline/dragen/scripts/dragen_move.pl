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


my $limit;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'step=s'=> \$step,
	'type=s' => \$type,
	#'low_calling=s' => \$low_calling,
);
my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
 my $ssh = Net::SSH::Perl->new("10.1.2.9");
$ssh->login("$username");

#my $user = system("whoami");
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );

my $tm = "/staging/tmp/";

#system ("mkdir -p $dir_dragen/".$project->name );

my $patient = $project->getPatient($patients_name);
my $dir_pipeline = $patient->getDragenDirName("pipeline");
my $cmd_dir = qq{test -d $dir_pipeline || mkdir -p $dir_pipeline};
my ($out, $err, $exit) = $ssh->cmd($cmd_dir);
#my $dir_pipeline = $patient->getDragenDirName("pipeline");
my $prefix = $patient->name;
my $bam_prod = $patient->getBamFileName("dragen-align");
#exit(0) if -e $bam_prod;
#warn "coucou";
my $bam_pipeline = $dir_pipeline."/".$prefix.".bam";
($out, $err, $exit)=  $ssh->cmd("test -f $bam_pipeline");

my $gvcf_pipeline = "$dir_pipeline/".$prefix.".hard-filtered.gvcf.gz";
my $target_pipeline  ="$dir_pipeline/".$prefix.".target.counts.gz";
my $target_pipeline_gc  = "$dir_pipeline/".$prefix.".target.counts.gc-corrected.gz";
($out, $err, $exit)=  $ssh->cmd("test -f $gvcf_pipeline");
die($gvcf_pipeline ." probleme no gvcf") unless  $exit ==0;

($out, $err, $exit)=  $ssh->cmd("test -f $target_pipeline_gc");
die($target_pipeline_gc ." probleme no target gc") unless  $exit ==0;
#die() unless -e $target_pipeline;
my $url = qq{$username\@10.200.27.109};
move_bam($bam_pipeline,$patient);
move_gvcf($gvcf_pipeline,$patient);
move_target($target_pipeline,$target_pipeline_gc,$patient);
if($project->isGenome){
	my $cnv_file  = "$dir_pipeline/".$prefix.".cnv.vcf.gz";
	move_cnv($cnv_file,$patient)
}
exit(0);






sub move_bam {
	my ($bam,$patient) = @_;
	my $prod = $patient->getBamFileName("dragen-align");
	system("rsync -rav  $url".":$bam $prod ");
	system("rsync -rav  $url".":$bam.bai $prod.bai ");
	
}

sub move_gvcf {
	my ($gvcf,$patient) = @_;
	my $prod = $patient->gvcfFileName("dragen-calling");
	system("rsync -rav  $url:$gvcf $prod");
	system("rsync -rav  $url:$gvcf.tbi $prod.tbi");
}

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



