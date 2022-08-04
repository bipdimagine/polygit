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
 my $ssh = Net::SSH::Perl->new("10.1.2.9");
$ssh->login("$username");
my $project_name;
my $patient_name;
GetOptions(
	'project=s' => \$project_name,
	'patients=s' => \$patient_name,
);
my $tmp = "/staging/tmp";
my $buffer = GBuffer->new();
my $project = $buffer->newProjectCache( -name 			=> $project_name );
my $patient = $project->getPatient($patient_name);
my $run = $patient->getRun();
my $hps =  $run->getAllPatientsInfos();
my %contr_projects;

map {$contr_projects{$_->{project}} ++} @$hps;
my $dir_pipeline_bam = $patient->getDragenDir("pipeline");
my $dir_pipeline = $patient->getDragenDir("target");
my $ref_dragen = $project->getGenomeIndex("dragen");

my $f1= $dir_pipeline."/".$patient_name.".sv.vcf.gz";
my $bam = $patient->getBamFile();
die() unless -e $bam;
my $bamin = "$dir_pipeline_bam/$patient_name.bam";
system("rsync -rav $bam*  $dir_pipeline_bam/") unless -e $f1;


my $capture_file  = $patient->getCapture->gzFileName();
my $prefix = $patient->name;

# --cnv-enable-self-normalization true

my $cmd = qq{dragen -f  -r $ref_dragen  --output-directory $dir_pipeline --output-file-prefix $prefix -b $bamin --intermediate-results-dir $tmp --enable-map-align false --enable-cnv true  --cnv-target-bed $capture_file};
if ($project->isGenome){
	$cmd = qq{dragen -f  -r $ref_dragen  --output-directory $dir_pipeline --output-file-prefix $prefix -b $bamin --intermediate-results-dir $tmp --enable-map-align false --enable-cnv true };
}
my $exit = system(qq{$Bin/../run_dragen.pl -cmd=\"$cmd\"}) ;
die if $exit != 0;

my $target_pipeline_file  = "$dir_pipeline/".$prefix.".target.counts.gz";
my $target_pipeline_gc_file  = "$dir_pipeline/".$prefix.".target.counts.gc-corrected.gz";

my $target_file = $patient->targetFile();
my $targetGC_file = $patient->targetGCFile();

die unless -e $target_pipeline_gc_file;
system("rsync -rav $target_pipeline_file $target_file 2>/dev/null");
system("rsync -rav $target_pipeline_gc_file $targetGC_file 2>/dev/null");

exit(0);





