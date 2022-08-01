#!/usr/bin/perl
 
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/";
use lib "$Bin/../../../GenBo/lib/GenBoDB";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../../../GenBo/lib/kyoto/";
use lib "$Bin/../../../GenBo/lib/GenBoDB/writeDB";
use lib "$Bin/../packages";
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

use File::Temp qw/ tempfile tempdir /;; 

use Term::Menus;
 use Proc::Simple;
 use Storable;
use JSON::XS;
use Net::SSH::Perl; 

 
 my $ssh = Net::SSH::Perl->new("10.200.27.109");
$ssh->login("pnitschk");

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
 my $myproc;
my $low_calling;
my $predef_type;
my $define_steps;
my $step;

my $dir_dragen = "/run/";


my $limit;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'step=s'=> \$step,
	'type=s' => \$type,
	#'low_calling=s' => \$low_calling,
);


my $user = system("whoami");

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );





#system ("mkdir -p $dir_dragen/".$project->name );

my $output = "/data-isilon/dragen-pipeline/".$project->name."/cnv1/";
system("mkdir -p $output");
warn $output;
my $patient = $project->getPatient($patients_name);
my $ref_dragen = "/data-isilon/public-data/genome/HG19_MT/dragen/all.fa.k_21.f_16.m_149 ";

my $prefix = $patient->name;
my $bam = $patient->getBamFile();
my $tmp = "/staging/intermediate";
my $capture_file = $patient->getCapture->gzFileName();
my $cmd = qq{dragen -f  -r $ref_dragen  --output-directory $output --output-file-prefix $prefix -b $bam --intermediate-results-dir $tmp --enable-map-align false --enable-cnv true --cnv-enable-self-normalization true --cnv-target-bed $capture_file};
my ($out,$err,$exit) = $ssh->cmd("$cmd");

warn "OK";






