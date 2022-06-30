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

my $buffer = GBuffer->new();
my $project = $buffer->newProjectCache( -name 			=> $project_name );
my $patient = $project->getPatient($patient_name);
my $run = $patient->getRun();
my $hps =  $run->getAllPatientsInfos();
my %contr_projects;

map {$contr_projects{$_->{project}} ++} @$hps;

my $dir_pipeline = $patient->getDragenDir("SV");
my $ref_dragen = $project->getGenomeIndex("dragen");
my $dir_out= $patient->getCallingPipelineDir("dragen-sv");

my $bam = $patient->getBamFile("dragen-align");
my $cmd = qq{dragen -f -r $ref_dragen --output-directory $dir_pipeline --bam-input $bam --output-file-prefix $patient_name --enable-map-align false  --enable-sv true  };
my $exit = system(qq{$Bin/../run_dragen.pl -cmd="$cmd"});
warn qq{$Bin/../run_dragen.pl -cmd="$cmd"};
#my ($out,$err,$exit) = $ssh->cmd($cmd);
die($cmd) unless $exit == 0;
my $f1= $dir_pipeline."/".$patient_name.".cnv.vcf.gz";
my $dir_prod = $project->getVariationsDir("dragen-pon");
system("rsync -rav $f1 $dir_prod 2>/dev/null");

exit(0);

