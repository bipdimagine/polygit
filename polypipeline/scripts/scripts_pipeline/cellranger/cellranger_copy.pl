#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/";
use lib "$Bin/../../../../GenBo/lib/GenBoDB";
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages";
use Logfile::Rotate;
 use Cwd;
use PBS::Client;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
use MooseX::Method::Signatures;
#use bds_steps;   
use file_util;
use Class::Inspector;
use Digest::MD5::File ;
use GBuffer;
use GenBoProject;
use colored; 
use Config::Std;
use Text::Table;

use File::Temp qw/ tempfile tempdir /;

use Term::Menus;
 use Proc::Simple;
 use Storable;
use JSON::XS;


my $projectName;
my $patients_name;
my $step;
my $bcl_dir;
my $feature_ref;
my $no_exec;
my $aggr_name;
my $lane;
my $limit;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
);
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
my $patients = $project->get_only_list_patients($patients_name);
my $run = $project->getRun();
my $dir = $project->getCallingPipelineDir("cellranger");
my $dir_prod = $project->getCellRangerRootDir();

foreach my $patient (@{$patients}) {
		my $name = $patient->name();
		my $run_name = $run->plateform_run_name;
		my $dirout = "/data-isilon/SingleCell/$projectName";
		my $cp_cmd = "mkdir $dirout ; rsync -rav  $dir/$name $dirout/ && rsync -rav  $dir/$name $dir_prod/ ; ";
		system(qq{echo "$cp_cmd" | run_cluster.pl -cpu=5 });
		#die ("archive $dir/$run_name.tar.gz already exists") if -e $dir."/".$run.".tar.gz";
		system ($cp_cmd)  unless $no_exec==1;
		#or die "impossible $tar_cmd";
		print "\t#########################################\n";
		print "\t  cp to /data-isilon/SingleCell/$projectName  \n";
		print "\t#########################################\n";
	}