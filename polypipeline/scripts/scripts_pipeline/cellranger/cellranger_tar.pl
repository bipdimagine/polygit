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
my $dir_pipeline = $project->getCallingPipelineDir("cellranger");
	my $run_name = $run->plateform_run_name;
my $dir_prod = $project->getCellRangerRootDir();	
#	my $tar_cmd = "tar -cvzf $dir_pipeline/".$run_name.".tar.gz $dir_pipeline/*/outs/web_summary.html $dir_pipeline/*/outs/cloupe.cloupe $dir_pipeline/*/outs/vloupe.vloupe $dir_pipeline/*/outs/filtered_*_bc_matrix/* $dir_pipeline/*/outs/raw_*_bc_matrix/*";
	my $tar_cmd = "tar -cvzf $dir_pipeline/".$run_name.".tar.gz $dir_pipeline/*/outs/web_summary.html $dir_pipeline/*/outs/cloupe.cloupe $dir_pipeline/*/outs/*.*loupe $dir_pipeline/*/outs/filtered_*_bc_matrix/* $dir_pipeline/*/outs/raw_*_bc_matrix/*";
	$tar_cmd .= " && rsync $dir_pipeline/".$run_name.".tar.gz $dir_prod ";
	
	die ("archive $dir_pipeline/$run_name.tar.gz already exists") if -e $dir_pipeline."/".$run_name.".tar.gz";
	warn qq{echo "$tar_cmd" | run_cluster.pl -cpu=5};
	system (qq{echo "$tar_cmd" | run_cluster.pl -cpu=5});
	#  unless $no_exec==1;
	#or die "impossible $tar_cmd";
	die ("archive not found $dir_pipeline/$run_name.tar.gz already exists") unless -e $dir_pipeline."/".$run_name.".tar.gz";
	print "\t#########################################\n";
	print "\t  link to send to the users : \n";
	print "\t www.polyweb.fr/NGS/$projectName/$run.tar.gz \n";
	print "\t#########################################\n";