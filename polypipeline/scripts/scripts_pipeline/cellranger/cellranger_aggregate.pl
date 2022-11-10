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
my $exec = "cellranger";
my $aggr_file = $dir_pipeline."/jobs_aggr".time.".txt";
	my $id = $projectName;
	$id = $aggr_name if $aggr_name;
	open (JOBS_AGGR, ">$aggr_file");
	print JOBS_AGGR "sample_id,molecule_h5\n";
	foreach my $patient (@{$patients}) {
		die($dir_pipeline."/".$patient->name()."/"."outs/molecule_info.h5") unless -e $dir_pipeline."/".$patient->name()."/"."outs/molecule_info.h5";
		print JOBS_AGGR $patient->name().",".$dir_pipeline."/".$patient->name()."/"."outs/molecule_info.h5\n";
	}
	close JOBS_AGGR;
	my $aggr_cmd = qq{ echo "cd $dir_pipeline; $exec aggr --id=$id --csv=$aggr_file" | run_cluster.pl -cpu=20};
	system ($aggr_cmd)  unless $no_exec==1;
