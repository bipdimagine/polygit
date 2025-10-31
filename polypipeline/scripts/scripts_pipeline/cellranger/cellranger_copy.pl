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
use Moo;
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
my $all_outs;
my $no_exec;

GetOptions(
	'project=s'		=> \$projectName,
	'patients=s'	=> \$patients_name,
	'all_outs!'		=> \$all_outs,
	'no_exec'		=> \$no_exec,
) || die("Error in command line arguments\n");

die("--project argument is mandatory") unless ($projectName);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
my $patients = $project->get_only_list_patients($patients_name);
my $run = $project->getRun();
my $type = $run->infosRun->{method};

my $dir = $project->getCountingDir('cellranger');
$dir = $project->getCountingDir('spaceranger') if ($type eq 'spatial');
#warn $dir;

my $dirout = "/data-pure/SingleCell/$projectName/";
unless (-d $dirout) {
	my $cp_cmd = "mkdir $dirout";
	warn $cp_cmd;
	system ($cp_cmd) unless ($no_exec);
}

foreach my $patient (@{$patients}) {
	my $name = $patient->name();
	next if ($name =~ /^ADT_/);
	my $cp_cmd;
	
	if (-e "$dir$name/web_summary.html") {
		# Copy ALL outs
		if ($all_outs){
			$cp_cmd = "cp -ru $dir/$name $dirout";
		}
		
		# Copy ONLY web_summary.html
		else {
			$cp_cmd = "mkdir $dirout$name ; " unless (-d $dirout.$name);
			$cp_cmd .= "cp $dir$name/web_summary.html $dirout$name/web_summary.html";
		}
		
		warn $cp_cmd;
		system ($cp_cmd) unless $no_exec;
	}
	else {
		warn ("$name web summary not found: '$dir$name/web_summary.html'\nAre you sure the counting finished successfully ?");
	}
}

unless ($no_exec){
	print "\t------------------------------------------\n";
	print "\tAll outs copied to $dirout \n" if ($all_outs);
	print "\tWeb summaries copied to $dirout \n" unless ($all_outs);
	print "\t------------------------------------------\n";
}
	
