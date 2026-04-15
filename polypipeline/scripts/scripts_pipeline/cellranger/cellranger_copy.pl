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
use File::Path qw(make_path);
use Carp;


my $projectName;
my $patients_name;
my $all_outs;
my $no_exec;
my $verbose;

GetOptions(
	'project=s'		=> \$projectName,
	'patients=s'	=> \$patients_name,
	'all_outs!'		=> \$all_outs,
	'no_exec'		=> \$no_exec,
	'verbose'		=> \$verbose,
) || confess("Error in command line arguments\n");

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
make_path($dirout, { mode => 0775 }) unless (-d $dirout);
my $error;
foreach my $patient (@{$patients}) {
	my $name = $patient->name();
	my $group = $patient->somatic_group;
	my $pool = $patient->family;
	next if ($group =~ /adt/i or $name =~ /^ADT_/);
	my $cp_cmd;
	
	if (-e "$dir$name/web_summary.html") {
		# Copy ALL outs
		if ($all_outs){
			$cp_cmd = "cp -ru $dir/$name $dirout";
		}
		
		# Copy ONLY web_summary.html
		else {
			make_path("$dirout$name", { mode => 0775 }) unless (-d $dirout.$name);
			$cp_cmd = "cp $dir$name/web_summary.html $dirout$name/web_summary.html";
		}
		
		warn $cp_cmd if ($verbose);
		my $exit = system ($cp_cmd) unless ($no_exec);
		$error->{$name} ++ if ($exit);
	}
	
	# multi outs (cellranger_multi)
	elsif (-e "$dir$pool/$name/web_summary.html" or -e "$dir$pool/qc_report.html") {
		next if ($group =~ /vdj/i);
		make_path("$dirout$pool/per_sample_outs/", { mode => 0775 }) unless (-d $dirout.$pool.'/per_sample_outs/');
		# Copy ALL outs
		if ($all_outs){
			$cp_cmd = "cp -ru $dir$pool/$name $dirout$pool/per_sample_outs/";
			$cp_cmd .= " && cp -u $dir$pool/qc_report.html $dirout$pool/" if (-f "$dir$pool/qc_report.html");
		}
		
		# Copy ONLY web_summary.html
		else {
			make_path("$dirout$pool/per_sample_outs/$name", { mode => 0775 }) unless (-d $dirout.$pool.'/per_sample_outs/$name');
			$cp_cmd = "cp $dir$pool/$name/web_summary.html $dirout$pool/per_sample_outs/$name/web_summary.html";
			$cp_cmd .= " && cp -u $dir$pool/qc_report.html $dirout$pool/" if (-f "$dir$pool/qc_report.html");
		}
		
		warn $cp_cmd if ($verbose);
		my $exit = system ($cp_cmd) unless ($no_exec);
		$error->{$name} ++ if ($exit);
	}

	# multi outs (cellranger_fixed2)
	elsif (-e "$dir$group/$name/web_summary.html" or -e "$dir$group/qc_report.html") {
		make_path("$dirout$group/per_sample_outs/", { mode => 0775 }) unless (-d $dirout.$group.'/per_sample_outs/');
		# Copy ALL outs
		if ($all_outs){
			$cp_cmd = "cp -ru $dir$group/$name $dirout$group/per_sample_outs/";
			$cp_cmd .= " && cp -u $dir$group/qc_report.html $dirout$group/" if (-f "$dir$group/qc_report.html");
		}
		
		# Copy ONLY web_summary.html
		else {
			make_path("$dirout$group/per_sample_outs/$name", { mode => 0775 }) unless (-d $dirout.$group.'/per_sample_outs/$name');
			$cp_cmd = "cp $dir$group/$name/web_summary.html $dirout$group/per_sample_outs/$name/web_summary.html";
			$cp_cmd .= " && cp -u $dir$group/qc_report.html $dirout$group/" if (-f "$dir$group/qc_report.html");
		}
		
		warn $cp_cmd if ($verbose);
		my $exit = system ($cp_cmd) unless ($no_exec);
		$error->{$name} ++ if ($exit);
	}
	
	else {
		$error->{$name} ++;
		warn ("$name web summary not found: '$dir$name/web_summary.html' or '$dir$group/$name/web_summary.html'\nAre you sure the counting finished successfully ?");
	}
}

if ($error){
	die("ERROR: Patient(s) not copied to $dirout:\n".Dumper $error);
}
unless ($no_exec or $error){
	print "\t------------------------------------------\n";
	print "\tAll outs copied to $dirout \n" if ($all_outs);
	print "\tWeb summaries copied to $dirout \n" unless ($all_outs);
	print "\t------------------------------------------\n";
}
	
