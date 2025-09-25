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
use Time::Local 'timelocal';
use File::Temp qw/ tempfile tempdir /;
use Term::Menus;
use Proc::Simple;
use Storable;
use JSON::XS;
use XML::Simple qw(:strict);
use Cwd qw(cwd abs_path);


my $projectName;
my $patients_name;
my $out_dir = cwd;
my $help;

GetOptions(
	'project=s'		=> \$projectName,
#	'patients=s'	=> \$patients_name,
	'out_dir=s'		=> \ $out_dir,
) || die("Error in command line arguments\n");

usage() if $help;
usage() unless ($projectName);
die("'$out_dir' does not exist or is not a directory") unless (-d $out_dir);
$out_dir .= '/' unless ($out_dir =~ /\/$/);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
#$patients_name = $project->get_only_list_patients($patients_name);
$patients_name = $project->getPatients;
die("No patient in project $projectName") unless ($patients_name);



my $dir = $project->getCountingDir('cellranger');

my $info_file = "$out_dir$projectName-infos.txt";
open(my $fh, '>', $info_file);
print $fh $projectName.','.$project->description."\n";
foreach my $pat (@$patients_name) {
	my $pname = $pat->name;
	next if ($pname =~ /^ADT_/);
	
	# count directory
	my $dir_cellranger; 
	my $dir1 = $dir.$pat->name.'/';
	my $dir2 = $project->getProjectRootPath().$pat->name.'/';
	if (-d $dir1) {
		$dir_cellranger = $dir1;
		$dir_cellranger .= 'outs/' if (-d $dir1.'outs/');
	}
	elsif (-d $dir2) {
		$dir_cellranger = $dir2
	}
	else {
		die("Can't find directory with cellranger data for patient $pname")
	}
	my $transcriptome = $project->getGenomeIndex($pat->alignmentMethod);
	
	# cellranger version
	my $cellranger_version = "";
	my $file_version = $dir_cellranger.'_versions';
	open (my $fv, '<', $file_version) || warn ("Can't open $file_version: $!");
	while (my $line = readline($fv)) {
		chomp $line;
		$cellranger_version = $1 if ($line =~ /^\s*"pipelines": "(\w+ranger(-\w*)?-[\d\.]+)"$/);
	}
	
	print $fh $pname.','.$pat->getSequencesDirectory.','.$dir_cellranger.','.$cellranger_version."\n";
}
close($fh);

warn("Output file: ".abs_path($info_file)."\n");


sub usage {
	print "
$0
-------------
Mandatory:
	project <s>                project name
Optional:
	out_dir <s>                output directory
	help                       display this message

";
	exit(1);
}

