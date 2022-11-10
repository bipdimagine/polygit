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


my $projectsName;
my $patients_name;
my $step;
my $bcl_dir;
my $feature_ref;
my $no_exec;
my $aggr_name;
my $lane;


my $limit;
GetOptions(
	'project=s' => \$projectsName,
	'lane=s' => \$lane,
	#'patients=s' => \$patients_name,
);
my @aprojects  = split(",",$projectsName);
die("lane ?") unless $lane;
my $lines;
my $tmp_dir;
my $run_name;# =  $run->plateform_run_name;
my $fastq_dir;# =  $run->fastq_dir();
foreach my $projectName (@aprojects){
	my $buffer = GBuffer->new();
	my $project = $buffer->newProject( -name => $projectName );
	my $run = $project->getRun();
	$run_name =  $run->plateform_run_name unless $run_name;
	$fastq_dir =  $run->fastq_dir() unless $fastq_dir;
	$bcl_dir = $run->bcl_dir unless $bcl_dir;
	if ($bcl_dir && $run->bcl_dir() ne $bcl_dir){
		die("problem different bcl dir : $bcl_dir ".$run->bcl_dir);
	}
	
	foreach my $patient (@{$project->getPatients}) {
	my $name=$patient->name();
	my $bc = $patient->barcode();
	my $group = $patient->somatic_group();
	
	push(@$lines,$name.",".$bc."\n");
	}
}
my $sample_sheet = "$bcl_dir/samplesheet.".time.".csv";
open (SAS,">$sample_sheet");
print SAS "Lane,Sample,Index\n";
#print SAS qq{[Settings],,
#,,,
#[Data],,
#Lane,Sample,Index
#};



for (my $i=1;  $i<=$lane;$i++){
		foreach my $l (@$lines) {
			print SAS $i.",".$l;
		}
	}	
close(SAS);
my $id = $run_name."_".time;
my $cmd1= "/software/distrib/cellranger/latest mkfastq --id=$id --run=$bcl_dir --csv=$sample_sheet --barcode-mismatches=0 --output-dir=$fastq_dir";
system(qq{echo "$cmd1" | run_cluster.pl -cpu=40});

