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

  
my $bin_cecile=qq{$Bin/scripts/scripts_db_polypipeline};
my $bin_script_pipeline = qq{$Bin/scripts/scripts_pipeline};


my $projectName;
my $patients_name;
my $step;
my $bcl_dir;
my $run;
my $feature_ref;
my $no_exec;
my $aggr_name;
my $lane;


my $limit;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'feature_ref=s' =>  \$feature_ref,
	'nb_lane=s' => \$lane,
	#'low_calling=s' => \$low_calling,
);





my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
my $run = $project->getRun();
warn $run->plateform_run_name;

die("can't find bcl directory : ".$run->bcl_dir()) unless -e $run->bcl_dir();
$bcl_dir = $run->bcl_dir();
my $fastq_dir =  $run->fastq_dir();

#$patients_name = "all" unless $patients_name;
my $patients = $project->get_only_list_patients($patients_name);

my $sampleSheet = $bcl_dir."/sampleSheet.csv";
my $exec = "cellranger";


open (SAMPLESHEET,">$sampleSheet");
my $dir = $project->getProjectRootPath();
print SAMPLESHEET "Lane,Sample,Index\n";

my $fastq;
my %hSamples;
my %test;
my $type;
my @group;

foreach my $patient (@{$patients}) {
	my $name=$patient->name();
	my $bc = $patient->barcode();
	my $group = $patient->somatic_group();
	push(@group,$group);
	for (my $i=1;  $i<=$lane;$i++){
		print SAMPLESHEET $i.",".$name.",".$bc."\n";
	}
}
close(SAMPLESHEET);


#my $prog =  $patient->alignmentMethod();
my $tmp = $project->getAlignmentPipelineDir("cellranger_count");

#my $cmd = "cd $tmp; /software/bin/demultiplex.pl -dir=$bcl_dir -run=$run -hiseq=10X -sample_sheet=$sampleSheet -cellranger_type=$exec";
warn $exec;
my $cmd = "cd $tmp; $Bin/../demultiplex/demultiplex.pl -dir=$bcl_dir -run=$run -hiseq=10X -sample_sheet=$sampleSheet -cellranger_type=$exec";
warn $cmd;
if ($step eq "demultiplex" or $step eq "all"){
	warn $cmd;
#	system $cmd or die "impossible $cmd" unless $no_exec==1;
}


my @jobs;
my $dir_pipeline =  $project->getCallingPipelineDir("cellranger");
foreach my $patient (@{$patients}) {
	my $name=$patient->name();
	my $bc = $patient->barcode();
	my $run = $patient->getRun();
	my $group = "EXP";
	$group = $patient->somatic_group() if ($patient->somatic_group());
	next() unless uc($group) eq "EXP" ;
	$type = $run->infosRun->{method};
	my $prog =  $patient->alignmentMethod();
	my $index = $project->getGenomeIndex($prog);
	$fastq = $patient->getSequencesDirectory();
	my $cmd = "cd $dir_pipeline; $exec count --id=$name --sample=$name --fastqs=$fastq --transcriptome=$index ";
	$cmd .= " --include-introns " if $type eq "nuclei";
	$cmd .= " :::40 ";
	push(@jobs,$cmd);
}	
my $file_jobs = "$dir/jobs_count.".time.".txt";
open (JOBS, ">$file_jobs");
print JOBS join("\n",@jobs);
close JOBS;

system("cat $file_jobs | run_cluster.pl -cpu=5 --limit=5 ");


