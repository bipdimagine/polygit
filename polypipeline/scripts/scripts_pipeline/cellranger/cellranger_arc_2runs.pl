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
my $feature_ref;
my $no_exec;
my $aggr_name;
my $cmo_ref;
my $lane;


my $limit;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'step=s'=> \$step,
#	'bcl=s' => \$bcl_dir,
#	'run=s' => \$run,
	'feature_ref=s' =>  \$feature_ref,
	'no_exec=s' => \$no_exec,
	'aggr_name=s' => \$aggr_name,
	'nb_lane=s' => \$lane,
	#'low_calling=s' => \$low_calling,
);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
#$patients_name = "all" unless $patients_name;
$patients_name = $project->get_only_list_patients($patients_name);
warn "toto";

my $exec= "cellranger";
my @group;
my $dir = $project->getProjectRootPath();

foreach my $patient (@{$patients_name}) {
	my $name=$patient->name();
	warn $name;
	my $bc = $patient->barcode();
	my $group = $patient->somatic_group();
	warn $group;
	push(@group,$group);
	for (my $i=1;  $i<=$lane;$i++){
		print SAMPLESHEET $i.",".$name.",".$bc."\n";
	}
}
close(SAMPLESHEET);

foreach my $patient (@{$patients_name}) {
	warn $patient;
#	my $run = $patient->getRun();
#	my $type = $run->infosRun->{method};
	my $type_arc = 1 if map {uc($_) =~ /ARC/ } @group;
	if ($type_arc ){
		open (JOBS_ARC, ">$dir/jobs_arc.txt");
		$exec = "cellranger-arc";
		my @exp = grep { uc($_->somatic_group()) eq "EXP"} @{$patients_name};
		foreach my $e(@exp){
			my $ename = $e->name(); 
			warn $ename;
			my $efam = $e->family();
			warn $efam;
			my $egroup = uc($e->somatic_group());
			my $lib_file = $dir."/".$ename."_library.csv";
			open(LIB,">$lib_file") or die "impossible $lib_file";
			my @atac = grep {$_->family() eq $efam && uc($_->somatic_group()) eq "ARC"} @{$patients_name};
			warn Dumper(@atac);
			next() if scalar(@atac == 0);
			my $atac_name = $atac[0]->name() ;
			my $prog =  $e->alignmentMethod();
			my $index = $project->getGenomeIndex($prog)."_arc";
			my $lib = "fastqs,sample,library_type\n".$e->getSequencesDirectory().",".$ename.",Gene Expression\n";
			$lib .= $atac[0]->getSequencesDirectory().",".$atac_name.",Chromatin Accessibility\n";
			print LIB $lib;
			close(LIB);
			print JOBS_ARC "cd $dir ; $exec count --id=$ename --reference=$index  --libraries=$lib_file\n"
		}
	}
}
close(JOBS_ARC);