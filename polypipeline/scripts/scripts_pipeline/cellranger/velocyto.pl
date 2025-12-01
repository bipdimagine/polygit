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

  
my $bin_cecile=qq{$Bin/scripts/scripts_db_polypipeline};
my $bin_script_pipeline = qq{$Bin/scripts/scripts_pipeline};


my $projectName;
my $patients_name;
my $no_exec;
my $cpu = 20;
my $force;
my $help;

GetOptions(
	'project=s'		=> \$projectName,
	'patients=s'	=> \$patients_name,
	'cpu=i'			=> \$cpu,
	'no_exec'		=> \$no_exec,
	'force'			=> \$force,
	'help'			=> \$help,
) || die("Error in command line arguments\n");

#usage() if $help;
#usage() unless ($projectName);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
$patients_name = $project->get_only_list_patients($patients_name);
die("No patient in project $projectName") unless $patients_name;

my $dir_proj = $project->getProjectRootPath;
my $dir_cellranger = $project->getCountingDir('cellranger');
my $tmp_cellranger = $project->getAlignmentPipelineDir("cellranger");
#my $tmp_cellranger = $project->getAlignmentPipelineDir("velocyto");
my $dir_velocyto = $project->getCountingDir('velocyto');

my $file_jobs = $dir_velocyto.'jobs_velocyto.txt';
open(my $jobs, ">$file_jobs") || die ("Can't open file '$file_jobs': $!");
foreach my $pat (@$patients_name) {
	my $pname = $pat->name;
	my $index = $project->getGenomeIndex($pat->alignmentMethod);
	my $gtf = $project->getGenomeIndex($pat->alignmentMethod).'/genes/genes.gtf';
	my $dir_pat;
	if ((-e "$dir_velocyto$pname.loom" or -e "$tmp_cellranger/$pname/velocyto/$pname.loom") and not $force) {
		warn "NEXT: '$pname.loom' already exists";
		next;
	}
	my $cmd;
	unless (-d "$tmp_cellranger$pname/outs/") {
		$dir_pat = "$dir_proj/$pname" if (-d "$dir_proj/$pname");
		$dir_pat = "$dir_cellranger/$pname" if (-d "$dir_cellranger/$pname");
		die ("No directory found for patient '$pname' in project ".$project->name) unless ($dir_pat);
		$cmd = "cp -r $dir_pat $tmp_cellranger && " if (-d "$dir_pat/outs/");
		$cmd = "mkdir $tmp_cellranger/$pname/outs && cp -r $dir_pat/* $tmp_cellranger/$pname/outs/ && " unless (-d "$dir_pat/outs/");
	}
	$cmd .= "singularity run -B $tmp_cellranger$pname -B $index /data-beegfs/software/sif/velocyto.sif velocyto run10x -v -@ $cpu $tmp_cellranger$pname $index/genes/genes.gtf";
	$cmd .= " && cp $tmp_cellranger/$pname/velocyto/$pname.loom $dir_velocyto";
	print {$jobs} $cmd."\n";
	warn $cmd if ($no_exec);
	# Execution time is ~3h for a typical sample but might vary significantly by sequencing depth and cpu power.
	# ~6h avec ce script cpu 20
}
warn ("cat $file_jobs | run_cluster.pl -cpu=$cpu");
system ("cat $file_jobs | run_cluster.pl -cpu=$cpu") unless ($no_exec);







sub usage {
	print "
$0
-------------
Obligatoires:
	project <s>                nom du projet
Optionels:
	patients <s>               noms de patients/échantillons, séparés par des virgules
	cpu <i>                    nombre de cpu à utiliser, défaut: 20
	no_exec                    ne pas exécuter les commandes
	force                      relance le pipeline même s'il a déjà tourné
	help                       affiche ce message

";
	exit(1);
}


