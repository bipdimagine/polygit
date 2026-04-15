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
use Carp;


my $projectName;
my $patients_name;
my $create_bam;
my $multi;
my $all_vdj;
my $no_exec;

GetOptions(
	'project=s'		=> \$projectName,
	'patients=s'	=> \$patients_name,
	'create_bam!'	=> \$create_bam,
	'multi|flex'	=> \$multi,
	'all_vdj'		=> \$multi,
	'no_exec'		=> \$no_exec,
) || confess("Error in command line arguments\n");

die("--project argument is mandatory") unless ($projectName);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
my $patients = $project->getPatients($patients_name);
my $run = $project->getRun();
my $type = $run->infosRun->{method};

my $dir = $project->getCountingDir('cellranger');
$dir = $project->getCountingDir('spaceranger') if ($type eq 'spatial');
warn $dir;

my @groups = map {$_->somatic_group()} @$patients;
my @pools = map {$_->family()} @$patients;;

my $tar_cmd = "cd $dir && tar -cvzf $projectName.tar.gz ";
if ($type eq 'spatial') {
	$tar_cmd .= join(' ',map {$_->name} @$patients).' ';
}
elsif ($multi) {
	$tar_cmd .= "*/*/web_summary.html  */qc_report.html ";
	$tar_cmd .= "*/*/sample_cloupe.cloupe */*/sample_*_bc_matrix.h5 " if (grep (! /vdj/i, @groups));
	$tar_cmd .= "*/*/vdj_*/vloupe.vloupe " if (grep (/vdj/i, @groups) and not $all_vdj);
	$tar_cmd .= "*/*/vdj_* " if ($all_vdj and grep (/vdj/i, @groups));
	$tar_cmd .= "*/*/possorted_bam.bam* " if ($create_bam and $create_bam ne 'false');
}
else {
	$tar_cmd .= "*/web_summary.html ";
	$tar_cmd .= "*/cloupe.cloupe */*_bc_matrix.h5 " if (grep (! /vdj/i, @groups));
	$tar_cmd .= "VDJ_*/* " if ($all_vdj and grep (/vdj/i, @groups));
	$tar_cmd .= "*/vloupe.vloupe " if grep (/vdj/i, @groups and not $all_vdj);
	$tar_cmd .= "*/possorted_bam.bam* " if ($create_bam and $create_bam ne 'false');
	$tar_cmd .= "*/fragments.tsv.gz* */peak_annotation.tsv */singlecell.csv " if (grep {/ATAC/i} @groups);
	# scRNAseq: GEX (+ADT): cloupe.cloupe, web_summary.html, filtered and raw_feature_bc_matrix
	#			VDJ: vloupe.vloupe, web_summary.html
	# scATAseq : cloupe.cloupe, web_summary.html, filtered and raw_peak_bc_matrix, possibly filtered_tf_bc_matrix, fragments.tsv.gz, peak_annotation.tsv
}

warn $tar_cmd;
if (-e "$dir$projectName.tar.gz" and !$no_exec) {
	system("ls -lh $dir$projectName.tar.gz");	
	die("archive '$projectName.tar.gz' already exists") unless (prompt("archive '$dir$projectName.tar.gz' already exists. Overwrite ?  (y/n) ", -yes_no));
}
unless ($no_exec){
	my $exit = system ($tar_cmd);
	die("Error while making the archive") if ($exit);
	print "\t------------------------------------------\n";
#	print "\tlink to send to the users : \n";
#	print "\twww.polyweb.fr/NGS/$projectName/$projectName.tar.gz \n";
	print "\tArchive to send to the users : \n";
	print "\t$dir$projectName.tar.gz\n" if (-e "$dir$projectName.tar.gz");
	print "\t------------------------------------------\n\n";
}

