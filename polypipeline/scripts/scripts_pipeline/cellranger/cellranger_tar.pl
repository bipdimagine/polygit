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
my $create_bam;
my $no_exec;

GetOptions(
	'project=s'		=> \$projectName,
	'patients=s'	=> \$patients_name,
	'create_bam!'	=> \$create_bam,
	'no_exec'		=> \$no_exec,
) || die("Error in command line arguments\n");

die("--project argument is mandatory") unless ($projectName);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
my $patients = $project->getPatients($patients_name);
my $run = $project->getRun();
my $type = $run->infosRun->{method};

my $dir = $project->getCountingDir('cellranger');
$dir = $project->getCountingDir('spaceranger') if ($type eq 'spatial');
warn $dir;

my @group;
foreach my $pat (@{$patients}) {
	my $group = $pat->somatic_group();
	push(@group,$group);
}

my $tar_cmd = "tar -cvzf $dir$projectName.tar.gz $dir*/web_summary.html $dir*/cloupe.cloupe ";
$tar_cmd .= "$dir*/vloupe.vloupe " if grep (/vdj/i, @group);
$tar_cmd .= "$dir*/*_bc_matrix/* ";
$tar_cmd .= "$dir*/*.bam* " if ($create_bam and $create_bam ne 'false');
$tar_cmd .= "$dir*/fragments.tsv.gz $dir*/peak_annotation.tsv " if (grep {/ATAC/i} @group);
warn grep {/ATAC/i} @group;
# scRNAseq: GEX (+ADT): cloupe.cloupe, web_summary.html, filtered and raw_feature_bc_matrix
#			VDJ: vloupe.vloupe, web_summary.html
# scATAseq : cloupe.cloupe, web_summary.html, filtered and raw_peak_bc_matrix, possibly filtered_tf_bc_matrix, fragments.tsv.gz, peak_annotation.tsv

warn $tar_cmd;
if (-e "$dir/$projectName.tar.gz" and !$no_exec) {
	my $overwrite = prompt("'archive $dir/$projectName.tar.gz' already exists. Overwrite ?  (y/n) ", -yes);
	die("archive '$dir/$projectName.tar.gz' already exists") unless ($overwrite);
}
unless ($no_exec){
	system ($tar_cmd);
	print "\t------------------------------------------\n";
#	print "\tlink to send to the users : \n";
#	print "\twww.polyweb.fr/NGS/$projectName/$projectName.tar.gz \n";
	print "\tArchive to send to the users : \n";
	print "\t$dir$projectName.tar.gz\n" if (-e "$dir$projectName.tar.gz");
	print "\t------------------------------------------\n\n";
}

