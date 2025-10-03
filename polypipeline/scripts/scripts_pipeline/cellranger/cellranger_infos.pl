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
	my $dir1 = $project->getProjectRootPath().$pat->name.'/';
	my $dir2 = $dir.$pat->name.'/';
	if (-d $dir1) {
		$dir_cellranger = $dir1;
		$dir_cellranger .= 'outs/' if (-d $dir_cellranger.'outs/');
	}
	elsif (-d $dir2) {
		$dir_cellranger = $dir2;
		$dir_cellranger .= 'outs/' if (-d $dir_cellranger.'outs/');
	}
	else {
		die("Can't find directory with cellranger data for patient $pname")
	}
	
	# Read from web_summary
	my $pipeline_version = "";
	my $transcriptome_version = "";
	my $transcriptome_path = "";
	die("Web summary for patient '$pname' not found: '$dir_cellranger/web_summary.html'") unless (-e $dir_cellranger.'web_summary.html') ;
#	warn $dir_cellranger.'web_summary.html';
	# ["Pipeline Version",;"cellranger-9.0.1"] ["Transcriptome","MM39-Notch3-rat-"] ["V(D)J Reference","vdj_GRCh38_alts_ensembl-5.0.0"]
	my $web_summary = $dir_cellranger.'web_summary.html';
	my $regex_pipeline = qr/\[\"Pipeline [Vv]ersion\",\"((cellranger(-\w*)?|spaceranger)-\d+\.\d+\.\d+)\"\]/;
	my $regex_transciptome = qr/\[\"(Transcriptome|V\(D\)J Reference)\",\"([a-zA-Z0-9_\-.]*)\"\]/;
#	my $regex_vdj_reference = qr/\[\"V\(D\)J Reference\",\"([\w-]*)\"\]/;
	my $regex_transciptome_atac = qr/\[\"Reference path\",\"([\w\/-]*)\"\]/;
	open(my $ws, '<', $web_summary) || warn ("Can't open file '$web_summary' : $!");
	while (my $line = <$ws>) {
		$pipeline_version = $1 if ($line =~ $regex_pipeline);
		$transcriptome_version = $2 if ($line =~ $regex_transciptome);
		$transcriptome_path = $1 if ($line =~ $regex_transciptome_atac);
#		$vdj_reference = $1 if ($line =~ $regex_vdj_reference);
		last if ($pipeline_version and ($transcriptome_version or $transcriptome_path));
	}
	close($ws);
	$transcriptome_version =~ s/-$//;
	
	
	# cellranger version
	my $file_version = $dir_cellranger.'_versions';
	if (-e $file_version) {
		my $cellranger_version;
		open (my $fv, '<', $file_version) || warn ("Can't open file '$file_version': $!");
		while (my $line = readline($fv)) {
			chomp $line;
			$cellranger_version = $1 if ($line =~ /^\s*"pipelines": "((cellranger(-\w*)?|spaceranger)-[\d\.]+)"$/);
		}
		close($fv);
		die("Different pipeline versions from _version ($cellranger_version) and web_summary ($pipeline_version) for patient '$pname'.") unless ($cellranger_version eq $pipeline_version);
	}
	
	# transcriptome
	my $release = $project->annotation_genome_version;
	my $index = $project->getGenomeIndex($pat->alignmentMethod()) =~ s/\/\//\//gr;
	#$index =~ s/\/\//\//;
	my $methSeq = $project->getRun->infosRun->{method};
	if($methSeq eq 'atac') {
		$index .= '_atac';
		warn("Different transcriptome path from current release ($release --> $index) and web_summary ($transcriptome_path) for patient '$pname'.") unless ($index eq $transcriptome_path);
	}
	$index .= '_vdj' if (lc($pat->somatic_group()) eq 'vdj');
	my $json_file = $index.'/reference.json';
#	warn $json_file;
	open(my $json, '<', $json_file) || warn ("Can't open file '$json_file' : $!");
	my $json_string = do { local $/; <$json> };
	my $reference_json = decode_json($json_string);
	close($json);
	my $transcriptome = join('-', @{$reference_json->{'genomes'}}) if (ref $reference_json->{'genomes'} eq 'ARRAY');
	$transcriptome = $reference_json->{'genomes'} unless (ref $reference_json->{'genomes'} eq 'ARRAY');
	$transcriptome .= '-'.$reference_json->{'version'} if (exists $reference_json->{'version'});
	die("Different transcriptomes from current reference.json ($release --> $transcriptome) and web_summary ($transcriptome_version) for patient '$pname'.") unless ($transcriptome eq $transcriptome_version or $methSeq eq 'atac');
	
			# sample,    fastq path,					 cellranger out dir, cellranger version	   transcriptome
	print $fh $pname.','.$pat->getSequencesDirectory.','.$dir_cellranger.','.$pipeline_version.','.$transcriptome_version."\n";
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

