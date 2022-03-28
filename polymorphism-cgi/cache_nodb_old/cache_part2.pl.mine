#!/usr/bin/perl
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
use Getopt::Long;
use GBuffer;
use GenBoProject;

my $fork = 1;
my ($project_name, $chr_name);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
);

unless ($project_name) { die("\n\nERROR: -project option missing... Die...\n\n"); }

my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name );
my $dir_log = $project->getCacheBitVectorDir()."/log/";

my $cmd1 = "$RealBin/scripts/cache_global_infos.pl -project=$project_name";
eval { system("$cmd1") and die "\n\nERROR: $project_name -> cache_global_infos. Die.\n\n"};
if ($@) {
	die;
};
if (not -e $project->getCacheBitVectorDir()."/global_infos.freeze") {
	warn "\n\nNOT FOUND: ".$project->getCacheBitVectorDir()."/global_infos.freeze\n";
	warn "\nERROR: $project_name -> cache_global_infos. Die.\n\n";
	die;
}

my $cmd2 = "$RealBin/scripts/cache_coverage_polyquery.pl -project=$project_name -fork=$fork";
eval { system("$cmd2") and die "\n\nERROR: $project_name -> coverage. Die.\n\n"};
if ($@) {
	die;
};
if (not -e $project->getCacheBitVectorDir()."/../coverage_lite/$project_name.lite") {
	warn "\n\nNOT FOUND: ".$project->getCacheBitVectorDir()."/../coverage_lite/$project_name.lite\n";
	warn "\nERROR: $project_name -> cache_global_infos. Die.\n\n";
	die;
}

my $cmd3 = "$RealBin/scripts/cache_lite_dejavu.pl -project=$project_name -fork=$fork";
eval { system("$cmd3") and die "\n\nERROR: $project_name -> lite_dejavu. Die.\n\n"};
if ($@) {
	die;
};
