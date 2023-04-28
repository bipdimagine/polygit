#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo/lib/obj-nodb/";
use lib "$Bin/packages";
use GBuffer;
#use GenBoProject;
#use PBS::Client;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
use MooseX::Method::Signatures;
use bds_cache_steps;    
#use file_util;
use Class::Inspector;
#use check_utils;
use Text::Table;
use colored; 
use Term::Menus;
use Config::Std;
 
my ($projectName, $filename, $name, $patients_name, $steps_name, $force, $type, $fastq_ext, $somatic, $method, $no_cluster, $stdout, $help, $annot_version);
my $fork = 1;
my $yes = 0;
my $nocluster = 0;
my $menu= 0;
my $filename_cfg;
my $secret;
my $analyse_type;
my $nobackup;
my $giab;
GetOptions(
	'project=s' => \$projectName,
	'force=s' => \$force,
);

print "\n\n##### CACHE Junctions #####\n\n";
my $cmd = "$Bin/bds_cache.pl -type=rna_junctions";
$cmd .= " -project=".$projectName;
if ($force) {
	$cmd .= " -force=1";
}
system ($cmd);



