#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/obj-nodb/";
use lib "$Bin/../packages";
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
 
my ($projectName);

my $nobackup;
GetOptions(
	'project=s' => \$projectName,
);


my @projects = split(",",$projectName);

my $tname = "/tmp/dude.".time.".".rand(time);
open (TOTO,">$tname");

foreach my $pr (@projects){
	print TOTO "$Bin/scripts/scripts_pipeline/dude/dude.pl -project=$pr -fork=20\n";
	
}

system("cat $tname | run_cluster.pl -cpu=20 ");

unlink $tname;

 exit(0);