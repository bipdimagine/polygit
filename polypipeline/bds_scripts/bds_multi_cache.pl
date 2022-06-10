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
my $force;
GetOptions(
	'project=s' => \$projectName,
	'force=s' => \$force,
);


my @projects = split(",",$projectName);



my $tname = "/data-beegfs/tmp/cache.";
my $nb =0;
my @log;
foreach my $pr (@projects){
	$nb ++;
	my $tname2 = $tname."$pr.log";
	my $cforce ="";
	 $cforce = "-force=1" if $force;
	# warn "running $pr $nb/".scalar(@$pr);
	push(@log,$tname2);
	system ("bds_cache.sh -project=$pr -yes=1 $cforce >$tname2");

	
}

print "LOG FILES \n";
foreach my $l (@log) {
	print "$l\n";
}

 exit(0);