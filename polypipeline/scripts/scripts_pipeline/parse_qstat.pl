#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use File::Basename;
use  File::Temp;
use YAML::Syck;
use colored; 
use Text::Table;
use POSIX qw{strftime};
use IO::Prompt;

 #my $yaml = YAML::Tiny->read();
my $pid = $ARGV[0];


my @res = `qstat -f `;
my $current_id;
foreach my $line (@res){
	if ($line =~/Job Id:/){
		my ($text,$jobid)  = split (":",$line); 
		$current_id = $jobid;
	}
	
	if ($line =~/$pid/ && $current_id){
		warn "stop $current_id";
		system("qdel $current_id");
		$current_id = undef;
	}
	
	
}

