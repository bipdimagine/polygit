#!/usr/bin/perl

use FindBin qw($Bin);
use FindBin qw($RealBin);
use lib "$RealBin/../GenBo";
use lib "$RealBin/../GenBo/lib/GenBoDB";
use lib "$RealBin/../GenBo/lib/obj-nodb";
use Carp;
use strict;
use GBuffer; 
use Getopt::Long;

my $cmd;
GetOptions(
	'cmd=s' => \$cmd,
);
unless ($cmd) { die("\n\nERROR !! No command !! Exit...\n\n"); }


my @projectToDo;
my @projectAll;
my $buffer1 = GBuffer->new();	
my @projectAll = @{$buffer1->listProjects()}; 

#my @chromosomes = (1..22,'X','Y','MT');
my $nb_project = scalar(@projectAll);
my $nbp =0;
foreach my $project_name (@projectAll){
	$nbp ++;
	next if $project_name !~/NGS/;
	#next if $project_name eq "NGS2015_0775" or $project_name eq "NGS2015_0629" or $project_name eq "NGS2015_0608";
	print "running $project_name : $nbp/".$nb_project."\n";
	my $cmd = qq{$cmd -project=}.$project_name;
	system ("$cmd");
	warn "end $project_name \n";

}

exit(0);