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


my @projectToDo;
my @projectAll;
my $buffer1 = GBuffer->new();	
my @projectAll = @{$buffer1->listProjects()}; 

#my @chromosomes = (1..22,'X','Y','MT');
my $nb_project = scalar(@projectAll);
my $nbp =0;
foreach my $project_name (@projectAll){
	print $project_name."\n";	

}
exit(0);