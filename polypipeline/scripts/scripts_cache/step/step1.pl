#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use GBuffer;
use GenBoProject;
use Getopt::Long;
use Data::Dumper;
#use illumina_util;
use IO::Prompt;
use Sys::Hostname;
 use File::Find::Rule ;
use Text::Table;
use Term::Twiddle;

my $project_name;
my $chr_name;
my $ppn;
my $set;
my $name;
my $solo;
GetOptions(
	'project=s' => \$project_name,
	'fork=s' => \$ppn,
	'set=s' => \$set,
	'name=s' => \$name,
	'solo=s' => \$solo,
);
my @list;
my $dir2 = "";
if($solo){
	$dir2 ="solo";
}
$set = "set".$set unless $set=~/set/;
if($name){
 @list = `cat ../../../../defidiag/project/$name/$set.txt`;
 chomp(@list);
}
elsif($project_name){
	push(@list,split(",",$project_name));
}
die() unless @list;

foreach my $project_name (@list){
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $cmd ="$Bin/launch_chr_cache.pl -project=$project_name -fork=20 -chr=";
foreach my $chr (@{$project->getChromosomes}){
	print $cmd.$chr->name." :::20\n";
}
}
