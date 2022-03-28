#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../../GenBo/lib/obj-nodb/";
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

GetOptions(
	'project=s' => \$project_name,
	'fork=s' => \$ppn,
	'set=s' => \$set,
	'name=s' => \$name,
	#'fork=s' => \$ppn,
);
$set = "set".$set unless $set=~/set/;
my @list;
if($name){
 @list = `cat ../../../../../defidiag/project/$set/$name.txt`;
 chomp(@list);
}
elsif($project_name){
	push(@list,$project_name);
}
die() unless @list;

foreach my $project_name (@list){
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $cmd ="$Bin/launch_chr_cache.pl -project=$project_name -fork=$ppn -chr=";
foreach my $chr (@{$project->getChromosomes}){
	print $cmd.$chr->name."\n";
}
}
