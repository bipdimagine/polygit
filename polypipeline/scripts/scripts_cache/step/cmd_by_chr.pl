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

my $project_name2;
my $chr_name;
my $ppn;
my $set;
my $name;

GetOptions(
	'project=s' => \$project_name2,
	'fork=s' => \$ppn,
	'set=s' => \$set,
	'name=s' => \$name,
	#'fork=s' => \$ppn,
);
$set = "set".$set unless $set=~/set/;
my @list = split(",",$project_name2);



foreach my $project_name (@list){
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
warn $project->name;
warn $project;
my $yes;
warn $Bin;
warn "$Bin/../../../../../GenBo/lib/obj-nodb/";
my $cache_directory_actual = $project->getCacheDir();
		warn "As you wish .... ".$cache_directory_actual;
		system ("rm -rf $cache_directory_actual/*") if -e $cache_directory_actual;
		$cache_directory_actual = $project->getCacheDir();
		warn $project->rocks_cache_2_root_dir();
		my $tr = $project->rocks_cache_2_root_dir()."/".$project->name;
		my $tr2 = $project->tiny_rocks_cache_dir();
		system ("rm -rf $tr/*") if -e $tr2;
		system("rmdir $cache_directory_actual");
		system("rmdir $tr");
		

my $cmd ="$Bin/launch_chr_cache.pl -project=$project_name -fork=$ppn -chr=";
foreach my $chr (@{$project->getChromosomes}){
	print $cmd.$chr->name."\n";
}
}
