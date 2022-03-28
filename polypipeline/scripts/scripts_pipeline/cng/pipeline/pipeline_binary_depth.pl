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
 @list = `cat ../../../../../defidiag/project/$set/$dir2/$name.txt`;
 chomp(@list);
}
elsif($project_name){
	push(@list,$project_name);
}
die() unless @list;
my $bin_dev = "$Bin/../../";
foreach my $project_name (@list){
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $cmd ="$bin_dev/coverage_genome.pl  -fork=$ppn  -project=$project_name -patient=";
foreach my $patient (@{$project->getPatients}){
	my $name = $patient->name();
	print $cmd.$patient->name." &&  perl $bin_dev/coverage_statistics_genome.pl -patient=$name  -fork=$ppn  -project=$project_name :::$ppn\n";
	
	#my $ppn2 = int($ppn/2);
	#my $cmd2 = qq{perl $bin_dev/wisecondor.pl -patient=$name -project=$project_name -fork=$ppn2 && perl $bin_dev/calling_wisecondor.pl -patient=$name -project=$project_name :::$ppn2};
	#print $cmd2 ."\n";
}
}
