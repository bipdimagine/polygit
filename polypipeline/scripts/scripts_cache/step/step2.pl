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
 @list = `cat $Bin/../../../../defidiag/project/$name/$set.txt`;
 chomp(@list);
}
elsif($project_name){
	push(@list,split",",$project_name);
}
die() unless @list;
my $p = $project_name;
my $bin_dev = "$Bin/../../scripts_pipeline";
my $cmd_perl = "$Bin/../../../../polymorphism-cgi/cache_nodb/scripts";
foreach my $p (@list) {
	#cnv,dejavu,quality_check,identito_vigilence
	my $cmd = "perl $bin_dev/manue_cnv/SV_global.pl -project=$p -fork=5 -nodejavu=1 :::5";
	my $cmd2 = "perl $cmd_perl/cache_lite_dejavu.pl -project=$p :::2";
	my $cmd3 = "perl $bin_dev/quality_check.pl -project=$p -fork=5 -cache=1 :::5";
	my $cmd4 =  qq{perl $bin_dev/identito_vigilence.pl  -fork=1  -project=$p :::2 };
	
	print "$cmd\n";
	print $cmd2."\n";
	print $cmd3."\n";
	print $cmd4."\n";
}