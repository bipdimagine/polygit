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
my $p = $project_name;
foreach my $p (@list) {
	#cnv,dejavu,quality_check,identito_vigilence
	my $cmd = "perl $Bin/../../manue_cnv/SV_global.pl -project=$p -fork=$ppn -nodejavu=1";
	my $cmd2 = "perl $Bin/../../../../../polymorphism-cgi/cache_nodb/scripts/cache_lite_dejavu.pl -project=$p ";
	my $cmd3 = "perl $Bin/../../quality_check.pl -project=$p -fork=$ppn -cache=1 ";
	my $cmd4 =  qq{perl $Bin/../../identito_vigilence.pl  -fork=$ppn  -project=$p };
	
	print "$cmd\n";
	print $cmd2."\n";
	print $cmd3."\n";
	print $cmd4."\n";
}