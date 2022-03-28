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
	'fork=s' => \$ppn,
);
$ppn= 20 unless $ppn;
my @list;
my $dir2 = "";
if($solo){
	$dir2 ="solo";
}
$set = "set".$set unless $set=~/set/;
if($name){
 @list = `cat ../../../../defidiag/project/$set/$dir2/$name.txt`;
 chomp(@list);
}
elsif($project_name){
	push(@list,split(",",$project_name));
}
die() unless @list;

foreach my $project_name (@list){
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $cmd_perl = "/usr/bin/perl $Bin/../../../../polymorphism-cgi/cache_nodb/scripts";
my $hash_cmd;
my $ppn = 10;
foreach my $chr (@{$project->getChromosomes}){
	my $chr_name = $chr->name;
	print "$cmd_perl/polyviewer.pl -project=$project_name -chr=$chr_name -fork=$ppn  && touch ".$chr->lmdb_cache_dir()."/lmdb.ok ; test -e ".$chr->lmdb_cache_dir()." :::$ppn\n";
	
}
}

