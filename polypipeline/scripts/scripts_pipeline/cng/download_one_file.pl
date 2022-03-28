#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use Getopt::Long;
use Data::Dumper;


my $name;
#my $password;
my $dir_in;
#my $type;
my $set;
my $bc;
GetOptions(
	'name=s' => \$name,
	'set=s' => \$set,
	'bc=s' => \$bc,
);

my @list = `cat /software/polyweb/poly-disk/poly-src/defidiag/project/download/list.txt`;
chomp(@list);
warn Dumper (@list);
my ($l) = grep{$_ =~ /$name/} @list;
die() unless $l;

my($n,$password,$type) = split(" ",$l);

die("type = devodecode or defidiag") unless $type;


my $dir = "/data-isilon/download/$type/";
system("mkdir $dir") unless -e $dir;
die() unless -e $dir;

my $cmd_bam="wget --no-clobber --recursive --no-parent --user $name  --password  $password https://www.cnrgh.fr/data/$name/alignements/$set/ -A \"*$bc*\" -R cram";  
my $cmd_var="wget --no-clobber --recursive --no-parent --user $name  --password  $password https://www.cnrgh.fr/data/$name/variants/$set/ -A \"*$bc*\" -R cram";  

warn "cd $dir && $cmd_bam";
system("cd $dir && $cmd_bam");
system("cd $dir && $cmd_var");