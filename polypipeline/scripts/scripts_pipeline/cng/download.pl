#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use Getopt::Long;



my $name;
#my $password;
my $dir_in;
#my $type;
my $set;

GetOptions(
	'name=s' => \$name,
	'set=s' => \$set,
);

my @list = `cat list.txt`;
chomp(@list);
my ($l) = grep{$_ =~ /$name/} @list;
die() unless $l;

my($n,$password,$type) = split(" ",$l);

die("type = devodecode or defidiag") unless $type;


my $dir = "/data-isilon/download/$type/";
die($type." ".$dir) unless -e $dir;

my $cmd_bam="wget --recursive --no-parent --user $name  --password  $password https://www.cnrgh.fr/data/$name/alignements/$set/ -R cram";  
my $cmd_var="wget --recursive --no-parent --user $name  --password  $password https://www.cnrgh.fr/data/$name/variants/$set/ -R cram";  

warn "cd $dir && $cmd_bam";
system("cd $dir && $cmd_bam");
system("cd $dir && $cmd_var");