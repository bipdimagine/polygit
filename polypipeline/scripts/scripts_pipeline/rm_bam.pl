#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";

use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use colored;
use File::Basename;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;


my $buffer = GBuffer->new();
 my $dir_trash =  $buffer->config->{project_data}->{root}."/TRASH/";
 system("mkdir -p $dir_trash && chmod a+rwx $dir_trash") unless -e $dir_trash;
 

my $bam_prod = $ARGV[0];

die("can't find $bam_prod") unless -e $bam_prod;


my $bai_prod = $bam_prod.".bai";
unless (-e  $bai_prod){
	$bai_prod = $bam_prod;
	$bai_prod =~s/bam/bai/;
}

if (-e $bai_prod){
	my $cmd2 = "chmod a+w $bai_prod && rm $bai_prod";
	system ($cmd2);
}



my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $ymd = sprintf("%04d_%02d_%02d_%02d_%02d",$year+1900,$mon+1,$mday,$hour,$min);
 my($trash_file, $dirs, $suffix) = fileparse($bam_prod);

$trash_file =~s/\.bam/\.$ymd\.bam/;

$trash_file = $dir_trash."/".$trash_file;


my $cmd = qq{chmod a+w $bam_prod;  mv  $bam_prod $trash_file};
system ($cmd);

die("problem with delete $bam_prod ") if -e $bam_prod;

exit(0);






