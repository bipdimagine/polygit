#!/usr/bin/perl
use strict;
use Getopt::Long;
use Data::Dumper;
my $file;
GetOptions(
	'file=s'  => \$file,
);
die($file) unless -e $file;
#unless (-e )
my $fai = $file.'.fai';
die() unless -e $fai;
my $file_ori =  $fai.".ori";

system("cp $fai $file_ori");
die("mv $file.fai $file_ori") unless -e $file_ori;
 
open(FAI,"$file_ori");
open(FAI_OUT,">$fai");
while(<FAI>){
	
	chomp();
	my @z = split(" ",$_);
	my (@id) = split(/\|/,$z[0]);
	my ($enst) = grep {$_=~/ENST/} @id;
	$enst =~ s/\..*//;
	my $y ="";
	$y="_Y"	if $z[0] =~ /PAR_Y/;
    $z[0] = $enst.$y;
    
    print FAI_OUT join("\t",@z)."\n";
  
}
close (FAI_OUT);