#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use Data::Dumper;




while (<STDIN>){
	if ($_ =~/^#/){
		print $_;
		next;
	}
	my $line = $_;
	my @all = split(" ",$_);
	my @infos = split(";",$all[7]);
		my $dp = 0;
		my $dp4 = 0;
		#warn $all[5];
		if ($all[5]>50){
			print $line;
			next;
		}
	
	foreach my $info (@infos){
		next unless $info=~/DP/;
	
		my ($name,$value) = split("=",$info);
		if ($name eq "DP"){
			$dp = $value;
		}
		if ($name eq "DP4"){
			map{$dp4+=$_} split(",",$value);
			
		}
	}
	next if $dp4 ==0;
	next if $dp ==0;
	my $p = int(($dp4/$dp) *100);
	next if ($p < 10) ;
	print $line;
	
	
}