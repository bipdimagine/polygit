#!/usr/bin/perl

use strict;
use Getopt::Long;
my %col;
my $col_format;
my $limit ;
GetOptions(
	"limit=s"          => \$limit,

);
die() unless $limit;
$| =1;
foreach my $line (<STDIN>){
	chomp($line);
	if ( $line =~ /^#CHROM/ ) {
		my @array = split( " ", $line );

		for ( my $i = 0 ; $i < @array ; $i++ ) {
			if ($array[$i] eq "FORMAT"){
				$col_format = $i;
			}
		}
	}
		
	if ($line =~ /^\#/){
		print $line."\n";
		next;
	}
die() unless $col_format;
my @parser = split(" ",$line);
my @split_format = split(":",$parser[$col_format]);
my $A0;
my $GT;
my $R0;
for ( my $i = 0 ; $i < @split_format ; $i++ ) {
	$A0 = $i if $split_format[$i] eq "AO";
	$GT = $i if $split_format[$i] eq "GT";
	$R0 = $i if $split_format[$i] eq "RO";
}
unless ($A0){
	print $line."\n";
	next;
}
unless ($R0){
	print $line."\n";
	next;
}
for ( my $i = $col_format+1 ; $i < @parser ; $i++ ) {
	my @split_patient = split(":",$parser[$i]);
	next if $split_patient[$GT] ne "0/0";
	my $a = $split_patient[$A0];
	my $r = $split_patient[$R0];
	die() if $A0 =~/,/;
	next if $r == 0;
	if ($a >0){
		warn $parser[$i] if $r ==0;
		warn $parser[$col_format] if $r ==0;
		warn $line if $r ==0;
	my $pourcent = $a/$r;
	if ($pourcent > $limit ){
		
		$split_patient[$GT] = "0/1";
	}
	$parser[$i] = join(":",@split_patient);
	}
}

print join("\t",@parser)."\n";
	
}