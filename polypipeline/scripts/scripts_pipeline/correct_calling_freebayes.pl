#!/usr/bin/perl

use strict;
use Getopt::Long;
use Data::Dumper;
my %col;
my $col_format;
my $bam ;
my $samtools;
GetOptions(
	"bam=s"          => \$bam,
	"samtools=s" =>\$samtools,
	
);
$| =1;
foreach my $line (<STDIN>){
	chomp($line);
	if ( $line =~ /^#/ ) {
		print $line."\n";
		next;
	}
	my @arg = split(" ",$line);	
	my $score= $arg[5];
	if ($arg[5] < 3){
		next;
	}
#	if ($arg[5] < 100){
		my ($DP) = grep {$_=~/DP/} split(";",$arg[7]);
		next unless $DP;
		my ($a,$cov) = split("=",$DP);
		my $cmd = qq{$samtools depth $bam -r $arg[0]:$arg[1]-$arg[1]};
		my ($s) = `$cmd | cut -f 3`;
		chomp($s);
	#	warn " score:$arg[5] $s=>$cov $arg[0]:$arg[1]-$arg[1]";
		
		next if $cov*30 < $s;
#	}
	print $line."\n";
}
	