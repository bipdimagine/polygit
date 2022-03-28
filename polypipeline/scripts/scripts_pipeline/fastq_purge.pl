#!/usr/bin/perl

use strict;
use PBS::Client;
use Getopt::Long;
use Data::Dumper;
my $projectName;
my $filename;
my $dir;
my $ffile;
my $lim;

GetOptions(
	"file=s" => \$ffile,
	"lim=n"  => \$lim);

open (fastq, $ffile) or die $!;
my $n = 0;
my $name;
my $seq;
my $sep;
my $qual;
#warn $lim;

while (<fastq>) {
	my $line = $_;
	chomp $line;
#	warn $n." ".$line;
	if ( $n == 0 ) {
		$name = $line;
		if ( !( $name =~ m/^\@/ ) ) {
			die "probleme start seq\n";
		}
		$n++;
	}
	else {
		if ( $n == 1 ) {
			$seq = $line;
			$n++;
		}
		else {
			if ( $n == 2 ) {
				$sep = $line;

				if ( !( $sep =~ m/^\+/ ) ) {
					die "probleme separateur seq\n";
				}
				$n++;
			}
			else {
				if ( $n == 3 ) {
					$qual = $line;
					#warn length($seq);
					if ( ( length($seq) ) > $lim ) {
						if ( ( length( $seq) ) == ( length( $qual) ) ) {
							#warn ".";
							print $name. "\n" 
							  . $seq . "\n" 
							  . $sep . "\n"
							  . $qual . "\n";
						}
					}
					$n = 0;
				}
			}
		}
	}
}
