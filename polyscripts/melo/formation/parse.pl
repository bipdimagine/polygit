#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
use Data::Dumper;

use GBuffer;




open(my $file, '<', "/home/mperin/test.ped") # '<' read-only
	or die "Can't open the file\n";
	
# '>' write, '>>' write in append mode
open(my $outputFile, ">", "/home/mperin/parseOutput.txt");

# readline, read, getc, sysread
while (my $line = readline($file)) {
	chomp $line;
	(my $family, my $name, my $father, my $mother, my $sex, my $status) = split("\t", $line);
	
	print $outputFile $name.", ";
	
	if ($sex==1){print $outputFile "garçon, "}
	else {print $outputFile "fille, "}
	
	if ($status==1){print $outputFile "sain"}
	else {print $outputFile "malade"}
	print $outputFile "\n"
}
close ($file);

close($outputFile);

print "Done\n"
