#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

my $file = $ARGV[0];

foreach my $file (@ARGV){
die("$file not found ") unless -e $file;


if ($file =~/\.bam/) {
	system("$Bin/rm_bam.pl $file");
}
elsif ($file =~/\.log/) {
	warn $file;
	unlink $file;
	system("$Bin/rm_log.pl $file");
}
elsif ($file =~/cache/) {
	unlink $file;
	system("$Bin/rm_cache.pl $file");
}
else {
	system("$Bin/rm_vcf.pl $file");
}
}