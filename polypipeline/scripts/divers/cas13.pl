#!/usr/bin/perl
use Data::Dumper;
use File::Find;
use Getopt::Long;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/";
use lib "$Bin/../../../GenBo/lib/GenBoDB";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../packages";
use GBuffer;

my $out;
my $file;

GetOptions(
	'file=s' => \$file,
	'out=s' => \$out,
 );

my @bases; 
open(SEQ, $file) or die ("$file fichier inconnu ");
while(<SEQ>){
		chomp($_);
		@bases = split("",$_);
}

open(OUT, ">$out") or die ("$out fichier inconnu ");

warn scalar(@bases);
for(my $i=0; $i < scalar(@bases)-24;$i+=3){
	my $start = $i;
	my $end = $i+24; 
	warn $start;
	warn $end;
	my $XXX = $bases[$start].$bases[$start+1].$bases[$start+2];
	my $YYY = $bases[$end-2].$bases[$end-1].$bases[$end];
	print OUT ">BCL11-".$XXX."-".$YYY."\n";
	my @seq;
	for (my $y=$start; $y <= $end; $y++){
		push(@seq,$bases[$y]);
		
	}
	my $chain= join("",@seq);
	print OUT $chain."\n";
}