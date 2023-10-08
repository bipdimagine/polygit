#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Getopt::Long;


my $fork = 1;
my ($project_name, $chr_name);
GetOptions(
	'project=s'    => \$project_name,
	'chr=s'        => \$chr_name,
	'fork=s'        => \$fork,
);

my $opt= " -project=$project_name -chr=$chr_name -fork=$fork";
my @run = ("cache_store_ids.pl","cache_store_annotations.pl","cache_strict_denovo.pl","update_vector_score.pl","update_hash_variant_chromosome.pl");
#@run = ("cache_store_annotations.pl","cache_strict_denovo.pl","update_hash_variant_chromosome.pl");


foreach my $r (@run){
	my $a = system("$RealBin/$r $opt");
	die() if $a != 0;
}
