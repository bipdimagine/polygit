#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use strict; 
use GenBoNoSqlRocks;
use Data::Printer;
use Getopt::Long;
use Carp;
use Data::Dumper;
use Getopt::Long; 
use GBuffer;
use GenBoProject;
use Parallel::ForkManager;
use GenBoNoSql;
 use Time::ETA;
use List::MoreUtils qw(any uniq natatime);
use Bio::DB::HTS::Tabix;
use POSIX;
my $allsnps;
use Getopt::Long;
use Date::Tiny;
use Bio::DB::HTS::Faidx;
use GenBoNoSqlRocksGenome;

my $chr_name;
my $version;
my $fork;
my $genome_version;
my $merge;
GetOptions(
	'chr=s' => \$chr_name,
	'version=s' => \$version,
	'fork=s' => \$fork,
	'genome=s' => \$genome_version,
	'merge=s' => \$merge,
);
die("genome") unless $genome_version; 
my $dir_public= "/data-isilon/public-data/repository/$genome_version/alphamissense/$version/";
my $file = $dir_public."tabix/AlphaMissense_isoforms_aa_substitutions.tsv.gz";
my $dir_out= $dir_public."/rocks/";
my $AA  = [
    'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
];
my $hAA;
for (my $i=0;$i<=@$AA;$i++){
	$hAA->{$AA->[$i]} = $i;
}
warn Dumper $hAA;

my $x =0;
my $finalrg = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out,mode=>"r",name=>"alphamissense");
 $finalrg->rocks;
	my $iter = $finalrg->rocks->new_iterator->seek_to_first;
		while (my ($var_id, $value) = $iter->each) {
			my $v = $finalrg->decode($value);
			warn Dumper $v;
			$x++;
		}
		
		
warn $x;		