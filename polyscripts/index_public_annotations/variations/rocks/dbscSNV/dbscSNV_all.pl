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
use Parallel::ForkManager;
 use Time::ETA;
use List::MoreUtils qw(any uniq natatime);
use Bio::DB::HTS::Tabix;
use POSIX;
my $allsnps;
use Getopt::Long;
use Date::Tiny;
use Bio::DB::HTS::Faidx;
use GenBoNoSqlRocksGenome;
use GenBoNoSqlRocksAnnotation;
my $chr_name;
my $version;
my $fork;
my $genome_version;
my $merge;
GetOptions(
	'chr=s' => \$chr_name,
	'version=s' => \$version,
	'genome=s' => \$genome_version,
);
die("genome") unless $genome_version; 
my $dir_public= "/data-isilon/public-data/repository/$genome_version/dbscSNV/$version/";
my $dir_in =  "$dir_public"."rocksdb_split/";
my $dir_out = "$dir_public"."rocksdb/";
if ($chr_name eq "MT"){
	my $finalrg = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out,mode=>"t",name=>$chr_name,pack=>[],version=>"",description=>{},pipeline=>1);	
	$finalrg->put_batch_raw("date",time);
	$finalrg->write_batch();
	$finalrg->close();
	exit(0);
}

my $rg = GenBoNoSqlRocksGenome->new(dir=>$dir_in,mode=>"r",index=>"genomic",chromosome=>$chr_name,genome=>$genome_version);	
warn Dumper $rg->chunks;
warn $rg->pack();
warn Dumper $rg->description();
warn $dir_out;

my $factor = ["0.001","0.001"];
#my $rg = GenBoNoSqlRocksGenome->new(dir=>$dir_out,mode=>"c",chromosome=>$chr_name,fasta=>$genome_fasta,compress=>1,index=>"vcf");
my $finalrg = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out,mode=>"t",name=>$chr_name,factor=>$factor,pack=>$rg->pack,version=>$version,description=>$rg->description,pipeline=>1);	

my $nb_regions = 0;
	foreach my $r (@{$rg->regions}) {
		my $no = $rg->_chunks_rocks($r->{id});
		$nb_regions ++;
		warn "chunk : ".$nb_regions."/".scalar(@{$rg->regions});
		my $iter = $no->rocks->new_iterator->seek_to_first;
		while (my ($key, $value) = $iter->each) {
    		printf "%s: %s\n", $key, $value;
    		#my $h = $no->get($key);
    		#foreach my $k (keys %$h){
    		#	$finalrg->put_batch_raw($key."!".$k,$h->{$k});
    		#}
    		#warn Dumper $h;
    		#die();
    		$finalrg->put_batch_raw($key,$value);
		}
		$finalrg->write_batch();
		
		$no = undef;
	
	}
$finalrg->close();

  exit(0);



