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
use GenBoNoSqlRocksAnnotation;

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
die() unless $fork;
die("genome") unless $genome_version; 
my $dir_public= "/data-isilon/public-data/repository/$genome_version/revel/$version/";
my $file = $dir_public."tabix/revel.$genome_version.txt.gz";

die($file) unless -e $file;

my $dir_out;

$dir_out= "/data-isilon/public-data/repository/$genome_version/revel/$version/"."rocksdb/";
my $pack = "C A";
my $factor = ["0.001",undef];
my $chromosomes = [1..22,'X','Y','MT'];
my $pm = new Parallel::ForkManager($fork);
foreach my $chr_name (@$chromosomes){
	my $pid = $pm->start and next;
	my $rg;

	$rg = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out,mode=>"c",name=>$chr_name,pack=>$pack,description=>["score"],pipeline=>1);


	run_chr_key($chr_name,$file,$rg);
	$rg->write_batch;
	$rg->close();
	$pm->finish(0,{});
	}
	
  $pm->wait_all_children;

  exit(0);


sub run_chr_key {
	my ($chr,$file_dbscSNV,$no) = @_;

	die($file_dbscSNV) unless -e $file_dbscSNV;

	#my $no = $rg->nosql_rocks($r);
	my $nb= 0;
	my $obj;
	my $vdebug;
	my $tabix  = Bio::DB::HTS::Tabix->new( filename => $file_dbscSNV );
	my $res = $tabix->query($chr.":1-10000000000");
	my $pos;
	my $xx = 1;
	my $vs ={};
	my $old_id;
	my $id;
	my $test =0;
	my $hash;
	 $no->put_batch("nothing",{rien=>"nada"}) unless $res;
	 return unless $res;
	 while ( my $line = $res->next ) {
				my(@array) = split(" ",$line);
				warn $chr." ".$test if $test%200_000 == 0;
				my $start = $array[1];
				my $revel =  $array[6];
				my @v;
				my $i=0;
				my $alternate_allele = $array[3];
				my $ref_allele = $array[2];
				my $aa_alternate = $array[5];
				$id = sprintf("%010d", $start);
				$old_id = $id unless $old_id;
				my $zid = compress_vcf_position($ref_allele,$alternate_allele);
				#$revel = int ($revel * 1000);
				push(@{$hash->{$id."!".$zid}}, [$revel,$aa_alternate]);
				if ($id ne $old_id){
					 foreach my $id (keys %$hash){
	 					$no->put_batch($id,$hash->{$id});
	  				}
	  				$hash = {};
	  				$old_id = $id;
				}
				$test ++; 	 

	  	 }
	  	 foreach my $id (keys %$hash){
	 			$no->put_batch($id,$hash->{$id});
	  	}
	  	
return ;
}



sub compress_vcf_position {
	my ($ref,$alt) = @_;
	my $l1 = length($ref);
	my $l2 = length($alt);
	if ($l1 ==1 && $l2 ==1){
		return $alt;
	}
	if ($l1 ==1 && $l2 > 1){
		return "+".substr($alt, 1);
	}
	if ($l1 >1 && $l2 == 1){
		$ref = substr($ref, 1);
		return ($l1 -1);
	}
	confess();
	
}
