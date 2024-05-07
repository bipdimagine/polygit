#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin/../../../../../GenBo/lib/obj-nodb/";
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
my $dir_public= "/data-isilon/public-data/repository/$genome_version/dbscSNV/$version/";
my $file = $dir_public."tabix/dbscSNV1.1.$genome_version.txt.gz";

die($file) unless -e $file;

my $dir_out;

$dir_out= "/data-isilon/public-data/repository/$genome_version/dbscSNV/$version/"."rocksdb/";
my $pack = "s2";
my $factor = ["0.001","0.001"];
my $chromosomes = [1..22,'X','Y','MT'];
my $pm = new Parallel::ForkManager($fork);
foreach my $chr_name (@$chromosomes){
	my $pid = $pm->start and next;
	my $rg;
	$rg = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out,factor=>$factor,mode=>"c",name=>$chr_name,genome=>$genome_version,pack=>$pack,description=>["ada","rf"],pipeline=>1);



	run_chr_key($chr_name,$file,$rg);
	$rg->write_batch;
	$rg->close();
	$pm->finish(0,{});
	}
	
  $pm->wait_all_children;

  exit(0);


sub run_chr_key {
	my ($chr,$file_dbscSNV,$rg) = @_;

	die($file_dbscSNV) unless -e $file_dbscSNV;

	my $nb= 0;
	my $obj;
	my $vdebug;
	my $tabix  = Bio::DB::HTS::Tabix->new( filename => $file_dbscSNV );
	my $res = $tabix->query($chr.":1-10000000000");
	 $rg->put_batch("nothing",{rien=>"nada"}) unless $res;
	 return unless $res;
	my $pos;
	my $xx = 1;
	my $vs ={};
	my $old_id;
	my $id;
	my $test =0;
	my %dj;
	my $pos_rf = 16;
	my $pos_ada = 17;
	if ($genome_version =~ /38/){
		$pos_rf=4;
		$pos_ada =5;
	}
	 while ( my $line = $res->next ) {
				my(@array) = split(" ",$line);
				my $start = $array[1];
				warn $test if $test %1000000 ==0 ;
				my $alternate_allele = $array[3];
				my $ref_allele = $array[2];
#				warn ($ref_allele." ".$alternate_allele) if length($ref_allele) >1 && length($alternate_allele)>1;
				#next if length($ref_allele) >1 && length($alternate_allele)>1;
				
				my $ada = int($array[$pos_rf]/$factor->[0]);
				my $rf = int($array[$pos_ada]/$factor->[1]);
				$id = sprintf("%010d", $start);
				# if $self->index eq "genomic";;
				my $zid = $rg->return_rocks_id($start,$ref_allele,$alternate_allele);
				my $cc= pack("$pack",($ada,$rf));
				confess($zid) if exists $dj{$zid};
				$dj{$zid." ".$start." ".$ref_allele." ".$alternate_allele} ++;
#				warn $zid." ".$cc;
				$rg->put_batch_raw($zid,$cc);
				$test ++;
				#last if $test > 10000;

	  	 }
return ;
}




