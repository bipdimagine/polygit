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
die() unless $fork;
die("genome") unless $genome_version; 
my $dir_public= "/data-isilon/public-data/repository/$genome_version/cadd/$version/";
my $file = $dir_public."tabix/whole_genome_SNVs.tsv.gz";

die($file) unless -e $file;
my $file_indel = $dir_public."tabix/InDels.tsv.gz";
my $hash_proc ={};
my $dir_out;

$dir_out= "/data-isilon/public-data/repository/$genome_version/cadd/$version/"."rocksdb_split/";
my $pack = "C";


my $rg;
$rg = GenBoNoSqlRocksGenome->new(pipeline=>1,dir=>$dir_out,mode=>"c",index=>"genomic",chromosome=>$chr_name,genome=>$genome_version,pack=>$pack,description=>["cadd_score"]);



my $pm = new Parallel::ForkManager($fork);
#my $lmdb_out  =  $save_out->get_lmdb_db("snps");
my $nbb =0;
my $nb_regions = scalar(@{$rg->regions});
 	my $eta = Time::ETA->new(
    milestones =>$nb_regions,
);
my $paths;
	$pm->run_on_finish (
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hres) = @_;
			$nbb ++;
			delete $hash_proc->{$hres->{process}};
			my $cc= $hres->{chr};
			my $name = $hres->{data};
			my $rid = $hres->{region_id};
			 $eta->pass_milestone();
 			my $remaining = $eta->get_remaining_time();
			my $elapsed = $eta->get_elapsed_time();
			my $p = int($eta->get_completed_percent());
			print $chr_name." complete $p % => elapsed :".$elapsed." remaining : ".$remaining."\n";
		});


	my $i;
	my $idp = time;
	my $cname = $chr_name;
	foreach my $r (@{$rg->regions}) {
	$idp ++;
	$hash_proc->{$idp} ++;
	my $pid = $pm->start and next;
	my $name = $chr_name.".".$r->{start}.".".$r->{end};
	my $no = $rg->nosql_rocks($r);
	run_chr_key($chr_name,$r,$file,$no);
	run_chr_key($chr_name,$r,$file_indel,$no,1);
	warn "compact";
	$no->rocks->compact_range();
	my $h = $no->rocks->get("0005264624!A");
	warn Dumper  $h;
	$no->close();
	
	#$rg->nosql_rocks($r)->close;
	my $hres;
	$hres->{chr} = $cname;
	$hres->{data} = $name;
	$hres->{process} = $idp;
	$hres->{region_id} = $r->{id};
	$pm->finish(0,$hres);
	}
	
  $pm->wait_all_children;
die(Dumper $hash_proc) if scalar(keys %$hash_proc);
  exit(0);


sub run_chr_key {
	my ($chr,$r,$file_cadd,$no,$add) = @_;

	die($file_cadd) unless -e $file_cadd;

	
	my $ACGT = {'A'=>0,'C'=>1,'G'=>2,'T'=>3};
	my $nb= 0;
	my $obj;
	my $variation;
	my $vdebug;
	my $tabix  = Bio::DB::HTS::Tabix->new( filename => $file_cadd );
	#warn $r->{start}." ".$r->{end};
	my $res = $tabix->query($r->{tabix});
	my $pos;
	my $xx = 1;
	my $vs ={};
	my $old_id;
	my $id;
	my $test =0;
	 while ( my $line = $res->next ) {
				my(@array) = split(" ",$line);
				my $start = $array[1];
				unless ($pos) {
					$pos = $start;
				}
				
				my $alternate_allele = $array[3];
				my $ref_allele = $array[2];
#				warn ($ref_allele." ".$alternate_allele) if length($ref_allele) >1 && length($alternate_allele)>1;
				next if length($ref_allele) >1 && length($alternate_allele)>1;
				unless (exists $variation->{$start}->{s}) {
					$variation->{$start}->{s} = [0,0,0,0];
				}
				my $v = int($array[5]+0.5);
				$id = sprintf("%010d", $start);
				# if $self->index eq "genomic";;
				my $zid = compress_vcf_position($ref_allele,$alternate_allele);
				my $cc= pack("$pack",$v);
				my $zid = $no->return_rocks_id($start,$ref_allele,$alternate_allele);
				#warn $zid;
				warn $id."!".$zid ." ".$v." $pack++=".$cc if $id eq "0005264624";
				$no->put_batch_raw($zid,$cc);
				$test ++;
		#		last if $id."!".$zid eq "0005264624!G";
				#last if $test > 10000;

	  	 }
	 warn "write";
	 sleep(5);
	 $no->write_batch();
	 warn $pack;
	 warn unpack("$pack",$no->get_raw("0005264624!G"));
	#$no->rocks->write($no->batch);
	
	#delete $no->{batch};
	warn "end";
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
	if ($l1 >1 ){
		$ref = substr($ref, 1);
		return ($l1 -1);
	}
	
}
