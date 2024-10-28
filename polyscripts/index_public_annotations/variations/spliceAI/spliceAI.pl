#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use strict; 
use Data::Printer;
use Getopt::Long;
use Carp;
use JSON::XS;
use Data::Dumper;
use Getopt::Long; 
use GBuffer;
use GenBoProject;
use Parallel::ForkManager;
 use Time::ETA;
use List::MoreUtils qw(any uniq natatime);
 use POSIX;
 use Date::Tiny;
 use GenBoNoSqlRocksGenome;
 use Compress::LZ4;
my $version;
my $fork;
my $chr_name;
my $genome_version;
my $merge;
GetOptions(
	'chr=s' => \$chr_name,
	'version=s' => \$version,
	'genome=s' => \$genome_version,
	'fork=s' => \$fork,
	'merge=s' => \$merge,
);
die("fork") unless $fork;
die("genome") unless $genome_version; 

#die("-type=genome or exome") unless $type_vcf;
die("-version= (/data-xfs/public-data/(genome_version)/vcf/gnomad/(version)/") unless $version;

my $dir_root = "/data-isilon/public-data/repository/$genome_version/spliceAI/$version";
my $file_snv = $dir_root."/tabix/spliceai_scores.masked.snv.$genome_version.vcf.gz";
my $file_indel = $dir_root."/tabix/spliceai_scores.masked.indel.$genome_version.vcf.gz";
#my $file = "/public-data/repository/HG19/spliceAI/$version/tabix/spliceai_scores.masked.snv.hg19.vcf.gz";
die($file_snv) unless -e $file_snv;



die($file_indel) unless -e $file_indel;
my $dir_out;
my $mode;

if($merge){
	$dir_out = "/data-isilon/public-data/repository/$genome_version/cadd-gnomad-spliceAI/rocksdb/$chr_name.g.rocks/";
	$mode = "w";
}
else {
	$dir_out = $dir_root."/rocksdb_split/";
	$mode ="c";
}


system("mkdir -p $dir_out") unless -e $dir_out;

my $description;
$description->{name} = "spliceAI";
$description->{version} = $version;;
$description->{file} = join(";",($file_snv,$file_indel));
 my $d = Date::Tiny->now;
$description->{date} =  $d->as_string;
$description->{code_unpack} =  "W4 C4";
open(JSON,">$dir_out/description.json");
print JSON encode_json $description;
close JSON;

#die("$dir_out/description.json");

my $pack = "C4 s4 A*";
my $factor = ["0.01","0.01","0.01","0.01","1","1","1","1",""];
my $description_a = ["AG","AL","DG","DL","PAG","PAL","PDG","PDL","gene"];
my $pm = new Parallel::ForkManager($fork);

#DS_AG	Delta score (acceptor gain)
#DS_AL	Delta score (acceptor loss)
#DS_DG	Delta score (donor gain)
#DS_DL	Delta score (donor loss)
#DP_AG	Delta position (acceptor gain)
#DP_AL	Delta position (acceptor loss)
#DP_DG	Delta position (donor gain)
#DP_DL	Delta position (donor loss)
my $rg = GenBoNoSqlRocksGenome->new(dir=>$dir_out,mode=>"t",chromosome=>$chr_name,index=>"genomic",genome=>$genome_version,factor=>$factor,pack=>$pack, description=>$description_a);#,pack=>"C4 s4");
#my $lmdb_out  =  $save_out->get_lmdb_db("snps");
my $nb_regions = scalar(@{$rg->regions});
warn $nb_regions;
 	my $eta = Time::ETA->new(
    milestones =>$nb_regions,
);

my $jobs;

$pm->run_on_finish (
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hres) = @_;
			warn Dumper $jobs;
			warn $hres->{job};
			die() unless exists $jobs->{$hres->{job}};
			delete $jobs->{$hres->{job}};
			 $eta->pass_milestone();
 			my $remaining = $eta->get_remaining_time();
			my $elapsed = $eta->get_elapsed_time();
			my $p = int($eta->get_completed_percent());
			print $chr_name." complete $p % => elapsed :".$elapsed." remaining : ".$remaining."\n";
		});

my $zzid = time;
foreach my $r (@{$rg->regions}){
	$zzid ++;
	$jobs->{$zzid} =1;
	my $pid = $pm->start and next;
		warn "start";
	my $hv = run_region2($r,$file_snv);

	 $hv = run_region2($r,$file_indel,1);
	 $rg->nosql_rocks($r)->close;
	 #last;
	my $hres;
	$hres->{variants} = $hv;
	$hres->{job} = $zzid;
	warn "end $zzid";
	$pm->finish(0,$hres);

	
}

$pm->wait_all_children;
confess(Dumper $jobs) if scalar (keys %$jobs) >0;
  warn "FINISH ";
exit(0);

sub run_region2 {
	my ($region,$vcf,$add) = @_;
	my $no = $rg->nosql_rocks($region);
	$no->intspan_keys;
	my $tabix  = Bio::DB::HTS::Tabix->new( filename => $vcf );
	# my $tabix  = new Tabix(-data =>$file);
	 my $res = $tabix->query($region->{tabix});
	  my $pos;
	  my $xx = 1;
	  my $old_id = 0;
	  my $vs = {};
	  my $hgenes = {};
	  my $gene_id = 0;
	  my $uid;
	  my $old_value;
	  my $genes;
	  my $array = [];
	 while ( my $line = $res->next ) {
	  	 		chomp($line);
				my(@array) = split("\t",$line);
				
				my $alternate_allele = $array[4];
				
				my $ref_allele = $array[3];
				my $start = $array[1];
				my @sAi = split(/\|/,$array[7]);
				my $gene = $sAi[1];
				my $id = $start;
				my $debug;
				$debug = 1 if $start == 178512799;
				my $zid = compress_vcf_position($ref_allele,$alternate_allele);
				$uid = $no->return_rocks_id($id,$ref_allele,$alternate_allele);
				#my $tab = [ceil($sAi[2]*100),ceil($sAi[3]*100),ceil($sAi[4]*100),ceil($sAi[5]*100),$sAi[6],$sAi[7],$sAi[8],$sAi[9]];
				if ($sAi[2] ==0 && $sAi[3] ==0 && $sAi[4] ==0 && $sAi[5] ==0){
					$no->intspan_keys->add($id);
					next;					
				}
				my $pp = pack($pack,ceil($sAi[2]*100),ceil($sAi[3]*100),ceil($sAi[4]*100),ceil($sAi[5]*100),$sAi[6],$sAi[7],$sAi[8],$sAi[9],uc($gene));
				if ($old_id ne $uid and $old_id ne 0){
    				$no->put_batch($old_id,$array);
    				$array=[];
    			}
    			$old_id = $uid;
    		
    			push(@$array,$pp);	
				
	}
	$no->put_batch($old_id,$array) if @$array; 
	$no->rocks->write($no->batch);
	delete $no->{batch};
	$no->rocks->compact_range();
	return ;

}



sub compress_vcf_position {
	my ($ref,$alt) = @_;
	my $l1 = length($ref);
	my $l2 = length($alt);
	if ($l1 ==1 && $l2 ==1){
		return $alt;
	}
	elsif ($l1 ==1 && $l2 > 1){
		return "+".substr($alt, 1);
	}
	elsif ($l1 >1 ){
		$ref = substr($ref, 1);
		return ($l1 -1);
	}
	confess();
	
}
