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
	$dir_out = $dir_root."/rocksdb/$chr_name.g.rocks/";
	$mode ="c";
}


system("mkdir -p $dir_out") unless -e $dir_out;


my $description;
$description->{name} = "spliceAI";
$description->{version} = $version;;
$description->{file} = join(";",@files);
 my $d = Date::Tiny->now;
$description->{date} =  $d->as_string;
$description->{code_unpack} =  "W4 C4";
open(JSON,">$dir_out/description.json");
print JSON encode_json $description;
close JSON;

#die("$dir_out/description.json");

	
my $pm = new Parallel::ForkManager($fork);



my $rg = GenBoNoSqlRocksGenome->new(dir=>$dir_out,mode=>$mode,chromosome=>$chr_name,fasta=>$genome_fasta,compress=>1,index=>"genomic",genome=>$genome_version);#,pack=>"C4 s4");
#my $lmdb_out  =  $save_out->get_lmdb_db("snps");
my $nb_regions = scalar(@{$rg->regions});
 	my $eta = Time::ETA->new(
    milestones =>$nb_regions,
);


$pm->run_on_finish (
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hres) = @_;
		
			 $eta->pass_milestone();
 			my $remaining = $eta->get_remaining_time();
			my $elapsed = $eta->get_elapsed_time();
			my $p = int($eta->get_completed_percent());
			print $chr_name." complete $p % => elapsed :".$elapsed." remaining : ".$remaining."\n";
		});


foreach my $r (@{$rg->regions}){
	my $pid = $pm->start and next;
	my $hv = run_region($r,$file_snv);
	$merge = 1;
	 $hv = run_region($r,$file_indel,1);
	 $rg->nosql_rocks($r)->close;
	 #last;
	my $hres;
	$hres->{variants} = $hv;
	$pm->finish(0,$hres);

	
}

$pm->wait_all_children;

warn "END";

  
exit(0);





sub run_region {
	my ($region,$vcf) = @_;
	my $no = $rg->nosql_rocks($region);
	my $tabix  = Bio::DB::HTS::Tabix->new( filename => $vcf );
	# my $tabix  = new Tabix(-data =>$file);
	 my $res = $tabix->query($region->{tabix});
	  my $pos;
	  my $xx = 1;
	  my $old_id;
	  my $vs = {};
	  my $hgenes = {};
	  my $gene_id = 0;
	 while ( my $line = $res->next ) {
	  	 		chomp($line);
				my(@array) = split("\t",$line);
				
				my $alternate_allele = $array[4];
				my $ref_allele = $array[3];
				my $start = $array[1];
				my @sAi = split(/\|/,$array[7]);
				my $gene = $sAi[1];
				my $id = $start;
				unless ($old_id){
					$old_id = $id  unless $old_id;
					$vs = {};
					if ($merge){
						$vs = $no->get($id);
						$vs = {} unless $vs;
					}
				}
				if ($id ne $old_id ){
					#warn $id." ".Dumper $vs;
					$no->put_batch($old_id,$vs);# if @$vs;
					$old_id = $id;
					$vs = {};
					if ($merge){
						$vs = $no->get($id);
						$vs = {} unless $vs;
					}
				}
				my $zid = compress_vcf_position($ref_allele,$alternate_allele);
				if ($sAi[2] ==0 && $sAi[3] ==0 && $sAi[4] ==0 && $sAi[5] ==0){
					push(@{$vs->{$zid}},undef);
					next;
				}
				else {
					#$vs->{$ref_allele."_".$alternate_allele}->{s} = pack(",pack=>"C4 s4"")
					unless (exists $hgenes->{$gene}){
						$hgenes->{$gene} = $gene_id;
						$gene_id ++;
					}
					my $zid = compress_vcf_position($ref_allele,$alternate_allele);
					my $tab = [$gene,ceil($sAi[2]*100),ceil($sAi[3]*100),ceil($sAi[4]*100),ceil($sAi[5]*100),$sAi[6],$sAi[7],$sAi[8],$sAi[9]];
					push(@{$vs->{$zid}},$tab);
				}
	}
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
