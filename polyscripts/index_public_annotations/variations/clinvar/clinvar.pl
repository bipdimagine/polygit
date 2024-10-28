#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../../../../GenBo/lib/obj-nodb";
use strict; 
use Bio::SearchIO;
use strict;

use Data::Printer;
use Getopt::Long;
use Carp;

use Data::Dumper;
use Getopt::Long; 
use Storable qw/retrieve freeze thaw nfreeze nstore_fd nstore/;
use GBuffer;
use Sys::Hostname;
use Parallel::ForkManager;
use Vcf;
use JSON::XS;
use Digest::MD5::File qw( file_md5_hex );
use Date::Tiny;
#use GenBoNoSql;
#require "$Bin/dbsnp.pm";
 require "$Bin/../packages/save.pm";
my $allsnps;

my $chromosomes = [1..22,'X','Y','MT'];
#my $chromosomes = [2];
#my $dir = "/data-xfs/public-data/HG19/snp/exac/latest/";
my $pm = new Parallel::ForkManager(3);
my $buffer = GBuffer->new();
my $version;
my $genome_version;
GetOptions(
	'version=s'   => \$version,
	'genome=s' => \$genome_version,
);
die("please add -version=") unless $version;
die("please add  genome") unless $genome_version; 
my $dir_public = "/data-isilon/public-data/repository/$genome_version/clinvar/$version";
my $dir_vcf  = $dir_public."/vcf/";
  my @files = `ls $dir_vcf/*.vcf.gz`;
  chomp(@files);
  my $file = $files[0];
  my $md5_1 = file_md5_hex( $file);
  warn $file." ".$md5_1;

die("clinar file doesn t exists " .$file ) unless -e $file;
my $dir_out= $dir_public."/rocksdb/";

system("mkdir $dir_out") unless -d $dir_out;
#run_chr("1");
my $hrun;
my $id = time;
my $nb;


my $no_patho = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out,mode=>"c",name=>"clinvar_pathogenic",version=>$version,pipeline=>1);
foreach my $chr (@$chromosomes){
	my $no = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out,mode=>"c",name=>"$chr",pack=>[1,1],version=>$version,description=>["clnsig","genes","clinvar_id","score"],pipeline=>1);
	warn "start $chr";
	run_chr($chr,$no);
	warn "end $chr";
	$no->write_batch;
	$no->close();
}#end chr
my $hversion;
$hversion->{name} = "clinvar";
$hversion->{version} = "$version";
$hversion->{file} = $file;
  my $d = Date::Tiny->now;

$hversion->{date} =  $d->as_string;
$hversion->{md5sum} = $md5_1;
open(JSON,">$dir_out/version.json") or die();
warn "$dir_out/clinvar.json";
print JSON encode_json $hversion;
close JSON;


if (keys %$hrun){
	warn Dumper $hrun;
	warn "ERROR ON YOUR LMDB DATABASE ";
	die();
}
warn Dumper $hrun;
warn "your database is here : ".$dir_out;

$no_patho->write_batch();
$no_patho->close();
exit(0);






sub run_chr {
	my ($chr,$no) = @_;



die() unless -e $file;

my $v = Bio::DB::HTS::VCF->new( filename => "$file" );
	
my $h = $v->header();
my $iter  = $v->query($chr.":1-".100_000_000_000); 
	my $b =0;
	return unless $iter;
	while (my $result = $iter->next) {
			$b ++;
			warn $b if $b%10000 ==0;
			my $alleles = $result->get_alleles();
		my $debug;
		# warn $result->position() ;
		my $varAllele = $alleles->[0];
		die() if scalar(@$alleles) >1;
		#if $debug;
		my $hfreq;
		my $variation;
		 my $ref_allele = $result->reference();

		my $type = $result->get_variant_type(1);
		foreach my $alternate_allele (@$alleles) {
		my $pos  = $result->position();
		 my ($score,$sig) = parse_clnsig(value_from_infos($result,$h, "CLNSIG"));
		 my $clinvar_id = value_from_infos($result,$h, "ALLELEID") ;
		 die() if $result->get_info($h, "PM") ne "ID_NOT_FOUND";
		 my $vcf_genes;
		 my $genes = value_from_infos($result,$h, "GENEINFO");
		 my @genes_infos = split (/\|/,$genes);
		 my @agenes;
		 foreach my $gl  (@genes_infos){
		 	my ($n,$n2) = split(":",$gl);
		 	push(@agenes,$no->get_id_dictionary($n));
		 }
		 my $id = sprintf("%010d", $pos);
		my $zids  = compress_vcf_position($ref_allele,$alternate_allele);
		foreach my $zid (@$zids){
			 $zid = $id."!".$zid;
			$no->put_batch($zid,[$no->get_id_dictionary($sig),\@agenes,$clinvar_id,$score]);
			 $no_patho->put_batch_raw($chr."!".$zid,1) if $score ==5;	
		}
		 
		}
	}
	$no->write_batch;
}

sub value_from_infos {
	my ($result,$h,$type) = @_;
	my $t = $result->get_info($h, $type);
	return "ID_NOT_FOUND" if $t eq "ID_NOT_FOUND";
	warn  $type." ".$result->id." ".Dumper $result->get_info($h, "ALLELEID")  if $t eq "ID_NOT_FOUND";
	die($type) if @$t > 1;
	return $t->[0];
}
				


sub compress_vcf_position {
	my ($ref,$alt) = @_;
	my $l1 = length($ref);
	my $l2 = length($alt);
	if ($l1 ==1 && $l2 ==1){
		return [$alt];
	}
	if ($l1 ==1 && $l2 > 1){
		return ["+".substr($alt, 1)];
	}
	if ($l1 >1 && $l2 == 1 ){
		$ref = substr($ref, 1);
		return [($l1 -1)];
	}
	if ($l1 == 2 && $l2 ==2){
		my @a =split("",$ref);
		my @b = split("",$alt);
		my $res;
		for(my $i=0;$i<@a;$i++ ){
			next if $a[$i] eq $b[$i];
			push(@$res,$alt);
		}
		return $res;
	}
	warn $ref." ".$alt;
	return [$ref.":".$alt];
	
}


sub parse_clnsig {
	my ($string) = @_;
	return (-9,"other") if $string eq "ID_NOT_FOUND";
	my $score ={
	"drug_response" =>1,
	"benign" => 2,
	"likely_benign" =>3,
	"likely_pathogenic" => 4,
	"pathogenic" => 5,
	"not_provided" => -1,
	"uncertain_significance" =>-2,
	"histocompatibility" => -3,
	"association" =>-4,
	"risk_factor" =>-5,
	"protective" =>-6,
	"affects" =>-7,
	"conflicting_interpretations_of_pathogenicity" => -8,
	"association_not_found" => -9,
	"uncertain_risk_allele" => -1,
	"likely_risk_allele" => 0,
	"other" => -9,
	
};
	my ($st1,$st2) = split(",",$string);
	my @values;
	push(@values,split( /[\|\/]/, $st1));
	
	my @res;
	my $maxv= -50;
	my $max_text ="other";
	foreach my $v (@values){
		
		$v = lc ($v);
		
		my $scorev = $score->{$v};
		
		unless (exists $score->{$v}){
			
			 $scorev = -9 ; 
		}
		if ($scorev > $maxv){
			$maxv = $scorev;
			$max_text = $v;
		}
		
	}
	return ($maxv,$max_text);
}



