#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../../../../../../GenBo/lib/obj-nodb";
use strict; 
use Bio::SearchIO;
use strict;
use Data::Printer;
use Getopt::Long;
use Carp;
use GenBoNoSqlRocksAnnotation;
use Data::Dumper;
use Getopt::Long; 
use GBuffer;
use Parallel::ForkManager;
use Vcf;
use JSON::XS;
use Digest::MD5::File qw( file_md5_hex );
use Date::Tiny;
 require "$Bin/../packages/rocks_util.pm";
my $allsnps;

my $chromosomes = [1..22,'X','Y','MT'];
#my $chromosomes = [2];
#my $dir = "/data-xfs/public-data/HG19/snp/exac/latest/";
my $pm = new Parallel::ForkManager(24);
my $buffer = GBuffer->new();
my $version;
my $genome;
GetOptions(
	'version=s'   => \$version,
	'genome=s'   => \$genome,
);
die("please add -version=") unless $version;
my $dir_vcf  = "/data-isilon/public-data/repository/$genome/hgmd/$version/vcf/";
  my @files = `ls $dir_vcf/*.vcf.gz`;
  chomp(@files);
  my $file = $files[0];
  my $md5_1 = file_md5_hex( $file);
  warn $file." ".$md5_1;

die("clinar file doesn t exists " .$file ) unless -e $file;
my $dir_out  = "/data-isilon/public-data/repository/$genome/hgmd/$version/rocksdb/";




my $hrun;
my $id = time;
my $nb;

$pm->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    	
		delete $hrun->{$data->{id}};		    
    	
    }
  );

foreach my $chr (@$chromosomes){
	$id ;
	$hrun->{$id."_".$chr}++;
	my $pid = $pm->start and next;
	warn "start $chr";
	run_chr($chr);
	$pm->finish(0,{id=>$id."_".$chr});
}#end chr
  $pm->wait_all_children;


if (keys %$hrun){
	warn Dumper $hrun;
	warn "ERROR ON YOUR LMDB DATABASE HGMD";
	die();
}
my $hversion;
$hversion->{name} = "hgmd";
$hversion->{version} = "$version";
$hversion->{file} = $file;
  my $d = Date::Tiny->now;

$hversion->{date} =  $d->as_string;
$hversion->{md5sum} = $md5_1;
open(JSON,">$dir_out/lmdb/version.json") or die();
warn "$dir_out/hgmd.json";
print JSON encode_json $hversion;
close JSON;






warn "your database is here : ".$dir_out;
exit(0);





sub run_chr {
	my ($chr) = @_;



die() unless -e $file;
	my $finalrg = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out,mode=>"t",name=>$chr,pack=>[],version=>"",description=>{},pipeline=>1);
	$db->dir_lmdb($dir_out);
	my @objs;
	 my $vcf = Vcf -> new (file=>$file, region=>$chr,tabix=>"tabix"); 
	  	eval {
	    	$vcf->parse_header();
			
	  	};
	 my $nb =0;
	 	my $hvariation;
	 	my %sig;
			 my $t = time;
		while (my $x = $vcf->next_data_hash()) {
				my $chrom = $$x{'CHROM'};
				my $varStart = $$x{'POS'};
				my $debug;
				my $varEnd;
				my $vcfRefAllele = $$x{'REF'};
				my $varAllele = $$x{'ALT'};
				my $alleles;
				
				my @af;
				
				my $nb_alt =0;
				my $halleles;		
				my $previous_allele;
				$x->{ID} =~ s/\.//;
		
			#	p $x if scalar(@$varAllele) >2 && $caf;
				#die()  if scalar(@$varAllele) >2 && $caf;
				my $nb_all = 0;
				my $stype="";
				foreach my $vcfVarAllele (@$varAllele) {
					$sig{$x->{INFO}->{CLNSIG}} ++;
				#	 warn $x->{INFO}->{CLNSIG} if $x->{INFO}->{CLNSIG} =~/\// ;
				
					 
					$nb_all ++;
					$vcfRefAllele = $$x{'REF'};
					
					my $variation;
					$variation->{ref_all} = $$x{'REF'};
					$variation->{alt_all} = $vcfVarAllele;
					$variation->{clinical} = 1 ;
					
					$variation->{rs} = $x->{ID};
					$nb_alt ++;
					my $type;
					my $len;
					($type,$len, $variation->{ht}) = $vcf->event_type($vcfRefAllele."", $vcfVarAllele);
					$previous_allele = $variation->{ht};
					my $lenVar = length($vcfVarAllele);
					my $id;
					$variation->{class} = $x->{INFO}->{CLASS};
					$variation->{rsname} = $x->{INFO}->{DB};
					$variation->{phen} = $x->{INFO}->{PHEN};
					$variation->{hgmd_score} = $x->{INFO}->{RANKSCORE};
					$variation->{hgmd_id} =  $$x{'ID'};
					$variation->{gene} =  $x->{INFO}->{GENE};
					my $rocks_id = rocks_util::rocks_id($$x{'POS'},$$x{'REF'},$vcfVarAllele);
					$no->put_batch($rocks_id,$variation);						
						
				}# for index_alt
				$nb++;		
			 	if ($nb%100000 == 0){
			 	
			 		
			 		  print $nb." => $chr \t"." ".abs(time-$t)."\n" if $nb%100000 == 0;
			 		  $t =time;
			 	}
			
		
    }# while parse vcf
    $no->close();
	#warn Dumper (keys %sig); 
	
	return \@objs;
}




sub parse_clnsig {
	my ($string) = @_;
	my $score ={

	"drug_response" =>1,
	"benign" => 2,
	"likely_benign" =>3,
	"likely_pathogenic" => 4,
	"pathogenic" => 5,
	"not_provided" => -1,
	"uncertain_significance" =>-2,
	" histocompatibility" => -3,
	"association" =>-4,
	"risk_factor" =>-5,
	"protective" =>-6,
	"affects" =>-7,
	"conflicting_interpretations_of_pathogenicity" => -8,
	"other" => -9,
	
};
	my ($st1,$st2) = split(",",$string);
	my @values = split("/",$st1);
	my @res;
	foreach my $v (@values){
		$v = lc ($v);
		push(@res,$score->{$v});
		die($v) unless exists $score->{$v};
	}
	return join(";",@res);
}



