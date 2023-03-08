#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../../../../../lib/obj-nodb";
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
my $pm = new Parallel::ForkManager(24);
my $buffer = GBuffer->new();
my $version;
GetOptions(
	'version=s'   => \$version,
);
die("please add -version=") unless $version;
my $dir_vcf  = "/public-data/repository/HG19/hgmd/$version/vcf/";
  my @files = `ls $dir_vcf/*.vcf.gz`;
  chomp(@files);
  my $file = $files[0];
  my $md5_1 = file_md5_hex( $file);
  warn $file." ".$md5_1;

die("clinar file doesn t exists " .$file ) unless -e $file;
my $dir_out  = "/public-data/repository/HG19/hgmd/$version/";

if (-e $dir_out."/lmdb"){
	die("hey the output directory already exists !!! $dir_out/lmdb");
}


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
	my $db = new save({name=>"lmdb",chromosome=>$chr,mode=>"c",integer=>1,db_type=>"lmdb"});
	$db->dir_lmdb($dir_out);
	#my $db = new save({name=>"dbsnp",chromosome=>$chr});
	my @objs;
	#   warn $file;
	 #  $file =~ s/dibayes/unifiedgenotyper/;
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
					$variation->{fields}->{class} = $x->{INFO}->{CLASS};
					$variation->{fields}->{rsname} = $x->{INFO}->{DB};
					$variation->{fields}->{phen} = $x->{INFO}->{PHEN};
					$variation->{fields}->{hgmd_score} = $x->{INFO}->{RANKSCORE};
					$variation->{fields}->{hgmd_id} =  $$x{'ID'};
					$variation->{start} =  $$x{'POS'};
					
					my $hrelation_variant_gene;
					$hrelation_variant_gene->{type_id} = 'hgmd_id';
					$hrelation_variant_gene->{hgmd_id} = $variation->{fields}->{hgmd_id};
					$hrelation_variant_gene->{gene} = $x->{INFO}->{GENE};
					$db->add_relation_variant_gene($hrelation_variant_gene);
	
					if ($type eq 's') {
						# SNP
						$db->add_snp($variation);
					}
					elsif ($type eq "i" && $len < 0 ){
						# deletions
					
						$db->add_deletion($variation);
						
					
					}
					elsif ($type eq "i" && $len >= 0 ){
						# insertion
						$db->add_insertion($variation);
					}
					elsif ($type eq 'r') {   $vcfVarAllele = $vcfRefAllele }
					elsif ($type eq 'o') { next;warn($$x{'REF'}." ". $vcfVarAllele." ".$variation->{ht}." ->$type");next;}
					else {warn ("\n\nERROR $type: strange vcf line...\nFile:".$file."\n$x\n"); next;}
					
				#	warn $id;
					
				
						
				}# for index_alt
				$nb++;		
			 	if ($nb%100000 == 0){
			 	
			 		
			 		  print $nb." => $chr \t"." ".abs(time-$t)."\n" if $nb%100000 == 0;
			 		  $t =time;
			 	}
			
		
    }# while parse vcf
	$vcf->close();
	#warn Dumper (keys %sig); 
	$db->save_intspan();
	#die() if $debug;
	$db->close();
	
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



