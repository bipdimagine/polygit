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
GetOptions(
	'version=s'   => \$version,
);
die("please add -version=") unless $version;
my $dir_vcf  = "/data-isilon/public-data/repository/HG19/clinvar/$version/vcf/";
  my @files = `ls $dir_vcf/*.vcf.gz`;
  chomp(@files);
  my $file = $files[0];
  my $md5_1 = file_md5_hex( $file);
  warn $file." ".$md5_1;

die("clinar file doesn t exists " .$file ) unless -e $file;
my $dir_out = "/tmp/clinvar/$version/";
system("mkdir -p /tmp/clinvar/$version/");

#run_chr("1");
#die();
my $hrun;
my $id = time;
my $nb;


my $max_gencode = $buffer->getQuery->getMaxGencodeVersion( );
my $gencode_dir = $buffer->public_data_root().'/HG19/'.$buffer->gencode->{$max_gencode}->{directory};
my $projname = $buffer->get_random_project_name_with_this_annotations_and_genecode();
my $omim_dir = $buffer->newProject(-name => $projname)->get_public_data_directory("omim");

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
	warn "end $chr";
	$pm->finish(0,{id=>$id."_".$chr});
}#end chr
  $pm->wait_all_children;

my $hversion;
$hversion->{name} = "clinvar";
$hversion->{version} = "$version";
$hversion->{file} = $file;
  my $d = Date::Tiny->now;

$hversion->{date} =  $d->as_string;
$hversion->{md5sum} = $md5_1;
open(JSON,">$dir_out/lmdb/version.json") or die();
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
exit(0);





sub run_chr {
	my ($chr) = @_;



die() unless -e $file;
	my $db = new save({name=>"lmdb",chromosome=>$chr,mode=>"c",integer=>1,db_type=>"lmdb"});
	$db->dir_lmdb($dir_out);



my $no_genes_syno = GenBoNoSqlAnnotation->new(
	dir  => $gencode_dir,
	mode => "r"
);

my $no_genes_obj = GenBoNoSqlLmdb->new(
	name        => "genbo",
	dir         => $gencode_dir,
	mode        => "r",
	is_compress => 1,
	vmtouch		=> $buffer->vmtouch
);

my $no_omim = GenBoNoSqlLmdb->new(
	name        => "omim",
	dir         => $omim_dir,
	mode        => "r",
	is_compress => 1,
	vmtouch=>$buffer->vmtouch
);

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
				$debug = 1 if $varStart == 110395675;
			
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
				
					warn "1" if $debug;
					
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
					$variation->{clinical} =2 if  exists $x->{INFO}->{PM};
					$variation->{id} =  $x->{INFO}->{ALLELEID};
					$variation->{clinvar_id} =  $x->{INFO}->{ALLELEID};
					#warn  $x->{INFO}->{ALLELEID};
					$variation->{clinvar_id};
					die() unless  $x->{INFO}->{ALLELEID};
					$variation->{id} =~s/\|/;/g;
					my $sig = parse_clnsig( $x->{INFO}->{CLNSIG});
					$variation->{sig} =  $sig;
					#$variation->{sig} =~s/\|/;/g;
					$variation->{start} =  $$x{'POS'};
					
					
					my @l_genes_infos = split('\|', $x->{INFO}->{GENEINFO});
					my $CLNVI_OMIM_gene_id;
					if (scalar(@l_genes_infos) > 1 and exists $x->{INFO}->{CLNVI}) {
						foreach my $clnvi (split('\|', $x->{INFO}->{CLNVI})) {
							my @lTmp_clnvi = split(':', $clnvi);
							if ($lTmp_clnvi[0] eq 'OMIM_Allelic_Variant') {
								my @lTmp_clnvi_omim = split('\.', $lTmp_clnvi[1]);
								$CLNVI_OMIM_gene_id = $lTmp_clnvi_omim[0];
							}
						}
					}
					foreach my $gene_info (@l_genes_infos) {
						my @lTmpG = split(':',$gene_info);
						my $gene_external_name = $lTmpG[0];
						my $gene_omim_id = undef;
						if ($CLNVI_OMIM_gene_id) {
							my $is_ok_check_CLNVI_OMIM_gene_id;
							my $h_synos = $no_genes_syno->get_key_values("synonyms", $gene_external_name);
							my $has_done;
							foreach my $gene_syno (values %{$h_synos}) {
								next if ($is_ok_check_CLNVI_OMIM_gene_id);
								my $gene;
								eval { $gene = $no_genes_obj->get($gene_syno); };
								if ($@) { next; }
								next unless ($gene);
								confess("\n\n".$gene->external_name()." ne $gene_external_name\n\n") unless ($gene->external_name() eq $gene_external_name);
								my $h_omim_gene = $no_omim->get($gene->id());
								my $gene_omim_id = $h_omim_gene->{omim_id};
								if ($gene_omim_id eq $CLNVI_OMIM_gene_id) {
									my $hrelation_variant_gene;
									$hrelation_variant_gene->{type_id} = 'clinvar_id';
									$hrelation_variant_gene->{clinvar_id} = $variation->{clinvar_id};
									$hrelation_variant_gene->{gene} = $gene_external_name;
									$db->add_relation_variant_gene($hrelation_variant_gene);
									$has_done = 1;
								}
							}
							unless ($has_done) {
								my $hrelation_variant_gene;
								$hrelation_variant_gene->{type_id} = 'clinvar_id';
								$hrelation_variant_gene->{clinvar_id} = $variation->{clinvar_id};
								$hrelation_variant_gene->{gene} = $gene_external_name;
								$db->add_relation_variant_gene($hrelation_variant_gene);
							}
						}
						else {
							my $hrelation_variant_gene;
							$hrelation_variant_gene->{type_id} = 'clinvar_id';
							$hrelation_variant_gene->{clinvar_id} = $variation->{clinvar_id};
							$hrelation_variant_gene->{gene} = $gene_external_name;
							$db->add_relation_variant_gene($hrelation_variant_gene);
						}
					}
					
					
	#				warn 
					warn "2" if $debug;
					if ($type eq 's') {
						# SNP
						warn "3" if $debug;
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
					
				
						warn "4" if $debug;
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
	$no_genes_syno->close();
	$no_genes_obj->close();
	$no_omim->close();
	
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
	"association_not_found" => -9,
	"other" => -9,
	
};
	my ($st1,$st2) = split(",",$string);
	my @values = split("/",$st1);
	
	my @res;
	foreach my $v (@values){
		
		$v = lc ($v);
		
		my $scorev = $score->{$v};
		unless (exists $score->{$v}){
			warn " !!!!!!!!".$v ;
			 $scorev = -9 ; 
		}
		
		
		push(@res,$scorev);
		#warn $v unless exists $score->{$v};
		
		#die($v) unless exists $score->{$v};
	}
	return join(";",@res);
}



