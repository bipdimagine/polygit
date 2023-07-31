#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use strict; 
use Set::IntSpan::Fast::XS ;
use Data::Dumper;
use Parallel::ForkManager;
use Vcf;
use GenBoNoSql;
use JSON::XS;
use Digest::MD5::File qw( file_md5_hex );
use Date::Tiny;
use GBuffer;
use Getopt::Long; 
 require "$Bin/../packages/save.pm";

my $allsnps;

my $chromosomes = [1..22,'X','Y','MT'];
#my $chromosomes = [21];

my $pm = new Parallel::ForkManager(8);
my $buffer = GBuffer->new();
my $version;
GetOptions(
	'version=s'   => \$version,
);
die("please add -version=") unless $version;
my $dir_out  = "/data-isilon/public-data/repository/HG19/cosmic/$version/";
my $dir_vcf = "$dir_out/vcf/";


 my @files = `ls $dir_vcf/*.vcf.gz`;
 chomp(@files);

 


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
	my $this_id = $id."_".$chr;
	$hrun->{$this_id}++;
	my $pid = $pm->start and next;
	warn "start $chr";
	run_chr($chr);
	my $hres;
	$hres->{id} = $this_id;
	$pm->finish(0, $hres);
}
$pm->wait_all_children;

if (keys %$hrun){
	warn Dumper $hrun;
	warn "ERROR ON YOUR LMDB DATABASE HGMD";
	die();
}
my $hversion;
$hversion->{name} = "cosmic";
$hversion->{version} = "$version";
$hversion->{file} = join(";",@files);
  my $d = Date::Tiny->now;

$hversion->{date} =  $d->as_string;
my @md5;
foreach my $f (@files){
	push(@md5,file_md5_hex( $f));
}

$hversion->{md5sum} = join(";",@md5);
warn "$dir_out/lmdb/";
open(JSON,">$dir_out/lmdb/version.json") or die();
print JSON encode_json $hversion;
close JSON;



warn "your database is here : ".$dir_out;


exit(0);




sub run_chr {
	my ($chr) = @_;

	#my $db = new save({name=>"cosmic",chromosome=>$chr,mode=>"c",integer=>1,compress=>1});
	my $db = new save({name=>"lmdb",chromosome=>$chr,mode=>"c",integer=>1,db_type=>"lmdb"});
	$db->dir_lmdb($dir_out);

foreach my $file (@files){
	my @objs;
	#   warn $file;
	 #  $file =~ s/dibayes/unifiedgenotyper/;
	 die() unless -e $file;
	    my $vcf = Vcf -> new (file=>$file, region=>$chr,tabix=>"tabix"); 
	  	eval {
	    	$vcf->parse_header();
			
	  	};
	 my $nb =0;
	 my $t =time;
		while (my $x = $vcf->next_data_hash()) {

				my $chrom = $$x{'CHROM'};
				my $varStart = $$x{'POS'};
				my $debug;
				my $varEnd;
				my $vcfRefAllele = $$x{'REF'};
				my $varAllele = $$x{'ALT'};
				my $alleles;
	
				my $nb_alt =0;
				my $halleles;		
				my $previous_allele;
				$x->{ID} =~ s/\.//;
				my $cnt =  $x->{'INFO'}->{'CNT'};
				foreach my $vcfVarAllele (@$varAllele) {
					$vcfRefAllele = $$x{'REF'};
					my $variation;
					$variation->{ref_all} = $$x{'REF'};
					$variation->{alt_all} = $vcfVarAllele;
					$variation->{frequence} = $cnt;
					$variation->{rs} = $x->{ID};
					$nb_alt ++;
					my $type;
					my $len;
					($type,$len, $variation->{ht}) = $vcf->event_type($vcfRefAllele."", $vcfVarAllele);
					$previous_allele = $variation->{ht};
					my $lenVar = length($vcfVarAllele);
					my $id;
							$variation->{len} = abs($len);
					my $id;
					$variation->{start} =  $$x{'POS'};
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
					elsif ($vcfVarAllele eq "." ){
							$variation->{ht} = $variation->{ref_all} ;
							$db->add_deletion($variation,1);
					}
					elsif ($type eq 'r') { warn "# File: ".$file."\n -> No reason to parse this variant because it's identical to reference )"; warn Dumper $x }
					elsif ($type eq 'o') { next;warn($$x{'REF'}." ". $vcfVarAllele." ".$variation->{ht}." ->$type");next;}
					else {next;warn ("\n\nERROR $type: strange vcf line...\nFile:".$file."\n$x\n"); next;}
					
				}# for index_alt
				$nb++;		
			 	if ($nb%10000 == 0){
			 	
			 		  print $nb." => $chr \t"." ".abs(time-$t)."\n" ;
			 		  $t =time;
			 	}
			
		
    }# while parse vcf
	$vcf->close();
}
	$db->save_intspan();
	$db->close();
	return ;
}

