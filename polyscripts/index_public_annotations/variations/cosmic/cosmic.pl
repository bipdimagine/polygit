#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use strict; 
use Set::IntSpan::Fast::XS ;
use Data::Dumper;
use Parallel::ForkManager;
use Vcf;
use GenBoNoSqlRocksAnnotation;
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

my $dir_rocks;
$dir_rocks = "$dir_out/rocksdb/";

unless (-e $dir_rocks){
	system("mkdir $dir_rocks;chmod a+rwx $dir_rocks");
}

 my @files = `ls $dir_vcf/*.vcf.gz`;
 chomp(@files);

 



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
	run_chr_rocks($chr);
	my $hres;
	$hres->{id} = $this_id;
	$pm->finish(0, $hres);
	last;
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




sub run_chr_rocks {
	my ($chr) = @_;

	#my $db = new save({name=>"cosmic",chromosome=>$chr,mode=>"c",integer=>1,compress=>1});
	my $nop = GenBoNoSqlRocksAnnotation->new(dir=>$dir_rocks,factor=>[1],pack=>"",mode=>"c",name=>$chr,version=>"110",pipeline=>1,description=>["rs","frequence"]);	
	

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
					$variation->{start} =  $$x{'POS'};
					$variation->{rs} = $x->{ID};
						my $gnomad_id = join("-",$chr,$variation->{start},$variation->{ref_all},$vcfVarAllele);
					$nb_alt ++;
					my $type;
					my $len;
					($type,$len, $variation->{ht}) = $vcf->event_type($vcfRefAllele."", $vcfVarAllele);
					$previous_allele = $variation->{ht};
					my $lenVar = length($vcfVarAllele);
					my $id;
							$variation->{len} = abs($len);
					my $id;
					if ($type eq 'r') { warn "# File: ".$file."\n -> No reason to parse this variant because it's identical to reference )"; warn Dumper $x }
					elsif ($type eq 'o') { next;warn($$x{'REF'}." ". $vcfVarAllele." ".$variation->{ht}." ->$type");next;}
					$variation->{start} =  $$x{'POS'};
					$nop->put_cosmic($gnomad_id,[$x->{ID},$cnt]);
					
				}# for index_alt
				$nb++;		
			 	if ($nb%10000 == 0){
			 	
			 		  print $nb." => $chr \t"." ".abs(time-$t)."\n" ;
			 		  $t =time;
			 	}
			
		
    }# while parse vcf
	$vcf->close();
}
	$nop->close();
	return ;
}

