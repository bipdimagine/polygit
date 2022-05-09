#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin/../../../../../lib/obj-nodb";
use strict; 
use Bio::SearchIO;
use strict;
use Tabix;
use Data::Printer;
use Getopt::Long;
use Carp;
use JSON::XS;
use Data::Dumper;
use Getopt::Long; 
use Storable qw/retrieve freeze thaw nfreeze nstore_fd nstore/;
use Set::IntSpan::Fast::XS ;
use Sys::Hostname;
use GBuffer;
use GenBoProject;
use  Set::IntSpan::Island;
use Parallel::ForkManager;
use Vcf;
use GenBoNoSql;
 use Time::ETA;
use List::MoreUtils qw(any uniq natatime);
use GenBoNoSqlLmdb; 
use GenBoNoSqlLmdbInteger;
 use POSIX;
 use Date::Tiny;
my $allsnps;
use Getopt::Long;
my $c;
my $version;
GetOptions(
	'chr=s' => \$c,
	'version=s' => \$version,
);

my @files =("/public-data/repository/HG19/spliceAI/$version/tabix/spliceai_scores.masked.snv.hg19.vcf.gz","/public-data/repository/HG19/spliceAI/$version/tabix/spliceai_scores.masked.indel.hg19.vcf.gz");
my $file = "/public-data/repository/HG19/spliceAI/$version/tabix/spliceai_scores.masked.snv.hg19.vcf.gz";
die() unless -e $file;

my $buffer = new GBuffer;
my $project_name= "NGS2022_5187";

my $dir_out = "/tmp/lmdb/spliceAI/$version/";

system("mkdir -p $dir_out") unless -e $dir_out;

my $description;
$description->{name} = "spliceAI";
$description->{version} = $version;;
$description->{file} = $file;
 my $d = Date::Tiny->now;
$description->{date} =  $d->as_string;
$description->{code_unpack} =  "W4 C4";
#$description->{order} =  {'A'=>0,'C'=>1,'G'=>2,'T'=>3};
$description->{scale} =  100;
open(JSON,">$dir_out/description.json");
#open(JSON,">$dir_out/description.json");
print JSON encode_json $description;
close JSON;
#die("$dir_out/description.json");
my $project = $buffer->newProject( -name => $project_name,);
my $nb;
my $size = 500000;

	my $chr = $project->getChromosome($c);
	my $no = GenBoNoSqlLmdb->new(dir=>$dir_out,mode=>"c",name=>$chr->name,is_compress=>0,is_integer=>1);
	$no->create();
	$no->close();

	my $length = $chr->length();
	my $region;
	
		
		my $from =1;
    		my $to = $length;
  
    	my $start;
    	my $end;
    	
    	while ($from < $to){
        $start = $from;
        $end = $from + $size;
        if ($end > $to){
            $end = $to;
        }
        push(@$region,{start=>$from,end=>$end});
      #  push(@$res,{start=>$from,end=>$end,intspan=>$intspan_region,ext_gvcf=>$self->ucsc_name.".$from.$end.g.vcf",chromosome=>$self->ucsc_name}) unless $intspan_region->is_empty;
      
        #print chrom_name + ":" + str(region_start) + "-" + str(end)
        $from = $end;
       }
	

	
my $pm = new Parallel::ForkManager(10);


#my $lmdb_out  =  $save_out->get_lmdb_db("snps");
my $nbb =0;
my $nb_regions = scalar(@$region);
 	my $eta = Time::ETA->new(
    milestones =>$nb_regions,
);
my $intspan =  Set::IntSpan::Fast::XS->new();
my $intspan2 =  Set::IntSpan::Fast::XS->new();
my $hb;
	$pm->run_on_finish (
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hres) = @_;
			$nbb ++;
			my $cc= $hres->{chr};
			
		#	
		#	my $lmdb_out3  =  $save_out->get_lmdb_db("snps3");
			my $data = $hres->{data};
			
			#warn Dumper $data;
		#	die();
			return unless $data;
			my $t = time;
			my $nb =0;
		
			my $no = GenBoNoSqlLmdb->new(dir=>$dir_out,mode=>"w",name=>$chr->name,is_compress=>1,is_integer=>1);
			foreach my $value (keys %{$data->{hb}}){
				foreach my $alt (keys %{$data->{hb}->{$value}}){
					$hb->{$value}->{$alt} = $data->{hb}->{$value}->Clone;# if exists $hb->{$value}->{$alt};
				}
			}	
			foreach my $k (sort {$a <=> $b} keys %$data){
				$nb ++;
				#warn $k;
			
#				for (my $i=0;$i<4;$i++){
#					#warn Dumper $data->{$k}->[$i];
#					foreach my $gene (keys %{$data->{$k}->[$i]}){
#						
#					#warn $gene;
#					}
#				}
#				
				#my $rec = pack("C$nb",@{$data->{$k}->{s}});
			
				if ($intspan2->contains($k)){
					my $h = $no->get($k);
					die() unless $h;
					
					my %hh =(%$h,%{$data->{$k}});
				
						$no->put($k,\%hh);
					
				}
				else {
					$no->put($k,$data->{$k});
				}
			#	$no->put($k,$data->{$k});
				$intspan2->add($k);
			}
			 $eta->pass_milestone();
 			my $remaining = $eta->get_remaining_time();
			my $elapsed = $eta->get_elapsed_time();
			my $p = int($eta->get_completed_percent());
			$intspan = $intspan->union($hres->{intspan});
			print $chr->name." complete $p % => elapsed :".$elapsed." remaining : ".$remaining."\n";
			$no->close;
		});


#	$size = int($length/10);
	my $i;
	my $cname = $chr->name();
	
	while (my $r = shift(@$region)){
	my $pid = $pm->start and next;
	my ($array,$intspan,$hb) = run_chr($cname,$r->{start},$r->{end},$chr->length);
	my $hres;
	$hres->{chr} = $cname;
	$hres->{data} = $array;
	$hres->{intspan} = $intspan;	
	$hres->{hb} = $hb;	
	$pm->finish(0,$hres);
	}
	
  $pm->wait_all_children;
  warn $i." ".$length;

my $no = GenBoNoSqlLmdb->new(dir=>$dir_out,mode=>"w",name=>$chr->name,is_compress=>1,is_integer=>1);
$no->put(0,$intspan);
$no->close();
#$lmdb_out->close();
exit(0);





sub run_chr {
	my ($chr,$start,$end,$length) = @_;

die() unless -e $file;
my $ACGT = {'A'=>0,'C'=>1,'G'=>2,'T'=>3};
my $nb= 0;
	my $obj;
my $variation ={};
 my $intspan =  Set::IntSpan::Fast::XS->new();
 my $hash_bivector;
foreach my $file (@files){
	 my $tabix  = new Tabix(-data =>$file);
	  my $res = $tabix->query($chr,$start,$end);
	
	  	 while(my $line = $tabix->read($res)) {
	  	 		chop($line);
				my(@array) = split("\t",$line);
				
				
				my $alternate_allele = $array[4];
				my $start = $array[1];
				my @sAi = split(/\|/,$array[7]);
				my $ref_allele = $array[3];
				if (length($ref_allele)<length($alternate_allele)){
					#$alternate_allele = 
					substr($alternate_allele, 0, 1) = "+";

					$start += 1;
				}
				elsif (length($ref_allele)>length($alternate_allele)){
					#$alternate_allele = 
					substr($ref_allele, 0, 1) = "-";
					$start += 1;
					$alternate_allele = $ref_allele;
				}
				
			#	my $id = $ACGT->{$alternate_allele};
				#push(@{$variation->{$start}->[$id]},pack('w4 c4 A*',ceil($sAi[2]*100),ceil($sAi[3]*100),ceil($sAi[4]*100),ceil($sAi[5]*100),$sAi[6],$sAi[7],$sAi[8],$sAi[9],$sAi[1]));
				#$variation->{$start}->{$alternate_allele}->{$sAi[1]} = pack('w4 c4',ceil($sAi[2]*100),ceil($sAi[3]*100),ceil($sAi[4]*100),ceil($sAi[5]*100),$sAi[6],$sAi[7],$sAi[8],$sAi[9]);#[ceil($sAi[2]*100),ceil($sAi[3]*100),ceil($sAi[4]*100),ceil($sAi[5]*100),$sAi[6],$sAi[7],$sAi[8],$sAi[9]];
				if ($sAi[2] ==0 && $sAi[3] ==0 && $sAi[4] ==0 && $sAi[5] ==0){
					
					#warn 'coucou';
				#	$variation->{$start}->{$alternate_allele}->{$sAi[1]} = undef;
				#delete $variation->{$start}->{$alternate_allele};
				
				}
				else {
					$variation->{$start}->{$alternate_allele}->{$sAi[1]} = pack('w4 c4',ceil($sAi[2]*100),ceil($sAi[3]*100),ceil($sAi[4]*100),ceil($sAi[5]*100),$sAi[6],$sAi[7],$sAi[8],$sAi[9]);#[ceil($sAi[2]*100),ceil($sAi[3]*100),ceil($sAi[4]*100),ceil($sAi[5]*100),$sAi[6],$sAi[7],$sAi[8],$sAi[9]];
				}
				my $max;
				my @values = (0.2,0.5,0.8);
				foreach my $value (@values){
					if ($max > $value){
						$hash_bivector->{$value}->{$alternate_allele} = Bit::Vector->new($length) unless exists $hash_bivector->{$alternate_allele};
						$hash_bivector->{$value}->{$alternate_allele}->Bit_On($start);
					}
				}

	  	 }
}

my @base = ["A","T","C","G"];
return ($variation,$intspan,$hash_bitvector);
}
