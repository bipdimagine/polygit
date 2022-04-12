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
my $allsnps;
use Getopt::Long;
use Date::Tiny;
my $c;
my $version;
GetOptions(
	'chr=s' => \$c,
	'version=s' => \$version,
);
my $file = "/public-data/repository/HG19/cadd/$version/tabix/whole_genome_SNVs.tsv.gz";
die($file) unless -e $file;
#my $chromosomes = [1..22,'X','Y','MT'];
#my $chromosomes = [22];
#my $dir = "/data-xfs/public-data/HG19/snp/exac/latest/";
my $buffer = new GBuffer;
my $project_name= "NGS2016_1186";
my $dir_out = "/data-beegfs/marc/cadd_1.6/$version/";
system("mkdir -p $dir_out");

my $description;
$description->{name} = "cadd";
$description->{version} = "1.6";
$description->{file} = $file;
 my $d = Date::Tiny->now;
$description->{date} =  $d->as_string;
$description->{code_unpack} =  "C4";
$description->{order} =  {'A'=>0,'C'=>1,'G'=>2,'T'=>3};
$description->{pack_string} =  "C4";
$description->{factor} = [1,1,1,1];
$description->{score_description} = ["A","C","G","T"];

system("mkdir -p $dir_out");
open(JSON,">$dir_out/description.json");
#open(JSON,">$dir_out/description.json");
print JSON encode_json $description;
close JSON;
die("write $dir_out/description.json");
my $project = $buffer->newProject( -name => $project_name,);
my $nb;
my $size = 500000;

	my $chr = $project->getChromosome($c);
	my $no = GenBoNoSqlLmdb->new(dir=>$dir_out,mode=>"c",name=>$chr->name,is_compress=>0,is_integer=>1);
	$no->put("-1","coucou");
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
	$pm->run_on_finish (
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hres) = @_;
			$nbb ++;
			my $cc= $hres->{chr};
			
		#	
		#	my $lmdb_out3  =  $save_out->get_lmdb_db("snps3");
			my $data = $hres->{data};
			return unless $data;
			my $t = time;
			my $nb =0;
			my $no = GenBoNoSqlLmdb->new(dir=>$dir_out,mode=>"w",name=>$chr->name,is_compress=>0,is_integer=>1);
			foreach my $k (sort {$a <=> $b} keys %$data){
				$nb ++;
			#	$lmdb_out3->put($k,$data->{$k}->{h});
				delete $data->{$k}->{h};
				#$lmdb_out->put($k,$data->{$k}->{s});
				my $nb = scalar(@{$data->{$k}->{s}});
				#warn join(";",@{$data->{$k}->{s}}) ." ".$k;
				my $rec = pack("C$nb",@{$data->{$k}->{s}});
				$no->put($k,$rec);
			
				
			}
			 $eta->pass_milestone();
 			my $remaining = $eta->get_remaining_time();
			my $elapsed = $eta->get_elapsed_time();
			my $p = int($eta->get_completed_percent());
			print $chr->name." complete $p % => elapsed :".$elapsed." remaining : ".$remaining."\n";
			$no->close;
		});


#	$size = int($length/10);
	my $i;
	my $cname = $chr->name();
	
	while (my $r = shift(@$region)){
	my $pid = $pm->start and next;
	my $array = run_chr($cname,$r->{start},$r->{end});
	my $hres;
	$hres->{chr} = $cname;
	$hres->{data} = $array;
			
	$pm->finish(0,$hres);
	}
	
  $pm->wait_all_children;
  warn $i." ".$length;

#$lmdb_out->close();
exit(0);





sub run_chr {
	my ($chr,$start,$end) = @_;

die() unless -e $file;
my $ACGT = {'A'=>0,'C'=>1,'G'=>2,'T'=>3};
#my $db = new save({name=>"caad",chromosome=>$chr});
my $nb= 0;
	my $obj;
my $variation;
	#   warn $file;
	 #  $file =~ s/dibayes/unifiedgenotyper/;
	 my $tabix  = new Tabix(-data =>$file);
	  my $res = $tabix->query($chr,$start,$end);
	 
	  	 while(my $line = $tabix->read($res)) {
				my(@array) = split(" ",$line);
				my $start = $array[1];
				
				my $alternate_allele = $array[3];
				#$variation->{$start}->{$alternate_allele}->{c} = [$array[114],$array[115]];
				unless (exists $variation->{$start}->{s}) {
					$variation->{$start}->{s} = [0,0,0,0];
				}
				my $id = $ACGT->{$alternate_allele};
				#$variation->{$start."_".$alternate_allele} = $array[5];
				my $v = ceil($array[5]);
				#$variation->{$start}->{h}->{$alternate_allele} = $v;
				#push(@{$variation->{$start}->{a}},$alternate_allele);	
				#push(@{$variation->{$start}->{s}},$v);			
				$variation->{$start}->{s}->[$id] = $v;
				#$variation->{$start}->[$id] =$v;
				

	  	 }
#	  	 $db->close();
#	  	my $db2 = $db->get_snp_db_lite("caad2"); 
#		warn Dumper $db2->stat();
#		 	my $db3 = $db->get_snp_db_lite("caad"); 
#		warn Dumper $db3->stat();
#warn Dumper $obj;
my @base = ["A","T","C","G"];

return $variation;
}
