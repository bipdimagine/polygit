#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use strict; 
use GenBoNoSqlRocks;
use Data::Printer;
use Getopt::Long;
use Carp;
use Data::Dumper;
use Getopt::Long; 
use Parallel::ForkManager;
 use Time::ETA;
use List::MoreUtils qw(any uniq natatime);
use Bio::DB::HTS::Tabix;
use POSIX;
my $allsnps;
use Getopt::Long;
use Date::Tiny;
use Bio::DB::HTS::Faidx;
use GenBoNoSqlRocksGenome;
use GenBoNoSqlRocksAnnotation;
use Sereal qw(sereal_encode_with_object sereal_decode_with_object write_sereal_file read_sereal_file);
my $chr_name;
my $version;
my $type;
my $genome_version;
my $merge;
GetOptions(
	'chr=s' => \$chr_name,
	'version=s' => \$version,
	'type=s' => \$type,
	'genome=s' => \$genome_version,
	'merge=s' => \$merge,
);
$type = "spliceAI";
die("spliceai,cadd,gnomad") unless $type;
die("genome") unless $genome_version; 
my $dir_public= "/data-isilon/public-data/repository/$genome_version/$type/$version/";



#my %bit_codes = (
#        T => 0b00,
#        A => 0b11,
#        C => 0b10,
#        G => 0b01,
#        );
#
## add the reverse mapping too
#@bit_codes{values %bit_codes} = keys %bit_codes;
#
#use constant WIDTH => 2;
#
#my $bits = '';
#my @bases = split //, 'CCGGAGAGATTAC';
#
#foreach my $i ( 0 .. $#bases ) {
#        vec( $bits, $i, WIDTH ) = $bit_codes{ $bases[$i] };
#        }
#
#print "Length of string is " . length( $bits ) . "\n";
#warn $bits;
#my $base = vec $bits, 2, WIDTH;
#printf "The third element is %s\n", $bit_codes{ $base };
#
#die($base);




die($dir_public) unless -e $dir_public;
my $dir_in =  "$dir_public"."rocksdb_split/";
my $dir_out = "$dir_public"."rocksdb/";
my $dir_out2 = "$dir_public"."rocksdb2/";
if ($chr_name eq "MT"){
	my $finalrg = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out,mode=>"t",name=>$chr_name,pipeline=>1);
	$finalrg->put_batch_raw("date",time);
	$finalrg->write_batch();
	$finalrg->close();
	exit(0);
}


my $rg = GenBoNoSqlRocksGenome->new(dir=>$dir_in,mode=>"r",index=>"genomic",chromosome=>1,genome=>$genome_version);	
$rg->{pack} .= " A*";
unless ($rg->description){
$rg->{description} = ["AG","AL","DG","DL","PAG","PAL","PDG","PDL","gene"];
}
$rg->{description} = ["AG","AL","DG","DL","PAG","PAL","PDG","PDL","gene"];
$rg->{factor} = ["0.01","0.01","0.01","0.01","1","1","1","1",""];

#my $rg = GenBoNoSqlRocksGenome->new(dir=>$dir_out,mode=>"c",chromosome=>$chr_name,fasta=>$genome_fasta,compress=>1,index=>"vcf");
my $final_all = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out,mode=>"t",name=>"spliceai",pack=>$rg->pack,version=>$version,description=>$rg->description,pipeline=>1,factor=>$rg->factor);	
undef $rg;
my @chromosomes =(1..22,"X","Y");
@chromosomes = ($chr_name) if $chr_name;
my $cyto_id  = 0;
foreach my $chr (@chromosomes){
	$cyto_id ++;
	my $rg = GenBoNoSqlRocksGenome->new(dir=>$dir_in,mode=>"r",index=>"genomic",chromosome=>$chr,genome=>$genome_version);	
	 warn Dumper $rg->chunks;
unless ($rg->description){
$rg->{description} = ["AG","AL","DG","DL","PAG","PAL","PDG","PDL","gene"];
}
$rg->{description} = ["AG","AL","DG","DL","PAG","PAL","PDG","PDL","gene"];
$rg->{factor} = ["0.01","0.01","0.01","0.01","1","1","1","1",undef];
my $finalrg = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out,mode=>"t",name=>$chr,pack=>$rg->pack,version=>$version,description=>$rg->description,pipeline=>1,factor=>$rg->factor);	
#my $finalrg2 = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out2,mode=>"t",name=>$chr_name,pack=>$rg->pack,version=>$version,description=>$rg->description,pipeline=>1);	
my $genes;
my $nb_regions = 0;
my $old_id = 0;
my $array = [];
my $gid = 0;
foreach my $r (@{$rg->regions}) {
		my $no = $rg->_chunks_rocks($r->{id});
		$nb_regions ++;
		warn "chunk : ".$nb_regions."/".scalar(@{$rg->regions});
		my $iter = $no->rocks->new_iterator->seek_to_first;
		while (my ($key, $value) = $iter->each) {
    		$finalrg->put_batch_raw($key,$value);
    	
    		#$finalrg->put_batch_raw($key,$value);
		}
		$finalrg->write_batch();
		$finalrg->intspan_keys($finalrg->intspan_keys->union($no->intspan_keys));
		warn "before undef";
		#warn $finalrg->intspan_keys->as_string();
		$no = undef;
		warn "end ";
	
	}
	
$finalrg->close();
undef $rg;

}
$final_all->close();
  exit(0);



