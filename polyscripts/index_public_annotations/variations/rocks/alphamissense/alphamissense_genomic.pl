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
use GBuffer;
use GenBoProject;
use Parallel::ForkManager;
use GenBoNoSql;
 use Time::ETA;
use List::MoreUtils qw(any uniq natatime);
use Bio::DB::HTS::Tabix;
use POSIX;
my $allsnps;
use Getopt::Long;
use Date::Tiny;
use Bio::DB::HTS::Faidx;
use GenBoNoSqlRocksGenome;

my $chr_name;
my $version;
my $fork;
my $genome_version;
my $merge;
GetOptions(
	'chr=s' => \$chr_name,
	'version=s' => \$version,
	'fork=s' => \$fork,
	'genome=s' => \$genome_version,
	'merge=s' => \$merge,
);
die("genome") unless $genome_version; 
my $dir_public= "/data-isilon/public-data/repository/$genome_version/alphamissense/$version/";
my $file = $dir_public."tabix/AlphaMissense_isoforms_aa_substitutions.tsv.gz";
my $file2 = $dir_public."tabix//AlphaMissense_".lc($genome_version).".tsv.gz";
my $dir_out= $dir_public."/rocks/";
my $AA  = [
    'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
];
my $hAA;
for (my $i=0;$i<=@$AA;$i++){
	$hAA->{$AA->[$i]} = $i;
}
warn Dumper $hAA;
my $no = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out,mode=>"c",name=>"alphamissense_genomic",pipeline=>1);	

open(TSV,"zcat $file2"."|");
my $array = [];
my $nb =0;
my %tt;
my $latest_uid;
my $value= [];
my $uid ="";
my $svalue="";
while (<TSV>) {
	chomp();
	my $line= $_;
	next if $line =~/^#/;
	my @z = split(" ",$line);
	my $chr = $z[0];
	$chr =~ s/chr//;
	$chr = "MT" if ($chr eq "M");
	
	my $pos = $z[1];
	
	my $alternate_allele = $z[3];
				
	my $ref_allele = $z[2];
	my ($enst,$v) =  split(/\./,$z[6]);
	
	$uid = $chr."!".$no->return_rocks_id($pos,$ref_allele,$alternate_allele);
	$z[7] =~ /([A-Z])(\d)([A-Z])/g;
	my $a = $1;
	my $b = $3;
	my $p = $2;
	my $pos_aa_array = $hAA->{$b};
	
	my $score = $z[8];
	my $cons = $z[9];
	
	my $d = $no->get_id_dictionary($enst);
#	if ($latest_uid && $uid ne $latest_uid ){
#		$no->put_batch_raw($uid, join(";",@$value));
#		$value = [];
#	}
	warn($line) if exists $tt{$uid};
	$tt{$uid} ++;
	my $d = $no->get_id_dictionary($enst);
	$uid = $uid."!".$d;
	$no->put_batch_raw($uid,$d.":".$score);
}
$no->put_batch_raw($uid, join(";",@$value));
$no->write_batch;
$no->close();




