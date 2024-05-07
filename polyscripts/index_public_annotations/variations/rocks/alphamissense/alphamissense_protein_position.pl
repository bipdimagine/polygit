#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin/../../../../../GenBo/lib/obj-nodb/";
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
my $finalrg = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out,mode=>"c",name=>"alphamissense",pack=>["1000"],version=>$version,description=>$AA,pipeline=>1);
open(TSV,"zcat $file2"."|");
my $ltest_enst;
my $array = [];
my $harray= {};
my $nb =0;
my $dj;
while (<TSV>) {
	chomp();
	my $line= $_;
	next if $line =~/^#/;
	my @z = split(" ",$line);
	my $chr = $z[0];
	my $pos = $z[1];
	my ($enst,$v) =  split(/\./,$z[6]);
	$z[7] =~ /([A-Z])(\d+)([A-Z])/g;
	my $a = $1;
	my $b = $3;
	my $p = $2;
	my $pos_aa_array = $hAA->{$b};
	die($b." ".$line) unless exists $hAA->{$b};
	my $score = $z[8];
	my $cons = $z[9];
	#my $d = $finalrg->get_id_dictionary($cons);
	$harray->{$enst}->[$p]->[$pos_aa_array] = $score;
}

foreach my $enst (keys %$harray){
	warn $nb;
	$nb++;
	$finalrg->put_batch($enst,$harray->{$enst});
}



$finalrg->write_batch;
 open(TSV,"zcat $file"."|");

 $ltest_enst = undef;
$array =[];
my $enst;
while (<TSV>) {
	chomp();
	my $line= $_;
	next if $line =~/^#/;
	next if $line =~/^transcr/;
	chomp($line);
	my @z = split(" ",$line);
	my $v;
	 ($enst,$v) =  split(/\./,$z[0]);
	$z[1] =~ /([A-Z])(\d+)([A-Z])/g;
	my $a = $1;
	my $b = $3;
	my $p = $2;
	
	if ($ltest_enst ne $enst) {
		warn $nb;
		$nb++;
		if ($ltest_enst) {
			warn scalar (@$array) ;#if $ltest_enst eq "ENST00000459813";
			$finalrg->put_batch($ltest_enst,$array);
		}
		$ltest_enst = $enst;
		$array =[];
	}
	my $pos_aa_array = $hAA->{$b};
	die($b." ".$line) unless exists $hAA->{$b};
	#my $d = $finalrg->get_id_dictionary($z[3]);
	$array->[$p]->[$pos_aa_array] = $z[2];
	my $x = $finalrg->get_id_dictionary($enst);
}

$finalrg->put_batch($enst,$array);


$finalrg->write_batch;



$finalrg->close();



