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

	'merge=s' => \$merge,
);

my $dir_public= "/data-isilon/public-data/repository/semantic/gnomad-constraint/$version/tsv";

my $file = $dir_public."/gnomad.v4.0.constraint_metrics.tsv.gz";

my $rocksdir =  "/data-isilon/public-data/repository/semantic/gnomad-constraint/$version/rocksdb/";

system("mkdir -p $rocksdir");


open (IN , "zcat $file |  cut -f 1,2,3,17,21 | ");
my (@h) = split("n ",<IN>);
my $factor = ["0.01","0.01"];
my $finalrg = GenBoNoSqlRocksAnnotation->new(dir=>$rocksdir,mode=>"c",name=>"pLI",pack=>"w2",factor=>$factor,version=>$version,description=>["pli","loeuf"]);
die() unless $h[4] ne "lof.oe_ci.upper";
die() unless $h[3] ne "lof.pLI";
my %hgenes;
my %htranscripts;
while (<IN>){
chomp();
	my @t = split(" ",$_);
	my $gene = $t[0];
	my $tr = $t[1];
	next unless $tr =~ /ENST/;
	my $mane = $t[2];
	my $pli = int($t[3]*100);
	my $loeuf = int($t[4]*100);
	
	$finalrg->put_raw($tr,pack("w2",$pli,$loeuf));
	if ($mane eq "true"){
		if ($gene eq "AURKA"){
				warn $gene;
	warn ($pli/100);
	warn " + ".($loeuf/100)." + ";
		warn $t[3];
		}
	
		$finalrg->put_raw($gene,pack("w2",$pli,$loeuf));
	}
	
}

$finalrg->close();

my $finalrg = GenBoNoSqlRocksAnnotation->new(dir=>$rocksdir,mode=>"r",name=>"pLI");
warn $finalrg->pli("AURKA");
die();

