#!/usr/bin/perl
use FindBin qw($Bin);
use lib $Bin;
use lib "$Bin/../../../../../lib/obj-nodb/";
use lib "$Bin/../packages/";
#use lib "/bip-d/perl/ensembl64/ensembl/modules";
use Bio::SearchIO;
use strict;
use Getopt::Long;
use Carp;
use JSON::XS;
use Getopt::Long; 
use List::MoreUtils qw(uniq);
#use ensembl_buffer;
use Storable qw/freeze thaw nfreeze nstore_fd nstore/;
use Set::IntSpan::Fast::XS ;
use Sys::Hostname;
use Storable qw/retrieve store/;
use  Set::IntSpan::Island;
use Parallel::ForkManager;
use Array::IntSpan;
use Data::Dumper;
use Array::IntSpan;
use GenBoNoSql;
use GenBoNoSqlAnnotation;
#use Bio::Ensembl::Registry;
use Storable qw(dclone);
use List::MoreUtils qw(uniq);
use Set::IntervalTree;
use Date::Tiny;
use Digest::MD5::File qw( file_md5_hex );
use GBuffer;
 use List::Util qw( max );
require "$Bin/../packages/parse_gff.pm";


require "$Bin/../packages/ensembl_buffer.pm";
#my $sqliteDir =  "/data-xfs/public-data/HG19/sqlite/75/annotation_test";

my $gene_code_version = "28";
#my $gff = "/data-xfs/public-data/HG19/gencode/v28/gencode.v28lift37.annotation.gff3.gz";

my $version;
GetOptions(
	'version=s' => \$version,
);
die("add version ") unless $version; 
my $sqliteDir =  "/tmp/lmdb/$gene_code_version/annotations";
system("mkdir -p $sqliteDir");
my $buffer = GBuffer->new();
my $project = $buffer->newProject(-name=>"NGS2019_2382");
$project->version("HG19.28");

my $dir_in = "/data-xfs/public-data/repository/HG19/gnomad-constraint/".$version."/";

my $file = $dir_in."txt/constraint.txt.bgz";

open (IN , "zcat $file | ");
my (@h) = split("n ",<IN>);
die() unless $h[21] ne "pLI";
my %hgenes;
my %htranscripts;
while (<IN>){
	chomp();
	my @t = split(" ",$_);
	my $tr = $t[1];
	my $pli = int($t[21]*100)/100;
	my $id =  $project->getGenBoId($tr);
	$htranscripts{$tr} = $pli;
	next unless  $id;
	my $transcript = $project->newTranscript($id);
	my $gene = $transcript->getGene();
	push(@{$hgenes{$gene->name}},$pli);
}

my $no  = GenBoNoSqlLmdb->new(name=>"pLI",dir=>$sqliteDir,mode=>"c",is_compress=>0);
foreach my $k (%hgenes){
	
	my $v = max(@{$hgenes{$k}});
	$v ="0.0" if $v ==0;
	warn $v if $k eq "ENSG00000177030";
	$no->put($k,$v);
	
}
foreach my $k (%htranscripts){
	
	my $v = $htranscripts{$k};
	
	$no->put($k,($v*1.0));
		$v ="0.0" if $v ==0;
}

warn $no->get("ENSG00000177030");
$no->close();