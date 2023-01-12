#!/usr/bin/perl
use FindBin qw($Bin);
use lib $Bin;
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
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
require "$Bin/../packages/parse_gff.pm";


require "$Bin/../packages/ensembl_buffer.pm";
#my $sqliteDir =  "/data-xfs/public-data/HG19/sqlite/75/annotation_test";
my $sqliteDir =  "/tmp/lmdb/annotation";
#my $gff = "/data-xfs/public-data/HG19/gencode/v28/gencode.v28lift37.annotation.gff3.gz";

my $dir;
my $version;
GetOptions(
	'version=s' => \$version,
);
my $dir =  "/tmp/lmdb/$version/annotations";

#my $gff = "$Bin/imagine.gff.gz";
my $gff = "$Bin/miss.gff3.gz";
my $translate_id ={};
my $gene_code_version = "IMG1";
 my $no2 = GenBoNoSqlAnnotation->new(dir=>$dir,mode=>"w");

my $genes = parse_gff::read_gff_genes($gff,$translate_id,$gene_code_version);
foreach my $id (keys %$genes){
	$genes->{$id}->{processed_transcript} =1;
}
 my ($transcripts,$proteins) = parse_gff::read_gff_transcripts($gff,$genes,$translate_id,$gene_code_version);
 foreach my $id (keys %$transcripts){
 	$transcripts->{$id}->{processed_transcript} =1;
 }

 parse_gff::save_sqlite ($no2,$genes,"gene",2);
 parse_gff::save_sqlite ($no2,$transcripts,"transcript",2);
 parse_gff::save_sqlite ($no2,$proteins,"protein",2);
 
 warn Dumper $no2->get("annotations","LCR");
 warn "end";
 
$no2->close();


exit(0);