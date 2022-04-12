#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../../../../../lib/obj-nodb";
use strict; 
use Bio::SearchIO;
use strict;

use Data::Printer;
use Getopt::Long;
use Carp;
use JSON::XS;
use Data::Dumper;
use Getopt::Long; 
use Storable qw/retrieve freeze thaw nfreeze nstore_fd nstore/;
use GBuffer;
use Sys::Hostname;
use Parallel::ForkManager;
use Vcf;
use Digest::MD5::File qw( file_md5_hex );
use Date::Tiny;
use Data::Compare;
#use GenBoNoSql;
#require "$Bin/dbsnp.pm";


my $buffer = GBuffer->new();

my $database = "clinvar";

my $dir_prod = "";
my $dir_temp = "";

my @chrs = (1..22);
my @types = ("snps","deletions","insertions");
foreach my $type (@types){
	my $lmdb_prod = $buffer->lmdb_public_dir."$database/".$type; 
		my $dir_lmdb_tmp ="/tmp/lmdb/$database/".$type; 
	foreach my $chr (@chrs){
		my $lmdb_prod =  GenBoNoSqlLmdb->new(dir=>$lmdb_prod,mode=>"r",name=>$chr,is_compress=>1,is_integer=>1);
		warn Dumper $lmdb_prod->get_with_sequence("949422","A");
		#my $pos =  $lmdb_prod->get_keys();
		
		my $lmdb_test =  GenBoNoSqlLmdb->new(dir=>$dir_lmdb_tmp,mode=>"r",name=>$chr,is_compress=>1,is_integer=>1);
		#my $keys = $lmdb_test->get_keys();
			while (my $k = $lmdb_prod->next_key){
				my $a = $lmdb_prod->get($k);
				my $b = $lmdb_test->get($k);
				if (Compare($a, $b)){
				} 
				else {
					warn $k;
					warn Dumper $a;
					warn Dumper $b;
				}
			}
			

	}
}