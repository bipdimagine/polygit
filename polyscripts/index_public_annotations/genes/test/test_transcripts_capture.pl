#!/usr/bin/perl
use FindBin qw($Bin);
use lib $Bin;
use lib "$Bin/../../../../../lib/obj-nodb/";
#use lib "/bip-d/perl/ensembl64/ensembl/modules";
use Bio::SearchIO;
use GBuffer;
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

my $sqliteDir =  "/tmp/lmdb/annotation";
my $annot =  GenBoNoSqlAnnotation->new(dir=>$sqliteDir,mode=>"r");
my $sqliteDir =  "/data-xfs/public-data/HG19/sqlite/75.1/annotation/";
my $annot2 =  GenBoNoSqlAnnotation->new(dir=>$sqliteDir,mode=>"r");
my $buffer = new GBuffer;
my $dbh = $buffer->dbh();

my $array = $dbh->selectall_arrayref("select ensembl_id,gene,description from PolyprojectNGS.transcripts as t,PolyprojectNGS.bundle_transcripts as bt , PolyprojectNGS.bundle b where  bt.transcript_id=t.id and b.bundle_id=bt.bundle_id;");
my $hh;
my $nb;
foreach my $a (@$array){
	#next unless $a->[1];
	next if $a->[2] eq "agilent_50mb_v5";
		next if $a->[2] =~ /omim/i;#"Omim_hg19";
		
		#	next if $a->[1] eq "";
	#next if $a->[0] =~/loop/;

				
	#eval{
	my $z = $annot->get("annotations",$a->[0]);
	$nb ++ if $z;
	unless ($z){
		$hh->{$a->[0]}->{gene} = $a->[1] ;
		$hh->{$a->[0]}->{capture} = $a->[2] ;
		#warn $a->[0]." ".$a->[1]." ".$a->[2]." ".$z->{ccds_id}
	}
	
	#die(Dumper ($a)) unless  $z;
	
#	}
	#die() unless $z;
	
	
	
}
warn "end ";
my @ts =keys %$hh;
warn scalar(@ts);
foreach my $t (@ts) {
	my $z = $annot2->get("annotations",$t);
	if ($z->{ccds_name}){
	$z->{ccds_name} =~s/\..*//;
	#warn $t." ".$z->{ccds_name};
	my $x = $annot->get("annotations",$z->{ccds_name});
	unless ($x){
		#warn $t ." ".Dumper $hh->{$t};
	}
	#warn $x->{ccds_name}." ".$t." ".$x->{name};
	}
	else {
		warn "$t ".$z->{external_name}." ".$hh->{$t}->{capture} ;
	}
	
}

warn scalar( keys %$hh)." ".$nb;

