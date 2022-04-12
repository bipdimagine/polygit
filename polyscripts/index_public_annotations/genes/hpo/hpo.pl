#!/usr/bin/perl
use FindBin qw($Bin);
use lib $Bin;
use lib "$Bin/../../../../../lib/obj-nodb/";
use lib "$Bin/../packages/";
use Bio::SearchIO;
use strict;
use Getopt::Long;
use Carp;
use JSON::XS;
use Getopt::Long; 
use List::MoreUtils qw(uniq);
use Storable qw/freeze thaw nfreeze nstore_fd nstore/;
use Set::IntSpan::Fast::XS ;
use Sys::Hostname;
use Storable qw/retrieve store/;
use Set::IntSpan::Island;
use Parallel::ForkManager;
use Array::IntSpan;
use Data::Dumper;
use Array::IntSpan;
use Storable qw(dclone);
use List::MoreUtils qw(uniq);
use Set::IntervalTree;
use Date::Tiny;
use Digest::MD5::File qw( file_md5_hex );
use GBuffer;



my $version;
my $dir_out =  "/tmp/lmdb/hpo.2020.08.11";
my $dir_root = "/public-data/repository/HG19/hpo/20200811/";
#my $dir_root ="/data-xfs/public-data/repository/OMIM/20.02.2019/";
my $dir_last_review = "/public-data/repository/HG19/hpo/20200811/lmdb/";

my $file_hpo = $dir_root.'/files/hp.obo';
my $file_genes_to_phenotypes = $dir_root.'/files/genes_to_phenotype.txt';

unless (-e $file_hpo) {
	warn "\n\nFILE $file_hpo doesn't exist. Die.\n\n";
	die();
}
unless (-e $file_genes_to_phenotypes) {
	warn "\n\nFILE $file_genes_to_phenotypes doesn't exist. Die.\n\n";
	die();
}

my $h_parse_values = {
          'remark' => 'value',
          'data-version' => 'value',
          'ontology' => 'value',
          'format-version' => 'value',
          'saved-by' => 'value',
          'logical-definition-view-relation' => 'value',
          'default-namespace' => 'value',
          'name' => 'value',
          'id' => 'value',
          'def' => 'value',
          'created_by' => 'value',
          'creation_date' => 'value',
        };

my ($h_hpo, $h_genes_to_hpo);
my $term_name = 'log';
my $term_id = $term_name;
my $h_term_infos;

#Parsing de tous les TERM Hpo
open (FILE, $file_hpo);
while (<FILE>) {
	my $line = $_;
	chomp($line);
	next unless ($line);
	if ($line eq '[Term]') {
		if ($term_name ne 'log' and $term_name ne 'All') {
			#warn Dumper $h_term_infos; die;
		}
#		$h_hpo->{$h_term_infos->{name}} = $h_term_infos;
		$h_hpo->{$h_term_infos->{id}} = $h_term_infos;
		$term_name = undef;
		$term_id = undef;
		$h_term_infos = undef;
	}
	else {
		my ($cat_id, $cat_value) = split(': ', $line);
		if ($cat_id eq 'name') { $term_name = $cat_value; }
		if ($cat_id eq 'id') { $term_id = $cat_value; }
		if (exists $h_parse_values->{$cat_id}) {
			if (exists $h_term_infos->{$cat_id}) {
				warn "\n\n";
				warn "Pb for term id:$term_id / name:$term_name\n";
				warn "$cat_id already exists...";
				warn "\n\n";
				warn Dumper $h_term_infos;
				warn 'This line: '.$line;
				warn "\n\n";
				die;
			}
			$h_term_infos->{$cat_id} = $cat_value;
		}
		else {
			push(@{$h_term_infos->{$cat_id}}, $cat_value);
		}
	}
}
close (FILE);


#Parsing de toutes les relations Genes / Phenotypes HPO
my $h_header_genes_to_pheno;
open (FILE, $file_genes_to_phenotypes);
while (<FILE>) {
	my $line = $_;
	chomp($line);
	next if ($line =~ /#Format/);
	my @lCol = split("\t", $line);
	my $gene_name = $lCol[1]; 
	my $hpo_id = $lCol[2]; # col du type HP:0000815
	my $hpo_name = $lCol[3];
	push(@{$h_hpo->{$hpo_id}->{genes}}, $gene_name);
	$h_genes_to_hpo->{$gene_name}->{$hpo_id} = $hpo_name;
}
close (FILE);



warn "# LMDB hpo_phenotypes put\n";
my $annot =  GenBoNoSqlLmdb->new(name=>"hpo_phenotypes",dir=>$dir_out,mode=>"c",is_compress=>1);
foreach my $hpo_id (keys %$h_hpo){
	next unless $hpo_id;
	$annot->put($hpo_id, $h_hpo->{$hpo_id});
}
$annot->close();

warn "# LMDB hpo_genes_to_phenotypes put\n";
my $annot2 =  GenBoNoSqlLmdb->new(name=>"hpo_genes_to_phenotypes",dir=>$dir_out,mode=>"c",is_compress=>1);
foreach my $gene_name (keys %$h_genes_to_hpo){
	next unless $gene_name;
	$annot2->put($gene_name, $h_genes_to_hpo->{$gene_name});
}
$annot2->close();

warn '# LMDB done: '.$dir_out."\n";



