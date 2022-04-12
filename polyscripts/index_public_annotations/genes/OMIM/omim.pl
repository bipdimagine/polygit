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
require "$Bin/../packages/parse_gff.pm";
use GBuffer;


require "$Bin/../packages/ensembl_buffer.pm";
#my $sqliteDir =  "/data-xfs/public-data/HG19/sqlite/75/annotation_test";
my $dir_out =  "/tmp/lmdb/omim.23.06.2020";
#my $gff = "/data-xfs/public-data/HG19/gencode/v28/gencode.v28lift37.annotation.gff3.gz";

my $version;
my $dir_root = "/public-data/repository/HG19/omim/23.06.2020/";
#my $dir_root ="/data-xfs/public-data/repository/OMIM/20.02.2019/";

my $dir_last_review = "/public-data/repository/HG19/omim/20.02.2019/lmdb/";

my $genes;
my $syno;

my $buffer = new GBuffer;

my $genemap_file = $dir_root."/txt/genemap2.txt";
warn "# PARSING $genemap_file\n";
my $h_genemap;
my @lHeaders;
open(FILE,"$genemap_file");
while(<FILE>){
	chomp();
	if ($_ =~ /# Chromosome/) {
		$_ =~ s/# //;
		my @t = split("\t", $_);
		foreach my $cat (@t) { push(@lHeaders, lc($cat)); }
	}
	elsif ($_ =~ /#/) { next; }
	$_ =~ s/, /,/g;
	my @t = split("\t",$_);
	$t[0] =~ s/chr//;
#	my $chr_name = $t[0];
#	my $gene_name = $t[8];
	my $gene_symbols = $t[6];
	my @lGeneNames = split(',', $gene_symbols);
	my $i = 0;
	foreach my $value (@t) {
		my $cat_name = $lHeaders[$i];
		$value =~ s/,/, /g;
#		if ($cat_name eq 'ensembl gene id') { $value .= '_'.$chr_name; }
		foreach my $gene_name (@lGeneNames) {
			$h_genemap->{$gene_name}->{omim_details}->{$cat_name} = $value;
			$h_genemap->{$gene_name}->{is_morbid} = undef;
		}
		$i++;
	}
	#TODO: rajouter tous les synonimes genes name en cle de hash aussi
}
close(FILE);

my $morbid_file = $dir_root."/txt/morbidmap.txt";
warn "# PARSING $morbid_file\n";
my @lHeaders2;
open(FILE,"$morbid_file");
while(<FILE>){
	chomp();
	if ($_ =~ /# Phenotype/) {
		$_ =~ s/# //;
		my @t = split("\t", $_);
		foreach my $cat (@t) { push(@lHeaders2, lc($cat)); }
	}
	elsif ($_ =~ /#/) { next; }
	$_ =~ s/, /,/g;
	my @t = split("\t",$_);
	$t[0] =~ s/chr//;
	my $gene_symbols = $t[1];
	my @lGeneNames = split(',', $gene_symbols);
	my $i = 0;
	foreach my $value (@t) {
		my $cat_name = $lHeaders2[$i];
		$value =~ s/,/, /g;
		foreach my $gene_name (@lGeneNames) {
			confess('ERROR '.$gene_name) unless (exists $h_genemap->{$gene_name});
			$h_genemap->{$gene_name}->{is_morbid} = 1;
			$h_genemap->{$gene_name}->{omim_morbid_details}->{$cat_name} = $value;
		}
		$i++;
	}
}
close(FILE);

my $mart_file = $dir_root."/txt/mart_omim_genes.txt";
warn "# PARSING $mart_file\n";
my @lHeaders3;
open(FILE,"$mart_file");
my $j = 0;
while(<FILE>){
	$j++;
	chomp();
	if ($j == 1){
		my @t = split("\t", $_);
		foreach my $cat (@t) { push(@lHeaders3, lc($cat)); }
		next;
	}
	my @t = split("\t",$_);
	my $gene_name = $t[0];
	my $i = 0;
	foreach my $value (@t) {
		my $cat_name = $lHeaders3[$i];
		$h_genemap->{$gene_name}->{mart_omim_details}->{$cat_name} = $value;
		$i++;
	}
	unless (exists $h_genemap->{$gene_name}->{is_morbid}) {
		$h_genemap->{$gene_name}->{is_morbid} = undef;
	}
}
close(FILE);

my $mart_morbid_file = $dir_root."/txt/mart_omim_genes_morbid.txt";
warn "# PARSING $mart_morbid_file\n";
my @lHeaders4;
open(FILE,"$mart_morbid_file");
$j = 0;
while(<FILE>){
	$j++;
	chomp();
	if ($j == 1){
		my @t = split("\t", $_);
		foreach my $cat (@t) { push(@lHeaders4, lc($cat)); }
		next;
	}
	my @t = split("\t",$_);
	my $gene_name = $t[0];
	my $i = 0;
	foreach my $value (@t) {
		my $cat_name = $lHeaders4[$i];
		$h_genemap->{$gene_name}->{mart_omim_morbid_details}->{$cat_name} = $value;
		$i++;
	}
	$h_genemap->{$gene_name}->{is_morbid} = 1;
}
close(FILE);

my $annot_last_review =  GenBoNoSqlLmdb->new(name=>"omim",dir=>$dir_last_review,mode=>"r",is_compress=>1);

foreach my $gene_name (keys %{$h_genemap}) {
	my ($gene_id, $chr_name, $omim_id, @l_phenotypes);
	if (exists $h_genemap->{$gene_name}->{omim_details}) {
		$chr_name = $h_genemap->{$gene_name}->{omim_details}->{'chromosome'};
		$omim_id = $h_genemap->{$gene_name}->{omim_details}->{'mim number'} if (not $omim_id and exists $h_genemap->{$gene_name}->{omim_details}->{'mim number'});
		$gene_id = $h_genemap->{$gene_name}->{omim_details}->{'ensembl gene id'};
		$h_genemap->{$gene_name}->{phenotypes}->{omim_details} = $h_genemap->{$gene_name}->{omim_details}->{'phenotypes'};
	}
	if (exists $h_genemap->{$gene_name}->{mart_omim_details}) {
		$omim_id = $h_genemap->{$gene_name}->{mart_omim_details}->{'mim number'} if (not $omim_id and exists $h_genemap->{$gene_name}->{mart_omim_details}->{'mim number'});
		$gene_id = $h_genemap->{$gene_name}->{mart_omim_details}->{'gene stable id'};
		$h_genemap->{$gene_name}->{phenotypes}->{mart_omim_details} = $h_genemap->{$gene_name}->{mart_omim_details}->{'phenotype description'};
	}
	if (exists $h_genemap->{$gene_name}->{omim_morbid_details}) {
		push(@l_phenotypes, $h_genemap->{$gene_name}->{omim_morbid_details}->{'phenotype'});
		$omim_id = $h_genemap->{$gene_name}->{omim_morbid_details}->{'mim gene accession'} if (not $omim_id and exists $h_genemap->{$gene_name}->{omim_morbid_details}->{'mim gene accession'});
	}
	if (exists $h_genemap->{$gene_name}->{mart_omim_morbid_details}) {
		$chr_name = $h_genemap->{$gene_name}->{mart_omim_morbid_details}->{'chromosome/scaffold name'};
		$omim_id = $h_genemap->{$gene_name}->{mart_omim_morbid_details}->{'mim gene accession'} if (not $omim_id and exists $h_genemap->{$gene_name}->{mart_omim_details}->{'mim number'});
		$gene_id = $h_genemap->{$gene_name}->{mart_omim_morbid_details}->{'gene stable id'};
		$h_genemap->{$gene_name}->{phenotypes}->{mart_omim_morbid_details} = $h_genemap->{$gene_name}->{omim_details}->{'phenotype description'};
	}
	$gene_id .= '_'.$chr_name;
	$h_genemap->{$gene_name}->{chr_name} = $chr_name;
	$h_genemap->{$gene_name}->{omim_id} = $omim_id;
	$h_genemap->{$gene_name}->{gene_id} = $gene_id;
	if (exists $h_genemap->{$gene_name}->{phenotypes}->{omim_details}) {
		$h_genemap->{$gene_name}->{phenotype}->{omim} = $h_genemap->{$gene_name}->{phenotypes}->{omim_details};
	}
	elsif (exists $h_genemap->{$gene_name}->{phenotypes}->{mart_omim_details}) {
		$h_genemap->{$gene_name}->{phenotype}->{omim} = $h_genemap->{$gene_name}->{phenotypes}->{mart_omim_details};
	}
	elsif (exists $h_genemap->{$gene_name}->{phenotypes}->{mart_omim_morbid_details}) {
		$h_genemap->{$gene_name}->{phenotype}->{omim} = $h_genemap->{$gene_name}->{phenotypes}->{mart_omim_morbid_details};
	}
	if ($h_genemap->{$gene_name}->{phenotype}->{omim} =~/recessive/i ){
		$h_genemap->{$gene_name}->{inheritance}->{omim}->{recessive} = 1;
	}
	if ($h_genemap->{$gene_name}->{phenotype}->{omim} =~/dominant/i ){
		$h_genemap->{$gene_name}->{inheritance}->{omim}->{dominant} = 1;
	}
	if ($chr_name eq "X" ){
		$h_genemap->{$gene_name}->{inheritance}->{omim}->{"X-linked"} = 1;
	}
	if ($chr_name eq "Y" ){
		$h_genemap->{$gene_name}->{inheritance}->{omim}->{"Y-linked"} = 1;
	}
	my $h_last = $annot_last_review->get($gene_id);
	if ($h_last) {
		$h_genemap->{$gene_name}->{'new'} = undef;
		if ($h_genemap->{$gene_name}->{is_morbid} == 1) {
			# TODO: VALABLE a partir de la version superieur a 30.06.2020
			if ($h_last->{$gene_id}->{is_morbid} == 1) { $h_genemap->{$gene_name}->{'new_morbid'} = undef; }
			elsif (exists $buffer->hash_genes_omim_morbid->{$gene_name}) { $h_genemap->{$gene_name}->{'new_morbid'} = undef; }
			else { $h_genemap->{$gene_name}->{'new_morbid'} = 1; }
		}
		else {
			$h_genemap->{$gene_name}->{'new_morbid'} = undef;
		}
	}
	else {
		$h_genemap->{$gene_name}->{'new'} = 1;
		if ($h_genemap->{$gene_name}->{is_morbid} == 1) {
			$h_genemap->{$gene_name}->{'new_morbid'} = 1;
		}
		else {
			$h_genemap->{$gene_name}->{'new_morbid'} = undef;
		}
	}
	$h_genemap->{$gene_id} = $h_genemap->{$gene_name};
}
$annot_last_review->close();


#warn Dumper $h_genemap->{ENSG00000186862_10}; die;

#my $mart_file = $dir_root."/txt/mart.ensg.mim.pheno.name.txt";
#if (-e $mart_file) {
#	open(MART,"$mart_file");
#	while(<MART>){
#		chomp();
#		my @t =split("\t",$_);
#		my $ens = $t[0];
#		my $oid = $t[1];
#		my $pheno = $t[2];
#		my $gene_name = $t[3];
#		my $chr = $t[4];
#		next if $oid eq "." && $pheno eq ".";
#		my $genboid = $ens."_".$chr;
#		$genes->{$genboid}->{id} = $ens;
#		$genes->{$genboid}->{chromosome} = $chr;
#		$genes->{$genboid}->{omim_id} = $oid;
#		if (exists $genes->{$genboid}->{phenotype}->{ensembl}){
#						$genes->{$genboid}->{phenotype}->{ensembl} =$genes->{$genboid}->{phenotype}->{ensembl}.";".lc($pheno);
#						warn $genes->{$genboid}->{phenotype}->{ensembl} ;
#		}
#		else {
#			$genes->{$genboid}->{phenotype}->{ensembl} = lc($pheno);
#		}
#		$genes->{$genboid}->{name} = $gene_name;
#		$syno->{$oid} = $genboid;
#		$syno->{$gene_name} = $genboid;
#	}
#	close MART;
#}
#my $sqliteDir = "/data-xfs/public-data/repository/HG19/annotations/lmdb/gencode.28/";
#
#my $no = GenBoNoSqlAnnotation->new(dir=>$sqliteDir,mode=>"r");
##warn Dumper $no->get("synonyms","UGT1A10");
##warn Dumper $no->get("synonyms","UGT1A10");
##die();
#my $omim = $dir_root."/txt/genemap2.txt ";
#warn $omim;
#open (OMIM," cat ".$omim." | cut -f 1,6,9,11,13 | grep -v '#'  | ") or die();
#while(<OMIM>){
#	chomp();
#	my @t =split("\t",$_);
#	my $chr = $t[0];
#	$chr =~ s/chr//;
#	$chr="MT" if $chr eq "M";
#	my $oid = $t[1];
#	my $name =$t[2];
#	my $ens =$t[3];
#	my $pheno = $t[4];
#	next if $pheno eq ".";
#	my $gid = $ens."_".$chr;
#	unless (exists $genes->{$gid}){
#		$gid = $syno->{$oid};
#		unless (exists $genes->{$gid}){
#			$gid = $syno->{$name};
#		}
#		unless (exists $genes->{$gid}){
#			$gid="toto";
#			eval {
#				#warn $name;
#			my $new  = $no->get("synonyms",$name) if $name;
#			if ($new){
#				my ($ens,$t) = split("_",$new);
#				$genes->{$new}->{id} = $ens;
#				$genes->{$new}->{chromosome} = $chr;
#				$genes->{$new}->{omim_id} = $oid;
#				
#					$genes->{$new}->{phenotype}->{omim} = $pheno;
#				$genes->{$new}->{name} = $name;
#			}
#			warn $new;
#			};
#		}
#		unless (exists $genes->{$gid}){
#			#warn $oid." ==> ".$name." ".$ens." $oid" ;
#			next;
#		}
#	}
#	 $genes->{$gid}->{phenotype}->{omim} = $pheno;
#}
#
#my $annot2 =  GenBoNoSqlLmdb->new(name=>"omim",dir=>$dir_out,mode=>"c",is_compress=>1);
#foreach my $gid (keys %$genes){
#	next unless exists $genes->{$gid}->{phenotype}->{omim};
#	my $pheno = $genes->{$gid}->{phenotype}->{omim};
#	my $chr = $genes->{$gid}->{chromosome};
#		if ($pheno =~/recessive/i ){
#		$genes->{$gid}->{inheritance}->{omim}->{recessive} =1;
#	}
#	if ($pheno =~/dominant/i ){
#		$genes->{$gid}->{inheritance}->{omim}->{dominant} =1;
#	}
#	if ($chr eq "X" ){
#		$genes->{$gid}->{inheritance}->{omim}->{"X-linked"} =1;
#	}
#	if ($chr eq "Y" ){
#		$genes->{$gid}->{inheritance}->{omim}->{"Y-linked"} =1;
#	}
#}

warn "# LMDB put\n";
my $annot2 =  GenBoNoSqlLmdb->new(name=>"omim",dir=>$dir_out,mode=>"c",is_compress=>1);
foreach my $gid (keys %$h_genemap){
	next unless $gid;
	$annot2->put($gid, $h_genemap->{$gid});
}
$annot2->close();
warn '# LMDB done: '.$dir_out."\n";