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
my $buffer = GBuffer->new();
my $query_panel = $buffer->queryPanel();
my $query_phenotype = $buffer->queryPhenotype();

my $phenotypes = $query_phenotype->getAllPhenotypes();
my $genes = $query_panel->getAllGenes();

my $project = $buffer->newProject( -name => "NGS2020_2921");

my $version;
GetOptions(
	'version=s' => \$version,
);
die("add version ") unless $version; 
my $sqliteDir =  "/tmp/lmdb/omim/$version/lmdb/";

my $dir = "/public-data/repository/HG19/polyScore/$version/lmdb";
system ("mkdir -p $dir") unless -e $dir;

my $no = GenBoNoSqlLmdb->new(name=>"polyScore",dir=>$dir,mode=>"c",is_compress=>1);

foreach my $gid (@$genes){
	my $gene = $project->newGene($gid);
	my $score_gene;
	foreach my $ph_id (@$phenotypes){
		
		my $phenotype = $project->newPhenotype($ph_id);
		#warn $phenotype->concept;
		my $res = score_phenotype($gene,$phenotype);
		next if $res->{score} == 0;
		$score_gene->{$phenotype->id}->{score} = $res->{score}; 
		#$score_gene->{$phenotype->id}->{hgmd} = $res->{hgmd}; 
		$score_gene->{$phenotype->id}->{name} = $phenotype->name; 
		#warn $gene->name." ".$score;
	} 
$no->put($gid,$score_gene);
}
$no->close();
sub raw_score {
	my ($self) = @_;
	my $debug;
	$debug = 1 if  $self->external_name eq "MED12";
		my $score = 0;
		$score += 0.5 if $self->pLI ne "-" && $self->pLI > 0.95;
		my $gpheno = $self->phenotypes();
		my $pscore =0;
		$score += 1 if $self->is_omim_morbid();
		warn $score if $debug;
		my $apheno = $self->getPhenotypes();

	#	warn $score if $debug;
		foreach my $pheno (@{$self->getProject->getPhenotypes}){
			foreach my $k  (@{$pheno->keywords()}) {
				my $lk = lc($k);
				if ($gpheno =~ /$lk/){
					$score += 1 ;
					last;
				}
			}
			warn $score if $debug;
			if (exists $pheno->statistic_genes->{$self->id}){
				$score += 1;
				my $p = ($pheno->statistic_genes->{$self->id} / $pheno->nb_panels) *100;
				$score += 0.5 if $p >= 25;
				$score += 0.5 if $p >= 50;
				$score += 0.5 if $p >= 75;
				$score += 0.5 if $p >= 90;
				$score += 0.5 if $p >= 100;
			}
			else {
				my @o = grep{$pheno->id ne $_->id;} @$apheno; 
				$score += 0.2 * (scalar(@o));
			}
			warn $score if $debug;
		} 
		warn $score if $debug;
		#die() if $debug;
		
		
		return $score;
}



sub score_phenotype {
	my ($gene,$phenotype) = @_;
		my $res;
		my $score = 0;
		my $debug;
		
		#$score += 0.5 if $gene->pLI ne "-" && $gene->pLI > 0.95;
		my $gpheno = $gene->phenotypes();
		
		my $pscore =0;
		#$score += 1 if $gene->is_omim_morbid();
		
		#my $apheno = $gene->getPhenotypes();
		#$score += 0.2 * (@$apheno - 1 );
		my $query_hgmd = $gene->buffer->queryHgmd;
		my $gname = $gene->external_name; 
		my $stat = $phenotype->statistic_genes->{$gene->id};
		$stat =0 unless $stat;
		my $nb_panels = $phenotype->nb_panels;
		$debug =1 if $gname eq "MED12";
		my $hgmd_res = $query_hgmd->get_score_concept($gene->external_name,$phenotype->concept);
		
		if (exists $hgmd_res->{$gname}){
			$nb_panels ++;
			
			$stat ++ if $hgmd_res->{$gname}->{num_matching} > 0; 
			#warn Dumper $res->{$gname};
			#$res->{hgmd} = $hgmd_res->{$gname};
			$score ++ if $hgmd_res->{$gname}->{num_matching} > 10 && $hgmd_res->{$gname}->{percentage} > 50;
		}
	#	warn $score if $debug;
			foreach my $k  (@{$phenotype->keywords()}) {
				my $lk = lc($k);
				if ($gpheno =~ /$lk/){
					$score += 1 ;
					last;
				}
			}
			if ( $stat > 0){
				$score += 1;
				my $p = ($stat / $nb_panels) *100;
				if ($nb_panels > 2){
					$score += 0.5 if $p >= 25;
					$score += 0.5 if $p >= 50;
					$score += 0.5 if $p >= 75;
					$score += 0.5 if $p >= 90;
					$score += 0.5 if $p >= 100;
				}
				else {
					$score += 0.5;
				}
				
			}
			
		
		$res->{score} = $score;
		#$res->{score} = $score;
		
		return $res;
}