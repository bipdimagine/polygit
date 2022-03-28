#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use Parallel::ForkManager;
use Storable qw(store retrieve freeze thaw);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Cwd 'abs_path';
use Digest::MD5 qw(md5_hex);
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
use lib "$RealBin/../../../polymorphism-cgi/cache_nodb/scripts/";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use Set::IntSpan::Fast::XS;
use Sys::Hostname;
use Cache_Commons;
 my $host = hostname();


#warn "*_*_*_*_*_".$host."*_*_*_*_*_";
#use Cache_Commons;

my $fork = 1;
my $force;
my ($project_name, $chr_name);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'chr=s'        => \$chr_name,
	'force=s'  => \$force,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
#unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }
my $t = `hostname`;

my $nbErrors = 0;
my $buffer = new GBuffer;

my $project = $buffer->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'familial' );
my $chromosomes = $project->getChromosomes;
my $categories = Cache_Commons::categories();

foreach my $chr (@{$chromosomes}){
	
	if ($chr_name){
		next if $chr->name ne "$chr_name";
	}
	run_update($chr->name);
}
exit(0) if ($chr_name);
#system("$RealBin/tree_cache.pl -project=$project_name -fork=$fork");


sub run_update {
	my ($chr_name) = @_;
	
	my $project = $buffer->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'familial' );
	my $chr = $project->getChromosome($chr_name);
	my $no = $chr->lmdb_variations("r");
	my $intspan_genes_categories = {};
	my $hGenesIds;
	
	my $ranges = $no->ranges($fork);
	#warn Dumper $ranges;
	$no->close();
	my $pm = new Parallel::ForkManager($fork);
	
	$pm->run_on_finish(
	    sub { 
	    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$hres)=@_;
	  
	    	unless (defined($hres) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
		
			#################
			#For genes annotation categorie intspan on gene
			#################
			foreach my $g ( keys %{ $hres->{genes} } ) {
				$hGenesIds->{$g}++;
				my $hgene = $hres->{genes}->{$g};
				unless ( exists $intspan_genes_categories->{$g} ) {
					$intspan_genes_categories->{$g} = init_genes_intspan();
					$intspan_genes_categories->{$g}->{all} = Set::IntSpan::Fast::XS->new();
				}
				foreach my $cat ( keys %{$hgene} ) {
					confess($cat) unless exists $intspan_genes_categories->{$g}->{$cat};
					$intspan_genes_categories->{$g}->{$cat} = $intspan_genes_categories->{$g}->{$cat}->union( $hgene->{$cat} );
					$intspan_genes_categories->{$g}->{all} = $intspan_genes_categories->{$g}->{all}->union( $hgene->{$cat} );
				}
			}
	    }
	 );
	    
	   
	$project->getPatients;
	$project->buffer->dbh_deconnect();
	my $vector_categories_genes;
	foreach my $r (@$ranges){
		my $pid = $pm->start and next;
		my $nb =0;
		
		$project->buffer->dbh_reconnect();
		my $vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>"idefix");
		my $no = $chr->lmdb_variations("r");
		my $h;
		my $intspan_genes_categories_fork;
		for (my $i=$r->[0];$i<$r->[1]+1;$i++){
			$nb ++;
#			warn $nb if $nb%10000 == 0;
			my $variation = $no->get_index($i);
			my $lmdb_index = $i;
			$variation->{project} =  $project;
			$variation->{buffer} = $buffer;
			#purge data
			if ($force){
				$variation->purge_deja_vu();
				$variation->purge_public_data();
			}
	
			my $debug;
			$debug =1 if  $variation->name eq "rs747758472";
			foreach my $g (@{$variation->getGenes()}) {
				unless ( exists $intspan_genes_categories_fork->{ $g->id } ) {
					$intspan_genes_categories_fork->{ $g->id } = init_genes_intspan();
				}
				my $h_spliceAI = $variation->spliceAI_score($g);
				if ($h_spliceAI) {
					foreach my $cat (keys %$h_spliceAI) {
						my $value = $h_spliceAI->{$cat};
						next if ($value eq '-');
						if (defined($value) and $value >= $buffer->config->{spliceAI}->{high}) {
							$intspan_genes_categories_fork->{ $g->id }->{spliceAI_high}->add($lmdb_index);
							$intspan_genes_categories_fork->{ $g->id }->{spliceAI_medium}->add($lmdb_index);
							$intspan_genes_categories_fork->{ $g->id }->{predicted_splice_site}->add($lmdb_index);
						}
						elsif (defined($value) and $value >= $buffer->config->{spliceAI}->{medium}) {
							$intspan_genes_categories_fork->{ $g->id }->{spliceAI_medium}->add($lmdb_index);
							$intspan_genes_categories_fork->{ $g->id }->{predicted_splice_site}->add($lmdb_index);
						}
					}
				}
				
				#$debug =1 if $g->id eq  "ENSG00000253317_8";
				warn $variation->name if $debug;
				warn $lmdb_index if $debug;
				my $cons_text = $variation->variationType($g);
				warn Dumper $cons_text if $debug;
				#die() if $variation->name eq "rs747758472";
				
				foreach my $c ( split( ",", $cons_text ) ) {
					warn $c." ".$lmdb_index if $debug;
					confess( $cons_text . " " . $c ) unless exists $intspan_genes_categories_fork->{ $g->id }->{$c};
					$intspan_genes_categories_fork->{ $g->id }->{$c}->add($lmdb_index);
					warn 	$intspan_genes_categories_fork->{ $g->id }->{$c}->as_string if $debug;
					my $prediction = $variation->categorie_frequency_predicion($g);
					$intspan_genes_categories_fork->{ $g->id }->{$prediction}->add($lmdb_index);
				}
			}
			
		
			######
			# END
			#######
			
		}
		$no->close();
		$no = undef;
		$h->{genes} = $intspan_genes_categories_fork;
		$pm->finish(0, $h);
	}
	$pm->wait_all_children();
	$project->buffer->dbh_reconnect();
	
	my $no4 =  $chr->_get_lmdb('c',"tmp_genes");
	construct_bitVector_for_gene( $no4, $intspan_genes_categories, 1 );
	$no4->close();
	
	my $dir_out = $chr->lmdb_cache_dir();
	my $f1_1 = $dir_out.'/tmp_genes';
	my $f2_1 = $dir_out.'/tmp_genes_index';
	my $f1_2 = $dir_out.'/genes';
	my $f2_2 = $dir_out.'/genes_index';
	`mv $f1_1 $f1_2`;
	`mv $f2_1 $f2_2`;
	
	return 1;
}

sub init_genes_intspan {
	my $hintspan;
	foreach my $c ( keys %{ $categories->{genes}->{annotations} } ) {
		$hintspan->{$c} = Set::IntSpan::Fast::XS->new();
	}
	return $hintspan;
}

sub construct_bitVector_for_gene {
	my ( $no, $intspan_genes_categories, $debug ) = @_;
	# construction d'une table gene du gene et surtout du bitvector pour un gene
	# table de hash du gene : {start} {end} , premier et dernier id du gene donc du intspan. {intspan} et  {size}  et bien sur j'ai rajouté {bitvector} apres tout le process et je sauvegarde
	return unless $intspan_genes_categories;
	foreach my $g ( keys %{$intspan_genes_categories} ) {
		my $hgene = $intspan_genes_categories->{$g};
		# let's start with all variants in genes
		my @array = $hgene->{all}->as_array;
		my $h;
		$h->{start} = $array[0];
		$h->{end}   = $array[-1];
		$h->{size}  = ( $h->{end} - $h->{start} ) + 1;
		$h->{name}  = $g;
		$h->{bitvector}->{all} = Bit::Vector->new_Enum( $h->{size}, join( ',', map { $_ - $h->{start} } @array ) );
		foreach my $cat ( keys %{$hgene} ) {
			next if ( $hgene->{$cat}->is_empty() );
			$h->{intspan}->{$cat} = $hgene->{$cat};
			my @array = map { $_ - $h->{start} } $hgene->{$cat}->as_array;
			$h->{bitvector}->{$cat} = Bit::Vector->new_Enum( $h->{size}, join( ',', @array ) );
		}
		$no->put( $g, $h );
	}
}