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
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime); 
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp; 
use GBuffer;
use Cache_Commons;
use Sys::Hostname;
 my $host = hostname();


warn "*_*_*_*_*_".$host."*_*_*_*_*_";


my $fork = 1;
my ($project_name, $chr_name, $annot_version);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'annot_version=s'    => \$annot_version,
	'chr=s'        => \$chr_name,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }



warn "\n### CACHE: store annotations step\n";
my $nbErrors = 0;
my $buffer  = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProject( -name => $project_name );
if ($annot_version) {
	$project->changeAnnotationVersion($annot_version);
}
my $chr     = $project->getChromosome($chr_name);
my $no      = $chr->get_lmdb_variations("r");

if ($no->size == 0 or -e $project->getCacheBitVectorDir()."/lmdb_cache/".$chr->id().".empty") {
	my $no = $chr->get_lmdb_categories("c");
	$no->create();
	warn $project->getCacheBitVectorDir()."/lmdb_cache/".$chr->id();
	my $cmd = "touch ".$project->getCacheBitVectorDir()."/lmdb_cache/".$chr->id()."/genes";
	`$cmd`;
	my $cmd2 = "touch ".$project->getCacheBitVectorDir()."/lmdb_cache/".$chr->id()."/genes_index";
	`$cmd2`;
	warn "empty ". $project->getCacheBitVectorDir()."/lmdb_cache/".$chr->id().".empty";
	warn "here";
	warn $no->size();
	$no->close();
	my $no5 = $chr->get_lmdb_patients("c");
	$no5->create();
	my $h;
	my $vnull =  Bit::Vector->new(0);
	$h->{bitvector}->{'ho'} = $vnull->Shadow();
	$h->{bitvector}->{'he'} = $vnull->Shadow();
	$h->{bitvector}->{'all'} = $vnull->Shadow();
	foreach my $p (@{$project->getPatients}){
		$no5->put($p->name,$h);
	}
	$no5->close();
	exit(0);
}
#my $no      = $chr->get_lmdb_variations("r");
my $hGenesIds;
my $ranges = $no->ranges($fork);
$no->close;
$no = undef;
my $intspan_genes_categories  = {};
my $intspan_global_categories = {};
my $categories = Cache_Commons::categories();
foreach my $c ( keys %{ $categories->{global}->{frequency} } ) {
	$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
}

foreach my $c ( keys %{ $categories->{global}->{variation_type} } ) {
	$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
}

my $process;
my $pm = new Parallel::ForkManager($fork);
#for intergenic variant construct region
my $intspan_region_intergenic = Set::IntSpan::Fast::XS->new();
#and store intergenic id;
my $id_intergenic = Set::IntSpan::Fast::XS->new();
my $hintspan_patients;

####################################
# one fork is finished i could "union"  all the intspan construct in the "sub annotation()"
#run on finish 	aggregate all result (intpsan) from child
#########################################
$pm->run_on_finish(
	sub {
		my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $hres ) = @_;
		unless ( defined($hres) or $exit_code > 0 ) {
			$nbErrors++;
			print qq|No message received from child process $exit_code $pid!\n|;
			warn Dumper $hres;
			return;
		}
		if ( exists $hres->{start} ) { print qq|1- No message received from child process $exit_code $pid!\n|;
			return;
		}
		delete $process->{ $hres->{idp} };
		#################
		#For frequencies  annotation categorie (intspan on chromosome)
		#################
		$intspan_region_intergenic = $intspan_region_intergenic->union( $hres->{intergenic}->{intspan_region} );
		$id_intergenic = $id_intergenic->union( $hres->{intergenic}->{id} );
		
		foreach my $cat ( keys %{ $hres->{frequency} } ) {
			next unless ( exists $intspan_global_categories->{$cat} );
			$intspan_global_categories->{$cat} = $intspan_global_categories->{$cat}->union( $hres->{frequency}->{$cat} );
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
#end run on finish

##############
#fork annotation
##############
my $project_name = $chr->project->name();
my $chr_name     = $chr->name();
my $true = 0;
my $idp  = 0;
while ( $true < 2 ) {   
	#check if all the region end correctly unless restart region failed
	my $cpt = 1;
	foreach my $r (@$ranges) {
		warn "$cpt/".scalar(@$ranges) if ( $project->cache_verbose() );
		$cpt++;
		$process->{ $r->[0] . "-" . $r->[1] } = $r;
		my $pid;
		sleep(1);
		$pid = $pm->start and next;
		my $hres;
		$hres->{start} = $r->[0] . "-" . $r->[1];
		my $hres = get_annotations( $project_name, $chr_name, $r, $annot_version );
		$hres->{idp} = $r->[0] . "-" . $r->[1];
		delete $hres->{start};
		$hres->{done} = 1;
		$pm->finish( 0, $hres );
	}    #end for range range
	$pm->wait_all_children();
	warn "end process" if ( $chr->project->cache_verbose() );
	
	if ($nbErrors > 0) {
		confess("\n\nERRORS: $nbErrors errors found... confess...\n\n");
	}
	
	#all region are ok
	unless ( keys %$process ) { last; }

	#at least one region failed restart only this region
	$true++;
	@$ranges = values %$process;
	#warn "second run";
}
confess("problem store process") if ( keys %$process );
my $no = $project->getChromosome($chr_name)->get_lmdb_variations("r");
my $size_variants = $no->size();
$no->close();

#########################
#GLOBAL CATEGORIES
# ok now all part is finished I can store the global intspan and construct bitvector for each global categories
##########################
warn 'store 1/3: complete lmdb categories';
my $no3 = $chr->get_lmdb_categories("c");    #open lmdb database
foreach my $k ( keys %{$intspan_global_categories} ) {
	my $h;
	$h->{name}    = $k;
	$h->{intspan} = $intspan_global_categories->{$k};
	my $bitv = Bit::Vector->new_Enum( $size_variants, join( ',', $intspan_global_categories->{$k}->as_array ) );
	$h->{bitvector} = $bitv;
	$no3->put( $k, $h );
}
$no3->close();
#########################
#Patients  CATEGORIES
# ok now all part is finished saved patients intspan , all , he,  ho
##########################
warn 'store 2/3: lmdb patients';
my $no5 = $chr->get_lmdb_patients("c"); #open lmdb database cretae new one or erase older
foreach my $pn ( keys %{$hintspan_patients} ) {
	my $h;
	$h->{name} = $pn;
	foreach my $type ( keys %{ $hintspan_patients->{$pn} } ) {
		warn $pn;
		$h->{intspan}->{$type} = $hintspan_patients->{$pn}->{$type};
		my $bitv = Bit::Vector->new_Enum( $size_variants, join( ',', $hintspan_patients->{$pn}->{$type}->as_array ) );
		$h->{bitvector}->{$type} = $bitv;
	}
	$no5->put( $pn, $h );
}
$no5->close();
$no5 = undef;
die();
#################
# GENES
#################
warn 'store 3/3: lmdb genes';
my @intergenic;
my $tree = Set::IntervalTree->new;
my $iter = $intspan_region_intergenic->iterate_runs();
my $i = 0;
delete $chr->{lmdb};

my $no3 = $chr->get_lmdb_variations("r");
while ( my ( $from, $to ) = $iter->() ) {
	my $id = 'intergenic_' . $chr_name . '_' . $from . '_' . $to;
	$tree->insert( $id, $from, $to + 1 );
}
my @array = $id_intergenic->as_array;
my $t = time;
foreach my $lmdb_id (@array) {
	my $v = $no3->get_index($lmdb_id);
	my $results = $tree->fetch( $v->{start}, $v->{end} + 1 );
	confess() unless @$results;
	confess() if scalar(@$results) > 1;
	my $interid = $results->[0];
	unless ( exists $intspan_genes_categories->{$interid} ) {
		$intspan_genes_categories->{$interid} = init_genes_intspan();
		$intspan_genes_categories->{$interid}->{all} = Set::IntSpan::Fast::XS->new();
		$intspan_genes_categories->{$interid}->{intergenic} = Set::IntSpan::Fast::XS->new();
	}
	$intspan_genes_categories->{$interid}->{all}->add($lmdb_id);
	$intspan_genes_categories->{$interid}->{intergenic}->add($lmdb_id);
}
my $no4 = $chr->get_lmdb_genes("c");    #open lmdb database
construct_bitVector_for_gene( $no4, $intspan_genes_categories, 1 );
$no4->close();


sub get_annotations {
	my ( $project_name, $chr_name, $region, $annot_version ) = @_;
	my $buffer = new GBuffer;
	$buffer->vmtouch(1);
	my $project = $buffer->newProject( -name => $project_name );
	$project->preload_patients();
	$project->buffer->disconnect();
	if ($annot_version) {
		$project->changeAnnotationVersion($annot_version);
	}
	# prepare intspan for global categories
	my $intspan_global_categories;
	foreach my $c ( keys %{ $categories->{global}->{frequency} } ) {
		$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
	}
	foreach my $c ( keys %{ $categories->{global}->{variation_type} } ) {
		$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
	}
	my $intspan_genes_categories = {};
	#initialisattion categorie patients
	my $index_patients = 0;
	my $hpatients;
	my @patient_names = sort { $a cmp $b } map { $_->name } @{ $project->getPatients };
	my @categorie_patient = ( "all", "he", "ho" );
	foreach my $pname (@patient_names) {
		foreach my $c (@categorie_patient) {
			$hpatients->{$pname}->{$c} = Set::IntSpan::Fast::XS->new();
		}
	}
	#for intergenic variant construct region
	my $intspan_region_intergenic = Set::IntSpan::Fast::XS->new();
	#and store intergenic id;
	my $id_intergenic = Set::IntSpan::Fast::XS->new();

	my $no = $project->getChromosome($chr_name)->get_lmdb_variations("r");
#	my $cursor = $no->cursor( $region->[0], $region->[1] );
	my $nb     = 0;
	my $l      = abs( $region->[0] - $region->[1] );
	my $cpt    = $region->[0];
	foreach (my $i= $region->[0];$i<=$region->[1];$i++){
	#while ( my $var_id = $cursor->next_key ) {
		my $variation;
		my $lmdb_index = $i;
		$variation = $no->get_index($i);
		if ( $lmdb_index < 0 ) {
			confess("index out of range <0");
		}
		#warn "************************************\n.".$var_id." ".Dumper($variation)."*--------------------*\n********" if $lmdb_index == -2;
		unless ($variation) {
			confess();
		}
		confess() unless $variation->id;
		$variation->{buffer}  = $buffer;
		$variation->{project} = $project;

		################
		# patient  Annotation
		################
		foreach my $pn ( keys %{ $variation->{patients_details} } ) {
			my $type = $variation->{patients_details}->{$pn}->{type};
			$hpatients->{$pn}->{all}->add($lmdb_index);
			$hpatients->{$pn}->{$type}->add($lmdb_index);
		}

		################
		#	get freqency categories and construct intspan for each
		###############
		my $cat_freq = $variation->categorie_frequency();
		
		$intspan_global_categories->{$cat_freq}->add($lmdb_index);
		my $cat_freq_ho = $variation->categorie_frequency_ho();
		$intspan_global_categories->{$cat_freq_ho}->add($lmdb_index);
		if ( $variation->isClinical() ) {
			$intspan_global_categories->{'pheno_snp'}->add($lmdb_index);
		}

		################
		#	get sub/ins/del/large_del categories and construct intspan for each
		###############
		if ( $variation->isVariation() ) {
			$intspan_global_categories->{'substitution'}->add($lmdb_index);
		}
		elsif ( $variation->isInsertion() ) {
			$intspan_global_categories->{'insertion'}->add($lmdb_index);
		}
		elsif ( $variation->isDeletion() ) {
			$intspan_global_categories->{'deletion'}->add($lmdb_index);
		}
		elsif ( $variation->isLargeDeletion() ) {
			$intspan_global_categories->{'large_deletion'}->add($lmdb_index);
		}
		elsif ( $variation->isLargeDuplication() ) {
			$intspan_global_categories->{'large_duplication'}->add($lmdb_index);
		}
		else {
			confess;
		}

		################
		#	get strong/medium/low confidence categories and construct intspan for each
		###############
		my $ngs_score = $variation->getNGSScore();
		if ( $ngs_score == 2 ) {
			$intspan_global_categories->{'ngs_score2'}->add($lmdb_index);
		}
		elsif ( $ngs_score == 1 ) {
			$intspan_global_categories->{'ngs_score1'}->add($lmdb_index);
		}
		elsif ( $ngs_score == 0 ) {
			$intspan_global_categories->{'ngs_score0'}->add($lmdb_index);
		}
		
		my $cadd_score = $variation->cadd_score();
		if ($cadd_score == -1)    { $intspan_global_categories->{cadd_not}->add($lmdb_index);  }
		elsif ($cadd_score >= 40) { $intspan_global_categories->{cadd_40}->add($lmdb_index); }
		elsif ($cadd_score >= 35) { $intspan_global_categories->{cadd_35}->add($lmdb_index); }
		elsif ($cadd_score >= 30) { $intspan_global_categories->{cadd_30}->add($lmdb_index); }
		elsif ($cadd_score >= 25) { $intspan_global_categories->{cadd_25}->add($lmdb_index); }
		elsif ($cadd_score >= 20) { $intspan_global_categories->{cadd_20}->add($lmdb_index); }
		elsif ($cadd_score >= 15) { $intspan_global_categories->{cadd_15}->add($lmdb_index); }
		elsif ($cadd_score >= 10) { $intspan_global_categories->{cadd_10}->add($lmdb_index); }
		elsif ($cadd_score >= 5)  { $intspan_global_categories->{cadd_5}->add($lmdb_index);  }
		else					  { $intspan_global_categories->{cadd_0}->add($lmdb_index);  }
 
		################
		#	get Genes annotations prediction and functional
		###############
		my $genes = $variation->getGenes();
		my $debug;
		$debug =1 if  $variation->name eq "rs747758472";
		foreach my $g (@$genes) {
			unless ( exists $intspan_genes_categories->{ $g->id } ) {
				$intspan_genes_categories->{ $g->id } = init_genes_intspan();
			}
			
			my $h_spliceAI = $variation->spliceAI_score($g);
			if ($h_spliceAI) {
				foreach my $cat (keys %$h_spliceAI) {
					my $value = $h_spliceAI->{$cat};
					next if ($value eq '-');
					if (defined($value) and $value >= $buffer->config->{spliceAI}->{high}) {
						$intspan_global_categories->{spliceAI_high}->add($lmdb_index);
						$intspan_global_categories->{spliceAI_medium}->add($lmdb_index);
						$intspan_global_categories->{predicted_splice_site}->add($lmdb_index);
						$intspan_genes_categories->{ $g->id }->{spliceAI_high}->add($lmdb_index);
						$intspan_genes_categories->{ $g->id }->{spliceAI_medium}->add($lmdb_index);
						$intspan_genes_categories->{ $g->id }->{predicted_splice_site}->add($lmdb_index);
					}
					elsif (defined($value) and $value >= $buffer->config->{spliceAI}->{medium}) {
						$intspan_global_categories->{spliceAI_medium}->add($lmdb_index);
						$intspan_genes_categories->{ $g->id }->{spliceAI_medium}->add($lmdb_index);
						$intspan_global_categories->{predicted_splice_site}->add($lmdb_index);
						$intspan_genes_categories->{ $g->id }->{predicted_splice_site}->add($lmdb_index);
					}
				}
			}
			
			#$debug =1 if $g->id eq  "ENSG00000253317_8";
			my $cons_text = $variation->variationType($g);
			#die() if $variation->name eq "rs747758472";
			
			foreach my $c ( split( ",", $cons_text ) ) {
				confess( $cons_text . " " . $c ) unless exists $intspan_genes_categories->{ $g->id }->{$c};
				$intspan_genes_categories->{ $g->id }->{$c}->add($lmdb_index);
				warn 	$intspan_genes_categories->{ $g->id }->{$c}->as_string if $debug;
				my $prediction = $variation->categorie_frequency_predicion($g);
				$intspan_genes_categories->{ $g->id }->{$prediction}->add($lmdb_index);
			}
		}

		unless (@$genes) {
			################
			# hmm this variant seems to be intergenic
			# construct one intspan for intergenic region : $intspan_region_intergenic later I will intersect this region with real genomic intergenic region
			# and of course I also get the index in another intspan
			###############
			confess("intergenic problem") unless ( $variation->getChromosome()->intergenic_intspan()->contains( $variation->start ) );
			$intspan_region_intergenic->add_range( $variation->start - 100_000, $variation->end + 100_000 );
			$id_intergenic->add($lmdb_index);
		}
		delete $variation->{buffer};
		delete $variation->{project};
		#end frequency
		$nb++;
	}
	my $hres;
	$hres->{frequency} = $intspan_global_categories;
	$hres->{genes}     = $intspan_genes_categories;
	my $ints = $project->getChromosome($chr_name)->intergenic_intspan()->intersection($intspan_region_intergenic);
	$hres->{intergenic}->{intspan_region} = $ints;
	$hres->{intergenic}->{id} = $id_intergenic;
	$hres->{patients} = $hpatients;
	$no->close();
	my $chr     = undef;
	my $project = undef;
	my $buffer  = undef;
	return $hres;
}


#sub get_annotations {
#	my ( $project_name, $chr_name, $region, $annot_version ) = @_;
#	my $buffer = new GBuffer;
#	$buffer->vmtouch(1);
#	my $project = $buffer->newProject( -name => $project_name );
#	$project->preload_patients();
#	$project->buffer->disconnect();
#	if ($annot_version) {
#		$project->changeAnnotationVersion($annot_version);
#	}
#	# prepare intspan for global categories
#	my $intspan_global_categories;
#	foreach my $c ( keys %{ $categories->{global}->{frequency} } ) {
#		$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
#	}
#	foreach my $c ( keys %{ $categories->{global}->{variation_type} } ) {
#		$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
#	}
#	my $intspan_genes_categories = {};
#	#initialisattion categorie patients
#	my $index_patients = 0;
#	my $hpatients;
#	my @patient_names = sort { $a cmp $b } map { $_->name } @{ $project->getPatients };
#	my @categorie_patient = ( "all", "he", "ho" );
#	foreach my $pname (@patient_names) {
#		foreach my $c (@categorie_patient) {
#			$hpatients->{$pname}->{$c} = Set::IntSpan::Fast::XS->new();
#		}
#	}
#	#for intergenic variant construct region
#	my $intspan_region_intergenic = Set::IntSpan::Fast::XS->new();
#	#and store intergenic id;
#	my $id_intergenic = Set::IntSpan::Fast::XS->new();
#	my $id_intergenic_cadd_not = Set::IntSpan::Fast::XS->new();
#	my $id_intergenic_cadd_0 = Set::IntSpan::Fast::XS->new();
#	my $id_intergenic_cadd_5 = Set::IntSpan::Fast::XS->new();
#	my $id_intergenic_cadd_10 = Set::IntSpan::Fast::XS->new();
#	my $id_intergenic_cadd_15 = Set::IntSpan::Fast::XS->new();
#	my $id_intergenic_cadd_20 = Set::IntSpan::Fast::XS->new();
#	my $id_intergenic_cadd_25 = Set::IntSpan::Fast::XS->new();
#	my $id_intergenic_cadd_30 = Set::IntSpan::Fast::XS->new();
#	my $id_intergenic_cadd_35 = Set::IntSpan::Fast::XS->new();
#	my $id_intergenic_cadd_40 = Set::IntSpan::Fast::XS->new();
#	my $no = $project->getChromosome($chr_name)->get_lmdb_variations("r");
#	my $cursor = $no->cursor( $region->[0], $region->[1] );
#	my $nb     = 0;
#	my $l      = abs( $region->[0] - $region->[1] );
#	my $cpt    = $region->[0];
#	while ( my $var_id = $cursor->next_key ) {
#		my $variation;
#		my $lmdb_index = $cursor->current_index();
#		$variation = $no->get($var_id);
#		if ( $lmdb_index < 0 ) {
#			confess("index out of range <0");
#		}
#		#warn "************************************\n.".$var_id." ".Dumper($variation)."*--------------------*\n********" if $lmdb_index == -2;
#		unless ($variation) {
#			confess();
#			$variation = $project->_newVariant($var_id);
#		}
#		confess() unless $variation->id;
#		$variation->{buffer}  = $buffer;
#		$variation->{project} = $project;
#
#		################
#		# patient  Annotation
#		################
#		foreach my $pn ( keys %{ $variation->{patients_details} } ) {
#			my $type = $variation->{patients_details}->{$pn}->{type};
#			$hpatients->{$pn}->{all}->add($lmdb_index);
#			$hpatients->{$pn}->{$type}->add($lmdb_index);
#		}
#
#		################
#		#	get freqency categories and construct intspan for each
#		###############
#		my $cat_freq = $variation->categorie_frequency();
#		
#		$intspan_global_categories->{$cat_freq}->add($lmdb_index);
#		my $cat_freq_ho = $variation->categorie_frequency_ho();
#		$intspan_global_categories->{$cat_freq_ho}->add($lmdb_index);
#		#my $cat_freq_he = $variation->categorie_frequency_he();
#		#$intspan_global_categories->{$cat_freq_he}->add($lmdb_index);
#		if ( $variation->isClinical() ) {
#			$intspan_global_categories->{'pheno_snp'}->add($lmdb_index);
#		}
#
#		################
#		#	get sub/ins/del/large_del categories and construct intspan for each
#		###############
#		if ( $variation->isVariation() ) {
#			$intspan_global_categories->{'substitution'}->add($lmdb_index);
#		}
#		elsif ( $variation->isInsertion() ) {
#			$intspan_global_categories->{'insertion'}->add($lmdb_index);
#		}
#		elsif ( $variation->isDeletion() ) {
#			$intspan_global_categories->{'deletion'}->add($lmdb_index);
#		}
#		elsif ( $variation->isLargeDeletion() ) {
#			$intspan_global_categories->{'large_deletion'}->add($lmdb_index);
#		}
#		elsif ( $variation->isLargeDuplication() ) {
#			$intspan_global_categories->{'large_duplication'}->add($lmdb_index);
#		}
#
#		################
#		#	get strong/medium/low confidence categories and construct intspan for each
#		###############
#		my $ngs_score = $variation->getNGSScore();
#		if ( $ngs_score == 2 ) {
#			$intspan_global_categories->{'ngs_score2'}->add($lmdb_index);
#		}
#		elsif ( $ngs_score == 1 ) {
#			$intspan_global_categories->{'ngs_score1'}->add($lmdb_index);
#		}
#		elsif ( $ngs_score == 0 ) {
#			$intspan_global_categories->{'ngs_score0'}->add($lmdb_index);
#		}
#		
#		my $cadd_score = $variation->cadd_score();
#		if ($cadd_score == -1)    { $intspan_global_categories->{cadd_not}->add($lmdb_index);  }
#		elsif ($cadd_score >= 40) { $intspan_global_categories->{cadd_40}->add($lmdb_index); }
#		elsif ($cadd_score >= 35) { $intspan_global_categories->{cadd_35}->add($lmdb_index); }
#		elsif ($cadd_score >= 30) { $intspan_global_categories->{cadd_30}->add($lmdb_index); }
#		elsif ($cadd_score >= 25) { $intspan_global_categories->{cadd_25}->add($lmdb_index); }
#		elsif ($cadd_score >= 20) { $intspan_global_categories->{cadd_20}->add($lmdb_index); }
#		elsif ($cadd_score >= 15) { $intspan_global_categories->{cadd_15}->add($lmdb_index); }
#		elsif ($cadd_score >= 10) { $intspan_global_categories->{cadd_10}->add($lmdb_index); }
#		elsif ($cadd_score >= 5)  { $intspan_global_categories->{cadd_5}->add($lmdb_index);  }
#		else					  { $intspan_global_categories->{cadd_0}->add($lmdb_index);  }
# 
#		################
#		#	get Genes annotations prediction and functional
#		###############
#		my $genes = $variation->getGenes();
#		my $debug;
#		$debug =1 if  $variation->name eq "rs747758472";
#		foreach my $g (@$genes) {
#			unless ( exists $intspan_genes_categories->{ $g->id } ) {
#				$intspan_genes_categories->{ $g->id } = init_genes_intspan();
#			}
#			
#			my $h_spliceAI = $variation->spliceAI_score($g);
#			if ($h_spliceAI) {
#				foreach my $cat (keys %$h_spliceAI) {
#					my $value = $h_spliceAI->{$cat};
#					next if ($value eq '-');
#					if (defined($value) and $value >= $buffer->config->{spliceAI}->{high}) {
#						$intspan_global_categories->{spliceAI_high}->add($lmdb_index);
#						$intspan_global_categories->{spliceAI_medium}->add($lmdb_index);
#						$intspan_global_categories->{predicted_splice_site}->add($lmdb_index);
#						$intspan_genes_categories->{ $g->id }->{spliceAI_high}->add($lmdb_index);
#						$intspan_genes_categories->{ $g->id }->{spliceAI_medium}->add($lmdb_index);
#						$intspan_genes_categories->{ $g->id }->{predicted_splice_site}->add($lmdb_index);
#					}
#					elsif (defined($value) and $value >= $buffer->config->{spliceAI}->{medium}) {
#						$intspan_global_categories->{spliceAI_medium}->add($lmdb_index);
#						$intspan_genes_categories->{ $g->id }->{spliceAI_medium}->add($lmdb_index);
#						$intspan_global_categories->{predicted_splice_site}->add($lmdb_index);
#						$intspan_genes_categories->{ $g->id }->{predicted_splice_site}->add($lmdb_index);
#					}
#				}
#			}
#			
#			#$debug =1 if $g->id eq  "ENSG00000253317_8";
#			warn $variation->name if $debug;
#			my $cons_text = $variation->variationType($g);
#			warn Dumper $cons_text if $debug;
#			#die() if $variation->name eq "rs747758472";
#			
#			foreach my $c ( split( ",", $cons_text ) ) {
#				warn $c." ".$lmdb_index if $debug;
#				confess( $cons_text . " " . $c ) unless exists $intspan_genes_categories->{ $g->id }->{$c};
#				$intspan_genes_categories->{ $g->id }->{$c}->add($lmdb_index);
#				warn 	$intspan_genes_categories->{ $g->id }->{$c}->as_string if $debug;
#				my $prediction = $variation->categorie_frequency_predicion($g);
#				$intspan_genes_categories->{ $g->id }->{$prediction}->add($lmdb_index);
#			}
#		}
#
#		unless (@$genes) {
#			################
#			# hmm this variant seems to be intergenic
#			# construct one intspan for intergenic region : $intspan_region_intergenic later I will intersect this region with real genomic intergenic region
#			# and of course I also get the index in another intspan
#			###############
#			confess("intergenic problem") unless ( $variation->getChromosome()->intergenic_intspan()->contains( $variation->start ) );
#			$intspan_region_intergenic->add_range( $variation->start - 100_000, $variation->end + 100_000 );
#			$id_intergenic->add($lmdb_index);
#		}
#		delete $variation->{buffer};
#		delete $variation->{project};
#		#end frequency
#		$nb++;
#	}
#	my $hres;
#	$hres->{frequency} = $intspan_global_categories;
#	$hres->{genes}     = $intspan_genes_categories;
#	my $ints = $project->getChromosome($chr_name)->intergenic_intspan()->intersection($intspan_region_intergenic);
#	$hres->{intergenic}->{intspan_region} = $ints;
#	$hres->{intergenic}->{id} = $id_intergenic;
#	$hres->{patients} = $hpatients;
#	$no->close();
#	my $chr     = undef;
#	my $project = undef;
#	my $buffer  = undef;
#	return $hres;
#}

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
