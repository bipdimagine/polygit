package Cache_PolyQuery;
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use Parallel::ForkManager;
use Bio::DB::Sam;
use JSON;
use CacheGenesData_bitvector;
use Storable qw(store retrieve freeze thaw);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Cwd 'abs_path';
use Digest::MD5 qw(md5_hex);
use lib "$RealBin/../../GenBo/lib/obj-nodb/packages/";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;


my $categories = {
	global => {
		frequency => {
			freq_none    => 1,
			freq_1       => 1,
			freq_05      => 1,
			freq_01      => 1,
			freq_001     => 1,
			freq_0001    => 1,
			freq_ho_none => 1,
			freq_ho_1    => 1,
			freq_ho_05   => 1,
			freq_ho_01   => 1,
			freq_ho_001  => 1,
			freq_ho_0001 => 1,
			freq_he_none => 1,
			freq_he_1    => 1,
			freq_he_05   => 1,
			freq_he_01   => 1,
			freq_he_001  => 1,
			freq_he_0001 => 1,
			pheno_snp    => 1,
			ngs_score2   => 1,
			ngs_score1   => 1,
			ngs_score0   => 1,
			cadd_not		   => 1,
			cadd_0			   => 1,
			cadd_5			   => 1,
			cadd_10			   => 1,
			cadd_15			   => 1,
			cadd_20			   => 1,
			cadd_25			   => 1,
			cadd_30			   => 1,
			cadd_35			   => 1,
			cadd_40			   => 1,

		},
		variation_type => {
			substitution   => 1,
			insertion      => 1,
			deletion       => 1,
			large_deletion => 1,
		},
	},
	genes => {
		annotations => {
			utr                => 1,
			splicing           => 1,
			pseudogene         => 1,
			coding             => 1,
			maturemirna        => 1,
			essential_splicing => 1,
			phase              => 1,
			silent             => 1,
			intergenic         => 1,
			stop               => 1,
			ncrna              => 1,
			frameshift         => 1,
			intronic           => 1,
			"non-frameshift"   => 1,
			prediction_0       => 1,
			prediction_1       => 1,
			prediction_2       => 1,
			prediction_3       => 1,
		},
	},
	patients => {
		all => 1,
		he  => 1,
		ho  => 1,
	  }
};

sub get_regions {
	my ( $chr, $fork ) = @_;
	if ( $fork == 1 ) {
		my $regions;
		my $hregions;
		$hregions->{chromosome} = $chr->name;
		$hregions->{start}      = 1;
		$hregions->{end}        = $chr->length;
		push( @$regions, $hregions );
		return \@$regions;
	}
	my $tabix = $chr->buffer->software("tabix");
	my %hv;
	foreach my $patient ( @{ $chr->project->getPatients() } ) {
		my $calling_files = $patient->callingFiles();
		my @files;
		foreach my $m ( keys %{$calling_files} ) {
			foreach my $v ( values %{ $calling_files->{$m} } ) {
				push( @files, $v );
			}
		}
		foreach my $file (@files) {
			my $file = $patient->getVariationsFiles->[0];
			next unless -e $file;
			open( TOTO, "$tabix $file " . $chr->ucsc_name() . " | cut -f 2 |" );
			while ( my $pos = <TOTO> ) {
				chomp($pos);
				$hv{$pos}++;
			}
			close TOTO;
		}
	}
	my @snps = sort { $a <=> $b } keys %hv;
	if ( scalar(@snps) == 0 ) {
		my $hregions;
		$hregions->{chromosome} = $chr->name;
		$hregions->{start}      = 1;
		$hregions->{end}        = $chr->length;
		return [$hregions];
	}
	warn "end snp" if ( $chr->project->cache_verbose() );
	my $nb;
	if ( $fork == 1 ) {
		$nb = scalar(@snps);
	}
	else {
		$nb = int( scalar(@snps) / ( $fork - 1 ) );
	}
	$nb = 7_000 if $nb > 7_000;
	my $regions;
	my $iter = natatime $nb, @snps;
	my $old_end;
	while ( my @tmp = $iter->() ) {
		my $hregions;
		$hregions->{chromosome} = $chr->name;
		if ($old_end) {
			$hregions->{start} = $old_end + 1;
		}
		else {
			$hregions->{start} = 1;
		}
		$hregions->{end} = $tmp[-1] + 100;
		$old_end = $hregions->{end};
		push( @$regions, $hregions );
	}
	$regions->[0]->{start} = 1;
	$regions->[-1]->{end}  = $chr->length + 1000;
	return $regions;
}

sub get_ids {
	my ( $project_name, $region ) = @_;
	my $buffer = new GBuffer;
	my $project = $buffer->newProject( -name => $project_name );
	my $reference = $project->getChromosome( $region->{chromosome} )->getReferences( $region->{start}, $region->{end} )->[0];
	my @all;
	my @patient_names = sort { $a cmp $b } map { $_->name } @{ $project->getPatients };
	my $hpatients;
	for ( my $i = 0 ; $i < @patient_names ; $i++ ) {
		$hpatients->{ $patient_names[$i] } = $i;
	}
	foreach my $variation ( @{ $reference->getStructuralVariations } ) {
		my $hv;
		my $array_patients;
		my $aho = [];
		my $ap  = [];
		foreach my $pat ( @{ $variation->getPatients() } ) {
			my $pn = $pat->name();
			my $hp;
			$hp->{name} = $pn;
			my $patient_id = $hpatients->{$pn};
			push( @$ap, $patient_id );
			push( @$aho, $patient_id ) if ( $variation->isHomozygote($pat) );
			$variation->{patients_details}->{$pn}->{vcf_infos} = $variation->{check_id};
			$variation->{patients_details}->{$pn}->{he} = $variation->{annex}->{ $pat->id() }->{he};
			$variation->{patients_details}->{$pn}->{ho} = $variation->{annex}->{ $pat->id() }->{ho};
			$variation->{patients_details}->{$pn}->{he_ho_details} = "he";
			$variation->{patients_details}->{$pn}->{he_ho_details} = "ho" if ( $variation->{annex}->{ $pat->id() }->{ho} eq '1' );
			$variation->{patients_details}->{$pn}->{he_ho_details} .= ':'.$variation->{annex}->{ $pat->id() }->{nb_all_ref}.':'.$variation->{annex}->{ $pat->id() }->{nb_all_mut};
			$variation->{patients_details}->{$pn}->{type} = "he";
			$variation->{patients_details}->{$pn}->{type} = "ho" if ( $variation->isHomozygote($pat) );
		}
		delete $variation->{buffer};
		delete $variation->{project};
		$variation->{heho_string} = join( ",", sort { $a <=> $b } @$ap );
		if ( scalar(@$aho) ) {
			$variation->{heho_string} = $variation->{heho_string}." ".join( ",", sort { $a <=> $b } @$aho ) . " HO";
		}
		# line to prepare dejavu global;
		$hv->{obj}   = compress( freeze($variation) );
		$hv->{start} = $variation->start;
		$hv->{end}   = $variation->end;
		$hv->{id}    = $variation->id;
		push( @all, $hv );
	}
	my @sort = sort { $a->{start} <=> $b->{start} or $a->{end} <=> $b->{end} } @all;
	return \@sort;
}

sub store_ids {
	my ( $chr, $fork ) = @_;
	my $patients = $chr->project->getPatients();
	my $regions;
	my @lVarRegion;
	my $new_pos       = 1;
	my $old_pos       = 1;
	my $nb_variants   = 0;
	my $nb_var_region = 0;
	$regions = get_regions( $chr, $fork );
	my $pm = new Parallel::ForkManager($fork);
	my $all;
	$pm->run_on_finish(
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;
			unless (defined($hRes) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
			#append in the txt file all the variations found in this region and return by getIds
			my ($region) = grep { $_->{start} eq $hRes->{region}->{start} and $_->{end} eq $hRes->{region}->{end} } @$regions;
			die() unless $region;
			my $freeze_file = $chr->lmdb_cache_dir."/". $region->{chromosome}.".".$region->{start}."-".$region->{end}.'.freeze';
			$region->{freeze} = $freeze_file;
			unlink $freeze_file if -e $freeze_file;
			store($hRes->{'variants'}, $freeze_file);
			die() unless -e $freeze_file;
			$nb_variants += scalar(@{$hRes->{'variants'}});
			foreach my $hv (@{$hRes->{'variants'}}) {
				$chr->{cache_hash_get_var_ids}->{$hv->{id}}++;
			}
		}
	);
	my $cpt = 1;
	foreach my $region (@$regions) {
		warn "$cpt/".scalar(@$regions) if ( $chr->project->cache_verbose() );
		$cpt++;
		my $pid = $pm->start and next;
		unless ( exists $region->{chromosome} ) {
			$region->{chromosome} = $chr->id();
		}
		my $time_start = time;
		my $all        = get_ids( $chr->getProject->name, $region );
		my $time_end   = time;
		my $hres;
		$hres->{variants} = $all;
		$hres->{region}   = $region;
		$hres->{ttime}    = abs( $time_start - $time_end );
		$pm->finish( 0, $hres );
	}
	$pm->wait_all_children();
	warn "end process" if ( $chr->project->cache_verbose() );

	#end fork now I have a file with all variation and  json  it's time to sort this file and store it in lmdb database
	#construct intspan 	for patient I will store lmdb_id in intspan
	my $project = $chr->project;
	my @patient_names = sort { $a cmp $b } map { $_->name } @{ $project->getPatients };
	my $index_patients = 0;
	my $hpatients;
	my @categorie_patient = ( "all", "he", "ho" );
	foreach my $pname (@patient_names) {
		$hpatients->{$pname}->{index} = $index_patients++;
		$hpatients->{$pname}->{name}  = $pname;
		foreach my $c (@categorie_patient) {
			$hpatients->{$pname}->{intspan}->{$c} = Set::IntSpan::Fast::XS->new();
		}
	}

	#initialisation global categorie
	my @global = ( "substitution", "insertion", "deletion", "large_deletion" );
	my $intspan_global_type;
	foreach my $g ( keys %{ $categories->{global}->{variation_type} } ) {
		$intspan_global_type->{$g} = Set::IntSpan::Fast::XS->new();
	}

	my $no2 = $chr->get_lmdb_variations("c");    #open lmdb database
	#ok sort and read the file
	my $hh;
	my $uniq;
	my $size_variants = 0;
	foreach my $region (@$regions) {
		my $freeze = $region->{freeze};
		warn Dumper $region  unless -e $freeze;
		die($freeze) unless -e $freeze;
		my $hall = retrieve $freeze;
		foreach my $hv ( @{$hall} ) {
			my $var_id = $hv->{id};
			next if exists $uniq->{$var_id};
			$uniq->{$var_id}++;
			my $variation = thaw( decompress( $hv->{obj} ) );
			my $index_lmdb = $no2->put( $var_id, $variation );
			$size_variants++;
			$hh->{$var_id} = $variation->{heho_string};
			foreach my $pn ( keys %{ $variation->{patients_details} } ) {
				my $type = $variation->{patients_details}->{$pn}->{type};
				$hpatients->{$pn}->{intspan}->{all}->add($index_lmdb);
				$hpatients->{$pn}->{intspan}->{$type}->add($index_lmdb);
			}
			unless ( exists $intspan_global_type->{ $variation->type() } ) {
				warn "\n\nERROR: doesn't exists \$intspan_global_type->{".$variation->type() . "}\n\n";
				warn Dumper keys %$intspan_global_type;
				die;
			}
			$intspan_global_type->{ $variation->type }->add($index_lmdb);
			unlink $freeze;
		}
	}
	$no2->close();

	my $chr_name = $chr->name();
	#store htable only fort dejavu by project purpose
	store( $hh, $project->lmdb_cache_dir . "/$chr_name.dv.freeze" ) if $hh;    
	my $no3 = $chr->get_lmdb_categories("c");
	foreach my $k ( keys %{$intspan_global_type} ) {
		my $h;
		$h->{name}    = $k;
		$h->{intspan} = $intspan_global_type->{$k};
		my $bitv = Bit::Vector->new_Enum( $size_variants, join( ',', $intspan_global_type->{$k}->as_array ) );
		$h->{bitvector} = $bitv;
		$no3->put( $k, $h );
	}
	$no3->close();
	my $no4 = $chr->get_lmdb_patients("c");
	foreach my $pname (@patient_names) {
		my $h;
		$h->{name} = $pname;
		foreach my $c (@categorie_patient) {
			my $intspan = $hpatients->{$pname}->{intspan}->{$c};
			$h->{intspan}->{$c} = $intspan;
			my $bitv = Bit::Vector->new_Enum( $size_variants, join( ',', $intspan->as_array ) );
			$h->{bitvector}->{$c} = $bitv;
		}
		$no4->put( $pname, $h );
	}
	$no4->close;
	return $nb_variants;
}

sub get_annotations {
	my ( $project_name, $chr_name, $region ) = @_;
	my $buffer = new GBuffer;
	my $project = $buffer->newProject( -name => $project_name );
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
	my $id_intergenic_cadd_not = Set::IntSpan::Fast::XS->new();
	my $id_intergenic_cadd_0 = Set::IntSpan::Fast::XS->new();
	my $id_intergenic_cadd_5 = Set::IntSpan::Fast::XS->new();
	my $id_intergenic_cadd_10 = Set::IntSpan::Fast::XS->new();
	my $id_intergenic_cadd_15 = Set::IntSpan::Fast::XS->new();
	my $id_intergenic_cadd_20 = Set::IntSpan::Fast::XS->new();
	my $id_intergenic_cadd_25 = Set::IntSpan::Fast::XS->new();
	my $id_intergenic_cadd_30 = Set::IntSpan::Fast::XS->new();
	my $id_intergenic_cadd_35 = Set::IntSpan::Fast::XS->new();
	my $id_intergenic_cadd_40 = Set::IntSpan::Fast::XS->new();
	my $no = $project->getChromosome($chr_name)->get_lmdb_variations("r");
	my $cursor = $no->cursor( $region->[0], $region->[1] );
	my $nb     = 0;
	my $l      = abs( $region->[0] - $region->[1] );
	my $cpt    = $region->[0];
	while ( my $var_id = $cursor->next_key ) {
		my $variation;
		my $lmdb_index = $cursor->current_index();
		$variation = $no->get($var_id);
		if ( $lmdb_index < 0 ) {
			die("index out of range <0");
		}
		#warn "************************************\n.".$var_id." ".Dumper($variation)."*--------------------*\n********" if $lmdb_index == -2;
		unless ($variation) {
			die();
			$variation = $project->_newVariant($var_id);
		}
		die() unless $variation->id;
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
		#my $cat_freq_he = $variation->categorie_frequency_he();
		#$intspan_global_categories->{$cat_freq_he}->add($lmdb_index);
		if ( $variation->isClinical() ) {
			$intspan_global_categories->{'pheno_snp'}->add($lmdb_index);
		}

		################
		#	get sub/ins/del/large_del categories and construct intspan for each
		###############
		if ( $variation->isVariation() ) {
			$intspan_global_categories->{'substitution'}->add($lmdb_index);
			
			#TODO: here
		}
		elsif ( $variation->isInsertion() ) {
			$intspan_global_categories->{'insertion'}->add($lmdb_index);
		}
		elsif ( $variation->isDeletion() ) {
			$intspan_global_categories->{'deletion'}->add($lmdb_index);
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
		foreach my $g (@$genes) {
			my $cons_text = $variation->variationType($g);
			unless ( exists $intspan_genes_categories->{ $g->id } ) {
				$intspan_genes_categories->{ $g->id } = init_genes_intspan();
			}
			foreach my $c ( split( ",", $cons_text ) ) {
				die( $cons_text . " " . $c ) unless exists $intspan_genes_categories->{ $g->id }->{$c};
				$intspan_genes_categories->{ $g->id }->{$c}->add($lmdb_index);
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
			die("intergenic problem") unless ( $variation->getChromosome()->intergenic_intspan()->contains( $variation->start ) );
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

sub store_annotations {
	my ( $project_name, $chr_name, $fork, $nb_part ) = @_;
	my $buffer  = new GBuffer;
	my $project = $buffer->newProject( -name => $project_name );
	my $chr     = $project->getChromosome($chr_name);
	my $no      = $chr->get_lmdb_variations("r");
	my $hGenesIds;
	unless ($nb_part) { $nb_part = $fork; }
	my $ranges = $no->ranges($nb_part);
	$no->close;
	$no = undef;
	my $intspan_genes_categories  = {};
	my $intspan_global_categories = {};
	foreach my $c ( keys %{ $categories->{global}->{frequency} } ) {
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
				print qq|No message received from child process $exit_code $pid!\n|;
				warn Dumper $hres;
				return;
			}
			if ( exists $hres->{start} ) { print qq|1- No message received from child process $exit_code $pid!\n|;
				warn Dumper $hres;
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
					die($cat) unless exists $intspan_genes_categories->{$g}->{$cat};
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
		foreach my $r (@$ranges) {
			$process->{ $r->[0] . "-" . $r->[1] } = $r;
			my $pid;
			sleep(1);
			$pid = $pm->start and next;
			my $hres;
			$hres->{start} = $r->[0] . "-" . $r->[1];
			my $hres = get_annotations( $project_name, $chr_name, $r );
			$hres->{idp} = $r->[0] . "-" . $r->[1];
			delete $hres->{start};
			$hres->{done} = 1;
			$pm->finish( 0, $hres );
		}    #end for range range
		$pm->wait_all_children();
		
		#all region are ok
		unless ( keys %$process ) { last; }

		#at least one region failed restart only this region
		$true++;
		@$ranges = values %$process;
		#warn "second run";
	}
	die("problem store process") if ( keys %$process );

	my $no = $project->getChromosome($chr_name)->get_lmdb_variations("r");
	my $size_variants = $no->size();
	$no->close();

	#########################
	#GLOBAL CATEGORIES
	# ok now all part is finished I can store the global intspan and construct bitvector for each global categories
	##########################
	my $no3 = $chr->get_lmdb_categories("w");    #open lmdb database

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
	my $no5 = $chr->get_lmdb_patients("c"); #open lmdb database cretae new one or erase older
	foreach my $pn ( keys %{$hintspan_patients} ) {
		my $h;
		$h->{name} = $pn;
		foreach my $type ( keys %{ $hintspan_patients->{$pn} } ) {
			$h->{intspan}->{$type} = $hintspan_patients->{$pn}->{$type};
			my $bitv = Bit::Vector->new_Enum( $size_variants, join( ',', $hintspan_patients->{$pn}->{$type}->as_array ) );
			$h->{bitvector}->{$type} = $bitv;
		}
		$no5->put( $pn, $h );
	}
	$no5->close();
	$no5 = undef;

	#################
	# GENES
	#################
	my @intergenic;
	my $tree = Set::IntervalTree->new;
	my $iter = $intspan_region_intergenic->iterate_runs();
	my $i = 0;
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
		die() unless @$results;
		die() if scalar(@$results) > 1;
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
	return $hGenesIds;
}

sub create_global_infos_freeze {
	my $project = shift;
	my $hashInfos = $project->infosProject();
	my ($hash, $hash_captures, $hash_align, $hash_calling);
	$hash->{global_infos}->{name} = $hashInfos->{name};
	$hash->{global_infos}->{id} = $hashInfos->{id};
	$hash->{global_infos}->{description} = $hashInfos->{description};
	$hash->{global_infos}->{creation_date} = $hashInfos->{creation_date};
	$hash->{global_infos}->{project_type} = $hashInfos->{projectType};
	$hash->{global_infos}->{project_type_id} = $hashInfos->{projectTypeId};
	$hash->{global_infos}->{dbname} = $hashInfos->{dbname};
	$hash->{analyse}->{build} = $project->version();
	$hash->{analyse}->{exome} = $project->isExome();
	$hash->{analyse}->{genome} = $project->isGenome();
	$hash->{analyse}->{diagnostic} = $project->isDiagnostic();
	foreach my $capture (@{$project->getCaptures()}) {
		$hash_captures->{$capture->name()} = undef;
	}
	my @lCapturesName = sort keys %$hash_captures;
	$hash->{analyse}->{capture} = \@lCapturesName;
	my @lTypes = ('evs', '1000genomes', 'dbsnp', 'prediction_matrix', 'cosmic', 'exac');
	foreach my $type (@lTypes) {
		my $root_dir = $project->buffer->config->{public_data}->{$project->version()};
		my $extend =  $project->buffer->config->{kyoto}->{$type};
		my $dir = $root_dir.$extend;
		my $dir_abs_path = abs_path($dir);
		my @lTmp = split('/', $dir_abs_path);
		$hash->{analyse}->{version}->{$type} = $lTmp[-1];
	}
	foreach my $patient (@{$project->getPatients()}) {
		foreach my $method (@{$patient->callingMethods()}) { $hash_calling->{$method} = undef; }
		foreach my $method (@{$patient->alignmentMethods()}) { $hash_align->{$method} = undef; }
		foreach my $file (@{$patient->getVariationsFiles()}) {
			$hash->{check}->{vcf}->{variations}->{$patient->name()}->{$file} = md5_hex($file);
		}
		foreach my $file (@{$patient->getIndelsFiles()}) {
			$hash->{check}->{vcf}->{indels}->{$patient->name()}->{$file} = md5_hex($file);
		}
	}
	my @lMethodsAlign = sort keys(%$hash_align);
	my @lMethodsCalling = sort keys(%$hash_calling);
	$hash->{analyse}->{alignment} = \@lMethodsAlign;
	$hash->{analyse}->{calling} = \@lMethodsCalling;
	$hash->{analyse}->{cache}->{dejavu} = strftime '%Y-%m-%d', localtime;
	unless (-d $project->getCacheBitVectorDir()) {
		my $cmd1 = "mkdir ".$project->getCacheBitVectorDir();
		my $cmd2 = "chmod 777 ".$project->getCacheBitVectorDir();
		`$cmd1`;
		`$cmd2`;
	}
	my $freeze_infos = $project->getCacheBitVectorDir().'/global_infos.freeze';
	`rm $freeze_infos` if (-e $freeze_infos);
	store($hash, $freeze_infos);
	`chmod 777 $freeze_infos`;
	return;
}

sub cache_strictdenovo {
	my ($project_cache, $fork, $chr_name) = @_;
	warn "\n### Cache Strict-Denovo\n";
	my $project_name = $project_cache->name();
	my $buffer = new GBuffer;
	my $project = $buffer->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'familial' );
	$project->getPatients();
	my @lChr;
	foreach my $chr (@{$project->getChromosomes()}){
		next if ($chr_name ne 'all' and $chr_name ne $chr->id());
		next if ($chr->not_used());
		push(@lChr, $chr);
	}
	mkdir ($project->getCacheBitVectorDir().'/strict-denovo') unless (-d $project->getCacheBitVectorDir().'/strict-denovo');
	return 'not_necessary' unless (@lChr);
	my $fasta_file = $project->getGenomeFasta();
	my $hFamInfos;
	foreach my $family (@{$project->getFamilies()}) {
		my $fam_name = $family->name();
		foreach my $patient_name (keys %{$family->children_ill()}) {
			$hFamInfos->{$fam_name}->{children_ill}->{$patient_name}->{list_bam} = $project->getPatient($patient_name)->getBamFiles();
		}
		foreach my $patient_name (keys %{$family->healthy()}) {
			$hFamInfos->{$fam_name}->{healthy}->{$patient_name}->{list_bam} = $project->getPatient($patient_name)->getBamFiles();
		}
	}
	my $hResults;
	foreach my $chr (@lChr) {
		if ($chr->not_used()) {
			my $hResults;
			foreach my $patient (@{$chr->getPatients()}) {
				$hResults->{$patient->name()} = $chr->getNewVector();
			}
			my $nosql = GenBoNoSql->new(dir => $chr->project->getCacheBitVectorDir().'/strict-denovo', mode => 'c');
			$nosql->put_bulk($chr->id(), $hResults);
			$nosql->close();
			next;
		}
		warn '-> Launching denovo model for chr'.$chr->id()."\n";
		my ($hVarPatNoHealthy, $hVarIds_ok, $hVarIds_del);
		$chr->delete_variants();
		my $nb_var_init = $chr->countThisVariants($chr->getVariantsVector());
		$chr->model('denovo');
		$chr->checkModel();
		my $nb_var_denovo = $chr->countThisVariants($chr->getVariantsVector());
		
		#je cree un nosql de suite
		my $hResInit;
		foreach my $patient (@{$chr->getPatients()}) {
			$hResInit->{$patient->name()} = $chr->getNewVector();
		}
		my $nosql = GenBoNoSql->new(dir => $chr->project->getCacheBitVectorDir().'/strict-denovo', mode => 'c');
		$nosql->put_bulk($chr->id(), $hResInit);
		$nosql->close();
		
		#regarde uniquement les variants denovo des familles ayant des patients sains
		my $var_fam_with_healthy = $chr->getNewVector();
		foreach my $family (@{$chr->getFamilies()}) {
			if (scalar(keys %{$family->healthy()} == 0)) {
				foreach my $patient (@{$family->getPatients()}) {
					$hVarPatNoHealthy->{$family->name()}->{$patient->name()} = $chr->getNewVector();
					$hVarPatNoHealthy->{$family->name()}->{$patient->name()} += $patient->getVariantsVector($chr);
					$patient->getVariantsVector($chr)->Empty();
				}
			}
			else {
				foreach my $patient (@{$family->getPatients()}) {
					$var_fam_with_healthy += $patient->getVariantsVector($chr);
				}
			} 
		}
		$chr->getVariantsVector->Intersection($chr->getVariantsVector(), $var_fam_with_healthy);
		my $nb_var_denovo_to_check = $chr->countThisVariants($chr->getVariantsVector());

		my $vector_size = $chr->getVariantsVector->Size();
		my $max_nb_var = 0;
		my $nb_variants = $chr->countThisVariants($chr->getVariantsVector());
		if ($nb_variants > 0) {
			if ($fork == 1) { $max_nb_var = $nb_variants; }
			else { $max_nb_var = int(($nb_variants / $fork) + 2); }
		}
		my $nb_errors = 0;
		my $nb_var = 0;
		my $max = 0;
		my $hRangeVar;
		my $nbRange = 1;
		foreach my $var (@{$chr->getStructuralVariations()}) {
			$nb_var++;
			$max++;
			push(@{$hRangeVar->{$nbRange}}, $var);
			if ($nb_var == $max_nb_var) {
				$nbRange++;
				$nb_var = 0;
			}
		}
		warn "-> Nb var init: $nb_var_init | Nb var denovo: $nb_var_denovo | Nb var to check: $nb_var_denovo_to_check\n";
		warn "-> Nb Var: $nb_variants | Nb Sub-list: ".scalar(keys %{$hRangeVar})." | Nb Var per Sub-list: $max_nb_var\n";
		my $i_progress = 0;
		if (scalar(keys %{$hRangeVar}) > 0) {
			print "-> Strict-denovo fork launch 0% ";
			foreach my $range (keys %{$hRangeVar}) { print "|"; }
			print " 100%\n";
			print "-> Strict-denovo in progress 0% ";
			my $pm = new Parallel::ForkManager($fork);
			$pm->run_on_finish (
				sub {
					my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hResGlobal) = @_;
					if (defined($hResGlobal)) {
						foreach my $v_id (keys %{$hResGlobal}) {
							my $hRes = $hResGlobal->{$v_id};
							my $var_id = $hRes->{id};
							my $vector_id = $hRes->{vector_id};
							foreach my $fam_name (keys %{$hRes->{families}}) {
								if ($hRes->{families}->{$fam_name} eq 'ok') {
									$hVarIds_ok->{$fam_name}->{$var_id} = $vector_id;
								}
								else {
									$hVarIds_del->{$fam_name}->{$var_id} = $hRes->{families}->{$fam_name};
								}
							}
						}
					}
					else {
						$nb_errors++;
						print qq|No message received from child process $pid!\n|;
					}
				}
			);
			foreach my $nb_range (keys %{$hRangeVar}) {
				$pm->start() and next;
				my $hSam;
				my $hTmp;
				foreach my $var (@{$hRangeVar->{$nb_range}}) {
					my $vector_id = $var->vector_id();
					$hTmp->{$vector_id}->{id} = $var->id();
					$hTmp->{$vector_id}->{vector_id} = $vector_id;
					foreach my $fam_name (keys %{$hFamInfos}) {
						$hTmp->{$vector_id}->{families}->{$fam_name} = 'ok';
						my $can_check_bam = 0;
						#verifie si au mois un des enfants malade de cette famille possede bien ce variant.
						foreach my $patient_name (keys %{$hFamInfos->{$fam_name}->{children_ill}}) {
							if ($project->getPatient($patient_name)->getVariantsVector($chr)->contains($vector_id)) {
								$can_check_bam = 1;
							}
						}
						next if ($can_check_bam == 0);
						#check si l un des patients sains passe le strict-denovo
						foreach my $patient_name (keys %{$hFamInfos->{$fam_name}->{healthy}}) {
							foreach my $bam_file (@{$hFamInfos->{$fam_name}->{healthy}->{$patient_name}->{list_bam}}) {
								next if ($hTmp->{$vector_id}->{families}->{$fam_name} eq 'to_del');
								unless (exists $hSam->{$bam_file}) {
									$hSam->{$bam_file} = Bio::DB::Sam->new(-bam=>$bam_file, -fasta=>$fasta_file);
								}
								my $sam = $hSam->{$bam_file};
								my $chr_name = $chr->id();
								my $start = $var->start();
								my $end = $var->end();
								my $all_var = $var->var_allele();
								$all_var = '+' if ($var->isInsertion());
								$all_var = '-' if ($var->isDeletion());
								my $count = 0;
								my $limit = 2;
								my $nb_mut = 0;
								my $callback = sub {
				         			my ($seqid, $pos1, $pileups) = @_;
				         			return if ($pos1 ne $start);
				         			if (scalar(@$pileups) < 3) {
				         				$nb_mut = 99;
				         				return;
				         			}
				         			foreach my $pileup (@$pileups){
										my $b     = $pileup->alignment;
										my $qbase  = substr($b->qseq,$pileup->qpos,1);
										$nb_mut ++ if ($qbase eq $all_var);
										if ($all_var eq "-"){
											$nb_mut ++ if $pileup->indel < 0;
										}
										elsif ($all_var eq "+"){
											$nb_mut ++ if $pileup->indel > 0;
										}
										last if ($nb_mut >= $limit);
				         			}
					 			};
								$sam->fast_pileup("chr$chr_name:$start-$end", $callback);
								$sam = undef;
							if ($nb_mut >= $limit) {
									$hTmp->{$vector_id}->{families}->{$fam_name} = 'to_del';
								}
							} 
						}
					}
				}
				$chr->purge();
				print "|";
				$pm->finish(0, $hTmp);
			}
			$pm->wait_all_children();
			if ($nb_errors > 0) {
				warn ("\n\nERROR: $nb_errors found in fork... DIE...\n\n");
				die();
			}
		}
		else{
			warn "-> Skipped: not necessary\n\n";
		}
		foreach my $family (@{$chr->getFamilies()}) {
			if (exists $hVarPatNoHealthy->{$family->name()}) {
				foreach my $patient (@{$family->getPatients()}) {
					$hResults->{$patient->name()} = $chr->getNewVector();
					$hResults->{$patient->name()} += $hVarPatNoHealthy->{$family->name()}->{$patient->name()};
				}
			}
			else {
				my $vector_ill = $chr->getNewVector();
				foreach my $patient_name (keys %{$hFamInfos->{$family->name()}->{children_ill}}) {
					$vector_ill += $project->getPatient($patient_name)->getVariantsVector($chr);
				}
				my @lVectorIds_fam;
				foreach my $var_id (keys %{$hVarIds_ok->{$family->name()}}) {
					push(@lVectorIds_fam, $hVarIds_ok->{$family->name()}->{$var_id});
				}
				my $vector_fam_ok = Bit::Vector->new_Enum($vector_size, join(',', sort @lVectorIds_fam));
				$vector_fam_ok->Intersection($vector_fam_ok, $vector_ill);
				foreach my $patient (@{$family->getPatients()}) {
					$hResults->{$patient->name()} = $chr->getNewVector();
					$hResults->{$patient->name()}->Intersection($vector_fam_ok, $patient->getVariantsVector($chr));
				}
			}
		}
		my $nosql = GenBoNoSql->new(dir => $chr->project->getCacheBitVectorDir().'/strict-denovo', mode => 'w');
		$nosql->put_bulk($chr->id(), $hResults);
		$nosql->close();
		$chr->purge();
		warn " Done\n\n";
	}
}

sub cache_somatic_loh {
	my ($project_cache, $fork, $chr_name) = @_;
	return unless ($project_cache->isSomaticStudy());
	warn "\n### Cache Somatic LOH\n";
	my $project_name = $project_cache->name();
	my $buffer = new GBuffer;
	my $project = $buffer->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'familial' );
	$project->getPatients();
	my @lChr;
	my $info = "[Somatic LOH Model] Checking CHR$chr_name  ";
	foreach my $chr (@{$project->getChromosomes()}){
		next if ($chr_name ne 'all' and $chr_name eq $chr->id());
		next if ($chr->not_used());
		push(@lChr, $chr);
	}
	mkdir ($project->getCacheBitVectorDir().'/somatic_loh') unless (-d $project->getCacheBitVectorDir().'/somatic_loh');
	return 'not_necessary' unless (@lChr);
	my $fasta_file = $project->getGenomeFasta();
	my $hResults;
	my $i_progress = 0;
	my $pr = String::ProgressBar->new( max => scalar(@lChr), info => $info );
	$pr->write();
	my $pm = new Parallel::ForkManager($fork);
	$pm->run_on_finish (
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;
			$i_progress++;
			$pr->update($i_progress);
			$pr->write();
			if (defined($hRes)) {
				my $nosql = GenBoNoSql->new(dir => $ project->getCacheBitVectorDir().'/somatic_loh', mode => 'c');
				$nosql->put_bulk($hRes->{chr_name}, $hRes->{results});
				$nosql->close();
			}
			else { print qq|No message received from child process $pid!\n|; }
		}
	);
	foreach my $chr (@lChr) {
		$pm->start() and next;
		my $hTmp;
		if ($chr->not_used()) {
			my $hResults;
			foreach my $patient (@{$chr->getPatients()}) {
				$hResults->{$patient->name()} = $chr->getNewVector();
			}
			$hTmp->{chr_name} = $chr->id();
			$hTmp->{results} = $hResults;
			$pm->finish(0, $hTmp) && next;
		}
		$chr->delete_variants();
		$chr->getGenes_no_contruct();
		my $var_tmp = $chr->getNewVector();
		my $var_somatic_ho  = $chr->getNewVector();
		my $var_somatic_he  = $chr->getNewVector();
		my $var_patient = $chr->getNewVector();
		my $var_group = $chr->getNewVector();
		my $var_global = $chr->getNewVector();
		my $no_var_infos = $chr->get_lmdb_variations();
		foreach my $group (@{$chr->getGroups()}) {
			$var_tmp->Empty();
			$var_group->Empty();
			$var_somatic_ho->Empty();
			$var_somatic_he->Empty();
			$group->used_model('loh');
			my @lPatSomatic;
			foreach my $patient (@{$group->getPatientsSomatic()}) {
				$patient->used_model('loh');
				push(@lPatSomatic, $patient);
				$var_somatic_he += $patient->getHe($chr);
				$var_somatic_ho += $patient->getHo($chr);
			}
			foreach my $patient (@{$group->getPatientsGerminal()}) {
				$patient->used_model('loh');
				$var_patient->Empty();
				$var_tmp->Intersection( $var_somatic_ho, $patient->getHe($chr) );
				$var_patient += $var_tmp;
				$var_patient -= $patient->getHo($chr);
				$var_tmp->Intersection( $var_somatic_he, $patient->getHe($chr) );
				foreach my $id (@{$chr->getListVarVectorIds($var_tmp)}) {
					my $var_id = $chr->getVarId($id);
					my $sample_var_1 = $no_var_infos->get($var_id)->{patients_details}->{$patient->name()}->{he_ho_details};
					unless ($sample_var_1) {
						warn "\n\nERROR: pb with $var_id and patient ".$patient->name()."... Die\n\n";
						die ();
					}
					my ($t,$a,$c)  = split(":", $sample_var_1);
					my ($t1,$b,$d);
					my $already;
					foreach my $pat_somatic (@lPatSomatic) {
						next if ($already);
						next unless ($pat_somatic->getHe($chr)->contains($id));
						my $sample_var_2 = $no_var_infos->get($var_id)->{patients_details}->{$pat_somatic->name()}->{he_ho_details};
						unless ($sample_var_2) {
							warn "\n\nERROR: pb with $var_id (somatic patients: ".join(', ', @lPatSomatic).")... Die\n\n";
							die ();
						}
						($t1,$b,$d) = split(":", $sample_var_2);
						$a = 0 unless ($a);
						$b = 0 unless ($b);
						$c = 0 unless ($c);
						$d = 0 unless ($d);
						my $p = fisher::fishers_exact( $a, $b, $c, $d, 1);
						next if $p > 0.01;
						$var_patient->Bit_On($id);
						$already = 1;
					}
				}
				$patient->getVariantsVector($chr)->Intersection( $patient->getVariantsVector($chr), $var_patient );
				$var_group  += $var_patient;
				$var_global += $var_patient;
			}
			foreach my $patient (@{$group->getPatientsSomatic()}) {
				$patient->getVariantsVector($chr)->Intersection( $patient->getVariantsVector($chr), $var_group );
				$hResults->{$patient->name()} = $patient->getVariantsVector($chr);
			}
		}
		$chr->purge();
		$no_var_infos->close();
		$hTmp->{chr_name} = $chr->id();
		$hTmp->{results} = $hResults;
		$pm->finish(0, $hTmp);
	}
	$pm->wait_all_children();
}

sub check_cache_done {
	my ($project, $type, $chr_name, $only_check, $fork) = @_;
	my $project_name = $project->name();
	# check if path project cache exists
	my $path_project_cache = $project->buffer->getDataDirectory("cache").'/'.$project->getVersion().'/'.$project_name;
	check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_project_cache not found...\n\n") unless (-d $path_project_cache);
	warn "Project path cache: OK !\n";
	check_cache_done_polyquery($project, $path_project_cache, $chr_name, $only_check, $fork) if ($type eq 'polyquery');
	check_cache_done_strict_denovo($project, $path_project_cache, $chr_name) if ($type eq 'strict-denovo');
	check_cache_done_somatic_loh($project, $path_project_cache, $chr_name) if ($type eq 'somatic_loh');
	check_cache_done_global_infos($project, $path_project_cache) 	if ($type eq 'global_infos');
	check_cache_done_coverage($project, $path_project_cache)		if ($type eq 'coverage');
	check_cache_done_polydiag($project, $path_project_cache) 		if ($type eq 'polydiag');
	return 1;
}

sub get_ids_for_check_cache {
	my $chr = shift;
	my $regions = get_regions($chr, 1);
	foreach my $region (@{$regions}) {
		unless ( exists $region->{chromosome} ) { $region->{chromosome} = $chr->id(); }
		my $all = get_ids( $chr->getProject->name, $region );
		foreach my $hv (@{$all}) {
			$chr->{cache_hash_get_var_ids}->{$hv->{id}}++;
		}
	}
}

sub print_log_step {
	my ($file, $statut, $message) = @_;
	if (-e $file.'.OK') { `rm $file.OK`; }
	if (-e $file.'.ERROR') { `rm $file.ERROR`; }
	open (LOG_CHR, ">$file.$statut");
	print LOG_CHR $message;
	close (LOG_CHR);
}

sub estim_nb_varids_from_vcf {
	my ($project, $fork) = @_;
	warn "Parsing VCF files [FORK=$fork]\n";
	my (@lVcf, $hPos);
	my @lPatients = @{$project->getPatients()};
	foreach my $patient (@lPatients) {
		foreach my $vcf (@{$patient->getVariationsFiles()}) { push (@lVcf, $vcf) if (-e $vcf); }
		foreach my $vcf (@{$patient->getIndelsFiles()}) { push (@lVcf, $vcf) if (-e $vcf); }
	}
	print '  -> VCF todo: ';
	foreach my $f (@lVcf) { print "|"; }
	print "\n  -> VCF done: ";
	my $pm = new Parallel::ForkManager($fork);
	$pm->run_on_finish(
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;
			if (not defined($hRes) or not exists $hRes->{msg}) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
			#warn $hRes->{msg};
			foreach my $chr (keys %{$hRes->{pos}}) {
				foreach my $id (keys %{$hRes->{pos}->{$chr}}) {
					$hPos->{$chr}->{$id}++;
				}
			}
			print "|";
		}
	);
	foreach my $vcf (@lVcf) {
		my $pid = $pm->start and next;
		my $hres;
		$vcf =~ s/\/\//\//g;
		if ($vcf =~ /.gz/) { open(VCF, "zcat $vcf |"); }
		else { open(VCF, $vcf); }
		while(<VCF>) {
			chomp ($_);
			my $ligne = $_;
			next if ($ligne =~ /#/);
			my @lCol = split(' ', $ligne);
			#only variants, not indels
			next if (length($lCol[3]) != 1);
			next if (length($lCol[4]) != 1);
			$lCol[0] =~ s/chr//;
			$hres->{pos}->{$lCol[0]}->{$lCol[1]}++;
		}
		close(VCF);
		my @lFields = split('/', $vcf);
		my $pat_name = $lFields[-1];
		$pat_name =~ s/\.vcf//;
		$pat_name =~ s/\.sam//;
		$pat_name =~ s/\.gz//;
		$hres->{msg} = '  -> VCF parsed '.$lFields[-2].' - '.$pat_name.' - '.$lFields[-3]."\n";
		$pm->finish( 0, $hres );
	}
	$pm->wait_all_children();
	print "\n";
	return $hPos;
}

sub check_cache_done_polyquery {
	my ($project, $path_project_cache, $chr_name, $only_check, $fork) = @_;
	my $project_name = $project->name();
	# list all used chr in each VCF patients 
	my ($hChr, @lChr);
	my $hPos;
	foreach my $chr (@{$project->getChromosomes()}) {
		next if ($chr_name ne 'all' and $chr_name ne $chr->id());
		if ($only_check == 1) {
			unless ($hPos) { ($hPos) = estim_nb_varids_from_vcf($project, $fork); }
			if (exists $hPos->{$chr->id()}) { $chr->{cache_hash_get_var_ids} = $hPos->{$chr->id()} };
		}
		if (scalar(keys %{$chr->cache_hash_get_var_ids()})) { push(@lChr, $chr); }
		else { warn '  -> WARN: no variation found in CHR'.$chr->id()."\n"; }
	}
	
	# check if path project vector cache exists
	my $path_vector = $path_project_cache.'/vector';
	check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_vector not found...\n\n")  unless (-d $path_vector);
	
	# check if path project lmdb_cache exists
	my $path_lmdb_cache = $path_vector.'/lmdb_cache';
	check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_lmdb_cache not found...\n\n")  unless (-d $path_lmdb_cache);
	
	# check if path project patients cache exists
	my $path_patients = $path_lmdb_cache.'/patients';
	check_cache_in_errors($project, "\n\nERROR CACHE $path_patients -> $path_patients not found...\n\n")  unless (-d $path_patients);
	
	foreach my $chr (@lChr) {
		# check if path project patients cache for each chr exists
		my $path_patients_chr = $path_patients.'/'.$chr->id();
		check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_patients_chr not found...\n\n")  unless (-e $path_patients_chr);
	}
	
	# check if path project variations cache exists
	my $path_variations = $path_lmdb_cache.'/variations';
	check_cache_in_errors($project, "\n\nERROR CACHE $path_patients -> $path_variations not found...\n\n")  unless (-d $path_variations);
	
	foreach my $chr (@lChr) {
		# check if path project variations for each chr cache exists
		my $path_variations_chr = $path_variations.'/'.$chr->id();
		check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_variations_chr not found...\n\n")  unless (-e $path_variations_chr);
		
		# check if path project variations index for each chr cache exists
		my $path_variations_index_chr = $path_variations.'/'.$chr->id().'_index';
		check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_variations_index_chr not found...\n\n")  unless (-e $path_variations_index_chr);
	}
	
	foreach my $chr (@lChr) {
		# check if path project dv freeze for each chr cache exists
		my $path_chr_dv_freeze = $path_lmdb_cache.'/'.$chr->id().'.dv.freeze';
		check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_chr_dv_freeze not found...\n\n")  unless (-e $path_chr_dv_freeze);
		
		# check if path dir for each chr cache exists
		my $path_chr = $path_lmdb_cache.'/'.$chr->id();
		check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_chr not found...\n\n")  unless (-d $path_chr);
		
		# check if path categories_annotations for each chr cache exists
		my $path_chr_categories_annotations = $path_lmdb_cache.'/'.$chr->id().'/categories_annotations';
		check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_chr_categories_annotations not found...\n\n")  unless (-e $path_chr_categories_annotations);
		
		# check if path categories_annotations_index for each chr cache exists
		my $path_chr_categories_annotations_index = $path_lmdb_cache.'/'.$chr->id().'/categories_annotations_index';
		check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_chr_categories_annotations_index not found...\n\n")  unless (-e $path_chr_categories_annotations_index);
		
		# check if path genes for each chr cache exists
		my $path_chr_genes = $path_lmdb_cache.'/'.$chr->id().'/genes';
		check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_chr_genes not found...\n\n")  unless (-e $path_chr_genes);
		
		# check if path genes_index for each chr cache exists
		my $path_chr_genes_index = $path_lmdb_cache.'/'.$chr->id().'/genes_index';
		check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_chr_genes_index not found...\n\n")  unless (-e $path_chr_genes_index);
	}
	
	# check (for each chr) if nb var found from Cache_PolyQuery::get_ids() is equals to nb var found from vectors
	
	foreach my $chr (@lChr) {
		my $path_vector_log = $path_vector.'/log/';
		mkdir ($path_vector_log) unless (-d $path_vector_log);
		my $log_file = $path_vector_log.'/check_cache.'.$chr->id().'.var.log';
		open (LOG, ">$log_file");
		print LOG "#CHR\tVAR_COUNT_GET_IDS\tVAR_COUNT_VECTORS\tSTATUS\n";
		my $nb_var_from_get_ids = scalar(keys %{$chr->cache_hash_get_var_ids()});
		my $buffer_cache  = new GBuffer;
		my $project_cache = $buffer_cache->newProjectCache( -name => $project->name(), -cache => '1', -typeFilters => 'individual' );
		$project_cache->getPatients();
		my $chr_cache = $project_cache->getChromosome($chr->id());
		if ($chr_cache->not_used()) {
			print LOG "chr".$chr->id()."\t$nb_var_from_get_ids\t0\tERROR\nDIE...\n";
			close (LOG);
			`mv $log_file $log_file.ERROR`;
			check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> In chr".$chr->id().", found $nb_var_from_get_ids variants from Cache_PolyQuery::get_ids() methods but 0 variant in vectors...\n\n");
		}
		if ($only_check == 1) {
			$chr_cache->delete_variants('insertion,deletion,large_deletion');
		}
		else {
			$chr_cache->delete_variants();
		}
		my $nb_var_from_cache = $chr_cache->countThisVariants($chr_cache->getVariantsVector());
		if ($only_check == 1) {
			if ($nb_var_from_get_ids > $nb_var_from_cache) {
				print LOG "chr".$chr->id()."\t$nb_var_from_get_ids\t$nb_var_from_cache\tERROR\nDIE...\n";
				close (LOG);
				`mv $log_file $log_file.ERROR`;
				check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> In chr".$chr->id().", found $nb_var_from_get_ids variants from Cache_PolyQuery::get_ids() methods but $nb_var_from_cache variants in vectors...\n\n");
			}
			else {
				print LOG "chr".$chr->id()."\t$nb_var_from_get_ids\t$nb_var_from_cache\tOK\n";
				warn '  -> [CHR'.$chr->id()."] get_ids $nb_var_from_get_ids - vectors: $nb_var_from_cache var - OK !\n";
			}
		}
		else {
			if ($nb_var_from_get_ids != $nb_var_from_cache) {
				print LOG "chr".$chr->id()."\t$nb_var_from_get_ids\t$nb_var_from_cache\tERROR\nDIE...\n";
				close (LOG);
				`mv $log_file $log_file.ERROR`;
				check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> In chr".$chr->id().", found $nb_var_from_get_ids variants from Cache_PolyQuery::get_ids() methods but $nb_var_from_cache variants in vectors...\n\n");
			}
			else {
				print LOG "chr".$chr->id()."\t$nb_var_from_get_ids\t$nb_var_from_cache\tOK\n";
				warn '  -> [CHR'.$chr->id()."] get_ids $nb_var_from_get_ids - vectors: $nb_var_from_cache var - OK !\n";
			}
		}
		`mv $log_file $log_file.OK`;
		close (LOG);
	}
	unless ($only_check == 1) {
		# check (for each chr) if nb genes (not include ''pseudogenes'' from intergenic variants) found from Cache_PolyQuery::get_ids() is equals to nb var found from vectors
		foreach my $chr (@lChr) {
			my $path_vector_log = $path_vector.'/log/';
			my $log_file = $path_vector_log.'/check_cache.'.$chr->id().'.genes.log';
			open (LOG, ">$log_file");
			print LOG "#CHR\tGENES_COUNT_GET_IDS\tGENES_COUNT_VECTORS\tSTATUS\n";
			my $nb_genes_from_get_ids = scalar(keys %{$chr->cache_hash_get_gene_ids()});
			my $buffer_cache  = new GBuffer;
			my $project_cache = $buffer_cache->newProjectCache( -name => $project->name(), -cache => '1', -typeFilters => 'individual' );
			$project_cache->getPatients();
			my $chr_cache = $project_cache->getChromosome($chr->id());
			if ($chr_cache->not_used()) {
				print LOG "chr".$chr->id()."\t$nb_genes_from_get_ids\t0\tERROR\nDIE...\n";
				close (LOG);
				`mv $log_file $log_file.ERROR`;
				check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> In chr".$chr->id().", found $nb_genes_from_get_ids genes from Cache_PolyQuery::get_ids() methods but 0 gene in vectors...\n\n");
			}
			$chr_cache->delete_variants('intergenic');
			$chr_cache->getGenes_no_contruct();
			my $nb_genes_from_cache = scalar(keys %{$chr_cache->{available_genes_ids}});
			if ($nb_genes_from_get_ids != $nb_genes_from_cache) {
				print LOG "chr".$chr->id()."\t$nb_genes_from_get_ids\t$nb_genes_from_cache\tERROR\nDIE...\n";
				check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> In chr".$chr->id().", found $nb_genes_from_get_ids genes from Cache_PolyQuery::get_ids() methods but $nb_genes_from_cache genes in vectors...\n\n");
			}
			else {
				print LOG "chr".$chr->id()."\t$nb_genes_from_get_ids\t$nb_genes_from_cache\tOK\n";
				warn '  -> [CHR'.$chr->id()."] get_ids $nb_genes_from_get_ids - vectors: $nb_genes_from_cache var - OK !\n";
			}
			`mv $log_file $log_file.OK`;
			close (LOG);
		}
	}
	warn "Project cache PolyQuery: OK !\n";
	return;
}

sub check_cache_done_global_infos {
	my ($project, $path_project_cache) = @_;
	my $project_name = $project->name();
	
	# check if path project vector cache exists
	my $path_vector = $path_project_cache.'/vector';
	check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_vector not found...\n\n")  unless (-d $path_vector);
	
	# check if path global_infos file exists
	my $global_infos_file = $path_vector.'/global_infos.freeze';
	check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $global_infos_file not found...\n\n")  unless (-e $global_infos_file);
	
	warn "Project cache Global Infos: OK !\n";
	return;
}

sub check_cache_done_polydiag {
	my ($project, $path_project_cache) = @_;
	my $project_name = $project->name();
	# check if path project polydiag cache exists (if diag)
	if ($project->isDiagnostic()) {
		my $path_polydiag_lite = $path_project_cache.'/polydiag_lite';
		check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_polydiag_lite not found...\n\n") unless (-d $path_polydiag_lite);
	}
	warn "Project cache PolyDiag: OK !\n";
	return;
}

sub check_cache_done_coverage {
	my ($project, $path_project_cache) = @_;
	my $project_name = $project->name();
	# check if path project coverage cache exists
	my $path_coverage = $path_project_cache.'/coverage_lite';
	check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_project_cache not found...\n\n") unless (-d $path_coverage);
	warn "Project cache Coverage: OK !\n";
	return;
}

sub check_cache_done_strict_denovo {
	my ($project, $path_project_cache, $chr_name) = @_;
	my $project_name = $project->name();
	# check if path project vector cache exists
	my $path_vector = $path_project_cache.'/vector';
	check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_vector not found...\n\n")  unless (-d $path_vector);
	
	# check if path project strict-denovo model cache exists (if ped)
	my $path_strict_denovo = $path_vector.'/strict-denovo';
	check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_strict_denovo not found...\n\n")  unless (-d $path_strict_denovo);
	if ($chr_name ne 'all') {
		my $path_strict_denovo_chr = $path_strict_denovo.'/'.$chr_name.'.lite';
		check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_strict_denovo_chr not found...\n\n")  unless (-e $path_strict_denovo_chr);
	}
	
	warn "Project cache Strict-Denovo model: OK !\n";
	return;
}

sub check_cache_done_somatic_loh {
	my ($project, $path_project_cache, $chr_name) = @_;
	my $project_name = $project->name();
	
	# check if path project vector cache exists
	my $path_vector = $path_project_cache.'/vector';
	check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_vector not found...\n\n")  unless (-d $path_vector);
	
	# check if path project somatic_loh model cache exists (is somatic)
	if ($project->isSomaticStudy()) {
		my $path_somatic_loh = $path_vector.'/somatic_loh';
		check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_somatic_loh not found...\n\n")  unless (-d $path_somatic_loh);
		if ($chr_name ne 'all') {
			my $path_somatic_loh_chr = $path_somatic_loh.'/'.$chr_name.'.lite';
			check_cache_in_errors($project, "\n\nERROR CACHE $project_name -> $path_somatic_loh_chr not found...\n\n")  unless (-e $path_somatic_loh_chr);
		}
	}
	
	warn "Project cache Somatic LOH model: OK !\n";
	return;
}

sub check_cache_in_errors {
	my ($project, $error_msg) = @_;
	warn $error_msg;
	my $path_project_cache = $project->buffer->getDataDirectory("cache").'/'.$project->getVersion().'/'.$project->name();
	check_cache_in_errors_move_dir($path_project_cache, 1);
	die;
}

sub check_cache_in_errors_move_dir {
	my ($path_project_cache, $i) = @_;
	my $path_project_cache_errors = $path_project_cache.'_ERROR_'.$i;
	if (-d $path_project_cache_errors) {
		$i++;
		check_cache_in_errors_move_dir($path_project_cache, $i);
	}
	else {
		my $cmd = "mv $path_project_cache $path_project_cache_errors";
		`$cmd`;
		die;
	}
}

sub cache_lite_for_dejavu {
	my ($project_name, $chr_name, $fork) = @_;
	warn "\n### Cache For Deja Vu\n";
	my $buffer1 = new GBuffer;
	my $project = $buffer1->newProject( -name => $project_name, -verbose => 1 );
	my $root_dir = $project->deja_vu_lite_dir() . "/projects/";
	mkdir $root_dir unless -e $root_dir;
	unlink $root_dir . "/" . $project_name . ".lite" if -e $root_dir . "/" . $project_name . ".lite";
	my @chr_names = map { $_->name } @{ $project->getChromosomes };
	my @patient_names =
	sort { $a cmp $b } map { $_->name } @{ $project->getPatients };
	my $dir_out = $project->getCacheBitVectorDir() . "/lmdb_cache";
	my $hpatients;
	for ( my $i = 0 ; $i < @patient_names ; $i++ ) {
		$hpatients->{ $patient_names[$i] } = $i;
	}
	my $no = GenBoNoSql->new( dir => $root_dir, mode => "c" );
	$no->put( $project_name, "patients", $hpatients );
	my $atleast;
	foreach my $chr ( @{ $project->getChromosomes } ) {
		next if ($chr_name ne 'all' and $chr_name ne $chr->id());
		my $fileout = $dir_out . "/" . $chr->name . ".dv.freeze";
		warn "miss $fileout " unless -e $fileout;
		next unless -e $fileout;
		$atleast++;
		my $h = retrieve $fileout;
		$no->put( $project_name, $chr->name, $h );
	}
	confess() unless $atleast;
}

1;
