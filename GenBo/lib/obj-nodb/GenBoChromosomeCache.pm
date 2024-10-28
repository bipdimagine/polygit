package GenBoChromosomeCache;
use strict;
use Storable qw(freeze thaw retrieve);
use Moo;
use Carp;
use Bit::Vector;
use Bit::Vector::Overload;
use Set::Intersection;
use Data::Dumper;
use Set::IntSpan::Fast::XS;
use Math::Combinatorics;
use List::MoreUtils qw(uniq);
use base 'Exporter';
use GenBoGeneCache;
use GenBoTranscriptCache;
use GenBoPatientCache;
use GenBoFamilyCache;
use GenBoPrimerCache;
use GenBoNoSql;
use packages::fisher;
use Storable qw(store retrieve freeze dclone thaw);
use Set::IntervalTree;
use List::MoreUtils qw(any uniq natatime lower_bound upper_bound);
extends 'GenBoChromosome', 'GenBoCache';



# objet GenBoChromosomeCache
has chromosome => (
	is      => 'ro',
	lazy    => 1,
	default => sub { return shift; }
);

# tag true if is GenBoChromosome / GenBoChromosomeCache
has isChromosome => (
	is       => 'rw',
	default	 => 1,
);





# tag true if is an empty chromosome
has not_used => (
	is       => 'rw',
	lazy	 => 1,
	default	 => sub {
		my $self = shift;
		confess();#to deleteX
		return 1 unless $self->get_lmdb_categories("r")->is_lmdb_exists('categories_annotations');
		eval { $self->get_lmdb_categories("r")->get_keys(); };
		if ($@) { return 1; }
		return 1 unless $self->get_lmdb_categories("r")->get_keys();
		if ($self->project->filter_chromosome()) {
			return 1 unless (exists $self->project->filter_chromosome->{$self->id()});
		}
		if ($self->getVariantsVector->is_empty()) { return 1; }
		return;
	},
);

# vector defined for all variants in this chromosome
has values_lmdb_genes => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my ($self) = @_;
		confess();
		my $no_genes = $self->get_lmdb_genes("r"); #intspan pour les genes
		 return [] unless $no_genes;
		my $genes = $no_genes->get_values();
		$no_genes->close();
		return $genes;
	},
);
	

# vector defined for all variants in this chromosome
has variants => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		$self->getVectorOrigin();
	},
);

sub getVectorOrigin {
	my ($self,$chr) = @_;
	my $size = $self->size_vector();
	
	$size = 0 unless $size;
	my $v1 = Bit::Vector->new($size);
	$v1->Fill;
	return $v1;
}

sub getNewVector {
	my ($self) = @_;
	my $size = $self->size_vector();
	
	$size = 0 unless $size;
	my $v1 = Bit::Vector->new($size);
	return $v1;
}
sub intersection {
	my ($self,$vector1)=@_;
	confess();
	my $v = $self->variants();
	$v &= $vector1;
 	$self->variants($v); 
}


sub  genes_panels_filters {
	my ($self) = @_;
	my $v = $self->variants;
	my $intspan = Set::IntSpan::Fast::XS->new();
	my $size = $v->Size();
	foreach my $g  (@{$self->project->genes_panels($self)}){
		next unless $g->categories_intspan();
		$intspan = $intspan->union($g->intspan->{all});
		
	}
	my $vector1 = Bit::Vector->new_Enum($size, $intspan->as_string);
	$self->intersection($vector1);
	
	#$self->intersection($vector1):
}

 
sub  patients_only_filters {
	my ($self) = @_;
	my $v1;
		foreach my $patient (@{$self->project->getPatients}){
				unless ($v1){
					$v1 =   $patient->getVectorOrigin($self)->Clone;
				}
				else {
					$v1 += $patient->getVectorOrigin($self);
				}
			}
			$self->intersection($v1);
}


# hash with each each patient variants vector (all, ho, he for each patient)
has patients_categories => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $h;
		unless ($self->size_vector()) {
			foreach my $patient (@{$self->getProject->getPatients()}) {
				$h->{$patient->name()} = Bit::Vector->new(0);
				$h->{$patient->name().'_he'} = Bit::Vector->new(0);
				$h->{$patient->name().'_ho'} = Bit::Vector->new(0);
			}
			return $h;
		}
		foreach my $patient (@{$self->getProject->getPatients()}) {
			my $name = $patient->{'name'};
			$h->{$name} = $patient->getVectorOrigin($self);
			$h->{$name.'_he'} = $patient->getVectorOriginHe($self);
			$h->{$name.'_ho'} = $patient->getVectorOriginHo($self);
		}
		return $h;
	},
);

# hash with each each global categories variants INTSPAN
has global_categories_intspan => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		confess(); # to deleteX
		
		my $h;
		my $no_categories = $self->get_lmdb_categories("r");
		my $categories = $no_categories->get_values();
		foreach my $cat (@$categories) {
			$h->{$cat->{'name'}} = $cat->{'intspan'};
		}
		$no_categories->close();
		return $h;
	},
);

# hash with each each global categories variants VECTOR
has global_categories => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $h = {};
		#HERE 
		confess();
		my $categories = $self->get_list_categories();
		
		my $hh = $self->get_vector_categories($categories);
		foreach my $cat (@$categories) {
			$h->{$cat} = $hh->{$cat};
		}
		
		#return $h;
		
		my ($h_2) = $self->correct_vectors_filters($h);
		return $h_2;
	},
);


# hash with each each categories variants INTSPAN
has categories_intspan => (
	is		=> 'rw',
	lazy    => 1,
	default => sub { {} },
);

# hash with each each global categories variants VECTOR
has categories => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		$self->getGenes();
		return $self->{categories};
	},
);

# excluded variants vector defined by each patient
has variants_excluded => (
	is      => 'rw',
	lazy    => 1,
	default => undef,
);

# intersected variants vector defined by each patient
has variants_intersected => (
	is      => 'rw',
	lazy    => 1,
	default => undef,
);

# hash with all patients informations from cache
has hash_freeze_file_patients => (
	is      => 'rw',
	lazy    => 1,
	default => undef,
);

# hash with all genes informations from cache
has hash_freeze_file_genes => ( 
	is => 'rw',
	lazy    => 1,
	default => undef,
);

# special hash used to select needed genes from annotations (size = nb genes)
has hash_freeze_file_all_genes => (
	is => 'rw',
	lazy    => 1,
	default => undef,
);

sub existsVariationsCacheDir {
	my $self = shift;
	return 1 if -e $self->project->getCacheBitVectorDir().'/lmdb_cache/variations/'.$self->id();
	return;
}

# NOSQL with vectors model somatic loh
has sqlite_loh => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return GenBoNoSql->new(dir => $self->project->getCacheBitVectorDir().'/somatic_loh', mode => 'r');
	},
);

# NOSQL with vectors model strict denovo
has dir_sqlite_strict_denovo => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->project->getCacheBitVectorDir().'/strict-denovo';
	},
);

has sqlite_strict_denovo => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return GenBoNoSql->new(dir => $self->dir_sqlite_strict_denovo(), mode => 'r');
	},
);

has hash_cache_strict_denovo => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $nosql = $self->sqlite_strict_denovo();
		my ($h) = $nosql->get_bulk($self->name());
		$nosql->close();
		return $h;
	},
);

# size of each chromosomes vectors
has size_vector => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		if ($self->project->isRocks){
			
			return $self->rocks_vector("r")->size;
		}
		
		my $no_categories = $self->get_lmdb_categories("r");
		return 0 unless $no_categories->is_lmdb_exists('categories_annotations');
		
		eval { $self->get_lmdb_categories("r")->get_keys(); };
		if ($@) { return 0 }
		return 0  unless $no_categories->get_keys();
		my $b =  $no_categories->get("substitution");
		return  0 unless $b;
		my $v = $b->{bitvector};
		return  $v->Size();
	},
);

# list filter region
has filter_region => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		if (exists $self->project->filter_region->{$self->id()}) { return $self->project->filter_region->{$self->id()}; }
		return;
	},
);

has genes_ids_regions => (
	is      => 'rw',
	lazy    => 1,
	default => undef,
);

has variants_regions_add => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->getNewVector();
	},
);

has variants_regions_exclude => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->getNewVector();
	},
);


# hash stats used for each group
has stats_groups => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $hash;
		return unless ($self->project->isSomaticStudy());
		foreach my $group ( @{ $self->getSomaticGroups() } ) {
			$self->project->print_dot(100);
			$hash->{ $group->name() } = $group->stats($self);
		}
		return $hash;
	},
);

# hash stats for regions ho / rec
has stats_regions_ho_rec => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		confess();
		return if ($self->project->filter_nbvar_regionho() == 0);
		my (@lRegions, $hRegionsDone, $hRegionsByPos);
		if ($self->project->typeFilters() eq 'individual') {
			
			my (@lPatIntersect, @lPatExcluded, $vector_intersect, $vector_excluded);
			foreach my $patient (@{$self->getPatients()}) {
				if ($patient->intersected()) { push(@lPatIntersect, $patient); }
				if ($patient->excluded() and $patient->excluded() eq 'ho_reg') { push(@lPatExcluded, $patient); }
			}
			if (scalar @lPatIntersect > 0) {
				$vector_intersect = $self->getVector_ho_regions_after_intersect(\@lPatIntersect);
				return if ($vector_intersect->is_empty());
			}
			if (scalar @lPatExcluded > 0)  { $vector_excluded = $self->getVector_ho_regions_after_exclude(\@lPatExcluded); }
			foreach my $patient (@{$self->getPatients()}) {
				my $vector_region_ho = $patient->getRegionHo($self)->Clone();
				if ($vector_intersect) { $vector_region_ho->Intersection($vector_region_ho, $vector_intersect); }
				if ($vector_excluded)  { $vector_region_ho -= $vector_excluded; }
				next if ($vector_region_ho->is_empty());
				my @lReg = split(',', $vector_region_ho->to_Enum());
				foreach my $region (@lReg) {
					my $vector_region = Bit::Vector->new_Enum($self->getVariantsVector->Size(), $region);
					next if ($vector_region->is_empty());
					my ($v_id_start, $v_id_end) = split('-', $region);
					my ($start, $end);
					if ($v_id_start == 0) { $start = $self->getVarObject(($v_id_start - 1))->start(); }
					else { $start = $self->getVarObject($v_id_start)->start(); }
					if ($v_id_end) {
						if ($v_id_end == ($self->size_vector() - 1)) { $end = $self->getVarObject($v_id_end)->end(); }
						else { $end = $self->getVarObject(($v_id_end + 1))->end(); }
					}
					else {
						if ($v_id_start == ($self->size_vector() - 1)) { $end = $self->getVarObject($v_id_start)->end(); }
						else { $end = $self->getVarObject(($v_id_start + 1))->end(); }
					}
					my $region_id = $patient->name().';chr'.$self->id().':'.$start.'-'.$end;
					next if (exists $hRegionsDone->{$region_id});
					my $length = $end - $start + 1;
					my $hThisRegion;
					$hThisRegion->{id} = $region_id;
					$hThisRegion->{name} = $patient->name();
					$hThisRegion->{chr} = 'chr'.$self->id();
					$hThisRegion->{start} = $start;
					$hThisRegion->{end} = $end;
					$hThisRegion->{length} = $length;
					$hThisRegion->{nb_genes} = 0;
					my @lGenes;
					foreach my $gene (@{$self->getGenes()}) {
						my $vector_gene = $gene->getVariantsVector->Clone();
						$vector_gene->Intersection($vector_gene, $vector_region);
						unless ($vector_gene->is_empty()) {
							push(@lGenes, $gene->external_name());
							$hThisRegion->{nb_genes}++;
						}
					}
					$hThisRegion->{genes} = join(',', @lGenes);
					$vector_region->Intersection($vector_region, $self->{saved_model}->{init}->{$patient->name().'_ho'});
					$hThisRegion->{nb_var_ho_before} = $self->countThisVariants($vector_region);
					$vector_region->Intersection($vector_region, $patient->getHo($self));
					$hThisRegion->{nb_var_ho_after} = $self->countThisVariants($vector_region);
					$hRegionsDone->{$region_id} = Bit::Vector->new_Enum($self->getVariantsVector->Size(), $region);
					$hRegionsByPos->{$start}->{$patient->name()} = $hThisRegion;
					#push (@lRegions, $hThisRegion);
				}
			}
		}
		elsif ($self->project->typeFilters() eq 'familial') {
			my (@lFamIntersect, @lFamExcluded, $vector_intersect, $vector_excluded);
			foreach my $family (@{$self->getFamilies()}) {
				if ($family->intersected()) { push(@lFamIntersect, $family); }
				if ($family->excluded() and $family->excluded() eq 'rec_reg') { push(@lFamExcluded, $family); }
			}
			if (scalar @lFamIntersect > 0) {
				$vector_intersect = $self->getVector_rec_regions_after_intersect(\@lFamIntersect);
				return if ($vector_intersect->is_empty());
			}
			if (scalar @lFamExcluded > 0)  { $vector_excluded = $self->getVector_rec_regions_after_exclude(\@lFamExcluded); }
			foreach my $family (@{$self->getFamilies()}) {
				my $vector_region_rec = $family->getVectorRegionRec($self);
				if ($vector_intersect) { $vector_region_rec->Intersection($vector_region_rec, $vector_intersect); }
				if ($vector_excluded)  { $vector_region_rec -= $vector_excluded; }
				next if ($vector_region_rec->is_empty());
				my @lReg = split(',', $vector_region_rec->to_Enum());
				foreach my $region (@lReg) {
					my $vector_region = Bit::Vector->new_Enum($self->getVariantsVector->Size(), $region);
					next if ($vector_region->is_empty());
					my ($v_id_start, $v_id_end) = split('-', $region);
					my ($start, $end);
					if ($v_id_start == 0) { $start = $self->getVarObject(($v_id_start - 1))->start(); }
					else { $start = $self->getVarObject($v_id_start)->start(); }
					if ($v_id_end) {
						if ($v_id_end == ($self->size_vector() - 1)) { $end = $self->getVarObject($v_id_end)->end(); }
						else { $end = $self->getVarObject(($v_id_end + 1))->end(); }
					}
					else {
						if ($v_id_start == ($self->size_vector() - 1)) { $end = $self->getVarObject($v_id_start)->end(); }
						else { $end = $self->getVarObject(($v_id_start + 1))->end(); }
					}
					my $region_id = $family->name().';chr'.$self->id().':'.$start.'-'.$end;
					my $length = $end - $start + 1;
					my $hThisRegion;
					$hThisRegion->{id} = $region_id;
					$hThisRegion->{name} = $family->name();
					$hThisRegion->{chr} = 'chr'.$self->id();
					$hThisRegion->{start} = $start;
					$hThisRegion->{end} = $end;
					$hThisRegion->{length} = $length;
					$hThisRegion->{nb_genes} = 0;
					my @lGenes;
					foreach my $gene (@{$self->getGenes()}) {
						my $vector_gene = $gene->getVariantsVector->Clone();
						$vector_gene->Intersection($vector_gene, $vector_region);
						unless ($vector_gene->is_empty()) {
							push(@lGenes, $gene->external_name());
							$hThisRegion->{nb_genes}++;
						}
					}
					$hThisRegion->{genes} = join(',', @lGenes);
					my $vector_fam_ho_init = $self->getNewVector();
					foreach my $patient (@{$family->getPatients()}) {
						$vector_fam_ho_init += $self->{saved_model}->{init}->{$patient->name().'_ho'}
					}
					$vector_region->Intersection($vector_region, $vector_fam_ho_init);
					$hThisRegion->{nb_var_ho_before} = $self->countThisVariants($vector_region);
					$vector_region->Intersection($vector_region, $family->getHo($self));
					$hThisRegion->{nb_var_ho_after} = $self->countThisVariants($vector_region);
					$hRegionsByPos->{$start}->{$family->name()} = $hThisRegion;
					$hRegionsDone->{$region_id} = Bit::Vector->new_Enum($self->getVariantsVector->Size(), $region);
					#push (@lRegions, $hThisRegion);
				}
			}
		}
		
		foreach my $start (sort {$a <=> $b} keys %$hRegionsByPos) {
			foreach my $pat_fam_name (sort keys %{$hRegionsByPos->{$start}}) {
				my @lCommonRegions;
				my $hThisRegion = $hRegionsByPos->{$start}->{$pat_fam_name};
				foreach my $region_id (keys %{$hRegionsDone}) {
					my $v = $hRegionsDone->{$region_id}->Clone();
					$v->Intersection($v, $hRegionsDone->{$hThisRegion->{id}});
					unless ($v->is_empty()) { push(@lCommonRegions, $region_id); }
				}
				$hThisRegion->{common_regions} = join(',', @lCommonRegions);
				push (@lRegions, $hThisRegion);
			}
		}
		
		return if (scalar @lRegions == 0);
		return \@lRegions;
	},
);

# hash avec le nom de tous les filtres qui ont ete supprimes
has hash_filters_deleted => (
	is      => 'rw',
	lazy    => 1,
	default => sub { {} },
);

# hash avec le nom de tous les filtres qui ont ete gardes
has hash_filters_keeped => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $h;
		foreach my $filter_name (keys %{$self->project->hash_ensembl_annotations()}) { $h->{$filter_name} = undef; }
		foreach my $filter_name (keys %{$self->global_categories}) { $h->{$filter_name} = undef; }
		foreach my $filter_name (keys %{$self->project->hash_prediction_filters()}) { $h->{$filter_name} = undef; }
		return $h;
	},
);

# hash avec les genes id des genes dont on garde l'info globale
has available_genes_ids => (
	is      => 'rw',
	lazy    => 1,
	default => sub { {} }
);

# no limit nb gene to construct
has no_limit_genes => (
	is      => 'rw',
	lazy    => 1,
	default => undef,
);

has genes_object_ids => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub { {}	},
);

has cgi_object => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> undef,
);

has hash_gene_id_to_name => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> undef,
);

sub cache_lmdb_variations {
	my $self = shift;
	return $self->get_lmdb_variations("r");
}

sub cache_variations {
	my $self = shift;
	#return $self->get_rocks_variations("r");
	if ($self->project->isRocks){
		return $self->get_rocks_variations("r");
	}
	return $self->get_lmdb_variations("r");;
	
}




has saved_model => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub { {} },
);

has stats_categories => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->project->buffer->config->{stats_chromosomes};
	},
);



##### SET / GET ####



# renvoie un hash avec le comptage de chaque categories (nb snps, ins, del, intronic, silent, utr, etc...)
sub getHashCountTypeVariants {
	my ($self, $useThisVar) = @_;
	my ($hash, $hashTypeVariants);
	if (ref($useThisVar) eq 'Bit::Vector') {
		$hashTypeVariants = $self->getHashTypeVariants($useThisVar);
	}
	else { $hashTypeVariants = $self->getHashTypeVariants($self->getVariantsVector()); }
	foreach my $type (keys %$hashTypeVariants) {
		$hash->{$type} = $self->countThisVariants( $hashTypeVariants->{$type} );
	}
	if (exists $hash->{'large_deletion'}) { $hash->{'cnv'} += $hash->{'large_deletion'}; }
	if (exists $hash->{'large_duplication'}) {  $hash->{'cnv'} += $hash->{'large_duplication'}}; 	
	my @lCat = ('all_variations', 'substitution', 'deletion', 'insertion', 'he', 'ho', 'stop', 'coding', 'cnv');
	foreach my $type (@lCat) {
		$hash->{$type} = 0 unless (exists $hash->{$type});
	}
	$hash->{'he'} = $hash->{'substitution'} + $hash->{'deletion'} + $hash->{'insertion'} - $hash->{'ho'};
	return $hash;
}

# renvoie un hash avec le vector de chaque annotation
sub getHashTypeVariants {
	my ($self, $init_variants) = @_;
	my $hash;
	my @categories = ('categories', 'global_categories');
	foreach my $type_category (@categories) {
		foreach my $category (sort keys %{$self->{$type_category}}) {
			
			my $var_tmp = $self->getNewVector();
			$var_tmp->Intersection( $self->getCategoryVariantsVector($category), $init_variants );
			next if ($var_tmp->is_empty());
			if (exists $hash->{$category}) {
				$hash->{$category} += $var_tmp;
			}
			else { $hash->{$category} = $var_tmp; }
			if (($category eq 'substitution') or ($category eq 'deletion') or ($category eq 'large_deletion')  or ($category eq 'large_duplication') or ($category eq 'insertion')) {
				if (exists $hash->{'all_variations'}) {
					$hash->{'all_variations'} += $var_tmp;
				}
				else { $hash->{'all_variations'} = $var_tmp; }
			}
		}
	}
	if (exists $hash->{'he'} and exists $hash->{'ho'}) {
		$hash->{'he'} -= $hash->{'ho'};
	}
	return $hash;
}

# donne le vector d'une region donnee 
sub getFilterRegionVector {
	my ($self, $filter) = @_;
	confess(warn "\n\nTODO getFilterRegionVector\n\n");
	my ($chrId_filter, $start_filter, $end_filter, $type_filter) = split(':', $filter);
	my @lIdsOn;
	my $found;
	my $still_ok = 1;
	foreach my $index (@{$self->getIdsBitOn( $self->getVariantsVector() )}) {
		last unless ($still_ok);
		#my $varId = $self->getVarId($index);
		my $varId = $self->getProject->returnVariants($self->name."!".$index)->id;
		my ($chrId_var, $pos_var, $ref_var, $alt_var) = split( '_', $varId );
		if (int($start_filter) <= int($pos_var) and int($pos_var) <= int($end_filter)) {
			push (@lIdsOn, $index);
			$found = 1;
		}
		else { $still_ok = undef if ($found); }
	}
	my $var = Bit::Vector->new_Enum($self->getVariantsVector->Size(), join(',', @lIdsOn));
	return $var;
}

sub return_genbo_variant {
	my ($self,$vid) = @_;
	my $no = $self->cache_variations();
	if ($self->project->isRocks){
		
		return  $no->get_variation($vid);
	}
	else {
		return $no->get($vid);
	}
	
	
}


# recupere l'objet variation stocke
sub flushObjectVariantCache {
	my ($self, $vid, $can_create_variant) = @_;
	my $var_obj = $self->return_genbo_variant($vid);
	return undef unless $var_obj;
	$var_obj->{buffer} = $self->project->buffer();
	$var_obj->{project} = $self->project();
	$var_obj->{vector_id} = $var_obj->{index_lmdb};
	my $ref = ref($var_obj);
	if ($ref eq 'GenBoVariation'){
		bless $var_obj , 'GenBoVariationCache';
	}
	elsif  ($ref eq 'GenBoInversion'){
		bless $var_obj , 'GenBoInversionCache';
	}
	elsif  ($ref eq 'GenBoBoundary'){
		bless $var_obj , 'GenBoBoundaryCache';
	}
	elsif  ($ref eq 'GenBoLargeDeletion'){
		bless $var_obj , 'GenBoLargeDeletionCache';
	}
	elsif  ($ref eq 'GenBoLargeDuplication'){
		bless $var_obj , 'GenBoLargeDuplicationCache';
	}
	elsif  ($ref eq 'GenBoDeletion'){
		bless $var_obj , 'GenBoDeletionCache'; 
	}
	elsif  ($ref eq 'GenBoLargeInsertion'){
		bless $var_obj , 'GenBoLargeInsertionCache';
	}
	elsif  ($ref eq 'GenBoInversion'){
		bless $var_obj , 'GenBoInversionnCache';
	}
	elsif  ($ref eq 'GenBoInsertion'){
		bless $var_obj , 'GenBoInsertionCache';
	}
	elsif  ($ref eq 'GenBoMei'){
		bless $var_obj , 'GenBoMeiCache';
	}
	elsif  ($ref ne 'GenBoVariationCache' &&  $ref ne 'GenBoInsertionCache' && $ref ne 'GenBoDeletionCache' && $ref ne 'GenBoLargeDuplicationCache' && $ref ne 'GenBoLargeDeletionCache') {
		confess("\n\nERROR: $vid not found in cache project. Die.\n\n") unless ($can_create_variant);
		#Si l'objet n a pas ete stocke correctement pendant le cache, je construit a la volee sa verion NON cache avec GenBoProject
		my $this_project = $self->getProject();
		$var_obj = $self->project->SUPER::_newVariant($vid);
	}
	my $method;
	if ($var_obj->isVariation()) { $method = 'variations'; }
	elsif ($var_obj->isLargeDeletion()) { $method = 'large_deletions'; }
	elsif ($var_obj->isLargeDuplication()) { $method = 'large_duplications'; }
	elsif ($var_obj->isDeletion()) { $method = 'deletions'; }
	elsif ($var_obj->isInsertion()) { $method = 'insertions'; }
	#warn $var_obj;
	#die $var_obj unless $method;
	die($var_obj) unless $method;
	$self->project->{objects}->{$method}->{$var_obj->id} = $var_obj;
	return 	$self->project->{objects}->{$method}->{$var_obj->id};
}

sub getVarObject {
	my ($self, $v_id) = @_;
	return $self->project->myflushobjects([$self->name."!".$v_id],"variants")->[0];
}

# recupere liste d'objets variations stocke du vecteur
sub getListVarObjects {
	my ($self, $vector) = @_;
	my @lVarObj;
	my $zids;
	foreach my $v_id (@{$self->getListVarVectorIds($vector)}) {
		push(@$zids,$self->name."!".$v_id);
	}
	return $self->project->myflushobjects($zids,"variants");;
}


# recupere l'ID d'un variant a partir d'un vector_id
sub getVarId {
	my ($self, $index) = @_;
	if($self->project->isRocks){
		confess;
		return $self->rocks_vector("r")->get_varid($index);
	}
	return $self->cache_lmdb_variations->get_key_index($index);
}

# recuperer les objets transcripts
sub getTranscripts {
	my $self = shift;
	my $trans = $self->getProject()->myflushobjects($self->setTranscripts(), "transcripts");
	return $trans;
}

sub setTranscripts {
	my $self = shift;
	my $hTrans;
	foreach my $g (@{$self->getGenes()}) {
		foreach my $t (@{$g->getTranscripts()}) {
			$hTrans->{$t->id()} = undef;
		}
	} 
	return $hTrans;
}

#sub setGenes {
#	my $self = shift;
#	my $hIds;
#	return $hIds if ($self->getVariantsVector->is_empty());
#	my @lValuesGenes = @{$self->values_lmdb_genes()};
#	if (scalar(@lValuesGenes) == 0) { return $hIds; }
#	foreach my $hashArgs (@lValuesGenes) {
#		next unless ($hashArgs);
#		$self->project->print_dot(100);
#		my $vector_id = $hashArgs->{index_lmdb};
#		my @lTmp = split('_', $hashArgs->{name});
#		my $obj;
#		if ($lTmp[0] eq 'intergenic') {
#			$obj = $self->getProject->newIntergenic($lTmp[1], $lTmp[2], $lTmp[3]);
#		}
#		else {
#			$obj = $self->project->newGene($hashArgs->{name});
#		}
#		$obj->{vector_id} = $vector_id;
#		$obj->{chromosome} = $self;
#		$obj->{project} = $self->project;			
#		$hIds->{$obj->id()} = undef;
#	}
#	$self->values_lmdb_genes([]);
#	return $hIds;
#}

sub getIntergenicExternalName {
	my ($self, $start, $end) = @_;
	return 'intergenic_'.$self->id().'_'.$start.'_'.$end;
}

# methode qui renvoie un hash d'informations pour un 'pseudo gene' intergenic
sub getHashForIntergenic {
	my ($self, $id, $hashArgs) = @_;
	my $selfEnsemblId;
	$selfEnsemblId = $id;
	my ($description, $chr_name, $start, $end) = split('_', $selfEnsemblId);
	$end = $start unless ($end);
	my $hash;
	$hash->{'genbo_id'} = $selfEnsemblId;
	$hash->{'name'} = $selfEnsemblId;
	$hash->{'external_name'} = $self->getIntergenicExternalName($start, $end);
	$hash->{'chromosome'} = $chr_name;
	$hash->{'start'} = $start;
	$hash->{'end'} = $end;
	$hash->{'strand'} = '1';
	unless ($description or $start or $end or $chr_name) {
		warn 'Chromosome '.$self->id();
		warn $id;
		warn "$description from $start to $end in chr $chr_name";
		die;
	}
	$hash->{'description'} = $description.' from '.$start.' to '.$end.' in chr'.$chr_name;
	$hash->{'transcripts'} = ();
	return $hash;
}


# reecupere un vector (de meme taille que le chromosome) pour un gene donne
sub get_variants_from_this_gene {
	my ($self, $geneId) = @_;
	confess();
}

# recupere l'ID d'un gene a partir d'un vector_id
sub getGeneEnsemblId {
	my ($self, $index) = @_;
	return $self->{hash_gene_id_to_name}->{$index};# if ($self->is_lmdb_cache());
}

sub getPatients {
	my $self = shift;
	return $self->project->getPatients();
}

sub getPatient {
	my ($self, $name) = @_;
	return $self->project->getPatient($name);
}

sub getFamilies {
	my $self = shift;
	return $self->project->getFamilies();
}

sub getFamily {
	my ($self, $name) = @_;
	return $self->project->getFamily($name);
}

sub getSomaticGroups {
	my $self = shift;
	return $self->project->getSomaticGroups();
}

sub getSomaticGroup {
	my ($self, $name) = @_;
	return $self->project->getSomaticGroup($name);
}

sub getNbGenes {
	my $self = shift;
	my $hash;
	return 0 if ($self->getVariantsVector->is_empty());
	my $nb = 0;
	foreach my $gene (@{$self->getGenes()}) {
		next unless $gene->getVectorOrigin();
		next if ($gene->is_intergenic());
		my $v = $self->getVariantsVector & $gene->getCurrentVector();
		$nb++ if (not $v->is_empty()) ;
	}
	return $nb;
}

has genes_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hRes = $self->setGenes();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;	
	}
);

has intervaltree_vector => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
	my $self = shift;
	my $array_tree = $self->rocks_vector("r")->get_vector_gene( "vector_intervaltree" );
	my $tree = Set::IntervalTree->new;
	return $tree unless $tree;
	foreach my $a (@$array_tree){
		next unless @$a;
		
		$tree->insert(@$a);
	}
		return $tree;
	},
);

sub getGenesIdFromVector {
	my ($self,$vector) = @_;
	my $tree = $self->intervaltree_vector;
	my $start = 0;
	my $lh;
	my $toto;
	while (($start < $vector->Size()) &&
    	(my ($min,$max) = $vector->Interval_Scan_inc($start)))
	{
		$max ++;
		 foreach my $g (@{$tree->fetch($min,$max)}){
		 	next if $g =~ /intergenic/;
		 	next if exists $lh->{$g};
		 	push(@$toto,"!".$g);
		 	
		 	$lh->{$g} ++;
		 }
    	$start = $max + 2;
	}
#	$self->project->rocksGenBo->prepare($toto);
#	foreach my $t (@$toto){
#		 $self->project->rocksGenBo->get($t);
#	}
#	warn "end";
	
	return [keys %$lh];
}

sub getGenesFromVector {
	my ($self,$vector) = @_;
	my $tree = $self->intervaltree_vector;
	my $start = 0;
	my $lh;
	while (($start < $vector->Size()) &&
    	(my ($min,$max) = $vector->Interval_Scan_inc($start)))
	{
		$max ++;
		 foreach my $g (@{$tree->fetch($min,$max)}){
		 	next if $g =~ /intergenic/;
		 	
		 	$lh->{$g} ++;
		 	
		 }
    	$start = $max + 2;
	}
	return $self->getProject()->myflushobjects($lh, "genes");
}

sub getNbGenesIntergenic {
	my $self = shift;
	my $hash;
	return 0 if ($self->getVariantsVector->is_empty());
	my $nb = 0;
	foreach my $gene (@{$self->getGenes()}) {
		$nb++ if ($gene->is_intergenic() and not $gene->getCurrentVector->is_empty()) ;
	}
	return $nb;
}


sub vector_global_categories {
	my ($self,$cat) = @_;
	if($self->project->isRocks){
			if ( $cat eq "all" ) {
				my $v = $self->getNewVector();
				$v->Fill();
			return $v;
			}
		return 	$self->getNewVector() if $self->size_vector == 0;
		my $v = $self->rocks_vector("r")->get_vector_chromosome($cat);
		return $v;
	}
	else {
		my $v = $self->hash_vector_global_categories()->{$cat};
		return $self->getNewVector() unless $v; 
		return $v;
	}
	
}
has hash_vector_global_categories => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $h ={};
		my $no_categories = $self->get_lmdb_categories("r");
		my $categories = $no_categories->get_values();
		my $toto = $no_categories->get("pheno_snp");
		foreach my $cat (@$categories) {
			$h->{$cat->{'name'}} = $cat->{'bitvector'};
		}
		$no_categories->close();
		return $h;
	},
);

has public_hgmd_dm_tree => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $tree = Set::IntervalTree->new;
		foreach my $v_id (@{$self->getListVarVectorIds($self->getVectorLmdbDm())}) {
			$tree->insert('DM', $v_id, $v_id+1);
		}
		return $tree;
	},
);
 
#has vector_lmdb_dm => (
#	is		=> 'rw',
#	lazy    => 1,
#	default => sub {
#		my $self = shift;
#		my $vector_all = $self->getNewVector();
#		$vector_all->Fill();
#		my $vector = $self->getNewVector();
#		foreach my $v_id (@{$self->getListVarVectorIds($vector_all)}) {
#			my $var_id = $self->getVarId($v_id);
#			$vector->Bit_On($v_id) if (exists $self->hash_hgmd_class_DM_var_ids->{$var_id});
#		}
#		return $vector;
#	},
#);

sub getVectorLmdbDm {
	my $self = shift;
	if ($self->project->isUpdate()) { 
		return $self->getVectorScore('dm');
	}
	my $vector_dm = $self->getNewVector();
	foreach my $var (@{$self->getListVarObjects($self->getVariantsVector())}) {
		next unless ($var->isDM());
		$vector_dm->Bit_On($var->vector_id()); 
	}
	return $vector_dm;
}

has getVectorLmdbDm_newVersion => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $vector_new_dm = $self->getNewVector();
		my @lVar;
		if ($self->project->isUpdate()) { @lVar = @{$self->getListVarObjects($self->getVectorLmdbDm())}; }
		else { @lVar = @{$self->getListVarObjects($self->getVariantsVector())}; }
		foreach my $var (@lVar) {
			next unless ($var->isDM());
			next unless ($var->isNewHgmd());
			$vector_new_dm->Bit_On($var->vector_id); 
		}
		return $vector_new_dm;
	},
);

sub getVectorConsequence {
 	my ($self,$consequences) = @_;
 	confess() unless $consequences;
 	my $s;
 	foreach my $c (@$consequences){
 		$self->project->getMaskCoding($c);
 		my $v = $self ->getVectorScore($c);
 		unless ($s) {
 			$s = $v;
 		}
 		else {
 			$s |= $v;
 		}
 	}
	return $s;
}

sub getVectorNoFreq {
	my $self = shift;
	return $self->vector_global_categorie("freq_none");
	return $self->getNewVector();
}

sub getVectorPublicVariations {
	my $self = shift;
	my $vector = $self->getVectorVariations->Clone();
	
	return $vector if $self->countThisVariants($vector) == 0;
	$vector -= $self->getVectorNoFreq();
	return $vector;
}

sub getVectorNcboost {
	my ($self, $cat) = @_;
	confess("\n\nERROR: category name mandatory in GenBoChromosomeCache::getVectorNcboost() method. Die.\n\n") unless ($cat);
	confess("\n\nERROR: category name $cat doesn't exists in genbo.cfg. Die.\n\n") unless (exists $self->project->buffer->config->{scaled_score_ncboost});
	return $self->lmdb_score_impact("r")->get('ncboost_'.$cat);
}

sub getVectorVariations {
	my $self = shift;
	my $v = $self->getVectorSubstitutions()->Clone();
	return $v;
	$v += $self->getVectorInsertions();
	$v += $self->getVectorDeletions();
	$v += $self->getVectorLargeDeletions();
	$v += $self->getVectorLargeDuplications();
	return $v;
}

sub getVectorSubstitutions {
	my $self = shift;
	return $self->vector_global_categories("substitution");
}

sub getVectorInsertions {
	my $self = shift;
	return $self->vector_global_categories("insertion");
}

sub getVectorDeletions {
	my $self = shift;
	return $self->vector_global_categories("deletion");
}

sub getVectorLargeDeletions {
	my $self = shift;
	return $self->vector_global_categories("large_deletion");
}

sub getVectorLargeDuplications {
	my $self = shift;
	return $self->vector_global_categories("large_duplication");
}

sub getVectorCnv {
	my $self = shift;
	my $vector = $self->getNewVector();
	$vector += $self->getVectorLargeDeletions() if ($self->getVectorLargeDeletions());
	$vector += $self->getVectorLargeDuplications() if ($self->getVectorLargeDuplications());
	return $vector;
}

sub setVariants {
	my ($self, $type) = @_;
	my $vector = $self->getNewVector();
	if ($type eq 'variations') {
		$vector = $self->getVectorSubstitutions();
	
	}
	elsif ($type eq 'insertions') {
		$vector = $self->getVectorInsertions();
	}
	elsif ($type eq 'deletions') {
			$vector = $self->getVectorDeletions();
	
	}
	elsif ($type eq 'large_deletions') {
			$vector = $self->getVectorLargeDeletions();
		$vector->Intersection( $self->getVariantsVector(), $self->global_categories->{large_deletion} ) if (exists $self->global_categories->{large_deletion});
	}
	elsif ($type eq 'large_duplications') {
		$vector = $self->getVectorDeletions();
		$vector->Intersection( $self->getVariantsVector(), $self->global_categories->{large_duplication} ) if (exists $self->global_categories->{large_duplication});
	}
	else {
		confess();
	}
	
	foreach my $var (@{$self->getListVarObjects($vector)}) {
		$self->{$var->type_object()}->{$var->id()} = undef;
		unless (exists $self->project->{objects}->{$type}->{$var->id()}) {
			$self->project->{objects}->{$type}->{$var->id()} = $var;
		}
	}
}
sub getTreeVariants {
	my ($self,$patient) = @_;
	 my $vector = $patient->getHe($self);
	 my $tree = Set::IntervalTree->new;
	 foreach my $v_id (@{$self->getListVarVectorIds($vector)}) {
	 	my $v = $self->cache_lmdb_variations->get_index($v_id);
	 	$tree->insert($v->id,$v->start,$v->end+1);
		#push(@lVarObj, $self->getVarObject($v_id));
	}
	return $tree;
}

sub getJunctionsVector {
	my $self = shift;
	my $vector = $self->getNewVector();
	$vector += $self->global_categories->{junction} if (exists $self->global_categories->{junction});
	return $vector;
}

sub setJunctions {
	my $self = shift;
	my $vector = $self->getJunctionsVector();
	foreach my $junction (@{$self->getListVarObjects($vector)}) {
		$self->{$junction->type_object()}->{$junction->id()} = undef;
		unless (exists $self->project->{objects}->{junctions}->{$junction->id()}) {
			$self->project->{objects}->{junctions}->{$junction->id()} = $junction;
		}
	}
	return $self->{junctions_object} ;
}
sub return_hash_from_vector{
	my ($self,$vector,$type) = @_;
	foreach my $var (@{$self->getListVarObjects($vector)}) {
		$self->{$var->type_object()}->{$var->id()} = undef;
		unless (exists $self->project->{objects}->{$type}->{$var->id}) {
			$self->project->{objects}->{$type}->{$var->id()} = $var;
		}
	}
}

sub setVariations {
	my $self = shift;
	$self->return_hash_from_vector($self->getVectorVariations,"variations");
	return $self->{variations_object} ;
}

sub setInsertions {
	my $self = shift;
	$self->setVariants('insertions');
	return $self->{insertions_object} ;
}

sub setDeletions {
	my $self = shift;
	$self->setVariants('deletions');
	return $self->{deletions_object} ;
}

sub setLargeDeletions {
	my $self = shift;
	$self->setVariants('large_deletions');
	return $self->{large_deletions_object} ;
}


sub setLargeDuplications {
	my $self = shift;
	$self->setVariants('large_duplications');
	return $self->{large_duplications_object} ;
}

sub getIntSpan {
	my $self = shift;
	my $intspan_chr = Set::IntSpan::Fast::XS->new( $self->getVariantsVector->to_Enum() );
	return $intspan_chr;
}



##### METHODS #####



# methode pour supprimer un type de filtre a GenBoChromosomeCache (interface_json.pl)
sub delete_variants {
	my ($self, $filter_names) = @_;
	confess("\n\nmethod $filter_names GenBoChromosome::construct_tree to supress ???\n\n");
	return unless ($filter_names);
	my $var = $self->getNewVector();
	my $var_tmp = $self->getNewVector();
	my $var_type_delete = $self->getNewVector();
	my $var_freq_delete = $self->getNewVector();
	my ($hFiltersKeeped, $hFilters);
	$self->hash_filters_keeped();
	foreach my $filter_name (split( ',', $filter_names )) {
		$hFilters->{$filter_name} = undef;
		$self->hash_filters_deleted->{$filter_name} = undef;
		if (exists $self->{hash_filters_keeped}->{$filter_name}) {
			delete $self->{hash_filters_keeped}->{$filter_name};
		}
	}
	my @lGlobalFiltersName;
	# keep variants categories NOT filtered
	$var +=  $self->get_vector_categories_not_filtered($hFilters);
	# intersect predictions variants if filtered by predictions
	$self->filter_predictions($hFilters);
	# variants avec type ou predictions a supprimer
	my ($supress_type_pat, $supress_freq_pat, $supress_confidence_pat, $supress_cadd_pat) = $self->get_vectors_type_freq_condidence_cadd_to_supress($hFilters);
	$var_type_delete += $supress_type_pat;
	$var_freq_delete += $supress_freq_pat;
	# recuperation des variants ayant plusieurs predictions dont au moins une a garder  
	my ($keep_type_pat, $keep_freq_pat) = $self->get_vectors_type_freq_to_keep($hFilters);
	$var_type_delete -= $keep_type_pat;
	$var_freq_delete -= $keep_freq_pat;
	$var -= $var_type_delete;
	$var -= $var_freq_delete;
	$var -= $supress_confidence_pat;
	$var -= $supress_cadd_pat;
	$self->setVariantsVector($var);
}

# return un vecteur avec tous les variants des categories NON filtrees
sub get_vector_categories_not_filtered {
	my ($self, $hFilters) = @_;
	my $var = $self->getNewVector();
	foreach my $filter_name (keys %{$self->global_categories}) {
		next if (exists $hFilters->{$filter_name});
		$self->hash_filters_keeped->{$filter_name} = undef;
		$var += $self->getCategoryVariantsVector($filter_name);
	}
	return $var;
}

# rajoute les categories de predictions NON filtrees dans la table de hash KEEP
sub filter_predictions {
	my ($self, $hFilters) = @_;
	foreach my $filter_name (keys %{$self->project->hash_prediction_filters()}) {
		
		if (exists $hFilters->{$filter_name}) {
			$self->hash_filters_deleted->{$filter_name} = undef;
			delete $self->hash_filters_keeped->{$filter_name} if (exists $self->{hash_filters_keeped}->{$filter_name});
		}
	}
	return;
}

# return trois vectors pour les variants de type (snp, ins, ...), les freq et les confidences (ngs_score) qui seront filtrees
sub get_vectors_type_freq_condidence_cadd_to_supress {
	my ($self, $hFilters) = @_;
	my $var_type_delete = $self->getNewVector();
	my $var_freq_delete = $self->getNewVector();
	my $var_confidence_delete = $self->getNewVector();
	my $var_cadd_delete = $self->getNewVector();
	foreach my $filter_name (keys %$hFilters) {
		if ( exists $self->{categories}->{$filter_name} ) {
			$self->hash_filters_deleted->{$filter_name} = undef;
			delete $self->{categories}->{$filter_name};
		}
		elsif ( exists $self->global_categories->{$filter_name} ) {
			if (exists $self->project->hash_var_type_filters->{$filter_name}) {
				$var_type_delete += $self->global_categories->{$filter_name};
			}
			elsif (exists $self->project->hash_frequence_filters->{$filter_name}) {
				$var_freq_delete += $self->global_categories->{$filter_name};
			}
			elsif (exists $self->project->hash_confidence_filters->{$filter_name}) {
				$var_confidence_delete += $self->global_categories->{$filter_name};
			}
			elsif (exists $self->project->hash_cadd_filters->{$filter_name}) {
				$var_cadd_delete += $self->global_categories->{$filter_name};
			}
			
			$self->hash_filters_deleted->{$filter_name} = 'global';
		}
	}
	return ($var_type_delete, $var_freq_delete, $var_confidence_delete, $var_cadd_delete);
}

# return deux vectors pour les variants de type (snp, ins, ...) et freq qui NE seront PAS filtrees
sub get_vectors_type_freq_to_keep {
	my ($self, $hFilters) = @_;
	my $var_type_delete = $self->getNewVector();
	my $var_freq_delete = $self->getNewVector();
	my $var_type_keep = $self->getNewVector();
	my $var_freq_keep = $self->getNewVector();
	foreach my $filter_name (keys %{$self->global_categories}) {
		if (exists $self->project->hash_var_type_filters->{$filter_name}) {
			next if (exists $hFilters->{$filter_name});
			$var_type_keep += $self->global_categories->{$filter_name};
		}
		elsif (exists $self->project->hash_frequence_filters->{$filter_name}) {
			next if (exists $hFilters->{$filter_name});
			$var_freq_keep += $self->global_categories->{$filter_name};
		}
	}
	return ($var_type_keep, $var_freq_keep);
}

# methode pour appliquer le filtre dejavu
sub get_vector_dejavu {
	my ($self, $vector, $max_dejavu, $is_only_ho) = @_;
	confess();
	return unless ($max_dejavu);
	if ($max_dejavu eq 'uniq') { $max_dejavu = 0; }
	my @lNOK;
	foreach my $var (@{$self->getListVarObjects($vector)}) {
		$self->project->print_dot(100);
		my $nb = $var->exome_projects();
		if ($nb > $max_dejavu) {
			push(@lNOK, $var->vector_id);
		}
	}
	return Bit::Vector->new_Enum($self->getVariantsVector->Size(), join(',', @lNOK));
}


# methode qui update le vector global du chromosome a partir des vector de chaque patient
sub update_from_patients {
	my $self = shift;
	my $var = $self->getNewVector();
	foreach my $patient (@{$self->getPatients()}) {
		next if ($patient->in_the_attic());
		$var += $patient->getVariantsVector($self);
	}
	$self->setVariantsVector($var);
}

# methode qui convertit un petit vector (celui d'un gene) a un vecteur de la taille du chromosome (afin de pouvoir les comparer) 
sub convertSubPartToCompleteVector {
	my ($self, $start_subpart, $len_subpart, $subpart) = @_;
	my $vector = $self->getNewVector();
	$vector->Interval_Copy($subpart, $start_subpart, 0, $len_subpart);
	return $vector;
}

# methode sur les filtres regions
sub check_each_var_filter_region {
	my ($self, $filter, $first_launch) = @_;
	my ($chrId, $start_filter, $end_filter, $include) = split(':', $filter);
	if ($include eq '0') {
		$self->check_each_var_filter_region_add($filter, $first_launch);
	}
	elsif ($include eq '-1') {
		$self->check_each_var_filter_region_exclude($filter, $first_launch);
	}
}

# recuperer les variants de cette region
sub check_each_var_filter_region_add {
	my ($self, $filter, $first_launch) = @_;
	my ($chrId, $start_filter, $end_filter, $include) = split(':', $filter);
	next if ($chrId ne $self->id());
	my $var_filter = $self->getFilterRegionVector($filter);
	if ($first_launch) {
		$self->variants_regions_add();
		$self->{variants_regions_add} += $var_filter;
		if ($var_filter->is_empty()) {
			my $hStats = {
				'id' => $filter,
				'chromosome' => $self->id(),
				'start' => int($start_filter),
				'end' => int($end_filter),
				'genes' => 0,
				'substitution' => 0,
				'insertion' => 0,
				'deletion' => 0,
				'include' => $include,
			};
			push(@{$self->project->stats_region()}, $hStats);
		}
		return;
	}
	foreach my $h (@{$self->project->stats_region()}) { return if ($h->{id} eq $filter); }
	my $nbGenes = 0;
	$self->{genes_ids_regions} = {} unless ($self->genes_ids_regions());
	foreach my $gene (@{$self->getGenes()}) {
		next if ($gene->end() < $start_filter);
		next if ($gene->start() > $end_filter);
		my $var_tmp = $self->getNewVector();
		$var_tmp->Intersection( $gene->getVariantsVector(), $var_filter );
		unless ($var_tmp->is_empty()) {
			$nbGenes++;
			$self->{genes_ids_regions}->{$gene->id()} = undef;
		}
	}
	my $v_sub = $self->getCategoryVariantsVector('substitution')->Clone();
	$v_sub->Intersection($v_sub, $var_filter);
	my $v_ins = $self->getCategoryVariantsVector('insertion')->Clone();
	$v_ins->Intersection($v_ins, $var_filter);
	my $v_del = $self->getCategoryVariantsVector('deletion')->Clone();
	$v_del->Intersection($v_del, $var_filter);
	my $hStats = {
		'id' => $filter,
		'chromosome' => $self->id(),
		'start' => int($start_filter),
		'end' => int($end_filter),
		'genes' => int($nbGenes),
		'substitution' => $self->countThisVariants($v_sub),
		'insertion' => $self->countThisVariants($v_ins),
		'deletion' => $self->countThisVariants($v_del),
		'include' => $include,
	};
	push(@{$self->project->stats_region()}, $hStats);
}

# exclue les variants de cette region
sub check_each_var_filter_region_exclude {
	my ($self, $filter, $first_launch) = @_;
	return if ($first_launch);
	my ($chrId, $start_filter, $end_filter, $include) = split(':', $filter);
	next if ($chrId ne $self->id());
	my $var_filter = $self->getFilterRegionVector($filter);
	return if ($var_filter->is_empty());
	$self->variants_regions_exclude();
	$self->{variants_regions_exclude} += $var_filter;
	my $nbGenes = 0;
	foreach my $gene (@{$self->getGenes()}) {
		next if ($gene->end() < $start_filter);
		next if ($gene->start() > $end_filter);
		my $var_tmp = $self->getNewVector();
		$var_tmp->Intersection( $gene->getVariantsVector(), $var_filter );
	}
	my $v_sub = $self->getCategoryVariantsVector('substitution')->Clone();
	$v_sub->Intersection($v_sub, $var_filter);
	my $v_ins = $self->getCategoryVariantsVector('insertion')->Clone();
	$v_ins->Intersection($v_ins, $var_filter);
	my $v_del = $self->getCategoryVariantsVector('deletion')->Clone();
	$v_del->Intersection($v_del, $var_filter);
	my $hStats = {
		'id' => $filter,
		'chromosome' => int($self->id()),
		'start' => int($start_filter),
		'end' => int($end_filter),
		'genes' => int($nbGenes),
		'substitution' => $self->countThisVariants($v_sub),
		'insertion' => $self->countThisVariants($v_ins),
		'deletion' => $self->countThisVariants($v_del),
		'include' => $include,
	};
	push(@{$self->project->stats_region()}, $hStats);
}

# on compare le vector du gene avec celui du chr et on regarde s'il reste des variants
sub convertGeneVectorToGlobalVector {
	my ($self, $hashArgs) = @_;
	my $variants = Bit::Vector->new_Enum($hashArgs->{len_subpart}, join(',',map{$_ - $hashArgs->{start_subpart}} $hashArgs->{getIntSpan}->as_array()));
	return $variants;
}

sub checkGeneExcludedRegion_byGene {
	my ($self, $var_chr_subpart, $hashArgs) = @_;
	my ($vector_global) = $self->convertSubPartToCompleteVector($hashArgs->{start_subpart}, $hashArgs->{len_subpart}, $var_chr_subpart);
	my $var_tmp = $self->getNewVector();
	$var_tmp->Intersection($vector_global, $self->{'variants_excluded'});
	unless ($var_tmp->is_empty()) {
		return $vector_global;
	}
	return;
}

# methode de filtre du quick search
sub check_filter_text {
	my ($self, $gene) = @_;
	return 1 unless ($self->project->filter_text());
	my $name = $gene->external_name();
	my ($description, $phenotypes);
	$description = $gene->description() if ($gene->description());
	$phenotypes = $gene->phenotypes() if ($gene->phenotypes()); 
	my $ok;
	my $search = lc($self->project->filter_text());
	$search =~ s/\+/ /;
	my @lSearch = split(' ', $search);
	foreach my $this_search (@lSearch) {
		if (lc($name) =~ /$this_search/) {
			$ok = 1;
		}
		elsif ($description and lc($description) =~ /$this_search/) {
			$ok = 1;
		}
		elsif ($phenotypes and lc($phenotypes) =~ /$this_search/) {
			$ok = 1;
		}
	}
#	unless ($ok) {
#		if ($self->project->gene_search_ok() and exists $self->project->gene_search_ok->{ $name }) {
#			$ok = 1;
#		}
#		if ($self->project->gene_search_not() and exists $self->project->gene_search_not->{ $name }) {
#			$ok = undef;
#		}
#	}
	return $ok;
}

# methode qui check si un gene passe le filtre region
sub check_gene_stats_region {
	my ($self, $start_gene, $end_gene) = @_;
	if ($self->filter_region()) {
		foreach my $filter (@{$self->filter_region()}) {
			my ($chrId, $start_filter, $end_filter, $type) = split(':', $filter);
			if ($start_gene <= $start_filter and $start_filter <= $end_gene) {
				$self->project->filter_region->{$filter}++;
			}
			elsif ($start_gene <= $end_filter and $end_filter <= $end_gene) {
				$self->project->filter_region->{$filter}++;
			}
			elsif ($start_filter <= $start_gene and $end_gene <= $end_filter) {
				$self->project->filter_region->{$filter}++;
			}
		}
	}
}

# methode qui check si un gene passe le filtre gene atlas
sub check_gene_atlas {
	my ($self, $hashPatGenes, $gene_name) = @_;
	if (exists $self->project->gene_atlas->{$gene_name}) {
		foreach my $this_id (keys %{$self->project->gene_atlas->{$gene_name}->{ids}}) {
			foreach my $patName (keys %$hashPatGenes) {
				$self->project->atlas_patients_cat->{$this_id}->{$patName}++;
			}
			$self->project->atlas_genes_cat->{$this_id}++;
		}
	}
}

# renvoie deux hashs pour donner les informations sur la presence des genes par patient
sub check_gene_count_nb_patient {
	my ($self, $gene_id, $var_gene, $start_subpart, $len_subpart, $listPatients, $hPatVar, $hashIdsByPat) = @_;
	$hashIdsByPat->{not_found} = {};
	my $var_tmp = $self->getNewVector();
	foreach my $patient ( @$listPatients ) {
		$var_tmp->Intersection( $var_gene, $hPatVar->{$patient->name()} );
		if ($var_tmp->is_empty()) {
			$hashIdsByPat->{not_found}->{$patient->name()} = undef;
			next;
		}
		$hashIdsByPat->{found}->{$patient->name()} = undef;
		push( @{$patient->hash_list_genes_ids->{$self->name()}}, $gene_id );
	}
	return $hashIdsByPat;
}

sub construct_tree {
	my ($self,$tree,$vector,$name) = @_;
	foreach my $id (@{$self->getIdsBitOn( $vector )}) {
		$tree->insert($name,$id,$id+1);
	}
}

has variations_hgmd_dm_tree => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $vector = $self->getVectorLmdbDm->Clone();
		$vector->Intersection($vector, $self->getVariantsVector());
		my $tree = Set::IntervalTree->new;
		$self->project->print_dot(1);
	 	$self->construct_tree($tree, $vector, 'HGMD_DM');
		return $tree;	
	}
);

has variations_patients_tree => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $tree = Set::IntervalTree->new;
		foreach my $p (@{$self->getPatients}){
			$self->project->print_dot(1);
		 	$self->construct_tree($tree,$p->getVariantsVector($self),$p->name);
		}	
		return $tree;	
	}
);

has variations_type_tree => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $tree = Set::IntervalTree->new;
		my @types = ("substitution", "insertion", "deletion", "large_deletion","large_duplication");
		foreach my $type (@types){
			$self->project->print_dot(1);
		 	$self->construct_tree($tree, $self->getCategoryVariantsVector($type), $type);
		 }	
		 return $tree;	
	}
);

has variations_genes_tree => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return Set::IntervalTree->new;
	}
);

has genes_tree => (
	is      => 'ro',
	lazy    => 1,
	default => sub { 
		my $self = shift;
		my $tree = Set::IntervalTree->new;
		my $genes = $self->getGenes();
		foreach my $gene (@$genes) {
			$self->project->print_dot(500);
			$self->construct_tree($tree,$gene->getVariantsVector,$gene->id);
		}
		return $tree;
	}
);

		
has fam_genetic_model_tree => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $tree = Set::IntervalTree->new;
		foreach my $model_name ('strict_denovo', 'denovo', 'dominant', 'recessif', 'uniparental_disomy') {
		#foreach my $model_name ('strict_denovo', 'denovo', 'dominant', 'recessif') {
		#foreach my $model_name ('denovo', 'dominant', 'recessif') {
		#foreach my $model_name ('denovo') {
			my $method_name = 'getModelVector_fam_'.$model_name;
			#$self->project->print_dot(1);
		 	$self->construct_tree($tree, $self->$method_name(), $model_name);
		}
		return $tree;
	}
);		
		
sub construct_tree_from_intspan {
	my ($self,$tree,$set,$name) = @_;
	my $iter = $set->iterate_runs();
    while (my ( $from, $to ) = $iter->()) {
       	$tree->insert($name,$from,$to+1);
    }
}

has transcripts_tree => (
	is      => 'ro',
	lazy    => 1,
	default => sub { 
	my $self = shift;
	 my $tree = Set::IntervalTree->new;
	  my $tree2 = Set::IntervalTree->new;
	 my $genes = $self->getGenes();
	my $t = time;
	my $project= $self->project;
	 foreach my $gene (@$genes) {
	 	#warn $gene->getVariantsVector;
	 	foreach my $tr (@{$gene->getTranscripts}){
	 		$tree->insert($tr->id,$tr->start-211,$tr->end+212);
	 		#$self->construct_tree($tree,$gene->getVariantsVector,$tr->id);
	 	}
	 }
	 return $tree;
		
		 }
);

sub getVector_ho_regions_after_exclude {
	my ( $self, $patients, $nbVar) = @_;
	my $vector_del_ho_regions = $self->getNewVector();
	foreach my $patient (@$patients) {
		$vector_del_ho_regions += $patient->getRegionHo($self, $nbVar);
		$patient->excluded('ho_reg');
	}
	return $vector_del_ho_regions;
}

sub getVector_rec_regions_after_exclude {
	my ( $self, $families, $nbVar) = @_;
	my $vector_del_rec_regions = $self->getNewVector();
	foreach my $family (@$families) {
		$vector_del_rec_regions += $family->getVectorRegionRec($self, $nbVar);
		$family->excluded('rec_reg');
	}
	return $vector_del_rec_regions;
}

sub getVector_ho_regions_after_intersect {
	my ( $self, $patients, $nbVar) = @_;
	my $vector_common_ho_regions = $patients->[0]->getRegionHo($self, $nbVar)->Clone();
	foreach my $patient (@$patients) {
		$vector_common_ho_regions->Intersection($vector_common_ho_regions, $patient->getRegionHo($self, $nbVar));
	}
	return $vector_common_ho_regions;
}

sub getVector_rec_regions_after_intersect {
	my ( $self, $families, $nbVar) = @_;
	my $vector_common_rec_regions = $families->[0]->getVectorRegionRec($self, $nbVar)->Clone();
	foreach my $family (@$families) {
		$vector_common_rec_regions->Intersection($vector_common_rec_regions, $family->getVectorRegionRec($self, $nbVar));
		$family->intersected(1);
	}
	return $vector_common_rec_regions;
}



##### MODELS #####



sub save_model_variants_all_patients {
	my ($self, $model_name) = @_;
	$self->{saved_model}->{$model_name}->{chromosome} = $self->getVariantsVector();
	foreach my $patient (@{$self->getPatients()}) {
		$self->{saved_model}->{$model_name}->{$patient->name()} = $self->getNewVector();
		$self->{saved_model}->{$model_name}->{$patient->name().'_he'} = $self->getNewVector();
		$self->{saved_model}->{$model_name}->{$patient->name().'_ho'} = $self->getNewVector();
		return unless ($self->size_vector());
		$self->{saved_model}->{$model_name}->{$patient->name()} += $patient->getVectorOrigin($self);
		$self->{saved_model}->{$model_name}->{$patient->name().'_he'} += $patient->getVectorOriginHe($self);
		$self->{saved_model}->{$model_name}->{$patient->name().'_ho'} += $patient->getVectorOriginHo($self);
	}
}

sub load_init_variants_all_patients {
	my ($self, $model_name) = @_;
	$self->setVariantsVector( $self->{saved_model}->{$model_name}->{chromosome} );
	foreach my $patient (@{$self->getPatients()}) {
		$self->patients_categories->{$patient->name()} = dclone $self->{saved_model}->{$model_name}->{$patient->name()};
		$self->patients_categories->{$patient->name().'_he'} =  dclone $self->{saved_model}->{$model_name}->{$patient->name().'_he'};
		$self->patients_categories->{$patient->name().'_ho'} =  dclone $self->{saved_model}->{$model_name}->{$patient->name().'_ho'};
		if ($model_name eq 'init') { $patient->used_model(''); }
		else { $patient->used_model($model_name); }
	}
	foreach my $family (@{$self->getFamilies()}) {
		if ($model_name eq 'init') { $family->used_model(''); }
		else { $family->used_model($model_name); }
	}
}

# applique le model dbl_evt (SOMATIC)
sub getModelGeneVector_som_dbl_evt {
	my ($self, $gene) = @_;
	my $ok;
	my $var_global = $self->getNewVector();
	foreach my $group (@{$self->getSomaticGroups()}) {
		$var_global += $group->checkModel_somatic_dbl_evt($gene, $self);
	}
	return ($var_global);
}

# applique le model strict_denovo (FAM)
sub checkModel_from_cache {
	my ($self, $var, $sqlite_name, $model_name) = @_;
	return 1 if ($var);
	confess unless ($self->$sqlite_name());
	my $nosql = $self->$sqlite_name();
	my ($h) = $nosql->get_bulk($self->name());
	$nosql->close();
	my $var_global = $self->getNewVector();
	if ($h) {
		foreach my $patient (@{$self->getPatients()}) {
			my $vector_filtre = $self->getNewVector();
			$vector_filtre += $h->{$patient->name()} if ($h and exists $h->{$patient->name()});
			if ($vector_filtre->is_empty()) {
				$self->patients_categories->{$patient->name()} = $self->getNewVector();
				$self->patients_categories->{$patient->name().'_he'} = $self->getNewVector();
				$self->patients_categories->{$patient->name().'_ho'} = $self->getNewVector();
			}
			else {
				$self->patients_categories->{$patient->name()}->Intersection($self->patients_categories->{$patient->name()}, $vector_filtre);
				$self->patients_categories->{$patient->name().'_he'}->Intersection($self->patients_categories->{$patient->name().'_he'}, $vector_filtre);
				$self->patients_categories->{$patient->name().'_ho'}->Intersection($self->patients_categories->{$patient->name().'_ho'}, $vector_filtre);
			}
			$var_global += $patient->getVariantsVector($self);
			$patient->used_model($model_name);
		}
		foreach my $family (@{$self->getFamilies()}) {
			$family->used_model($model_name);
		}
	}
	$self->setVariantsVector($var_global);
}

# applique le model loh (SOMATIC)
sub getModelVector_som_loh {
	my ($self, $var) = @_;
	return $self->checkModel_from_cache($var, 'sqlite_loh', 'loh');
}

sub getModelVector_fam_mosaique {
	my $self = shift;
	return $self->{vector_transmission}->{mosaique} if exists $self->{vector_transmission}->{mosaique};
	my $var_global  = $self->getNewVector();
	foreach my $family (@{$self->getFamilies()}) {
		$var_global += $family->getModelVector_fam_mosaique($self);
	}
	$self->get_lmdb_variations()->close();
	$self->{vector_transmission}->{mosaique} = $var_global;
	return $var_global;
}
	
# applique le model recessif (FAM)
sub getModelVector_fam_recessif {
	my ($self, $is_recessif_compound) = @_;
	return $self->{vector_transmission}->{recessif} if exists $self->{vector_transmission}->{recessif};
	my $var_global = $self->getNewVector();
	foreach my $family (@{$self->getFamilies()}) {
		$var_global += $family->getModelVector_fam_recessif($self);
	}
	$self->{vector_transmission}->{recessif} = $var_global;
	return $var_global;
}

# applique le model denovo (FAM)
sub getModelVector_fam_denovo {
	my $self = shift;
	return $self->{vector_transmission}->{denovo} if exists $self->{vector_transmission}->{denovo};
	my $var_global = $self->getNewVector();
	foreach my $family (@{$self->getFamilies()}) {
		$var_global += $family->getModelVector_fam_denovo($self);
	}
	$self->{vector_transmission}->{denovo} = $var_global;
	return $var_global;
}

# applique le model strict_denovo (FAM)
sub getModelVector_fam_strict_denovo {
	my ($self, $var) = @_;
	return 1 if ($var);
	return $self->{vector_transmission}->{strict_denovo} if exists $self->{vector_transmission}->{strict_denovo};
	my $var_global = $self->getNewVector();
	foreach my $family (@{$self->getFamilies()}) {
		$var_global += $family->getModelVector_fam_strict_denovo($self);
	}
	$self->{vector_transmission}->{strict_denovo} = $var_global;
	return $var_global;
}

# applique le model dominant (FAM)
sub getModelVector_fam_dominant {
	my ($self, $var) = @_;
	return 1 if ($var);
	return $self->{vector_transmission}->{dominant} if exists $self->{vector_transmission}->{dominant};
	my $var_global = $self->getNewVector();
	foreach my $family (@{$self->getFamilies()}) {
		$var_global += $family->getModelVector_fam_dominant($self);
	}
	$self->{vector_transmission}->{dominant} = $var_global;
	return $var_global;
}

sub getModelVector_fam_uniparental_disomy {
	my $self = shift;
	return $self->{vector_transmission}->{uniparental_disomy} if exists $self->{vector_transmission}->{uniparental_disomy};
	my $var_global = $self->getNewVector();
	foreach my $family (@{$self->getFamilies()}) {
		$var_global += $family->getModelVector_fam_uniparental_disomy($self);
	}
	$self->{vector_transmission}->{uniparental_disomy} = $var_global;
	return $var_global;
}

#TODO: a verifier;
# applique le model no_model (variants presents chez les deux parents) (FAM)
sub getModelVector_fam_both_parents {
	my $self = shift;
	my $vector_both_parents = $self->getNewVector();
	foreach my $family (@{$self->getFamilies()}) {
		my $vector_fam_ok = $self->getNewVector();
		foreach my $patient (@{$family->getChildren()}) {
			next if ($patient->in_the_attic());
			$vector_fam_ok += $patient->getVariantsVector($self);
		}
		next if ($vector_fam_ok->is_empty());
		my @lParents = @{$family->getParents()};
		next if (scalar(@lParents) < 2);
		my $vector_parents = $lParents[0]->getVariantsVector($self)->Clone();
		$vector_parents->Intersection($vector_parents, $lParents[1]->getVariantsVector($self));
		$vector_both_parents += $vector_parents;
	}
	return $vector_both_parents;
}


# applique le model compound (IND)
sub getModelGeneVector_indiv_compound {
	my ($self, $gene, $var_gene_annot_var1, $var_gene_annot_var2) = @_;
	my $ok;
	my $var_tmp = $self->getNewVector();
	my $var_global_patients = $self->getNewVector();
	foreach my $patient (@{$self->getPatients()}) {
		$patient->{model_indiv_compound}->{$self->id()} = $self->getNewVector() unless (exists $patient->{model_indiv_compound}->{$self->id()});
		$var_tmp->Empty();
		$patient->used_model('compound');
		$var_tmp->Intersection( $patient->getVariantsVector($self), $gene->getCurrentVector() );
		next if ($var_tmp->is_empty());
		unless ($self->countThisVariants($var_tmp) == 1) {
#			if ($var_gene_annot_var1 and $var_gene_annot_var2) {
#				$var_gene_annot_var1->Intersection($var_gene_annot_var1, $var_tmp);
#				$var_gene_annot_var2->Intersection($var_gene_annot_var2, $var_tmp);
#				if ($self->countThisVariants($var_gene_annot_var1) >= 1 and $self->countThisVariants($var_gene_annot_var2) >= 1) {
#					$ok = 1;
#					$patient->{model_vector_var_ok}->{$self->id()} = $self->getNewVector() unless (exists $patient->{model_vector_var_ok}->{$self->id()});
#					$patient->{model_vector_var_ok}->{$self->id()} += $var_gene_annot_var1;
#					$patient->{model_vector_var_ok}->{$self->id()} += $var_gene_annot_var2;
#					$patient->{model_indiv_compound}->{$self->id()} += $var_gene_annot_var1;
#					$patient->{model_indiv_compound}->{$self->id()} += $var_gene_annot_var2;
#					$var_global_patients += $patient->{model_vector_var_ok}->{$self->id()};
#				}
#			}
#			else {
				$ok = 1;
				$patient->{model_vector_var_ok}->{$self->id()} = $self->getNewVector() unless (exists $patient->{model_vector_var_ok}->{$self->id()});
				$patient->{model_vector_var_ok}->{$self->id()} += $var_tmp;
				$patient->{model_indiv_compound}->{$self->id()} += $var_tmp;
				$var_global_patients += $patient->{model_vector_var_ok}->{$self->id()};
#			}
		}
	}
	return $var_global_patients;
}

# applique le model recessif (IND)
sub getModelGeneVector_indiv_recessif {
	my ($self, $gene) = @_;
	my $ok;
	my $var_tmp = $self->getNewVector();
	my $var_tmp_2 = $self->getNewVector();
	my $var_global_patients = $self->getNewVector();
	foreach my $patient (@{$self->getPatients()}) {
		$patient->{model_indiv_recessif}->{$self->id()} = $self->getNewVector() unless (exists $patient->{model_indiv_recessif}->{$self->id()});
		$var_tmp->Empty();
		$patient->used_model('recessif');
		$var_tmp->Intersection( $patient->getVariantsVector($self), $gene->getCurrentVector() );
		next if ($var_tmp->is_empty());
		$patient->{model_vector_var_ok}->{$self->id()} = $self->getNewVector() unless (exists $patient->{model_vector_var_ok}->{$self->id()});
		if ($self->countThisVariants($var_tmp) == 1) {
			$var_tmp_2->Empty();
			$var_tmp_2->Intersection( $var_tmp, $patient->getVectorOriginHo($self, 'ho') );
			unless ($self->countThisVariants($var_tmp_2) == 0) {
				$ok = 1;
				$patient->{model_vector_var_ok}->{$self->id()} += $var_tmp_2;
				$patient->{model_indiv_recessif}->{$self->id()} += $var_tmp_2;
				$var_global_patients += $var_tmp_2;
			}
		} 
		else {
			$ok = 1;
			$patient->{model_vector_var_ok}->{$self->id()} += $var_tmp;
			$patient->{model_indiv_recessif}->{$self->id()} += $var_tmp;
			$var_global_patients += $var_tmp;
		}
	}
	return $var_global_patients;
}

sub get_vector_predicted_splice_region {
	my ($self, $hToFilter) = @_;
	my $vector;
	if ($self->lmdb_score_impact("r")->exists_db){
		$vector = $self->lmdb_score_impact->get("predicted_splice_region");	
	}
	return $vector if defined $vector;
	return $self->getNewVector();
}

sub get_vector_predicted_splice_region_high {
	my ($self, $hToFilter) = @_;
	return $self->get_vector_spliceAI_high();
}

sub get_vector_spliceAI_high {
	my ($self, $hToFilter) = @_;
	my $vector;
	if ($self->lmdb_score_impact->exists_db){
		$vector = $self->lmdb_score_impact->get("spliceAI_high");	
	}
	return $vector if defined $vector;
	return $self->getNewVector();
}

sub get_vector_spliceAI_medium {
	my ($self, $hToFilter) = @_;
	my $vector;
	if ($self->lmdb_score_impact->exists_db){
		$vector = $self->lmdb_score_impact->get("spliceAI_medium");	
	}
	return $vector if defined $vector;
	return $self->getNewVector();
}

sub get_vector_type_variants {
	my ($self, $hToFilter) = @_;
	return unless ($hToFilter);
	my $var = $self->getNewVector();
	foreach my $filter_name (keys %{$self->project->hash_var_type_filters()}) {
		if (exists $hToFilter->{$filter_name}) {
			next unless (exists $self->global_categories->{$filter_name});
			$var += $self->global_categories->{$filter_name};
		}
	}
	return $var;
}

sub get_vector_freq_variants {
	my ($self, $hToFilter) = @_;
	return unless ($hToFilter);
	my $var_freq = $self->getNewVector();
	foreach my $filter_name (keys %{$self->project->hash_frequence_filters()}) {
		next if exists $hToFilter->{$filter_name};
		$var_freq += $self->vector_global_categories($filter_name);
	}
	return $var_freq;
}

sub get_vector_confidence_variants {
	my ($self, $hToFilter) = @_;
	return unless ($hToFilter);
	my $var = $self->getNewVector();
	foreach my $filter_name (keys %{$self->project->hash_confidence_filters()}) {
		if (exists $hToFilter->{$filter_name}) {
			next unless (exists $self->global_categories->{$filter_name});
			$var += $self->global_categories->{$filter_name};
		}
	}
	return $var;
}

sub get_vector_cadd_variants {
	my ($self, $hToFilter) = @_;
	return unless ($hToFilter);
	my $var = $self->getNewVector();
	foreach my $filter_name (keys %{$self->project->hash_cadd_filters()}) {
		if (exists $hToFilter->{$filter_name}) {
			next unless (exists $self->global_categories->{$filter_name});
			$var += $self->global_categories->{$filter_name};
		}
	}
	return $var;
}



has tree_primers =>(
is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift; 
		my $tabix = $self->project->tabix_primers;
		my $res = $tabix->query_full($self->ucsc_name,$self->start,$self->end);
		my $tree = Set::IntervalTree->new;
		return $tree unless $res;
		while(my $line = $res->next){
		my($a,$b,$c,$pid) = split(" ",$line);
		
			$tree->insert($pid,$b,$c);
		}
		return $tree;
		
	}
);


sub set_ordered_primers {
	my $self = shift;
	$self->project->getPrimersByObjects($self);
	return $self->ordered_primers_ids();
}


sub ordered_primers_ids {
		my ($self,$patient) = @_;
		
		return $self->{ordered_primers}  if exists $self->{ordered_primers};
	
	$self->{ordered_primers}  = [];
	my $tabix = $self->project->tabix_primers;
	return unless $tabix;
	my $res = $tabix->query_full($self->ucsc_name,$self->start,$self->end);

	return $self->{ordered_primers} unless $res;
	 while(my $line = $res->next){
		my($a,$b,$c,$pid) = split(" ",$line);
		push (@{$self->{ordered_primers}},$pid);
	}
	
		return $self->{ordered_primers};
}

sub start_flush_primers {
	my ($self) = @_;
	die();
	 $self->{current_list} = [@{$self->ordered_primers_ids}];
	 
}

sub get_natatime{
	my ($self) = @_;
	return $self->{natatime} if exists $self->{natatime};
	$self->{natatime} = [];
	return $self->{natatime} unless $self->ordered_primers_ids;
	my $it = natatime 2000, @{$self->ordered_primers_ids};
		 my @array;
		 while (my @vals = $it->())
  		{
    		push(@array,[@vals]);
  		}
	$self->{natatime} = \@array;
	$self->{current_list} = [];
	return $self->{natatime};
}

sub flush_primer{
	my ($self,$objid) = @_;
	my $no = $self->get_lmdb_cnvs("r");
	return $no->get($objid);
	
}
sub next_primer {
		my ($self,$run) = @_;
		if ($self->{current_primer}){
			delete $self->{current_primer}->{buffer}; delete $self->{current_primer}->{project}; 
			delete $self->project->{objects}->{primers}->{$self->{current_primer}->id};
			$self->{current_primer} = undef;
		}
		
		if  ($self->{current_list} && scalar(@{$self->{current_list}})){
			my $o =shift @{$self->{current_list}};
			bless $o , 'GenBoPrimerCache';
			$o->{project} =  $self->project;
			$o->{buffer} = $self->buffer;
			$self->{primers_object}->{$o->id} = 0;
			$self->project->{objects}->{primers}->{$o->{id}} = $o;
			$self->{current_primer} = $o;
			return $o;
		}
		
		return undef unless scalar ( @{$self->get_natatime});
		my $no = $self->get_lmdb_cnvs("r");
			my $array = shift @{$self->get_natatime};
			while (my $pid = shift(@{$array})) {
				my $o = $no->get($pid);
				if($run){
					next unless exists $o->{runs_object}->{$run->id};
				}
				bless $o , 'GenBoPrimerCache';
				$o->{project} =  $self->project;
				$o->{buffer} = $self->buffer;
				$self->{primers_object}->{$o->id} = 0;
				$self->project->{objects}->{primers}->{$o->{id}} = $o;
				$self->{current_primer} = $o;
				push(@{$self->{current_list}},$o);
			}
		
			return shift @{$self->{current_list}};
    
		

}

sub getPrimers {
	my ($self) = @_;
	return $self->project->getPrimersByObjects($self);
}

#TODO: TEST pour voir ce que ca peut donner. A voir aussi avec de gros projets
has hash_genomic_coord_vector_id => (
        is      => 'rw',
        lazy    => 1,
        default => sub { {} },
);

has variants_intspan_genomic_coord => (
        is      => 'rw',
        lazy    => 1,
        default => sub {
                my $self = shift;
                my $vector = $self->getNewVector();
                $vector->Fill();
                my @lPos;
                foreach my $vector_id (@{$self->getIdsBitOn($vector)}) {
                        my $varId = $self->getVarId($vector_id);
                        my ($chr_id, $pos, $ref, $all) = split('_', $varId);
                        $self->hash_genomic_coord_vector_id->{$pos} = $vector_id;
                        push(@lPos, $pos);
                }
                my $intspan = Set::IntSpan::Fast::XS->new(@lPos);
                return $intspan;
        },
);

sub getVariantsVector_from_coord {
        my ($self, $start, $end) = @_;
        confess("\n\nERROR: need start argument in method GenBoChromosomeCache::getVariantsVector_from_coord(). Die...\n\n") unless ($start);
        confess("\n\nERROR: need end argument in method GenBoChromosomeCache::getVariantsVector_from_coord(). Die...\n\n") unless ($end);
        my $intspan_wanted = Set::IntSpan::Fast::XS->new($start.'-'.$end);
        my $intspan_found = $intspan_wanted->intersection($self->variants_intspan_genomic_coord());
        my @lIds = $intspan_found->as_array();
        my $nb_ids = scalar(@lIds);
        return $self->getNewVector() if ($nb_ids == 0);
        my $start_v_id = $self->hash_genomic_coord_vector_id->{$lIds[0]};
        my $end_v_id = $start_v_id;
        $end_v_id = $self->hash_genomic_coord_vector_id->{$lIds[-1]} if ($nb_ids > 1);
        return Bit::Vector->new_Enum($self->getVariantsVector->Size(), $start_v_id.'-'.$end_v_id);
}


has vectorClinvarPathogenic => (
        is      => 'rw',
        lazy    => 1,
        default => sub {
        	my $self = shift;
        	my $vector ;
        	if ($self->project->isRocks){
        		return $self->getNewVector() if $self->size_vector == 0;
        		return $self->rocks_vector->get_vector_chromosome("dm");
        	}
        	if ($self->lmdb_score_impact->exists_db){
        		
        		$vector = $self->lmdb_score_impact->get("clinvar_pathogenic");	
        		 $vector = $self->getNewVector() unless defined $vector;
        		 return $vector;
        	}
        	confess();
  			#return $vector if defined $vector;
  			#die();
  			#confess("problem ".$self->name);
 			$vector =  $self->getNewVector();
  			foreach my $v (@{$self->getStructuralVariations}){
  				if ($v->isClinvarPathogenic){
  					$vector->Bit_On($v->vector_id); 
  				}
  			}
  			return $vector;
        },
);

has vectorDM => (
        is      => 'rw',
        lazy    => 1,
        default => sub {
        	my $self = shift;
        	if ($self->project->isRocks){
        		return $self->getNewVector() if $self->size_vector == 0;
        		return $self->rocks_vector->get_vector_chromosome("dm");
        	}
        	my $vector ;
        	if ($self->lmdb_score_impact->exists_db){
        		$vector = $self->lmdb_score_impact("r")->get("dm");	
        	}
  			return $vector if defined $vector;
  			confess();
 			$vector =  $self->getNewVector();
  			foreach my $v (@{$self->getStructuralVariations}){
  				if ($v->isDM){
  					$vector->Bit_On($v->vector_id); 
  				}
  			}
  			return $vector;
        },
);



has tree_vector_transcripts => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->lmdb_score_impact("r")->get_interval_tree("intervaltree_transcripts".$self->name);
		#my $tree = $self->project->interval_tree_vector_transcripts->get("transcripts_vector",$self->name);
		#return $self->project->interval_tree_vector_transcripts->get("transcripts_vector",$self->name);
	},
);

sub get_transcripts_from_vector {
	my ($self,$vector) =@_;
	my @t = split(",",$vector->to_Enum);
	my $return ={};

	foreach my $p (@t){
		my ($s,$e) = split("-",$p);
		$e =$s unless $e;
		$e ++;
		foreach my $ttr (@{$self->tree_vector_transcripts->fetch($s,$e)}){
			$return->{$ttr} ++;# if $self->project->lmdbMainTranscripts->exists($ttr);
		}
		
	}
	return [keys %$return];#$return;
}

#validations 
##validation 

sub getVectorLocalValidations {
	my ($self) =@_;
	return $self->{vector_deja_vu_local_validation} if exists $self->{vector_deja_vu_local_validation};
	my $v = $self->project->lite_deja_vu_local_validation()->get($self->project->name."_".$self->name);
	if ($v and $v->Size() eq $self->getNewVector->Size()) {
		$self->{vector_deja_vu_local_validation} = $v;
	}
	else {
		$self->{vector_deja_vu_local_validation} = $self->getVectorLocalValidations_on_live();
	}
	return $self->{vector_deja_vu_local_validation};
}

sub getVectorLocalValidations_on_live {
	my ($self) = @_;
	my $search = $self->name()."_";
	my @z = grep { $_ =~ /!$search/ && $self->getProject->validations->{$_}->[-1]->{idacmg}>=1 } keys %{ $self->getProject->validations() };
	my $vLocal = $self->getNewVector();
	my @lIdsVal;
	foreach my $zid (@z) {
		my ( $g, $vid ) = split( "!", $zid );
		my ($a,$s,$r,$gr) = split("_",$vid);
		my $k = $self->cache_lmdb_variations->get($vid);
		next unless $k;
		push(@lIdsVal, $k->{index_lmdb});
		#$vLocal->Bit_On( $k->{index_lmdb} );
	}
	$vLocal->from_Enum(join(',', @lIdsVal));
	return $vLocal;
}

sub getVectorByIntspan {
	my ($self,$intspan) =@_;
		my $iter = $intspan->iterate_runs();
		my @t;
	my $vector_pos = $self->getNewVector();
	while ( my ( $from, $to ) = $iter->() ) {
	
		$vector_pos += $self->getVectorByPosition($from,$to);
	}
	return $vector_pos;
	
}
sub getVariantIndex {
	my ($self,$vid) =@_;
	my ($a,$start,$r,$gr) = split("_",$vid);
	my $end += $start+1;
	my $nb = $self->getVariantsVector->Size();
	my @ids =(0..$nb-1);
	my $lb = lower_bound {$self->getVariantPosition($_,"start") <=> $start}  @ids; # returns 2
	#my $ub = upper_bound {$self->getVariantPosition($_,"end") <=> $start}  @ids;
	#for  (my $i=$lb;$i<= $ub  ;$i++){
	#	my $no = $self->cache_lmdb_variations();
	#	my $id = $no->get_varid($i);
	#	return $i if $id eq $vid;
	#}
	return undef;
	#warn $lb." ".$ub;
}

sub getVectorByPosition {
	my ($self,$start,$end,$debug) =@_;
	$end ++ if $end == $start;

#	my @ids = (1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6, 7, 7, 7, 8, 8, 9, 9, 9, 9, 9, 11, 13, 13, 13, 17);
#	my $lb = upper_bound { $_ <=> 4 } @ids; # returns 10
#	warn $lb." ".scalar(@ids);
#	 $lb = upper_bound { $_ <=> -1 } @ids; # returns 10
#	 warn $lb." ".scalar(@ids);
#	  $lb = upper_bound { $_ <=> 100 } @ids; # returns 10
#	 warn $lb." ".scalar(@ids);
#	 die();
	my $nb = $self->getVariantsVector->Size();
	return $self->getNewVector() if $nb ==0;
	my $ilimit = $nb -1;

#	warn $start." ".$end if  $end <= 236557770;
	$debug =1 if $end==145878861; 
#	warn $start." ".$end if  $debug;
	#confess() if $debug;
		
	#->get_varid($vid);
	my @ids =(0..$ilimit);
	my $lb = lower_bound {$self->getVariantPosition($_,"start") <=> $start}  @ids; # returns 2
	
	if ($start > $self->getVariantPosition($ilimit,"start")){
		return $self->getNewVector();
	#	warn " *** LB : ".$lb." lb nb:$nb $start-$end ".scalar(@ids);
	#	warn $self->getVariantPosition($ilimit,"start");
	#	die($lb);
	}
	if ($end < $self->getVariantPosition(0,"start") ){
		return $self->getNewVector();
		
	}
#	warn "LB : ".$lb." lb nb:$nb $start-$end ".scalar(@ids) if $debug;
#	warn $ids[$lb];
#	warn "pos ".$self->getVariantPosition($lb-1,"start") if $debug;
#	warn ($start - $self->getVariantPosition($lb-1,"start")) if $debug;
#	warn "pos ".$self->getVariantPosition($lb,"start") if $debug;
	
	confess() if $debug;
	#die("coucou $lb") if $start > $self->getVariantPosition($lb,"start");
	my $ub ;
#	warn $lb." ".$ilimit if $debug;;
	
	unless ($lb == $ilimit){
	 $ub = upper_bound {$self->getVariantPosition($_,"end") <=> $end}  @ids; # returns 2
	 $ub = $ilimit if $ub>$ilimit;
#	 warn "RAW UB: $lb-".$ub." $ilimit" ;#if $debug;
	}
	else {
		$ub = $lb;
		$lb --;
	}
	#warn "UB: $lb-".$ub." " if $debug;
	if ($ub == $lb){
		# on est à la limite 
		$ub = 1;
	#	warn "pos ".$self->getVariantPosition($ub,"start")." ".$end ;
	#	die();
	}
	#warn "UB: $lb-".$ub." $start $end" if $debug;
	#warn $ub." ub" if $debug;
	#die() if $debug;
	
	#die($ub." ".$lb." $nb")  if $lb == $ub;
	#$ub +=2 if $lb == $ub;
#	warn $ub." ub $nb nb" if $debug;
	my $no = $self->cache_lmdb_variations();
	my $id = $no->get_varid($lb);
	my($chr,$pos,$a,$b) = split("_",$id);
	my $vector_pos = $self->getNewVector();
	
	$ub = $nb -1 if $ub > $nb;
	#if $lb>= $nb;
#	warn $self->name;
#	warn $lb."-$ub  $nb "." $start $end";# if $debug;
	die() if $debug;
	if ($lb+1 > $ub -1){
		#warn $lb." $nb ".($lb+1)." $start $end";# if $debug;
		#confess();# if $debug;
		$vector_pos->Interval_Fill($lb,$lb+1);
	}
	else {
		#warn "==>".$lb."-$ub:$nb";
		$vector_pos->Interval_Fill($lb,$ub);
	}
	if (length($a) > length($b)){
		for  (my $i = $lb+1;$i<$ub+1;$i++) {
			my $pos1 = $self->getVariantPosition($i,"start");
			last if ($pos1 >= $start);
			$vector_pos->Bit_Off($i);
		}
	}	
	return $vector_pos;
}
sub getVariantPosition {
	my ($self,$v,$startorend) = @_;
	my $no = $self->cache_lmdb_variations();
	my $id = $no->get_varid($v);
	confess() if $v == -1;
	warn $v unless $id;
	my($chr,$pos,$a,$b) = split("_",$id);
	
	if($startorend eq "start"){
		$pos += length($a)-1;
	}
	return $pos;	
}
#######################
#  methods CNV
######################

has large_duplication_interval_tree => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $tree = Set::IntervalTree->new;
		my $ld = $self->getListVarObjects($self->getVectorLargeDuplications);
		foreach my $v (@{$ld}) {
			my $start = $v->start();
			my $end = $v->start() + $v->length();
			$end++ if ($start == $end);
			$tree->insert([$start, $end, $v->vector_id], $start, $end);
		}
		return $tree;
	},
);

has large_deletion_interval_tree => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $tree = Set::IntervalTree->new;
		my $ld = $self->getListVarObjects($self->getVectorLargeDeletions);
		foreach my $v (@{$ld}) {
			$tree->insert([$v->start, $v->start+$v->length,$v->vector_id], $v->start, $v->start+$v->length);
		}
		return $tree;
	},
);


###################
## Vector
###################

#Patient

sub get_vector_patient {
	my ($self,$patient,$queries) = @_;
	if($self->project->isRocks){
		my $h = {};
		foreach my $cat (@$queries){
			$h->{$cat} = $self->rocks_vector("r")->get_vector_patient($patient,$cat);
			#$h->{$cat} = $self->rocks_vector("r")->get($patient->name."+".$cat);
		}
		return $h;
	}
	else {
		my $h =  $self->get_lmdb_patients("r")->get($patient->name);
		my $z ;
		foreach my $q (@$queries) {
			$z->{$q} = $h->{bitvector}->{$q};
		}
		return $z;
	}
}


#has patients_categories => (
#	is		=> 'rw',
#	lazy    => 1,
#	default => sub {
#		my $self = shift;
#		my $h;
#		#HERE
#		unless ($self->size_vector()) {
#			foreach my $patient (@{$self->getPatients()}) {
#				$h->{$patient->name()} = Bit::Vector->new(0);
#				$h->{$patient->name().'_he'} = Bit::Vector->new(0);
#				$h->{$patient->name().'_ho'} = Bit::Vector->new(0);
#				$h->{$patient->name().'_RI'} = Bit::Vector->new(0);
#				$h->{$patient->name().'_SE'} = Bit::Vector->new(0);
#			}
#			return $h;
#		}
#		my $no_patients = $self->get_lmdb_patients("r");
#		my $patients = $no_patients->get_values();
#		foreach my $patient (@$patients) {
#			my $name = $patient->{'name'};
#			$h->{$name} = $patient->{'bitvector'}->{'all'};
#			$h->{$name.'_he'} = $patient->{'bitvector'}->{'he'};
#			$h->{$name.'_ho'} = $patient->{'bitvector'}->{'ho'};
#			$h->{$name.'_RI'} = $patient->{'bitvector'}->{'RI'} if (exists $patient->{'bitvector'}->{'RI'});
#			$h->{$name.'_SE'} = $patient->{'bitvector'}->{'SE'} if (exists $patient->{'bitvector'}->{'SE'});
#		}
#		$no_patients->close();
#		return $h;
#	},
#);





sub get_vector_category {
	my ( $self, $category ) = @_;
	if($self->project->isRocks){
					if ( $category eq "all" ) {
				my $v = $self->getNewVector();
				$v->Fill();
			return $v;
			}
	 return $self->getNewVector() if $self->size_vector == 0;	
	return $self->rocks_vector("r")->get($category);
	}
	else {
		my $no_categories = $self->get_lmdb_categories("r");
		return $no_categories->get($category)->{bitvector};
	}

}
 sub get_vector_categories {
 	my ( $self,$categories) = @_;
 	if($self->project->isRocks){
 		my $h = {};
 			 return $self->getNewVector() if $self->size_vector == 0;	
		foreach my $cat (@$categories) {
			$h->{$cat} = $self->rocks_vector("r")->get($cat);
		}
		return $h;
	}
	else {
		my $no_categories = $self->get_lmdb_categories("r");
		my $h ={};
		foreach my $cat (@$categories) {
			$h->{$cat} = $no_categories->get($cat)->{bitvector};
		}
		return $h;
	}
 }

1;
