package GenBoGeneCache;
use strict;
use Storable qw(retrieve thaw);
use Moo;
use Bit::Vector;
use Bit::Vector::Overload;
use Data::Dumper;
use GenBoGene;
use Storable qw(store retrieve freeze dclone thaw);
use Carp;
extends 'GenBoGene', 'GenBoCache';


has chromosome => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {}
);

sub getChromosome {
	my $self = shift;
	return $self->chromosome() if (ref($self->chromosome()) eq 'GenBoChromosomeCache');
	return $self->SUPER::getChromosome();
}

# id de ce gene
has id => (
	is		 => 'ro',
	required => 1,
);

# vector id de ce gene

# ceci est la longueur du vector gene
has kyotoId => (
	is		 => 'ro',
	required => 1,
);


## renvoie l'ID ENSG du gene
has variants => (
	is		 => 'ro',
	lazy	 => 1,
	default	 => sub {
		my $self = shift;
		my $hPos;
		unless ($self->categories_intspan()) {
			return $self->getChromosome->getNewVector();
		}
		my $debug;
		$debug =1 if $self->id eq "ENSG00000250979_8";
		foreach my $cat (keys %{$self->categories_intspan()}) {
#			warn $cat ." ".$self->categories_intspan->{$cat}->as_string if $debug;
			foreach my $pos ($self->categories_intspan->{$cat}->as_array()) {
				$hPos->{$pos} = undef;
				my $name = $self->id.";".$cat;
				$self->getChromosome->variations_genes_tree->insert($name, $pos, $pos+1);
			}
		}
		$self->cache_specific_atrributes() unless exists $self->{origin_vector_start};
		die() unless exists $self->{origin_vector_start};
		my $len = abs($self->{origin_vector_start} - $self->{origin_vector_end}) +1;
		$self->start_subpart( $self->{origin_vector_start});
		$self->end_subpart( $self->{origin_vector_end});
		$self->len_subpart($len);
		return   $self->{origin_vector};# Bit::Vector->new_Enum($self->getChromosome->size_vector(), join(',', @lNewCoord));
	}
);


# les vector d'un gene est sous vector de celui du chr. Ceci est sa position start par rapport au vector chr.
has start_subpart => (
	is	 => 'rw',
	lazy 	=> 1,
	default	=> undef,
);

# les vector d'un gene est sous vector de celui du chr. Ceci est sa position end par rapport au vector chr.
has end_subpart => (
	is	 => 'rw',
	lazy 	=> 1,
	default	=> undef,
);

# ceci est la longueur du vector gene
has len_subpart => (
	is		 => 'rw',
	lazy 	=> 1,
	default	=> undef,
);

## hash avec le nom des patients n ayant pas ce gene des la creation de ce GenBoGeneCache. Il n aura donc JAMAIS ce patient (raccourci pour le comptage des genes par patient)
has patients_never_found => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub { {} },
);

has vector_id => (
	is		 => 'ro',
	lazy =>1,
		default	=> sub {
			my $self = shift;
			$self->cache_specific_atrributes();
		},
);

has bitvector => (
	is		 => 'ro',
	lazy =>1,
		default	=> sub {
			my $self = shift;
			$self->cache_specific_atrributes();
			return $self->{bitvector};
		},
);

has intspan => (
	is		 => 'ro',
	lazy =>1,
		default	=> sub {
			my $self = shift;
			confess();
			$self->cache_specific_atrributes();
			return $self->{intspan};
		},
);



sub cache_specific_atrributes{
	my $self = shift;
	confess();
	my $chr = $self->getChromosome();
	my $lmdb = $chr->get_lmdb_genes("r");
	
	$self->{intspan} = {};
	return unless $lmdb;
	my $hashArgs = $lmdb->get($self->id);
	
	#$hashArgs->{id} = $hashArgs->{name};;#$self->id().'_'.$vector_id;
	$self->{vector_id} =  $hashArgs->{index_lmdb};;
	
	
	$self->{kyotoId} = $hashArgs->{name};
	$self->{bitvector} =$hashArgs->{bitvector};
	#$self->{intspan} =$hashArgs->{intspan};
	$self->{intspan} =$hashArgs->{intspan};
	unless (exists $hashArgs->{start}) {
		$self->{origin_vector_start} = 0;
		$self->{origin_vector_end} = 0;
		$self->{origin_vector} =  Bit::Vector->new($chr->size_vector()) ;#Bit::Vector->new_Enum($chr->size_vector(), $hashArgs->{start}."-".$hashArgs->{end});
	}
	else {

		$self->{origin_vector_start} =  $hashArgs->{start};
		$self->{origin_vector_end} =  $hashArgs->{end};
		$self->{origin_vector} =  Bit::Vector->new_Enum($chr->size_vector(), $hashArgs->{start}."-".$hashArgs->{end});

	}
	#$hashArgs->{project} = $self->project();
	#$hashArgs->{chromosome} = $self;
	
}

sub getVectorPatient {
	my ($self,$patient) = shift;
	my $vector = dclone $self->getCurrentVector();
	$vector &= $patient->getVectorOrigin($self->getChromosome);# if $patient;
	return $vector;
}

sub _getVectorOrigin {
	my $self = shift;
	if ($self->project->isRocks){
		return $self->getChromosome->rocks_vector("r")->get_vector_gene($self->id);
	}
	else {
		my $chr = $self->getChromosome();
	 	my $lmdb = $chr->get_lmdb_genes("r");
	 	my $hashArgs = $lmdb->get($self->id);
	 	unless (exists $hashArgs->{start}) {
		$self->{origin_vector_start} = 0;
		$self->{origin_vector_end} = 0;
		return  Bit::Vector->new($chr->size_vector()) ;#Bit::Vector->new_Enum($chr->size_vector(), $hashArgs->{start}."-".$hashArgs->{end});
	}
	else {

		$self->{origin_vector_start} =  $hashArgs->{start};
		$self->{origin_vector_end} =  $hashArgs->{end};
		return   Bit::Vector->new_Enum($chr->size_vector(), $hashArgs->{start}."-".$hashArgs->{end});

	}
	 	
	}
	
	
}
sub getVectorOrigin {
	my $self = shift;
	return $self->{origin_vector} if exists $self->{origin_vector};
	$self->{origin_vector} = $self->_getVectorOrigin();
	return  $self->{origin_vector};
}

has global_categories => (
	is		=> 'rw',
	lazy 	=> 1,
	default => sub {
		my $self = shift;
		my $hash;
		confess();
		return {} unless $self->intspan;
		foreach my $cat (keys %{$self->intspan()}) {
			$hash->{$cat} = $self->convert_intspan_to_vector($self->intspan->{$cat}, $self->getChromosome());
		}
		foreach my $cat (keys %{$self->global_categories_intspan()}) {
			$hash->{$cat} = $self->convert_intspan_to_vector($self->intspan->{$cat}, $self->getChromosome());

		}
		return $hash;
	},
);

has categories => (
	is		=> 'rw',
	lazy 	=> 1,
	default => sub {
		my $self = shift;
		my $hash;
		foreach my $cat (keys %{$self->categories_intspan()}) {
			$hash->{$cat} = $self->convert_intspan_to_vector($self->intspan->{$cat}, $self->getChromosome());
		}
		return $hash;
	},
);

has global_categories_intspan => (
	is		=> 'rw',
	lazy 	=> 1,
	default => sub {
		my $self = shift;
		my $hash;
		#confess();
		foreach my $cat (keys %{$self->intspan()}) {
			next if (exists $self->project->hash_ensembl_annotations->{$cat});
			#confess();
			$hash->{$cat} = $self->intspan->{$cat};
		}
		return $hash;
	},
);

has categories_intspan => (
	is		=> 'rw',
	lazy 	=> 1,
	default => sub {
		my $self = shift;
		my $hash;
		return $hash unless ($self->intspan());
		foreach my $cat (keys %{$self->intspan()}) {
			next unless (exists $self->project->hash_ensembl_annotations->{$cat});
			$hash->{$cat} = $self->intspan->{$cat};
		}
		return $hash;
	},
);

has is_intergenic => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return 1 if (exists $self->categories->{intergenic});
		return;
	},
);

has patients_found => (
	is		=> 'ro',
	lazy    => 1,
	default => undef,
);

has families_found => (
	is		=> 'ro',
	lazy    => 1,
	default => undef,
);



##### METHODS #####




sub setTranscripts {
	my $self = shift;
	my $hTranscriptsId = {};
	my $lTranscriptsId = $self->transcripts();
	foreach my $id (@$lTranscriptsId) {
		unless ($id =~/_/){
			$id.="_".$self->getChromosome->name();
		}
		$hTranscriptsId->{$id} = undef;
	}
	return $hTranscriptsId;
}

# Used for PolyQuery (huge projects)
sub getCategoriesVariantsVector {
	my ($self, $hCat) = @_;
	my $vector = $self->getChromosome->getNewVector();
	return $vector if ($self->getVariantsVector->is_empty());
	foreach my $cat (keys %{$self->categories_intspan()}) {
		next if (exists $hCat->{$cat});
		unless (exists $self->{categories}->{$cat}) {
			$self->{categories}->{$cat} = $self->convert_intspan_to_vector($self->intspan->{$cat}, $self->getChromosome());
		}
		$vector += $self->{categories}->{$cat};
	}
	# FILTRE les variants ayant une prediction Sift / Polyphen que le n on veut pas
	foreach my $cat (%{$self->getChromosome->project->hash_prediction_filters()}) {
		next unless (exists $hCat->{$cat});
		if (exists $self->intspan()->{$cat} and not $self->{global_categories}->{$cat}) {
			$self->{global_categories}->{$cat} = $self->convert_intspan_to_vector($self->intspan->{$cat}, $self->getChromosome());
		}
		elsif (exists $self->global_categories_intspan()->{$cat} and not $self->{global_categories}->{$cat}) {
			$self->{global_categories}->{$cat} = $self->convert_intspan_to_vector($self->global_categories_intspan->{$cat}, $self->getChromosome());
		}
		else {
			next;
		}
		next unless (exists $self->{global_categories}->{$cat});
		$vector -= $self->{global_categories}->{$cat};
	}
	return $vector;
}

sub setPatients {
	my $self = shift;
	my $h;
	foreach my $patient (@{$self->getChromosome->project->getPatients()}) {
		if ($self->getVariantsVector->subset($patient->getVariantsVector($self->getChromosome()))) {
			$h->{$patient->id()} = undef;
		}
	}
	return $h;
}

sub getPatients {
	my $self = shift;
	my $patients = $self->patients_object();
	if ($patients) { 
		return $self->getProject()->myflushobjects($patients, "patients");
	}
	return;
}

sub setFamilies {
	my $self = shift;
	my $h;
	foreach my $pat (@{$self->getPatients()}){
		$h->{$pat->getFamily()->name()} = undef;
	}
	return $h;
}


sub setVariants {
	my ($self, $type) = @_;
	#die();
	my $chr = $self->getChromosome();
	my $vector = $chr->getNewVector();
	if ($type eq 'variations') {
		$vector->Intersection( $self->getVariantsVector(), $chr->{global_categories}->{substitution} ) if (exists $chr->{global_categories}->{substitution});
	}
	elsif ($type eq 'insertions') {
		$vector->Intersection( $self->getVariantsVector(), $chr->{global_categories}->{insertion} ) if (exists $chr->{global_categories}->{insertion});
	}
	elsif ($type eq 'deletions') {
		my $vector_del = $chr->getNewVector();
		$vector_del += $chr->{global_categories}->{deletion} if (exists $chr->{global_categories}->{deletion});
		$vector_del += $chr->{global_categories}->{large_deletion} if (exists $chr->{global_categories}->{large_deletion});
		$vector->Intersection( $self->getVariantsVector(), $vector_del );
	}
	elsif ($type eq 'large_deletions') {
		$vector->Intersection( $self->getVariantsVector(), $chr->{global_categories}->{large_deletion} ) if (exists $chr->{global_categories}->{large_deletion});
	}
	foreach my $var (@{$chr->getListVarObjects($vector)}) {
		$self->{$var->type_object()}->{$var->id()} = undef;
		unless (exists $self->project->{objects}->{$type}->{$var->id()}) {
			$self->project->{objects}->{$type}->{$var->id()} = $var;
		}
	}
}

sub setVariations {
	my $self = shift;
	$self->setVariants('variations');
	return $self->{variations_object} ;
}

sub getVariations {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->variations_object(), "variations");
}

sub setInsertions {
	my $self = shift;
	$self->setVariants('insertions');
	return $self->{insertions_object} ;
}

sub getInsertions {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->insertions_object(), "insertions");
}

sub setDeletions {
	my $self = shift;
	$self->setVariants('deletions');
	return $self->{deletions_object} ;
}

sub getDeletions {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->deletions_object(), "deletions");
}

sub setLargeDeletions {
	my $self = shift;
	$self->setVariants('large_deletions');
	return $self->{large_deletions_object} ;
}

sub getLargeDeletions {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->large_deletions_object(), "deletions");
}

sub getIndels{
	my $self = shift;
	my @lRes;
	push(@lRes,@{$self->getInsertions()});
	push(@lRes,@{$self->getDeletions()});
    return \@lRes;
}

sub getStructuralVariations {
	my $self = shift;
	my @lRes;
	push(@lRes,@{$self->getVariations()});
	push(@lRes,@{$self->getIndels()});
    return \@lRes;
}

sub hasVariantsForAllPatients {
	my ($self, $patients) = @_;
	my $nb_ok = 0;
	$self->getVariantsVector()->Intersection($self->getVariantsVector(), $self->getChromosome->getVariantsVector());
	my $var_tmp = $self->getChromosome->getNewVector();
	foreach my $patient (@$patients) {
		$var_tmp->Empty();
		#warn "\n\n";
		#warn 'BEFORE: '.$self->getChromosome->countThisVariants($self->getVariantsVector());
		$var_tmp->Intersection($self->getVariantsVector(), $patient->getVariantsVector($self->getChromosome()));
		#warn 'AFTER: '.$self->getChromosome->countThisVariants($var_tmp);
		if ($var_tmp->is_empty()) { return; }
	}
	return 1;
}

sub getFilteredVariants {
	my ($self,$patient) = @_;
	my $vector = dclone $self->getVariantsVector();
	 $vector->Intersection($patient->getVariantsVector($self->getChromosome),$vector);
	 
	 return $self->getChromosome->getListVarObjects($vector);
}

sub getCurrentVector {
	my $self = shift;
	unless (exists $self->{current}) {
		$self->{current} = dclone $self->variants();
	}
	
	$self->{current} &= $self->getChromosome->getVariantsVector();
	return $self->{current};
}


sub getVariants {
	my ($self,$patient,$filters) = @_;
	#confess();
	my $vector = dclone $self->getCurrentVector();
	#warn $vector if $self->getChromosome() eq "MT";
	$vector &= $patient->getVectorOrigin($self->getChromosome) if $patient;
	#warn $self->countThisVariants($v3)." ".$self->countThisVariants($v)." ".$self->countThisVariants($v2);
	my $vf;
	my $filter = $filters->{frequency};
	if ($filter){
		foreach my $f (keys %{$self->buffer->config->{frequence_filters}}){
			my $value = $self->buffer->config->{frequence_filters}->{$f};
			if ($value <= $filter){
				next unless  $self->getChromosome->vector_global_categories($f);
				$vf = $self->getChromosome->vector_global_categories($f) unless $vf;
				$vf += $self->getChromosome->vector_global_categories($f);
			}
			
		}
		return [] unless $vf;
		$vector &= $vf;
	}
	return $self->getChromosome->getListVarObjects($vector);
}

has transcripts_tree => (
	is      => 'ro',
	lazy    => 1,
	default => sub { 
	my $self = shift;
	 my $tree = Set::IntervalTree->new;
	 my $transcripts = $self->getTranscripts();
	 foreach my $tr (@$transcripts) {
		$tree->insert($tr->id,$tr->start-50,$tr->end+52);
	 }
	 return $tree;
		
		 }
);

has variants_tree => (
	is      => 'ro',
	lazy    => 1,
	default => sub { 
	my $self = shift;
		my $tree = Set::IntervalTree->new;
		my $vs = $self->getVariants();
		foreach my $v (@$vs) {
			$tree->insert($v->vector_id,$v->start,$v->end+1);
		}
	return $tree;
	}
);

# applique le model compound (FAM)
sub getModelGeneVector_fam_compound {
	my $self = shift;
#	my ($self, $self, $hVectorDenovoEachFam) = @_;
	my $var_res = $self->getChromosome->getNewVector();
	my $var = $self->getVariantsVector->Clone();
	$var->Intersection($var, $self->getChromosome->getVariantsVector());
	
#	#TODO: here;
#	if ($hVectorDenovoEachFam) {
#		warn Dumper $hVectorDenovoEachFam;
#	}
#	die;
	
	return $var_res if ($self->getChromosome->countThisVariants($var) < 2);
	my ($hVar, @keep_var);
	foreach my $family (@{$self->getChromosome->getFamilies}) {
		# CAS fam compound sans parents (= gene avec au moins 2 He pour un patient malade dont les var sont jamais vus chez sains)
		
		if (scalar(@{$family->getParents()}) == 0) {
			my $var_ok = $self->getChromosome->getNewVector();
			my $var_fam = $self->getChromosome->getNewVector();
			foreach my $patient (@{$family->getChildrenIll()}) {
				$patient->getVariantsVector($self->getChromosome);
				$var_fam += $patient->getHe($self->getChromosome);
			}
			foreach my $patient (@{$family->getChildrenHealthy()}) {
				$var_fam -= $patient->getVariantsVector($self->getChromosome);
			}
			$var_fam->Intersection($var, $var_fam);
			if ($self->getChromosome->countThisVariants($var_fam) >= 2) {
				foreach my $patient (@{$family->getChildrenIll()}) {
					$patient->getVariantsVector($self->getChromosome);
					my $var_pat = $self->getChromosome->getNewVector();
					$var_pat += $patient->getHe($self->getChromosome);
					$var_pat->Intersection($var_pat, $var_fam);
					if ($self->getChromosome->countThisVariants($var_pat) >= 2) {
						$var_ok += $var_pat;
					}
				}
				foreach my $v (@{$self->getChromosome->getListVarVectorIds( $var_ok )}) {
					push(@{$family->keep_var_compound->{$self->getChromosome->id()}}, $v);
				}
				$var_res += $var_ok;
			}
			$family->used_model('compound no parent');
			next;
		}
		
		# get common variants from all ill child 
		my $firs_pat_ill_done;
		my $var_pat_ill = $self->getChromosome->getNewVector();
		foreach my $patient (@{$family->getChildrenIll()}) {
			my $var_pat = $patient->getHe($self->getChromosome)->Clone();
			$var_pat->Intersection( $var, $var_pat );
			if ($firs_pat_ill_done) { $var_pat_ill->Intersection($var_pat_ill, $var_pat); }
			else {
				$var_pat_ill = $var_pat;
				$firs_pat_ill_done = 1;
			}
		}
		
		# delete variant seen HO in all ill child
		foreach my $patient (@{$family->getChildrenIll()}) {
			$var_pat_ill -= $patient->getHo($self->getChromosome);
		}
		
		# delete variants seen HO in each parent
		foreach my $patient (@{$family->getParents()}) {
			$var_pat_ill -= $patient->getHo($self->getChromosome);
		}
		
		#next if ($self->getChromosome->countThisVariants($var_pat_ill) < 2);
		
		# get variants seen HE in each parent and seen in ill children too
		my @lVarParents;
		my $var_all_parents = $self->getChromosome->getNewVector();
		foreach my $patient (@{$family->getParents()}) {
			my $var_this_parent = $patient->getHe($self->getChromosome)->Clone();
			$var_this_parent->Intersection($var_this_parent, $var_pat_ill);
			push(@lVarParents, $var_this_parent->Clone());
			$var_all_parents += $var_this_parent;
		}

		# delete variants seen HO in each parent
		my $nb_parents = 0;
		foreach my $patient (@{$family->getParents()}) {
			$var_all_parents -= $patient->getHo($self->getChromosome);
			$nb_parents++;
		}
		next unless ($nb_parents == 2);

		# delete common variants in parents
#		if ($nb_parents == 2) {
			my $var_both_parents = $self->getChromosome->getNewVector();
			$var_both_parents->Intersection($lVarParents[0], $lVarParents[1]);
			$lVarParents[0] -= $var_both_parents;
			$lVarParents[1] -= $var_both_parents;
			$var_all_parents -= $var_both_parents;
#		}
#		elsif ($nb_parents == 1) {
#			my $vecTmp = $self->getVariantsVector->Clone();
#			$vecTmp -= $family->getParents->[0]->getVariantsVector($self->getChromosome);
#			push(@lVarParents, $vecTmp);
#		}
		next if ($lVarParents[0]->is_empty());
		next if ($lVarParents[1]->is_empty());
		
		# combinaison var
		my $nbc = 0;
		my $combinaisons = {};
		for my $v1 (@{$self->getChromosome->getListVarVectorIds($lVarParents[0])}) {
			for my $v2 (@{$self->getChromosome->getListVarVectorIds($lVarParents[1])}) {
				$combinaisons->{$nbc} = $self->getChromosome->getNewVector();
				$combinaisons->{$nbc}->Bit_On($v1);
				$combinaisons->{$nbc}->Bit_On($v2);
				$nbc++;
			}
		}
		
		# cas ou on considere un parent et un denovo
#		if ($hVectorDenovoEachFam) {
#			for my $v1 (@{$self->getChromosome->getListVarVectorIds($lVarParents[0])}) {
#				for my $v2 (@{$self->getChromosome->getListVarVectorIds($hVectorDenovoEachFam->{$family->name()})}) {
#					$combinaisons->{$nbc} = $self->getChromosome->getNewVector();
#					$combinaisons->{$nbc}->Bit_On($v1);
#					$combinaisons->{$nbc}->Bit_On($v2);
#					$nbc++;
#				}
#			}
#			for my $v1 (@{$self->getChromosome->getListVarVectorIds($lVarParents[1])}) {
#				for my $v2 (@{$self->getChromosome->getListVarVectorIds($hVectorDenovoEachFam->{$family->name()})}) {
#					$combinaisons->{$nbc} = $self->getChromosome->getNewVector();
#					$combinaisons->{$nbc}->Bit_On($v1);
#					$combinaisons->{$nbc}->Bit_On($v2);
#					$nbc++;
#				}
#			}
#		}
		# TODO: on obtient un cas tordu avec des couples present aussi chez les enfants sains...
		
		# delete combinaison seen in all healthy child
		foreach my $patient (@{$family->getChildrenHealthy()}) {
			foreach my $nbc (keys %$combinaisons) {
				my $v_comb = $combinaisons->{$nbc}->Clone();
				$v_comb->Intersection($v_comb, $patient->getVariantsVector($self->getChromosome));
				if ($self->getChromosome->countThisVariants($v_comb) == 2) {
					delete $combinaisons->{$nbc};
				}
			}
		}
		
		foreach my $nbc (keys %$combinaisons) {
			my $v_comb = $combinaisons->{$nbc}->Clone();
			$v_comb->Intersection($v_comb, $self->getVariantsVector());
			$v_comb->Intersection($v_comb, $family->getVariantsVector($self->getChromosome));
			$v_comb->Intersection($v_comb, $self->getChromosome->getVariantsVector());
			$var_res += $combinaisons->{$nbc};
			foreach my $v (@{$self->getChromosome->getListVarVectorIds($combinaisons->{$nbc})}) {
				push(@{$family->keep_var_compound->{$self->getChromosome->id()}}, $v);
			}
		}
	}
	return $var_res;
}

1;
