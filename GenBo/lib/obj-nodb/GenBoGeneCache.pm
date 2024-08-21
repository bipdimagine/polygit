package GenBoGeneCache;
use strict;
use Storable qw(retrieve thaw);
use Moo;
use Bit::Vector;
use Bit::Vector::Overload;
use Data::Dumper;
use Storable qw(store retrieve freeze dclone thaw);
use Carp;
extends 'GenBoGene', 'GenBoCache';

has chromosome => (
	is      => 'rw',
	lazy    => 1,
	default => sub { }
);

sub getChromosome {
	my $self = shift;
	return $self->chromosome()
	  if ( ref( $self->chromosome() ) eq 'GenBoChromosomeCache' );
	return $self->SUPER::getChromosome();
}

# id de ce gene
has id => (
	is       => 'ro',
	required => 1,
);
has hash_annotation => (
	is		=> 'ro',
	lazy => 1,
	default => sub { 
		my $self = shift;
		my $hash  = $self->getChromosome->rocks_vector("r")->get_vector_gene( $self->id."_annotations");
	 },
);

has hgmd => (
	is		=> 'ro',
	lazy => 1,
	default => sub { 
		my $self = shift;
		
		return $self->hash_annotation->{hgmd};
	 }
);
has polyquery_phenotypes => (
	is		=> 'ro',
	lazy => 1,
	default => sub { 
		my $self = shift;
		
		return $self->hash_annotation->{polyquery_phenotypes};
	 }
);


# vector id de ce gene

# ceci est la longueur du vector gene
has kyotoId => (
	is       => 'ro',
	required => 1,
);

## renvoie l'ID ENSG du gene
has variants => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->getVectorOrigin();
	}
);
# les vector d'un gene est sous vector de celui du chr. Ceci est sa position start par rapport au vector chr.
has start_subpart => (
	is      => 'rw',
	lazy    => 1,
	default => undef,
);

# les vector d'un gene est sous vector de celui du chr. Ceci est sa position end par rapport au vector chr.
has end_subpart => (
	is      => 'rw',
	lazy    => 1,
	default => undef,
);

# ceci est la longueur du vector gene
has len_subpart => (
	is      => 'rw',
	lazy    => 1,
	default => undef,
);

## hash avec le nom des patients n ayant pas ce gene des la creation de ce GenBoGeneCache. Il n aura donc JAMAIS ce patient (raccourci pour le comptage des genes par patient)
has patients_never_found => (
	is      => 'rw',
	lazy    => 1,
	default => sub { {} },
);

has vector_id => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		$self->cache_specific_atrributes();
	},
);

has bitvector => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		$self->cache_specific_atrributes();
		return $self->{bitvector};
	},
);

has intspan => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		confess();
		$self->cache_specific_atrributes();
		return $self->{intspan};
	},
);

has enum => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->getChromosome->rocks_vector("r")->get_vector_gene( $self->id."_enum" );
	},
);

sub cache_specific_atrributes {
	my $self = shift;
	return;
	my $chr  = $self->getChromosome();
	
	confess();
	
	my $lmdb = $chr->get_lmdb_genes("r");

	$self->{intspan} = {};
	return unless $lmdb;
	my $hashArgs = $lmdb->get( $self->id );

	#$hashArgs->{id} = $hashArgs->{name};;#$self->id().'_'.$vector_id;
	$self->{vector_id} = $hashArgs->{index_lmdb};

	$self->{kyotoId}   = $hashArgs->{name};
	$self->{bitvector} = $hashArgs->{bitvector};

	#$self->{intspan} =$hashArgs->{intspan};
	$self->{intspan} = $hashArgs->{intspan};
	unless ( exists $hashArgs->{start} ) {
		$self->{origin_vector_start} = 0;
		$self->{origin_vector_end}   = 0;
		$self->{origin_vector}       = Bit::Vector->new( $chr->size_vector() )
		  ; #Bit::Vector->new_Enum($chr->size_vector(), $hashArgs->{start}."-".$hashArgs->{end});
	}
	else {

		$self->{origin_vector_start} = $hashArgs->{start};
		$self->{origin_vector_end}   = $hashArgs->{end};
		$self->{origin_vector} = Bit::Vector->new_Enum( $chr->size_vector(),
			$hashArgs->{start} . "-" . $hashArgs->{end} );

	}

	#$hashArgs->{project} = $self->project();
	#$hashArgs->{chromosome} = $self;
}


sub getVectorPatient {
	my ( $self, $patient ) = shift;
	my $vector = dclone $self->getCurrentVector();
	$vector &= $patient->getVectorOrigin( $self->getChromosome ); # if $patient;
	return $vector;
}





sub init_gene_vector {
	my ($self) = @_;
	my $numbers = $self->getChromosome->rocks_vector("r")->get_vector_gene( $self->id."_vector_characteristic" );
	#my $enum = $self->enum->{all};
		
	#my @numbers = $enum =~ /(\d+)/g;
	$self->{compact_vector_start} = $numbers->[0];
	$self->{compact_vector_length} = $numbers->[1];
		
}

has compact_vector_start => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		$self->init_gene_vector;
		return $self->{compact_vector_start};
	},
);


has compact_vector_length => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		$self->init_gene_vector;
		return $self->{compact_vector_length};
	},
);
has chromosome_vector_length  => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->getChromosome->getNewVector->Size;
	},
);

sub enlarge_compact_vector{
	my ($self,$small) = @_;
	my  $vector = $self->getChromosome->getNewVector;
	$vector->Interval_Substitute($small,$self->compact_vector_start,$self->compact_vector_length,0,$self->compact_vector_length);
	return $vector;
	
}

sub return_compact_vector {
	my ($self,$vector) = @_;
	 my $vsmall = Bit::Vector->new($self->compact_vector_length);
	 
	$vsmall->Interval_Copy($vector,0,$self->compact_vector_start,$self->compact_vector_length);
	return $vsmall;
}

has compact_vector => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $t  =  $self->getChromosome->rocks_vector("r")->get_vector_gene( $self->id."_compact_vector" );
		return $t;
	},
);





sub getNewCompactVector {
my ($self) = @_;
 return  Bit::Vector->new($self->compact_vector_length); 
}

sub getCompactVectorPatient{
	my ($self,$patient) = @_;
	return $self->getCompactVectorOriginCategory($patient->id."_all");
}

sub getVectorOriginCategory {
	my ($self, $cat) = @_; 
	return $self->enlarge_compact_vector($self->getCompactVectorOriginCategory($cat));
}
sub getCompactVectorOriginCategory {
	my ($self, $cat) = @_;
	return $self->compact_vector->{$cat} if exists $self->compact_vector->{$cat};
	
	$self->compact_vector->{$cat} = $self->getNewCompactVector();
	return $self->compact_vector->{$cat};
}

sub getCompactVector {
		my ($self,$cat) = @_;
		return $self->compact_vector->{all} unless $cat;
		return $self->getCompactVectorOriginCategory($cat);
}

sub getVectorOrigin {
		my ($self) = @_;
		return $self->enlarge_compact_vector($self->getCompactVector);
}

sub getCompactVectorOriginCategories {
		my ($self, $cats,$debug) = @_;
		my $small =  $self->getNewCompactVector();
		foreach my $cat (@$cats) {
 			$small += $self->getCompactVectorOriginCategory($cat);
 			warn "\t\t\t ".$cat." ".$small->Norm if $debug;
	}
	return $small;
}


sub getVectorOriginCategories {
	my ($self, $cats) = @_;
	return $self->enlarge_compact_filter($self->getCompactVectorOriginCategories($cats));
}

has is_intergenic => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $v_intergenic = $self->getCompactVectorOriginCategory('intergenic');
		return if $v_intergenic->is_empty();
		return 1;
	},
);

has patients_found => (
	is      => 'ro',
	lazy    => 1,
	default => undef,
);

has families_found => (
	is      => 'ro',
	lazy    => 1,
	default => undef,
);

##### METHODS #####

sub setTranscripts {
	my $self           = shift;
	my $hTranscriptsId = {};
	my $lTranscriptsId = $self->transcripts();
	foreach my $id (@$lTranscriptsId) {
		unless ( $id =~ /_/ ) {
			$id .= "_" . $self->getChromosome->name();
		}
		$hTranscriptsId->{$id} = undef;
	}
	return $hTranscriptsId;
}



sub setPatients {
	my $self = shift;
	my $h;
	foreach my $patient ( @{ $self->project->getPatients() } ) {
		next if $self->getCompactVectorPatient($patient)->is_empty;
			$h->{ $patient->id() } = undef;
	}
	return $h;
}

sub getPatients {
	my $self     = shift;
	my $patients = $self->patients_object();
	if ($patients) {
		return $self->getProject()->myflushobjects( $patients, "patients" );
	}
	return;
}

sub setFamilies {
	my $self = shift;
	my $h;
	foreach my $pat ( @{ $self->getPatients() } ) {
		$h->{ $pat->getFamily()->name() } = undef;
	}
	return $h;
}

sub setVariants {
	my ( $self, $type ) = @_;
	confess();
	my $chr    = $self->getChromosome();
	my $vsmall = $self->getCompactVectorOriginCategory($type);
	my $vector = $self->enlarge_comapct_vector($vsmall);
	foreach my $var ( @{ $chr->getListVarObjects($vector) } ) {
		$self->{ $var->type_object() }->{ $var->id() } = undef;
		unless ( exists $self->project->{objects}->{$type}->{ $var->id() } ) {
			$self->project->{objects}->{$type}->{ $var->id() } = $var;
		}
	}
}

sub setVariations {
	my $self = shift;
	$self->setVariants('variations');
	return $self->{variations_object};
}

sub getVariations {
	my $self = shift;
	return $self->getProject()
	  ->myflushobjects( $self->variations_object(), "variations" );
}

sub setInsertions {
	my $self = shift;
	$self->setVariants('insertions');
	return $self->{insertions_object};
}

sub getInsertions {
	my $self = shift;
	return $self->getProject()
	  ->myflushobjects( $self->insertions_object(), "insertions" );
}

sub setDeletions {
	my $self = shift;
	$self->setVariants('deletions');
	return $self->{deletions_object};
}

sub getDeletions {
	my $self = shift;
	return $self->getProject()
	  ->myflushobjects( $self->deletions_object(), "deletions" );
}

sub setLargeDeletions {
	my $self = shift;
	$self->setVariants('large_deletions');
	return $self->{large_deletions_object};
}

sub getLargeDeletions {
	my $self = shift;
	return $self->getProject()->myflushobjects( $self->large_deletions_object(), "deletions" );
}

sub getIndels {
	my $self = shift;
	my @lRes;
	push( @lRes, @{ $self->getInsertions() } );
	push( @lRes, @{ $self->getDeletions() } );
	return \@lRes;
}

sub getStructuralVariations {
	my $self = shift;
	my @lRes;
	push( @lRes, @{ $self->getVariations() } );
	push( @lRes, @{ $self->getIndels() } );
	return \@lRes;
}


#sub getFilteredVariants {
#	my ( $self, $patient ) = @_;
#	my $vector = dclone $self->getVariantsVector();
#	$vector->Intersection( $patient->getVariantsVector( $self->getChromosome ),
#		$vector );
#
#	return $self->getChromosome->getListVarObjects($vector);
#}


sub getCurrentCompactVector {
	my $self = shift;
	unless ( exists $self->{current_compact} ) {
		return $self->compact_vector->{all};
	}
	
	return $self->{current_compact};
}

sub setCurrentVector {
	my ($self,$vector) = @_;
	confess() unless defined  $vector;
	my $size =  $vector->Size;
	if ($size == $self->compact_vector_length){
		$self->{current_compact}  = $vector;
	}
	elsif ($size == $self->chromosome_vector_length) {
		$self->{current_compact}  = $self->return_compact_vector($vector);
	}
	else{
		confess($size." ". $self->compact_vector_length." ". $self->chromosome_vector_length);
	}	
	return 1;
}

sub getCurrentVector {
	my $self = shift;
	my $vc = $self->getCurrentCompactVector();
	return $self->enlarge_compact_vector($vc);
}



sub getVariants {
	my ( $self, $patient, $filters ) = @_;
	confess() if $filters;
	
	
	my $vector =  $self->getCurrentCompactVector() & $self->getCompactVectorPatient($patient);
		return $self->getChromosome->getListVarObjects($vector);

#	#warn $vector if $self->getChromosome() eq "MT";
#	
#
##warn $self->countThisVariants($v3)." ".$self->countThisVariants($v)." ".$self->countThisVariants($v2);
#	my $vf;
#	my $filter = $filters->{frequency};
#	if ($filter) {
#		foreach my $f ( keys %{ $self->buffer->config->{frequence_filters} } ) {
#			my $value = $self->buffer->config->{frequence_filters}->{$f};
#			if ( $value <= $filter ) {
#				next unless $self->getChromosome->vector_global_categories($f);
#				$vf = $self->getChromosome->vector_global_categories($f)
#				  unless $vf;
#				$vf += $self->getChromosome->vector_global_categories($f);
#			}
#
#		}
#		return [] unless $vf;
#		$vector &= $vf;
#	}

}

has transcripts_tree => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self        = shift;
		my $tree        = Set::IntervalTree->new;
		my $transcripts = $self->getTranscripts();
		foreach my $tr (@$transcripts) {
			$tree->insert( $tr->id, $tr->start - 50, $tr->end + 52 );
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
		my $vs   = $self->getVariants();
		foreach my $v (@$vs) {
			$tree->insert( $v->vector_id, $v->start, $v->end + 1 );
		}
		return $tree;
	}
);

# applique le model compound (FAM)
sub getModelGeneVector_fam_compound {
	my $self = shift;

	#	my ($self, $self, $hVectorDenovoEachFam) = @_;
	my $var_res = $self->getChromosome->getNewVector();
	my $var     = $self->getCurrentVector->Clone();
	$var->Intersection( $var, $self->getChromosome->getVariantsVector() );

	#	#TODO: here;
	#	if ($hVectorDenovoEachFam) {
	#		warn Dumper $hVectorDenovoEachFam;
	#	}
	#	die;

	return $var_res if ( $self->getChromosome->countThisVariants($var) < 2 );
	my ( $hVar, @keep_var );
	foreach my $family ( @{ $self->getChromosome->getFamilies } ) {

# CAS fam compound sans parents (= gene avec au moins 2 He pour un patient malade dont les var sont jamais vus chez sains)

		if ( scalar( @{ $family->getParents() } ) == 0 ) {
			my $var_ok  = $self->getChromosome->getNewVector();
			my $var_fam = $self->getChromosome->getNewVector();
			foreach my $patient ( @{ $family->getChildrenIll() } ) {
				$patient->getVariantsVector( $self->getChromosome );
				$var_fam += $patient->getHe( $self->getChromosome );
			}
			foreach my $patient ( @{ $family->getChildrenHealthy() } ) {
				$var_fam -= $patient->getVariantsVector( $self->getChromosome );
			}
			$var_fam->Intersection( $var, $var_fam );
			if ( $self->getChromosome->countThisVariants($var_fam) >= 2 ) {
				foreach my $patient ( @{ $family->getChildrenIll() } ) {
					$patient->getVariantsVector( $self->getChromosome );
					my $var_pat = $self->getChromosome->getNewVector();
					$var_pat += $patient->getHe( $self->getChromosome );
					$var_pat->Intersection( $var_pat, $var_fam );
					if (
						$self->getChromosome->countThisVariants($var_pat) >= 2 )
					{
						$var_ok += $var_pat;
					}
				}
				foreach my $v (
					@{ $self->getChromosome->getListVarVectorIds($var_ok) } )
				{
					push(
						@{
							$family->keep_var_compound
							  ->{ $self->getChromosome->id() }
						},
						$v
					);
				}
				$var_res += $var_ok;
			}
			$family->used_model('compound no parent');
			next;
		}

		# get common variants from all ill child
		my $firs_pat_ill_done;
		my $var_pat_ill = $self->getChromosome->getNewVector();
		foreach my $patient ( @{ $family->getChildrenIll() } ) {
			my $var_pat = $patient->getHe( $self->getChromosome )->Clone();
			$var_pat->Intersection( $var, $var_pat );
			if ($firs_pat_ill_done) {
				$var_pat_ill->Intersection( $var_pat_ill, $var_pat );
			}
			else {
				$var_pat_ill       = $var_pat;
				$firs_pat_ill_done = 1;
			}
		}

		# delete variant seen HO in all ill child
		foreach my $patient ( @{ $family->getChildrenIll() } ) {
			$var_pat_ill -= $patient->getHo( $self->getChromosome );
		}

		# delete variants seen HO in each parent
		foreach my $patient ( @{ $family->getParents() } ) {
			$var_pat_ill -= $patient->getHo( $self->getChromosome );
		}

		#next if ($self->getChromosome->countThisVariants($var_pat_ill) < 2);

		# get variants seen HE in each parent and seen in ill children too
		my @lVarParents;
		my $var_all_parents = $self->getChromosome->getNewVector();
		foreach my $patient ( @{ $family->getParents() } ) {
			my $var_this_parent =
			  $patient->getHe( $self->getChromosome )->Clone();
			$var_this_parent->Intersection( $var_this_parent, $var_pat_ill );
			push( @lVarParents, $var_this_parent->Clone() );
			$var_all_parents += $var_this_parent;
		}

		# delete variants seen HO in each parent
		my $nb_parents = 0;
		foreach my $patient ( @{ $family->getParents() } ) {
			$var_all_parents -= $patient->getHo( $self->getChromosome );
			$nb_parents++;
		}
		next unless ( $nb_parents == 2 );

		# delete common variants in parents
		#		if ($nb_parents == 2) {
		my $var_both_parents = $self->getChromosome->getNewVector();
		$var_both_parents->Intersection( $lVarParents[0], $lVarParents[1] );
		$lVarParents[0]  -= $var_both_parents;
		$lVarParents[1]  -= $var_both_parents;
		$var_all_parents -= $var_both_parents;

#		}
#		elsif ($nb_parents == 1) {
#			my $vecTmp = $self->getVariantsVector->Clone();
#			$vecTmp -= $family->getParents->[0]->getVariantsVector($self->getChromosome);
#			push(@lVarParents, $vecTmp);
#		}
		next if ( $lVarParents[0]->is_empty() );
		next if ( $lVarParents[1]->is_empty() );

		# combinaison var
		my $nbc          = 0;
		my $combinaisons = {};
		for my $v1 (
			@{ $self->getChromosome->getListVarVectorIds( $lVarParents[0] ) } )
		{
			for my $v2 (
				@{
					$self->getChromosome->getListVarVectorIds(
						$lVarParents[1]
					)
				}
			  )
			{
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
		foreach my $patient ( @{ $family->getChildrenHealthy() } ) {
			foreach my $nbc ( keys %$combinaisons ) {
				my $v_comb = $combinaisons->{$nbc}->Clone();
				$v_comb->Intersection( $v_comb,
					$patient->getVariantsVector( $self->getChromosome ) );
				if ( $self->getChromosome->countThisVariants($v_comb) == 2 ) {
					delete $combinaisons->{$nbc};
				}
			}
		}

		foreach my $nbc ( keys %$combinaisons ) {
			my $v_comb = $combinaisons->{$nbc}->Clone();
			$v_comb->Intersection( $v_comb, $self->getVariantsVector() );
			$v_comb->Intersection( $v_comb,
				$family->getVariantsVector( $self->getChromosome ) );
			$v_comb->Intersection( $v_comb,
				$self->getChromosome->getVariantsVector() );
			$var_res += $combinaisons->{$nbc};
			foreach my $v (
				@{
					$self->getChromosome->getListVarVectorIds(
						$combinaisons->{$nbc}
					)
				}
			  )
			{
				push(
					@{
						$family->keep_var_compound->{ $self->getChromosome->id()
						}
					},
					$v
				);
			}
		}
	}
	return $var_res;
}

1;
