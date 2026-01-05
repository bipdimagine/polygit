package GenBoSomaticGroupCache;
use strict;
use Storable qw(thaw retrieve);
use Moo;
use Bit::Vector;
use Bit::Vector::Overload;
use Bio::DB::Sam;
use Data::Dumper;
use Carp;
extends 'GenBoSomaticGroup', 'GenBoCache';



# hash avec le details du fichier somatic	
has somatic_details => (
	is		=> 'ro',
	required => 1,
);

# hash avec les infos du patients
has patients => (
	is		=> 'ro',
	required => 1,
);

# hash des avec les vector de chaque annotation
has categories => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub { {}	},
);

# tag si group in the attic
has in_the_attic => (
	is		=> 'rw',
	lazy	=> 1,
	default => undef,
);

# tag si group exclue (all, ho, he)
has excluded => (
	is		=> 'rw',
	lazy	=> 1,
	default => 0,
);

# tag si group intersect
has intersected => (
	is		=> 'rw',
	default => 0,
);

# nom du model utilise sur group
has used_model => (
	is		=> 'rw',
	lazy	=> 1,
	default => undef,
);


has list_stats_patients => (
	is		=> 'rw',
	default => sub { [] },
);



##### METHODS #####





sub getVariantsVector {
	my ($self, $chr_obj) = @_;
	confess("\n\nERROR: GenBoFamilyCache->getVariantsVector() method need a GenBoChromosomeCache object in argument. Die.\n\n") unless($chr_obj);
	my $var = $chr_obj->getNewVector();
	foreach my $patient (@{$self->getPatients()}) {
		$var += $patient->getVariantsVector($chr_obj);
	}
	$var->Intersection( $var, $chr_obj->getVariantsVector() );
	$self->{variants}->{$chr_obj->id()} = $var;
	return $self->{variants}->{$chr_obj->id()};
}

# for polyquery
sub setCurrentVariantsVector {
	my ($self, $chr_obj, $vector) = @_;
	$self->{current_variants}->{$chr_obj->id()} = $vector;
}

# for polyquery
sub addCurrentVariantsVector {
	my ($self, $chr_obj, $vector) = @_;
	$self->{current_variants}->{$chr_obj->id()} += $vector;
}

# for polyquery
sub getCurrentVariantsVector {
	my ($self, $chr_obj) = @_;
	if (not exists $self->{current_variants} or not exists $self->{current_variants}->{$chr_obj->id()}) {
		$self->{current_variants}->{$chr_obj->id()} = $self->getVariantsVector($chr_obj);
	}
	return $self->{current_variants}->{$chr_obj->id()};
}

sub countVariantsByType {
	my ($self, $type) = @_;
	my $var = $self->getNewVector();
	foreach my $patient (@{$self->getPatients()}) {
		my $category;
		$category = 'categories'  if (exists $patient->categories->{$type});
		$category = 'global_categories' if (exists $patient->global_categories->{$type});
		next unless ($category);
		my $var_tmp = $self->getNewVector();
		$var_tmp->Intersection( $patient->getVariantsVector(), $patient->$category->{$type} );
		$var += $var_tmp; 
	}
	return $self->countThisVariants( $var );
}

sub countAllVariations {
	my $self = shift;
	my $nb = 0;
	foreach my $patient (@{$self->getPatients()}) {
		$nb += $patient->countVariants();
	}
	return $nb;
} 

sub countTypeVariants {
	my ($self, $type) = @_;
	return $self->countThisVariants( $self->getTypeVariants($type) );
}

sub getTypeVariants {
	my ($self, $type) = @_;
	return $self->getHe() if ($type eq 'he');
	my $var = $self->getNewVector();
	foreach my $patient (@{$self->getPatients()}) {
		my $category;
		$category = 'categories'  if (exists $patient->categories->{$type});
		$category = 'global_categories' if (exists $patient->global_categories->{$type});
		next unless ($category);
		$var += $patient->$category->{$type};
	}
	$var->Intersection( $var, $self->getVariantsVector() );
	return $var;
}

sub getSubstitutions {
	my $self = shift;
	return $self->getTypeVariants('substitution');
}

sub getDeletions {
	my $self = shift;
	return $self->getTypeVariants('deletion');
}

sub getInsertions {
	my $self = shift;
	return $self->getTypeVariants('insertion');
}

sub getHo {
	my $self = shift;
	return $self->getTypeVariants('ho');
}

sub getHe {
	my $self = shift;
	my $var_ho = $self->getHo();
	my $var_he = $self->getNewVector();
	$var_he = $self->getVariantsVector() - $var_ho;
	return $var_he;
}


sub getVector_individual_loh {
	my ($self, $chr,$child,$compute) = @_;
	my $key = "som_loh_".$child->name;
	return $self->{vector_transmission}->{$key}->{$chr->id} if exists $self->{vector_transmission}->{$key}->{$chr->id};
	if ($self->project->isRocks && (! defined $compute)){
		$self->{vector_transmission}->{$key}->{$chr->id} = $chr->rocks_vector->get_vector_transmission($child,"som_loh");
		return $self->{vector_transmission}->{$key}->{$chr->id};
	}
	confess();
}

# methode pour intersecter les patients en mode SOMATIC
sub setIntersectPatients {
	my ($self, $patients, $chr) = @_;
	my $var = $self->getNewVector();
	unless (exists $self->{variants_intersected}->{$chr->id()}) {
		$self->{variants_intersected}->{$chr->id()} = $self->getNewVector();
	}
	foreach my $patient (@$patients) {
		$patient->intersected(1);
		$self->{variants_intersected}->{$chr->id()} += $patient->getVariantsVector($chr);
	}
	foreach my $patient (@$patients) {
		$self->{variants_intersected}->{$chr->id()}->Intersection($patient->getVariantsVector($chr), $self->{variants_intersected}->{$chr->id()});
	}
	$self->{variants_intersected}->{$chr->id()} -= $self->variants_excluded->{$chr->id()} if (exists $self->variants_excluded->{$chr->id()});
	foreach my $patient (@{$self->getPatients()}) {
		$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $self->{variants_intersected}->{$chr->id()});
		$var += $patient->getVariantsVector($chr);
	}
	return $var;
}

# methode pour exclure les patients en mode SOMATIC
sub setExcludePatients {
	my ($self, $lPatients, $status, $chr) = @_;
	my $var_excl = $chr->getNewVector();
	foreach my $patient (@$lPatients) {
		next if ($patient->in_the_attic());
		if ($status eq 'all') {
			$var_excl += $patient->getVariantsVector($chr);
			$patient->variants_excluded->{$chr->id()} = $patient->getVariantsVector($chr);
		}
		elsif ($status eq 'he')  {
			$var_excl += $patient->getHe($chr);
			$patient->variants_excluded->{$chr->id()} = $patient->getHe($chr);
		}
		elsif ($status eq 'ho')  {
			$var_excl += $patient->getHo($chr);
			$patient->variants_excluded->{$chr->id()} = $patient->getHo($chr);
		}
		$patient->excluded($status);
	}
	unless (exists $self->variants_excluded->{$chr->id()}) {
		$self->variants_excluded->{$chr->id()} = $var_excl;
	}
	else { $self->variants_excluded->{$chr->id()} += $var_excl; }
	my $var = $self->getNewVector();
	foreach my $patient (@{$self->getPatients()}) {
		next if ($patient->in_the_attic());
#		warn 'delete var '.$patient->name();
		$patient->delete_variants($chr, $var_excl);
		$var += $patient->getVariantsVector($chr);
	}
	return $var;
}

#sub setIntersectPatients {
#	my ($self, $patients_names, $status, $separator) = @_;
#	return unless ($status eq '1');
#	my $var_inter = $self->chromosome->getNewVector();
#	my @lPatInter = split($separator, $patients_names);
#	foreach my $patName (@lPatInter) {
#		my $patient = $self->getPatient($patName);
#		next if ($patient->in_the_attic());
#		$patient->intersected($status);
#		if ($var_inter->is_empty()) { $var_inter += $patient->getVariantsVector(); }
#		else { $var_inter->Intersection( $var_inter, $patient->getVariantsVector() ); }
#	}
#	foreach my $patient (@{$self->getPatients()}) {
#		next if ($patient->in_the_attic());
#		$patient->getVariantsVector->Intersection( $patient->getVariantsVector() , $var_inter );
#		$patient->update_categories();
#	}
#}
#
#sub setExcludePatients {
#	my ($self, $patients_names, $status, $separator) = @_;
#	my $var_excl = $self->chromosome->getNewVector();
#	my @lPatInter = split($separator, $patients_names);
#	foreach my $patName (@lPatInter) {
#		my $patient = $self->getPatient($patName);
#		next if ($patient->in_the_attic());
#		$patient->excluded($status);
#		if    ($status eq 'all') { $var_excl += $patient->getVariantsVector(); }
#		elsif ($status eq 'he')  { $var_excl += $patient->getHe(); }
#		elsif ($status eq 'ho')  { $var_excl += $patient->getHo(); }
#	}
#	foreach my $patient (@{$self->getPatients()}) {
#		next if ($patient->in_the_attic());
#		$patient->delete_variants($var_excl);
#		$patient->update_categories();
#	}
#}

#sub stats {
#	my ($self, $chr) = @_;
#	my $hashCount = $self->getHashCountTypeVariants($chr);
#	my $typeFilters = $chr->project->typeFilters();
#	my $hash;
#	$hash->{name} 			= $self->name();
#	$hash->{id} 			= $self->name();
#	$hash->{nb} 			= scalar(keys %{$self->patients()});
#	$hash->{model} 			= $self->used_model();
#	$hash->{model} 			= undef if ($typeFilters eq 'individual' or $typeFilters eq 'familial');
#	$hash->{include} 		= 1;
#	$hash->{include} 		= 0 if ($self->intersected());
#	$hash->{include} 		= 2 if ($self->in_the_attic());
#	$hash->{include} 		= -1 if ($self->excluded());
#	$hash->{substitution} 	= $hashCount->{substitution};
#	$hash->{insertion} 		= $hashCount->{insertion};
#	$hash->{deletion} 		= $hashCount->{deletion};
#	$hash->{homozygote} 	= $hashCount->{ho};
#	$hash->{heterozygote} 	= $hashCount->{he};
#	$hash->{genes} 			= 1;
#	#$hash->{genes} 			= $self->getNbGenes($chr) if ($hashCount->{total} > 0);
#	#warn Dumper $hash; die;
#	return $hash;
#}
#
#sub stats_null {
#	my ($self, $chr) = @_;
#	my $hashCount = $self->getHashCountTypeVariants($chr);
#	my $typeFilters = $chr->project->typeFilters();
#	my $hash;
#	$hash->{name} 			= $self->name();
#	$hash->{id} 			= $self->name();
#	$hash->{nb} 			= scalar(keys %{$self->patients()});
#	$hash->{model} 			= undef;
#	$hash->{include} 		= 1;
#	$hash->{include} 		= 0 if ($self->intersected());
#	$hash->{include} 		= 2 if ($self->in_the_attic());
#	$hash->{include} 		= -1 if ($self->excluded());
#	$hash->{substitution} 	= 0;
#	$hash->{insertion} 		= 0;
#	$hash->{deletion} 		= 0;
#	$hash->{homozygote} 	= 0;
#	$hash->{heterozygote} 	= 0;
#	$hash->{genes} 			= 0;
#	return $hash;
#}

# renvoie un hash avec le comptage ds variants par annotation
sub getHashCountTypeVariants {
	my ($self, $chr) = @_;
	my $hash;
	$hash->{'substitution'} = 0;
	$hash->{'insertion'} = 0;
	$hash->{'deletion'} = 0;
	$hash->{'ho'} = 0;
	$hash->{'he'} = 0;
	foreach my $patient (@{$self->getPatients()}) {
		my $v_pat = $patient->getVariantsVector($chr)->Clone();
		my $v_pat_he = $patient->getHe($chr)->Clone();
		my $v_pat_ho = $patient->getHo($chr)->Clone();
		
		my $v_sub = $chr->getVectorSubstitutions();
		$v_sub->Intersection($v_sub, $v_pat);
		
		my $v_ins = $chr->getVectorInsertions();
		$v_ins->Intersection($v_ins, $v_pat);
		
		my $v_del = $chr->getVectorDeletions();
		$v_del->Intersection($v_del, $v_pat);
		
		$hash->{'substitution'} += $v_sub->Norm();
		$hash->{'insertion'} += $v_ins->Norm();
		$hash->{'deletion'} += $v_del->Norm();
		$hash->{'ho'} += $v_pat_ho->Norm();
		$hash->{'he'} += $v_pat_he->Norm();
		$hash->{'total'} += $v_pat->Norm();
	}
	return $hash;
}

# renvoie un hash avec un vector pour chaque annotation, en mergant les differents patients de la famiile
sub getHashTypeVariants {
	my ($self, $chr) = @_;
	my $hash;
	foreach my $patient (@{$self->getPatients()}) {
		my $hashPat = $patient->getHashTypeVariants($chr);
		foreach my $category (keys %$hashPat) {
			if (exists $hash->{$category}) {
				$hash->{$category} += $hashPat->{$category};
			}
			else {
				$hash->{$category} = $hashPat->{$category};
			}
		}
		if (exists $hashPat->{'large_deletion'}) { $hash->{'deletion'} += $hashPat->{'large_deletion'}; }
	}
	return $hash;
}

# renvoies le nombre de genes de ce group
sub getNbGenes {
	my ($self, $chr) = @_;
	#return $self->nb_genes() if ($self->nb_genes());
	return 0 if ($self->excluded() eq '1');
	my $hash;
	my $nb = 0;
	my $vector_fam = $chr->getNewVector();
	my $vector_tmp = $chr->getNewVector();
	foreach my $patient (@{$self->getPatients()}) {
		$vector_fam += $patient->getVariantsVector($chr);
	}
	return 0 if ($vector_fam->is_empty());
	foreach my $gene (@{$chr->getGenes()}) {
		next if ($gene->is_intergenic());
		$vector_tmp->Empty();
		$vector_tmp->Intersection($vector_fam, $gene->getVariantsVector());
		unless ($vector_tmp->is_empty()) {
			$nb++;
		}
	}
	return $nb;
}

sub checkModel_somatic_dbl_evt {
	my ($self, $gene, $chr) = @_;
	my $var_tmp = $chr->getNewVector();
	my $var_global = $chr->getNewVector();
	my $var_sain = $chr->getNewVector();
	my $var_malade = $chr->getNewVector();
	my $var_common = $chr->getNewVector();
	my $var_global_patients = $chr->getNewVector();
	$self->used_model('dbl_evt');
	foreach my $patient (@{$self->getPatients()}) {
		$var_tmp->Empty();
		$patient->used_model('dbl_evt');
		$patient->{model_vector_var_ok}->{$chr->id()} = $chr->getNewVector() unless (exists $patient->{model_vector_var_ok}->{$chr->id()});
		$patient->{somatic_dbl_evt}->{$chr->id()} = $chr->getNewVector() unless (exists $patient->{somatic_dbl_evt}->{$chr->id()});
		$var_tmp->Intersection( $patient->getVariantsVector($chr), $gene->getVariantsVector() );
		next if ($var_tmp->is_empty());
		if ($patient->tissue() eq 'C') {
			$var_sain += $var_tmp;
		}
		elsif ($patient->tissue() eq 'T') {
			unless ($patient->getHo($chr)->is_empty()) {
				$patient->getVariantsVector($chr)->Intersection( $patient->getVariantsVector($chr), $patient->getHe($chr) );
			}
			$var_tmp->Intersection( $var_tmp, $patient->getHe($chr) );
			$var_malade += $var_tmp;
		}
	}
	$var_common->Intersection($var_sain, $var_malade);
	$var_sain   -= $var_malade;
	$var_sain   -= $var_common;
	$var_malade -= $var_sain;
	$var_malade -= $var_common;
	$var_tmp->Empty();
	$var_tmp = $var_common + $var_malade;
	unless ($var_malade->is_empty() || $var_common->is_empty()) {
		foreach my $patient (@{$self->getPatients()}) {
			if ($patient->tissue() eq 'C') {
				$patient->{model_vector_var_ok}->{$chr->id()} += $var_common;
				$patient->{somatic_dbl_evt}->{$chr->id()} += $var_common;
			}
			elsif ($patient->tissue() eq 'T') {
				$patient->{model_vector_var_ok}->{$chr->id()} += $var_common;
				$patient->{model_vector_var_ok}->{$chr->id()} += $var_malade;
				$patient->{somatic_dbl_evt}->{$chr->id()} += $var_common;
				$patient->{somatic_dbl_evt}->{$chr->id()} += $var_malade;
			}
			$var_global_patients += $patient->{model_vector_var_ok}->{$chr->id()};
		}
		$var_global += $var_global_patients;
	}
	return $var_global;
}

1;