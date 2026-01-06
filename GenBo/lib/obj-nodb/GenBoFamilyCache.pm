package GenBoFamilyCache;
use strict;
use Storable qw(thaw retrieve);
use Moo;
use Bit::Vector;
use Bit::Vector::Overload;
use Bio::DB::Sam;
use Data::Dumper;
extends 'GenBoFamily', 'GenBoCache';



# tag si famille in the attic
has in_the_attic => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		foreach my $patient (@{$self->getPatients()}) {
			next if ($patient->in_the_attic());
			return;
		}
		return 1
	},
);

# tag si famille exclue (all, ho, he)
has excluded => (
	is		=> 'rw',
	lazy	=> 1,
	default => 0,
);

# tag si famille intersect
has intersected => (
	is		=> 'rw',
	default => 0,
);

has hash_models_genetics_used => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub { {} },
);

# nom du model genetique utilise sur cette famille
sub used_model {
	my $self = shift;
	return if (scalar keys %{$self->hash_models_genetics_used()} == 0);
	my @lModels;
	foreach my $model (sort keys (%{$self->hash_models_genetics_used()})) {
		if ($model eq 'recessif') {
			my $has_parent;
			foreach my $parent (@{$self->getParents()}) {
				next if ($parent->in_the_attic());
				$has_parent = 1;
			}
			$model = 'recessif no parent' unless ($has_parent);
		}
		push(@lModels, $model);
	}
	return join(' + ', @lModels);
}

# variants issus du model compound a garder
has keep_var_compound => (
	is		=> 'rw',
	default => sub { {} },
);

# vector des variants exclus par les differents GenBoPatientCache
has variants_excluded => (
	is      => 'rw',
	lazy    => 1,
	default => sub { {} },
);

# vector des variants intersectes par les differents GenBoPatientCache
has variants_intersected => (
	is      => 'rw',
	lazy    => 1,
	default => undef,
);

has stats_categories => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->project->buffer->config->{stats_families};
	},
);



##### METHODS #####



sub get_vector_keep_var_compound {
	my ($self, $chr) = @_;
	my $vector = $chr->getNewVector();
	if (exists $self->keep_var_compound->{$chr->id()}) {
		$vector = Bit::Vector->new_Enum($chr->getVariantsVector->Size(), join(',', @{$self->keep_var_compound->{$chr->id()}}));
	}
	return $vector;
}

# renvoies le nombre de genes dans cette famille
sub getNbGenes {
	my ($self, $chr) = @_;
	#return $self->nb_genes() if ($self->nb_genes());
	return 0 if ($self->excluded() eq '1');
	return 0 if ($chr->getVariantsVector->is_empty());
	my $hash;
	my $nb = 0;
	my $vector_fam = $chr->getNewVector();
	foreach my $patient (@{$self->getPatients()}) {
		$vector_fam += $patient->getVariantsVector($chr);
	}
	$vector_fam->Intersection($vector_fam, $chr->getVariantsVector());
	foreach my $gene (@{$chr->getGenes()}) {
		next if ($gene->is_intergenic());
		my $v_tmp = $vector_fam->Clone();
		$v_tmp->Intersection($v_tmp, $gene->getVariantsVector());
		unless ($v_tmp->is_empty()) {
			$nb++;
			#warn '-> '.$gene->external_name().': '.$nb;
		}
	}
	return $nb;
}

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

# renvoie le vector MAJ pour l'annotation demandee
sub getTypeVariants {
	my ($self, $chr, $type) = @_;
	confess("\n\nERROR: need a GenBoChromosomeCache. Die...\n\n") unless ($chr);
	return $self->getHe($chr) if ($type eq 'he');
	my $var = $chr->getNewVector();
	foreach my $patient (@{$self->getPatients()}) {
		$var += $patient->getCategoryVariantsVector($chr, $type);
	}
	$var->Intersection( $var, $self->getVariantsVector($chr) );
	return $var;
}

# renvoie le vector des SUB
sub getSubstitutions {
	my ($self, $chr) = @_;
	confess("\n\nERROR: need a GenBoChromosomeCache. Die...\n\n") unless ($chr);
	return $self->getTypeVariants($chr, 'substitution');
}

# renvoie le vector des DEL
sub getDeletions {
	my ($self, $chr) = @_;
	confess("\n\nERROR: need a GenBoChromosomeCache. Die...\n\n") unless ($chr);
	return $self->getTypeVariants($chr, 'deletion');
}

# renvoie le vector des INS
sub getInsertions {
	my ($self, $chr) = @_;
	confess("\n\nERROR: need a GenBoChromosomeCache. Die...\n\n") unless ($chr);
	return $self->getTypeVariants($chr, 'insertion');
}

# renvoie le vector des HO
sub getHo {
	my ($self, $chr) = @_;
	confess("\n\nERROR: need a GenBoChromosomeCache. Die...\n\n") unless ($chr);
	return $self->getTypeVariants($chr, 'ho');
}

# renvoie le vector des HE
sub getHe {
	my ($self, $chr) = @_;
	confess("\n\nERROR: need a GenBoChromosomeCache. Die...\n\n") unless ($chr);
	my $var_ho = $self->getHo($chr);
	my $var_he = $chr->getNewVector();
	$var_he = $self->getVariantsVector($chr) - $var_ho;
	return $var_he;
}

# initialise les infos issues du PED
sub initPedigree {
	my $self = shift;
	$self->patients();
	$self->parents();
	$self->children();
}

# methode pour intersecter les patients en mode FAM
sub setIntersectPatients {
	my ($self, $patients, $chr) = @_;
	my $var = $chr->getNewVector();
	unless (exists $self->{variants_intersected}->{$chr->id()}) {
		$self->{variants_intersected}->{$chr->id()} = $chr->getNewVector();
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

# methode pour exclure les patients en mode FAM
sub setExcludePatients {
	my ($self, $lPatients, $status, $chr) = @_;
	my $var_excl = $chr->getNewVector();
	foreach my $patient (@$lPatients) {
		next if ($patient->in_the_attic());
		if ($status eq 'all') {
			$var_excl += $patient->getVectorOrigin($chr);
			$patient->variants_excluded->{$chr->id()} = $patient->getVectorOrigin($chr);
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
	my $var = $chr->getNewVector();
	foreach my $patient (@{$self->getPatients()}) {
		next if ($patient->in_the_attic());
#		warn 'delete var '.$patient->name();
		$patient->getVariantsVector($chr)->AndNot($patient->getVariantsVector($chr), $var_excl);
		$var += $patient->getVariantsVector($chr);
	}
	return $var;
}

# methode pour intersecter les patients en mode FAM
sub getGenes_intersect {
	my ($self, $patients, $chr) = @_;
	my ($hGenes, $hGenesToDel);
	foreach my $patient (@$patients) {
		$patient->intersected(1);
		foreach my $gene_id (@{$patient->hash_list_genes_ids->{$chr->name()}}) {
			$hGenes->{$gene_id}->{$patient->name()} = undef;
		}
	}
	foreach my $gene_id (keys %$hGenes) {
		unless (scalar(keys %{$hGenes->{$gene_id}}) == scalar(@$patients)) {
			delete $hGenes->{$gene_id};
		}
	}
	foreach my $patient (@{$self->getPatients()}) {
		foreach my $gene (@{$chr->getGenes()}) {
			unless (exists $hGenes->{$gene->id()}) {
				$patient->delete_variants_from_gene($chr, $gene);
			}
		}
		$patient->update_list_genes($chr);
	}
	return keys %$hGenes;
}

# applique le model denovo incoherence HE enfant malade mais HO chez les deux parents saints (FAM)
sub getModelVector_fam_pb_ho_he {
	my ($self, $chr) = @_;
	return $self->{vector_transmission}->{pb_ho_he}->{$chr->id} if exists $self->{vector_transmission}->{pb_ho_he}->{$chr->id};
	my $var_global = $chr->getNewVector();
	my $var_model_incoherence = $chr->getNewVector();
	$var_model_incoherence->Empty();
	foreach my $patient (@{$self->getPatientsIll()}) {
		$var_model_incoherence += $patient->getHe($chr);
	}
	my @lParents = @{$self->getParentsHealthy()};
	if (scalar(@lParents) == 2) {
		my $var_ho_parents = $chr->getNewVector();
		$var_ho_parents->Intersection( $lParents[0]->getHo($chr), $lParents[1]->getHo($chr) );
		$var_model_incoherence->Intersection( $var_model_incoherence, $var_ho_parents );
	}
	else { $var_model_incoherence->Empty(); }
	$var_global += $var_model_incoherence;
	$self->{vector_transmission}->{pb_ho_he}->{$chr->id} = $var_global;
	return $var_global;
}

# variant Ho enfant malade et present uniquement chez un parent et en plus HE. Jamais HO chez enfant sain
sub getModelVector_fam_uniparental_disomy {
	my ($self, $chr) = @_;
	return $self->{vector_transmission}->{uniparental_disomy}->{$chr->id} if exists $self->{vector_transmission}->{uniparental_disomy}->{$chr->id};
#	$self->used_model('uniparental_disomy');
	$self->getModelVector_fam_uniparental_disomy_chr_x($chr) if ($chr->id() eq 'X');
	my $var_global = $chr->getNewVector();
	my @lParents = @{$self->getParentsHealthy()};
	my $var_he_only_one_parent = $chr->getNewVector();
	if (scalar(@lParents) == 2) {
		my $var_par_1_he = $lParents[0]->getVectorHe($chr)->Clone();
		my $var_par_2_he = $lParents[1]->getVectorHe($chr)->Clone();
		my $var_par_commons_he = $chr->getNewVector();
		$var_par_commons_he->Intersection($var_par_1_he, $var_par_2_he);
		$var_par_1_he -= $var_par_commons_he;
		$var_par_2_he -= $var_par_commons_he;
		$var_he_only_one_parent += $var_par_1_he;
		$var_he_only_one_parent += $var_par_2_he;
		$var_he_only_one_parent -= $lParents[0]->getVectorHo($chr);
		$var_he_only_one_parent -= $lParents[1]->getVectorHo($chr);
#		$self->used_model('uniparental_disomy');
	}
	else { return $chr->getNewVector(); }
	
	my $var_all_pat_healthy_ho = $chr->getNewVector();
	foreach my $patient (@{$self->getPatientsHealthy()}) {
		$var_all_pat_healthy_ho += $patient->getVectorHo($chr);
	}
	
	foreach my $patient (@{$self->getPatientsIll()}) {
		my $var_child_ho = $patient->getVectorHo($chr)->Clone();
		$var_child_ho -= $var_all_pat_healthy_ho;
		$var_child_ho->Intersection($var_child_ho, $var_he_only_one_parent);
		$var_global += $var_child_ho;
	}
	
	#check strict uniparental disomy
	my $var_ok_from_parent_1 = $var_he_only_one_parent->Clone();
	$var_ok_from_parent_1->Intersection($var_ok_from_parent_1, $var_global);
	foreach my $var (@{$chr->getListVarObjects($var_ok_from_parent_1)}) {
		$self->project->print_dot(20);
		unless ($var->check_no_var_allele_in_bam($lParents[0], 4)) {
			$var_global->Bit_Off($var->vector_id());
		}
	}
	
	my $var_ok_from_parent_2 = $var_he_only_one_parent->Clone();
	$var_ok_from_parent_2->Intersection($var_ok_from_parent_2, $var_global);
	foreach my $var (@{$chr->getListVarObjects($var_ok_from_parent_2)}) {
		$self->project->print_dot(20);
		unless ($var->check_no_var_allele_in_bam($lParents[1], 4)) {
			$var_global->Bit_Off($var->vector_id());
		}
	}
	
#	foreach my $patient (@{$self->getPatients()}) {
#		$patient->used_model('uniparental_disomy');
#	}
	
	$self->{vector_transmission}->{uniparental_disomy}->{$chr->id} = $var_global;
	return $var_global;
}

sub getModelVector_fam_uniparental_disomy_chr_x {
	my ($self, $chr) = @_;
	return $self->{vector_transmission}->{uniparental_disomy}->{$chr->id} if exists $self->{vector_transmission}->{uniparental_disomy}->{$chr->id};
	my $var_global = $chr->getNewVector();
	my @lParents = @{$self->getParentsHealthy()};
	my $var_he_only_one_parent = $chr->getNewVector();
	if (scalar(@lParents) == 2) {
		my $var_par_1_he = $lParents[0]->getVectorHe($chr)->Clone();
		my $var_par_2_he = $lParents[1]->getVectorHe($chr)->Clone();
		my $var_par_commons_he = $chr->getNewVector();
		$var_par_commons_he->Intersection($var_par_1_he, $var_par_2_he);
		$var_par_1_he -= $var_par_commons_he;
		$var_par_2_he -= $var_par_commons_he;
		$var_he_only_one_parent += $var_par_1_he;
		$var_he_only_one_parent += $var_par_2_he;
		$var_he_only_one_parent -= $lParents[0]->getVectorHo($chr);
		$var_he_only_one_parent -= $lParents[1]->getVectorHo($chr);
	}
	else { return $chr->getNewVector(); }
	
	my $var_all_pat_healthy_ho = $chr->getNewVector();
	foreach my $patient (@{$self->getPatientsHealthy()}) {
		$var_all_pat_healthy_ho += $patient->getVectorHo($chr);
	}
	
	foreach my $patient (@{$self->getPatientsIll()}) {
#		$patient->used_model('uniparental_disomy');
#		warn $self->name().' - '.$patient->name().' -> '.$patient->sex();
		next if ($patient->sex() eq '1');
		my $var_child_ho = $patient->getVectorHo($chr)->Clone();
		$var_child_ho -= $var_all_pat_healthy_ho;
		$var_child_ho->Intersection($var_child_ho, $var_he_only_one_parent);
		$var_global += $var_child_ho;
	}
	
	#check strict uniparental disomy
	my $var_ok_from_parent_1 = $var_he_only_one_parent->Clone();
	$var_ok_from_parent_1->Intersection($var_ok_from_parent_1, $var_global);
	foreach my $var (@{$chr->getListVarObjects($var_ok_from_parent_1)}) {
		$self->project->print_dot(20);
		unless ($var->check_no_var_allele_in_bam($lParents[0], 3)) {
			$var_global->Bit_Off($var->vector_id());
		}
	}
	
	my $var_ok_from_parent_2 = $var_he_only_one_parent->Clone();
	$var_ok_from_parent_2->Intersection($var_ok_from_parent_2, $var_global);
	foreach my $var (@{$chr->getListVarObjects($var_ok_from_parent_2)}) {
		$self->project->print_dot(20);
		unless ($var->check_no_var_allele_in_bam($lParents[1], 3)) {
			$var_global->Bit_Off($var->vector_id());
		}
	}
	$self->{vector_transmission}->{uniparental_disomy}->{$chr->id} = $var_global;
	return $var_global;
}

sub getModelVector_fam_denovo() {
	my ($self, $chr) = @_;
	return $self->{vector_transmission}->{denovo}->{$chr->id} if exists $self->{vector_transmission}->{denovo}->{$chr->id};
	my $var_excluded = $chr->getNewVector();
	my $var_keeped = $chr->getNewVector();
	my $var_ok = $chr->getNewVector();
	#my $var_pb_ho_he = $self->getModelVector_fam_pb_ho_he($chr);
	foreach my $parent (@{$self->getParents()}) { $var_excluded += $parent->getVariantsVector($chr); }
	foreach my $child (@{$self->getChildrenHealthy()}) {
		$var_excluded += $child->getVariantsVector($chr);
	}
	foreach my $child (@{$self->getChildrenIll()}) {
		$var_keeped += $child->getVariantsVector($chr);
	}
	$var_keeped -= $var_excluded;
	foreach my $child (@{$self->getChildrenIll()}) {
		my $v_child = $child->getVariantsVector($chr)->Clone();
		$v_child->Intersection( $v_child, $var_keeped );
		$var_ok += $v_child;
	}
	$self->{vector_transmission}->{denovo}->{$chr->id} = $var_ok;
	return $var_ok;
}

sub getVectorParent {
	my ($self, $chr) = @_;
	return $self->{vector_transmission}->{parent}->{$chr->id} if exists $self->{vector_transmission}->{parent}->{$chr->id};
	 $self->{vector_transmission}->{parent}->{$chr->id} = $chr->getNewVector();
	foreach my $parent (@{$self->getParents()}) { $self->{vector_transmission}->{parent}->{$chr->id} += $parent->getVariantsVector($chr); }
	return  $self->{vector_transmission}->{parent}->{$chr->id};
}



sub getVectorChildrenHealthy {
	my ($self, $chr) = @_;
	return $self->{vector_transmission}->{children_healthy}->{$chr->id} if exists $self->{vector_transmission}->{children_healthy}->{$chr->id};
	$self->{vector_transmission}->{children_healthy}->{$chr->id} = $chr->getNewVector();
	
	foreach my $parent (@{$self->getChildrenHealthy()}) { $self->{vector_transmission}->{children_healthy}->{$chr->id} += $parent->getVariantsVector($chr); }
		
	return  $self->{vector_transmission}->{children_healthy}->{$chr->id};
}


sub getVector_individual_denovo() {
	my ($self, $chr,$child) = @_;
	my $key = "denovo_".$child->id;
	return $self->{vector_transmission}->{$key}->{$chr->id} if exists $self->{vector_transmission}->{$key}->{$chr->id};
	unless ($self->isTrio) {
		 $self->{vector_transmission}->{$key}->{$chr->id} = $chr->getNewVector();
		return  $self->{vector_transmission}->{$key}->{$chr->id};
	}
	$self->{vector_transmission}->{$key}->{$chr->id} = $child->getVariantsVector($chr);
	$self->{vector_transmission}->{$key}->{$chr->id} -= $self->getVectorParent($chr);
	return $self->{vector_transmission}->{$key}->{$chr->id};
}

sub getVector_family_denovo {
	my ($self, $chr) = @_;
	confess() unless $chr;
	return $self->{vector_transmission}->{denovo_fam}->{$chr->id} if exists $self->{vector_transmission}->{denovo_fam}->{$chr->id};
	unless ($self->isTrio) {
		$self->{vector_transmission}->{denovo_fam}->{$chr->id} = $chr->getNewVector();
		return $self->{vector_transmission}->{denovo_fam}->{$chr->id}
	}
	foreach my $child (@{$self->getChildrenIll()}) {
		unless ($self->{vector_transmission}->{denovo_fam}->{$chr->id}){
			$self->{vector_transmission}->{denovo_fam}->{$chr->id} = $self->getVector_individual_denovo($chr,$child)->Clone;
		} 
		else {
			$self->{vector_transmission}->{denovo_fam}->{$chr->id} &= $self->getVector_individual_denovo($chr,$child);
		}
	}
	if (exists $self->{vector_transmission}->{denovo_fam}->{$chr->id}) {
		$self->{vector_transmission}->{denovo_fam}->{$chr->id} -= $self->getVectorChildrenHealthy($chr);
	}
	else {
		$self->{vector_transmission}->{denovo_fam}->{$chr->id} = $chr->getNewVector();
	}
	return $self->{vector_transmission}->{denovo_fam}->{$chr->id};
}

sub getVector_individual_strict_denovo() {
	my ($self, $chr,$child) = @_;
	my $key = "strict_denovo_".$child->id;
	return $self->{vector_transmission}->{$key}->{$chr->id} if exists $self->{vector_transmission}->{$key}->{$chr->id};
	$self->{vector_transmission}->{$key}->{$chr->id} = $self->getVector_individual_denovo($chr,$child);
	$self->{vector_transmission}->{$key}->{$chr->id} &= $chr->hash_cache_strict_denovo->{$child->name()};
	return $self->{vector_transmission}->{$key}->{$chr->id};
}

sub getVector_family_strict_denovo() {
	my ($self, $chr) = @_;
	my $key = "strict_denovo_fam";
	return $self->{vector_transmission}->{$key}->{$chr->id} if exists $self->{vector_transmission}->{$key}->{$chr->id};
	 if (@{$self->getChildrenIll()}  == 1 ){
	 	$self->{vector_transmission}->{$key}->{$chr->id} = $self->getVector_individual_strict_denovo($chr,$self->getChildrenIll()->[0]);
	 }
	 else {
	 	foreach my $child (@{$self->getChildrenIll()}) {
		unless ($self->{vector_transmission}->{$key}->{$chr->id}){
			$self->{vector_transmission}->{$key}->{$chr->id} = $self->getVector_individual_denovo($chr,$child)->Clone;
		} 
		else {
			$self->{vector_transmission}->{$key}->{$chr->id} &= $self->getVector_individual_denovo($chr,$child);
		}
	}
	$self->{vector_transmission}->{$key}->{$chr->id} -= $self->getVectorChildrenHealthy($chr);
	return $self->{vector_transmission}->{$key}->{$chr->id};
	 	
	 }
}
sub getVector_individual_dominant {
	my ($self,$chr,$child) = @_;
	my $key = "dominant_".$child->id;
	return $self->{vector_transmission}->{$key}->{$chr->id} if exists $self->{vector_transmission}->{$key}->{$chr->id};
	
	$self->{vector_transmission}->{$key}->{$chr->id} = $child->getVectorOrigin($chr)->Clone;
	foreach my $parent (@{$self->getParents}){
		if ($parent->isIll){
			$self->{vector_transmission}->{$key}->{$chr->id} &= $parent->getVectorOrigin($chr);
		}
		else {
			$self->{vector_transmission}->{$key}->{$chr->id} -= $parent->getVectorOrigin($chr);
		}
	}
	return $self->{vector_transmission}->{$key}->{$chr->id};
	
}



sub getVector_family_dominant() {
	my ($self, $chr) = @_;
	my $key = "family_dominant";
	return $self->{vector_transmission}->{$key}->{$chr->id} if exists $self->{vector_transmission}->{$key}->{$chr->id};
	my @lPatIll = @{$self->getPatientsIll()};
	
	
	foreach my $child (@{$self->getChildrenIll()}) {
		unless ($self->{vector_transmission}->{$key}->{$chr->id}){
			$self->{vector_transmission}->{$key}->{$chr->id} = $self->getVector_individual_dominant($chr,$child)->Clone;
		} 
		else {
			$self->{vector_transmission}->{$key}->{$chr->id} &= $self->getVector_individual_dominant($chr,$child);
		}
	}
	
	foreach my $patient (@{$self->getChildrenHealthy()}) {
		unless ($self->{vector_transmission}->{$key}->{$chr->id}){
			$self->{vector_transmission}->{$key}->{$chr->id} = $chr->getNewVector();
		} 
		$self->{vector_transmission}->{$key}->{$chr->id} -= $patient->getVectorOrigin($chr);
	}
	
	return $self->{vector_transmission}->{$key}->{$chr->id};
}

sub getModelVector_som_only_tissues_somatic {
	my ($self, $chr) = @_;
	my $key = "only_tissues_somatic";
	return $self->{vector_transmission}->{$key}->{$chr->id} if exists $self->{vector_transmission}->{$key}->{$chr->id};
	$self->{vector_transmission}->{$key}->{$chr->id} = $chr->getNewVector();
	foreach my $pat_somatic (@{$self->getSomatics()}) {
		$self->{vector_transmission}->{$key}->{$chr->id} += $pat_somatic->getVariantsVector($chr);
	}
	foreach my $pat_germinal (@{$self->getGerminals()}) {
		$self->{vector_transmission}->{$key}->{$chr->id} -= $pat_germinal->getVariantsVector($chr);
	}
	return $self->{vector_transmission}->{$key}->{$chr->id};
}
 
sub getVector_individual_recessive() {
	my ($self,$chr,$child) = @_;
	my $key = "recessive_".$child->id;
	return $self->{vector_transmission}->{$key}->{$chr->id} if exists $self->{vector_transmission}->{$key}->{$chr->id};
	$self->{vector_transmission}->{$key}->{$chr->id} = $child->getHo($chr);
	#my $vparent = $chr->getNewVector();
	$self->{vector_transmission}->{$key}->{$chr->id} &= $self->getMother->getHe($chr) if ($self->getMother());
	unless ($chr->name eq "X"){
		$self->{vector_transmission}->{$key}->{$chr->id} &= $self->getFather->getHe($chr) if ($self->getFather());
	}
	
	#$self->{vector_transmission}->{$key}->{$chr->id} &= $vparent;
	return $self->{vector_transmission}->{$key}->{$chr->id};
}




sub getVector_family_recessive() {
	my ($self, $chr) = @_;
	my $key = "family_recessive";
	return $self->{vector_transmission}->{$key}->{$chr->id} if exists $self->{vector_transmission}->{$key}->{$chr->id};
	foreach my $child (@{$self->getChildrenIll()}) {
		unless ($self->{vector_transmission}->{$key}->{$chr->id}){
			$self->{vector_transmission}->{$key}->{$chr->id} = $self->getVector_individual_recessive($chr,$child)->Clone;
		} 
		else {
			$self->{vector_transmission}->{$key}->{$chr->id} &= $self->getVector_individual_recessive($chr,$child);
		}
	}
	unless (exists $self->{vector_transmission}->{$key}->{$chr->id}) {
		$self->{vector_transmission}->{$key}->{$chr->id} = $chr->getNewVector();
		return $self->{vector_transmission}->{$key}->{$chr->id};
	}
	foreach my $patient (@{$self->getPatientsHealthy()}) {
		$self->{vector_transmission}->{$key}->{$chr->id} -= $patient->getHo($chr);
	}
	return $self->{vector_transmission}->{$key}->{$chr->id};
	
}



sub getVector_individual_mother  {
	my ($self, $chr,$child,$debug) = @_;
	die() unless $child;
	my $key = "mother_".$child->id;
	return $chr->getNewVector()  unless $self->getMother;
	return $self->{vector_transmission}->{$key}->{$chr->id} if exists $self->{vector_transmission}->{$key}->{$chr->id};
	$self->{vector_transmission}->{$key}->{$chr->id} = $child->getVariantsVector($chr);
	$self->{vector_transmission}->{$key}->{$chr->id} &= $self->getMotherVector($chr) if $self->getMother;
	$self->{vector_transmission}->{$key}->{$chr->id} -= $self->getFatherVector($chr) if $self->getFather;
	return $self->{vector_transmission}->{$key}->{$chr->id};
}




sub getVector_family_mother() {
	my ($self, $chr) = @_;
	my $key = "mother_family";
	return $self->{vector_transmission}->{$key}->{$chr->id} if exists $self->{vector_transmission}->{$key}->{$chr->id};
	foreach my $child (@{$self->getChildrenIll()}) {
		unless ($self->{vector_transmission}->{$key}->{$chr->id}){
			$self->{vector_transmission}->{$key}->{$chr->id} = $self->getVector_individual_mother($chr,$child)->Clone;
		} 
		else {
			$self->{vector_transmission}->{$key}->{$chr->id} &= $self->getVector_individual_mother($chr,$child);
		}
	}
	unless (exists $self->{vector_transmission}->{$key}->{$chr->id}) {
		$self->{vector_transmission}->{$key}->{$chr->id} = $chr->getNewVector();
		return $self->{vector_transmission}->{$key}->{$chr->id};
	}
	foreach my $child (@{$self->getChildrenHealthy()}) {
		$self->{vector_transmission}->{$key}->{$chr->id} -= $self->getVector_individual_mother($chr,$child);
	}
	return $self->{vector_transmission}->{$key}->{$chr->id};
}

sub getVector_individual_father() {
	my ($self, $chr,$child) = @_;
	my $key = "father_".$child->id;
	return $chr->getNewVector()  unless $self->getFather;
	return $self->{vector_transmission}->{$key}->{$chr->id} if exists $self->{vector_transmission}->{$key}->{$chr->id};
	$self->{vector_transmission}->{$key}->{$chr->id} = $child->getVariantsVector($chr) ;
	$self->{vector_transmission}->{$key}->{$chr->id} &= $self->getFatherVector($chr) if $self->getFather;
	$self->{vector_transmission}->{$key}->{$chr->id} -= $self->getMotherVector($chr) if $self->getMother;
	return $self->{vector_transmission}->{$key}->{$chr->id};
}

sub getVector_family_father() {
	my ($self, $chr) = @_;
	my $key = "father_family";
	return $self->{vector_transmission}->{$key}->{$chr->id} if exists $self->{vector_transmission}->{$key}->{$chr->id};
	foreach my $child (@{$self->getChildrenIll()}) {
		unless ($self->{vector_transmission}->{$key}->{$chr->id}){
			$self->{vector_transmission}->{$key}->{$chr->id} = $self->getVector_individual_father($chr,$child)->Clone;
		} 
		else {
			$self->{vector_transmission}->{$key}->{$chr->id} &= $self->getVector_individual_father($chr,$child);
		}
	}
	unless (exists $self->{vector_transmission}->{$key}->{$chr->id}) {
		$self->{vector_transmission}->{$key}->{$chr->id} = $chr->getNewVector();
		return $self->{vector_transmission}->{$key}->{$chr->id};
	}
	foreach my $child (@{$self->getChildrenHealthy()}) {
		$self->{vector_transmission}->{$key}->{$chr->id} -= $self->getVector_individual_father($chr,$child);
	}
	return $self->{vector_transmission}->{$key}->{$chr->id};
}





sub getModelVector_fam_denovo_in_all_children() {
	my ($self, $chr) = @_;
	return $self->{vector_transmission}->{denovo_fam}->{$chr->id} if exists $self->{vector_transmission}->{denovo_fam}->{$chr->id};
	my $var_excluded = $chr->getNewVector();
	
	my $var_ok = $chr->getNewVector();
	#my $var_pb_ho_he = $self->getModelVector_fam_pb_ho_he($chr);
	foreach my $parent (@{$self->getParents()}) { $var_excluded += $parent->getVariantsVector($chr); }
	foreach my $child (@{$self->getChildrenHealthy()}) {
		$var_excluded += $child->getVariantsVector($chr);
	}
	my $var_denovo;#= $chr->getNewVector();
	foreach my $child (@{$self->getChildrenIll()}) {
		unless ($var_denovo){
			$var_denovo = $child->getVariantsVector($chr)->Clone;
		} 
		else {
		$var_denovo &= $child->getVariantsVector($chr);
		}
	}
	
	$var_denovo -= $var_excluded;
	
	
	$self->{vector_transmission}->{denovo_fam}->{$chr->id} = $var_denovo;
	return $var_ok;
}

#sub getModelVector_fam_incomplete_penetrance {
#	my ($self, $chr) = @_;
#	return $self->{vector_transmission}->{incomplete_penetrance}->{$chr->id} if exists $self->{vector_transmission}->{incomplete_penetrance}->{$chr->id};
#	my $var_model = $chr->getNewVector();
#	foreach my $patient (@{$self->getChildrenIll()}) {
#		my @lPat = @{$self->getParentsHealthy()};
#		next unless (scalar @lPat == 2);
#		my $var_child_1 = $patient->getVariantsVector($chr)->Clone();
#		$var_child_1->Intersection($var_child_1, $lPat[0]->getVariantsVector($chr));
#		$var_child_1 -= $lPat[1]->getVariantsVector($chr);
#		my $var_child_2 = $patient->getVariantsVector($chr)->Clone();
#		$var_child_2->Intersection($var_child_2, $lPat[1]->getVariantsVector($chr));
#		$var_child_2 -= $lPat[0]->getVariantsVector($chr);
#		$var_model += $var_child_1;
#		$var_model += $var_child_2;
#	}
#	$self->{vector_transmission}->{incomplete_penetrance}->{$chr->id} = $var_model;
#	return $var_model;
#}

sub getModelVector_fam_strict_denovo() {
	my ($self, $chr) = @_;
	
	return $self->{vector_transmission}->{strict_denovo}->{$chr->id} if exists $self->{vector_transmission}->{strict_denovo}->{$chr->id};
	my $var_model = $chr->getNewVector();
	return $var_model unless (-d $chr->dir_sqlite_strict_denovo);
	return $var_model unless (-e $chr->dir_sqlite_strict_denovo.'/'.$chr->name().'.lite');
	
	foreach my $patient (@{$self->getChildrenIll()}) {
		$var_model += $chr->hash_cache_strict_denovo->{$patient->name()};
	}
	foreach my $patient (@{$self->getPatientsHealthy()}) {
		$var_model -= $patient->getVariantsVector($chr);
	}
	$self->{vector_transmission}->{strict_denovo}->{$chr->id} = $var_model;
	return $var_model;
}

sub getModelVector_fam_dominant() {
	my ($self, $chr) = @_;
	return $self->{vector_transmission}->{dominant}->{$chr->id} if exists $self->{vector_transmission}->{dominant}->{$chr->id};
	my $var_model_classic = $chr->getNewVector();
	my $firstPatName;
	my @lPatIll = @{$self->getPatientsIll()};
	return $var_model_classic unless (@lPatIll);
	$var_model_classic += $lPatIll[0]->getVariantsVector($chr);
	shift(@lPatIll);
	foreach my $patient (@lPatIll) {
		$var_model_classic->Intersection( $var_model_classic, $patient->getVariantsVector($chr));
	}
	foreach my $patient (@{$self->getPatientsHealthy()}) {
		$var_model_classic -= $patient->getVariantsVector($chr);
	}
	$self->{vector_transmission}->{dominant}->{$chr->id} = $var_model_classic;
	return $var_model_classic;
}

sub getModelVector_fam_recessif() {
	my ($self, $chr, $is_recessif_compound) = @_;

	return $self->{vector_transmission}->{recessif}->{$chr->id} if exists $self->{vector_transmission}->{recessif}->{$chr->id};
	if ($chr->name() eq "X") {
		return $self->getModelVector_fam_recessif_chrX($chr,$is_recessif_compound);
	}

	my $var_sain = $chr->getNewVector();
	my $var_malade = $chr->getNewVector();
	my $var_he_parents = $chr->getNewVector();
	my $var_global = $chr->getNewVector();

	foreach my $child (@{$self->getChildrenHealthy()}) {
		next if ($child->in_the_attic());
		$var_sain += $child->getHo($chr);
	}
#	foreach my $child (@{$self->getChildrenHealthy()}) {
#		next if ($child->in_the_attic());
#		$chr->patients_categories->{$child->name()} -= $var_sain;
#	}
	foreach my $child (@{$self->getChildrenIll()}) {
		next if ($child->in_the_attic());
		$var_malade += $child->getHo($chr);
	}
	foreach my $child (@{$self->getChildrenIll()}) {
		next if ($child->in_the_attic());
		$var_malade->Intersection( $var_malade, $child->getHo($chr) );
	}
	$var_malade -= $var_sain;
	my $has_parents;
	foreach my $parent (@{$self->getParents()}) {
		next if ($parent->in_the_attic());
		if ($var_he_parents->is_empty()) { $var_he_parents += $parent->getHe($chr); }
		else { $var_he_parents->Intersection( $var_he_parents, $parent->getHe($chr) ); }
		$has_parents = 1;
	}
	$var_malade->Intersection( $var_malade, $var_he_parents ) if ($has_parents);
	my $has_patient;
	foreach my $patient (@{$self->getPatients()}) {
		next if ($patient->in_the_attic());
		my $v_pat = $patient->getVariantsVector($chr)->Clone();
		$v_pat->Intersection( $v_pat, $var_malade );
		$var_global += $v_pat;
		$has_patient = 1;
	}
#	unless ($is_recessif_compound) {
#		if ($has_parents) { $self->used_model('recessif'); }
#		elsif ($has_patient) { $self->used_model('recessif no parent'); }
#	}
	$self->{vector_transmission}->{recessif}->{$chr->id} = $var_global;
	return $var_global;
}

sub getModelVector_fam_recessif_chrX() {
	my ($self, $chr, $is_recessif_compound) = @_;
	return $self->{vector_transmission}->{recessif}->{$chr->id} if exists $self->{vector_transmission}->{recessif}->{$chr->id};
	my $var_sain = $chr->getNewVector();
	my $var_malade = $chr->getNewVector();
	my $var_he_parents = $chr->getNewVector();
	my $var_global = $chr->getNewVector();
	
	foreach my $child (@{$self->getChildrenHealthy()}) {
		next if ($child->in_the_attic());
		if ($child->sex() == 1) { $var_sain += $child->getVariantsVector($chr); }
		else { $var_sain += $child->getHo($chr); }
	}
#	foreach my $child (@{$self->getChildrenHealthy()}) {
#		next if ($child->in_the_attic());
#		$chr->patients_categories->{$child->name()} -= $var_sain;
#	}
	foreach my $parent (@{$self->getParents()}) {
		next if ($parent->in_the_attic());
		if ($parent->sex() == 1) {
			$var_sain += $parent->getVariantsVector($chr);
		}
	}
	foreach my $child (@{$self->getChildrenIll()}) {
		next if ($child->in_the_attic());
		if ($child->sex() == 1) { $var_malade += $child->getVariantsVector($chr); }
		else { $var_malade += $child->getHo($chr); }
	}
	foreach my $child (@{$self->getChildrenIll()}) {
		next if ($child->in_the_attic());
		if ($child->sex() == 1) { $var_malade->Intersection( $var_malade, $child->getVariantsVector($chr) ); }
		else { $var_malade->Intersection( $var_malade, $child->getHo($chr) ); }
	}
	$var_malade -= $var_sain;
	my $has_parents;
	foreach my $parent (@{$self->getParents()}) {
		next if ($parent->in_the_attic());
		if ($parent->sex() == 2) { $var_he_parents += $parent->getHe($chr); }
		$has_parents = 1;
	}
	$var_malade->Intersection( $var_malade, $var_he_parents ) if ($has_parents);
	my $has_patient;
	foreach my $patient (@{$self->getPatients()}) {
		next if ($patient->in_the_attic());
		my $v_pat = $patient->getVariantsVector($chr)->Clone();
		$v_pat->Intersection( $v_pat, $var_malade );
		$var_global += $v_pat;
		$has_patient = 1;
	}
#	unless ($is_recessif_compound) {
#		if ($has_parents) { $self->used_model('recessif'); }
#		elsif ($has_patient) { $self->used_model('recessif no parent'); }
#	}
	$self->{vector_transmission}->{recessif}->{$chr->id} = $var_global;
	return $var_global;
}

sub getModelVector_fam_compound() {
	my ($self, $gene) = @_;
	return $self->{vector_transmission}->{compound}->{$gene->id} if (exists $self->{vector_transmission}->{compound}->{$gene->id});
	my $var_res = $gene->getChromosome->getNewVector();
	my $var = $gene->getVariantsVector->Clone();
	$var->Intersection($var, $gene->getChromosome->getVariantsVector());
	return $var_res if ($gene->getChromosome->countThisVariants($var) < 2);
	
	# CAS fam compound sans parents (= gene avec au moins 2 He pour un patient malade dont les var sont jamais vus chez sains)
	if (scalar(@{$self->getParents()}) == 0) {
		my $var_ok = $gene->getChromosome->getNewVector();
		my $var_fam = $gene->getChromosome->getNewVector();
		foreach my $patient (@{$self->getChildrenIll()}) {
			$patient->getVariantsVector($gene->getChromosome);
			$var_fam += $patient->getHe($gene->getChromosome);
		}
		foreach my $patient (@{$self->getChildrenHealthy()}) {
			$var_fam -= $patient->getVariantsVector($gene->getChromosome);
		}
		$var_fam->Intersection($var, $var_fam);
		if ($gene->getChromosome->countThisVariants($var_fam) >= 2) {
			foreach my $patient (@{$self->getChildrenIll()}) {
				$patient->getVariantsVector($gene->getChromosome);
				my $var_pat = $gene->getChromosome->getNewVector();
				$var_pat += $patient->getHe($gene->getChromosome);
				$var_pat->Intersection($var_pat, $var_fam);
				if ($gene->getChromosome->countThisVariants($var_pat) >= 2) {
					$var_ok += $var_pat;
				}
			}
			foreach my $v (@{$gene->getChromosome->getListVarVectorIds( $var_ok )}) {
				push(@{$self->keep_var_compound->{$gene->getChromosome->id()}}, $v);
			}
			$var_res += $var_ok;
		}
		$self->used_model('compound no parent');
		return $var_res;
	}
	
	# get common variants from all ill child 
	my $firs_pat_ill_done;
	my $var_pat_ill = $gene->getChromosome->getNewVector();
	foreach my $patient (@{$self->getChildrenIll()}) {
		my $var_pat = $patient->getHe($gene->getChromosome)->Clone();
		$var_pat->Intersection( $var, $var_pat );
		if ($firs_pat_ill_done) { $var_pat_ill->Intersection($var_pat_ill, $var_pat); }
		else {
			$var_pat_ill = $var_pat;
			$firs_pat_ill_done = 1;
		}
	}
	
	# delete variant seen HO in all ill child
	foreach my $patient (@{$self->getChildrenIll()}) {
		$var_pat_ill -= $patient->getHo($gene->getChromosome);
	}
	
	# delete variants seen HO in each parent
	foreach my $patient (@{$self->getParents()}) {
		$var_pat_ill -= $patient->getHo($gene->getChromosome);
	}
	
	#next if ($gene->getChromosome->countThisVariants($var_pat_ill) < 2);
	
	# get variants seen HE in each parent and seen in ill children too
	my @lVarParents;
	my $var_all_parents = $gene->getChromosome->getNewVector();
	foreach my $patient (@{$self->getParents()}) {
		my $var_this_parent = $patient->getHe($gene->getChromosome)->Clone();
		$var_this_parent->Intersection($var_this_parent, $var_pat_ill);
		push(@lVarParents, $var_this_parent->Clone());
		$var_all_parents += $var_this_parent;
	}

	# delete variants seen HO in each parent
	my $nb_parents = 0;
	foreach my $patient (@{$self->getParents()}) {
		$var_all_parents -= $patient->getHo($gene->getChromosome);
		$nb_parents++;
	}
	return $var_res unless ($nb_parents == 2);

	# delete common variants in parents
#	if ($nb_parents == 2) {
		my $var_both_parents = $gene->getChromosome->getNewVector();
		$var_both_parents->Intersection($lVarParents[0], $lVarParents[1]);
		$lVarParents[0] -= $var_both_parents;
		$lVarParents[1] -= $var_both_parents;
		$var_all_parents -= $var_both_parents;
#	}
#	elsif ($nb_parents == 1) {
#		my $vecTmp = $gene->getVariantsVector->Clone();
#		$vecTmp -= $self->getParents->[0]->getVariantsVector($gene->getChromosome);
#		push(@lVarParents, $vecTmp);
#	}
	return $var_res if ($lVarParents[0]->is_empty());
	return $var_res if ($lVarParents[1]->is_empty());
		
	# combinaison var
	my $nbc = 0;
	my $combinaisons = {};
	for my $v1 (@{$gene->getChromosome->getListVarVectorIds($lVarParents[0])}) {
		for my $v2 (@{$gene->getChromosome->getListVarVectorIds($lVarParents[1])}) {
			$combinaisons->{$nbc} = $gene->getChromosome->getNewVector();
			$combinaisons->{$nbc}->Bit_On($v1);
			$combinaisons->{$nbc}->Bit_On($v2);
			$nbc++;
		}
	}
		
	# cas ou on considere un parent et un denovo
#	if ($hVectorDenovoEachFam) {
#		for my $v1 (@{$gene->getChromosome->getListVarVectorIds($lVarParents[0])}) {
#			for my $v2 (@{$gene->getChromosome->getListVarVectorIds($hVectorDenovoEachFam->{$self->name()})}) {
#				$combinaisons->{$nbc} = $gene->getChromosome->getNewVector();
#				$combinaisons->{$nbc}->Bit_On($v1);
#				$combinaisons->{$nbc}->Bit_On($v2);
#				$nbc++;
#			}
#		}
#		for my $v1 (@{$gene->getChromosome->getListVarVectorIds($lVarParents[1])}) {
#			for my $v2 (@{$gene->getChromosome->getListVarVectorIds($hVectorDenovoEachFam->{$self->name()})}) {
#				$combinaisons->{$nbc} = $gene->getChromosome->getNewVector();
#				$combinaisons->{$nbc}->Bit_On($v1);
#				$combinaisons->{$nbc}->Bit_On($v2);
#				$nbc++;
#			}
#		}
#	}
	# TODO: on obtient un cas tordu avec des couples present aussi chez les enfants sains...
	
	# delete combinaison seen in all healthy child
	foreach my $patient (@{$self->getChildrenHealthy()}) {
		foreach my $nbc (keys %$combinaisons) {
			my $v_comb = $combinaisons->{$nbc}->Clone();
			$v_comb->Intersection($v_comb, $patient->getVariantsVector($gene->getChromosome));
			if ($gene->getChromosome->countThisVariants($v_comb) == 2) {
				delete $combinaisons->{$nbc};
			}
		}
	}
		
	foreach my $nbc (keys %$combinaisons) {
		my $v_comb = $combinaisons->{$nbc}->Clone();
		$v_comb->Intersection($v_comb, $gene->getVariantsVector());
		$v_comb->Intersection($v_comb, $self->getVariantsVector($gene->getChromosome));
		$v_comb->Intersection($v_comb, $gene->getChromosome->getVariantsVector());
		$var_res += $combinaisons->{$nbc};
		foreach my $v (@{$gene->getChromosome->getListVarVectorIds($combinaisons->{$nbc})}) {
			push(@{$self->keep_var_compound->{$gene->getChromosome->id()}}, $v);
		}
	}
	$self->{vector_transmission}->{compound}->{$gene->id} = $var_res;
	return $var_res;
}

sub getModelVector_fam_mosaique() {
	my ($self, $chr) = @_;
	return $self->{vector_transmission}->{mosaic}->{$chr->id} if exists $self->{vector_transmission}->{mosaic}->{$chr->id};
	my $var_mosaique  = $chr->getNewVector();
	my $var_father  = $chr->getNewVector();
	my $var_mother  = $chr->getNewVector();
	my @lParents = @{$self->getParents()};
	my @lChildren_ill = @{$self->getChildrenIll()};
	if (scalar(@lParents) == 0 or scalar(@lChildren_ill) == 0) {
#		foreach my $patient (@{$self->getPatients()}) {
#			$patient->getVariantsVector($chr)->Empty();
#		}
#		$self->used_model('mosaic');
		return $var_mosaique;
	}
	foreach my $patient (@lParents) {
		if ($patient->sex() == '1') { $var_father += $patient->getHe($chr); }
		elsif ($patient->sex() == '2') { $var_mother += $patient->getHe($chr); }
	}
	my $var_to_del  = $chr->getNewVector();
	$var_to_del->Intersection( $var_father, $var_mother );
	foreach my $patient (@lParents) { $var_to_del += $patient->getHo($chr); }
	$var_father -= $var_to_del;
	$var_mother -= $var_to_del;
	foreach my $child (@lChildren_ill) {
		my $var_child = $chr->getNewVector();
		$var_child += $child->getHe($chr);
		foreach my $parent (@lParents) {
			my $var_inter_child_parent = $chr->getNewVector();
			if ($parent->sex() == '1') {
				$var_inter_child_parent->Intersection( $var_child, $var_father );
			}
			elsif ($parent->sex() == '2') {
				$var_inter_child_parent->Intersection( $var_child, $var_mother );
			}
			foreach my $id (@{$self->getListVarVectorIds($var_inter_child_parent)}) {
				$self->project->print_dot(100);
				my $var_id = $chr->getVarId($id);
				my $sample_var_1 = $chr->get_lmdb_variations()->get($var_id)->{patients_details}->{$child->name()}->{he_ho_details};
				unless ($sample_var_1) {
					warn "\n\nERROR: pb with $var_id and patient ".$child->name()."... Die\n\n";
					die ();
				}
				next unless ($parent->getHe($chr)->contains($id));
				my $sample_var_2 = $chr->get_lmdb_variations()->get($var_id)->{patients_details}->{$parent->name()}->{he_ho_details};
				unless ($sample_var_2) {
					warn "\n\nERROR: pb with $var_id (parent: ".$parent->name().")... Die\n\n";
					die ();
				}
				my ($p_fisher, $res_fisher) = $chr->project->buffer->test_fisher($sample_var_2, $sample_var_1, $self->project->mosaique_min_p());
				if ($res_fisher == 1) {
					$var_mosaique->Bit_On($id);
					my ($t,$a,$c)  = split(":", $sample_var_1);
					my ($t1,$b,$d) = split(":", $sample_var_2);
					my $perc1 = ($c / ($a+$c))*100;
					my $perc2 = ($d / ($b+$d))*100;
				}
			}
		}
	}
#	foreach my $patient (@{$self->getPatients()}) {
#		$patient->getVariantsVector($chr)->Intersection( $patient->getVariantsVector($chr), $var_mosaique );
#	}
#	my $model_name = 'mosaic';
#	$self->used_model($model_name);
	$self->{vector_transmission}->{mosaic}->{$chr->id} = $var_mosaique;
	return $var_mosaique;
}

#sub checkModel_compound() {
#	my ($self, $chr, $gene) = @_;
#	
##	warn "\n\n";
##	warn '### GENE '.$gene->external_name();
##	warn "\n\n";
##	warn "Family ".$self->name();
##	warn "\n";
#	
#	my $var_children_healthy = $chr->getNewVector();
#	foreach my $child (@{$self->getChildrenHealthy()}) {
#		$var_children_healthy += $child->getVariantsVector($chr);
#	}
#	
#	my $var_common_parents = $chr->getNewVector();
#	my $var_father = $chr->getNewVector();
#	my $var_mother = $chr->getNewVector();
#	foreach my $parent (@{$self->getParentsHealthy()}) {
#		if ($parent->isFather) { $var_father += $parent->getVariantsVector($chr); }
#		if ($parent->isMother) { $var_mother += $parent->getVariantsVector($chr); }
#	}
#	$var_common_parents->Intersection($var_father, $var_mother);
#	
#	
##	warn "Father: ".$chr->countThisVariants($var_father)."\n";
##	warn "Mother: ".$chr->countThisVariants($var_mother)."\n";
##	warn "Common Parents: ".$chr->countThisVariants($var_common_parents)."\n";
##	warn "Children Healthy: ".$chr->countThisVariants($var_children_healthy)."\n";
#	
##	warn "\n";
#	
#	my $var_ok = $chr->getNewVector();
#	foreach my $child (@{$self->getChildrenIll}) {
#		my $this_var = $chr->getNewVector();
##		warn '# '.$child->name();
#		$this_var += $child->getVectorHe($chr);
##		warn ' -> variants HE: '.$self->countThisVariants($this_var)."\n";
#		$this_var -= $var_children_healthy;
##		warn ' -> variants HE - children healthy: '.$self->countThisVariants($this_var)."\n";
#		$this_var -= $var_common_parents;
##		warn ' -> variants HE - children healthy - common parents: '.$self->countThisVariants($this_var)."\n";
#		
#		$this_var->Intersection($this_var, $gene->getVariantsVector());
##		warn ' -> variants intersect gene: '.$self->countThisVariants($this_var)."\n";
#		
#		next if ($this_var->is_empty());
#		next if ($self->countThisVariants($this_var) < 2);
#		
#		
#		my $common_child_father = $this_var->Clone();
#		$common_child_father->Intersection($common_child_father, $var_father);
#		next if ($common_child_father->is_empty());
#		
#		my $common_child_mother = $this_var->Clone();
#		$common_child_mother->Intersection($common_child_mother, $var_mother);
#		next if ($common_child_mother->is_empty());
#		
#		$var_ok += $this_var;
#		
##		die;
#	}
##	warn "\n\n";
#	return $var_ok;
#}
sub getChildVector() {
	my ($self,$chr) = @_;
	return  $self->getChild()->getVectorOrigin($chr);#  if $self->getFather();;
	return $chr->getNewVector();
}
sub getFatherVector() {
	my ($self,$chr) = @_;
	return  $self->getFather()->getVectorOrigin($chr)  if $self->getFather();;
	return $chr->getNewVector();
}
sub getMotherVector() {
	my ($self,$chr) = @_;
	return  $self->getMother()->getVectorOrigin($chr)  if $self->getMother();;
	return $chr->getNewVector();
}

sub getFatherVectorHe() {
	my ($self,$chr) = @_;
	return  $self->getFather()->getVectorOriginHe($chr)  if $self->getFather();;
	return $chr->getNewVector();
}
sub getMotherVectorHe() {
	my ($self,$chr) = @_;
	return  $self->getMother()->getVectorOriginHe($chr)  if $self->getMother();;
	return $chr->getNewVector();
}


sub getVectorRecessiveTransmission {
	my ($self,$chr,$child) = @_;
#	$chr->getModelVector_fam_strict_denovo;
#	warn Dumper $self->{vector};
#	die();
	unless ($child){
			confess();
			$child =  $self->getChild();
		}
	return 	$chr->getNewVector() unless $self->getMother() and $self->getFather();
	return $self->{vector_transmission_child}->{recessive}->{$chr->id}->{$child->id} if exists $self->{vector_transmission_child}->{recessive}->{$chr->id}->{$child->id};
	my $vector = $child->getVectorHo($chr)->Clone;
	if ($chr->name eq "X"){
		$vector->And($vector, $self->getMother()->getVectorHe($chr))  if $self->getMother();;
		$vector -= $self->getFather()->getVectorOrigin($chr)  if $self->getFather();;
	}
	else {
	$vector->And($vector, $self->getMother()->getVectorHe($chr)) if $self->getMother();
	$vector->And($vector, $self->getFather()->getVectorHe($chr)) if $self->getFather();
	
	}
	 $self->{vector_transmission_child}->{recessive}->{$chr->id}->{$child->id} = $vector;
	return $vector;
}



sub getVectorFatherTransmission {
	my ($self,$chr,$child) = @_;
		return $self->{vector_transmission_child}->{father}->{$chr->id} if exists $self->{vector_transmission_child}->{father}->{$chr->id};
		return 	$chr->getNewVector() unless $self->getFather();
		my $vector_father = $self->getFatherVector($chr);
		my $vector  =  $vector_father->Clone;
		my $vector_mother = $self->getMotherVector($chr);
		$vector -= $vector_mother;
		 $self->{vector_transmission_child}->{father}->{$chr->id} = $vector;
		return $vector;
}

sub getVectorMotherTransmission {
		my ($self,$chr,$child) = @_;
		return $self->{vector_transmission_child}->{mother}->{$chr->id} if exists $self->{vector_transmission_child}->{mother}->{$chr->id};
		return 	$chr->getNewVector() unless $self->getMother();
		my $vector_father = $self->getFatherVector($chr);
		my $vector_mother = $self->getMotherVector($chr);
		my $vector  =  $vector_mother->Clone;;
		$vector -= $vector_father;
		
		 $self->{vector_transmission_child}->{mother}->{$chr->id} = $vector;
		return $vector;
}

sub getVectorBothTransmission {
		my ($self,$chr,$child,$v,$debug) = @_;
		#return $self->{vector_transmission_child}->{both}->{$chr->id}->{$child->id} if exists $self->{vector_transmission_child}->{both}->{$chr->id}->{$child->id};
	
		return 	$chr->getNewVector()  unless $self->getMother() and $self->getFather();
		#my $vector_father = $self->getFatherVector($chr);
		#my $vector_mother = $self->getMother->getVariantsVector($chr);
		my $vector  =  $self->getMotherVector($chr)->Clone;
		
		$vector &= $self->getFatherVector($chr);
		
		my $vector2 = $self->getVectorRecessiveTransmission($chr,$child);
		$vector -= $vector2;
		 $self->{vector_transmission_child}->{both}->{$chr->id}->{$child->id} = $vector;
		#die() if $debug;
		return $self->{vector_transmission_child}->{both}->{$chr->id}->{$child->id};
}


sub getVectorDenovoTransmission {
		my ($self,$chr,$child) = @_;
		return $self->{vector_transmission_child}->{denovo}->{$chr->id}->{$child->id} if exists $self->{vector_transmission_child}->{denovo}->{$chr->id}->{$child->id};
		unless ($child){
			die();
			$child =  $self->getChild();
		}
### je ne sais pas quoi faire
			my $vector = $child->getVectorOrigin($chr)->Clone;

		if ($child->isHealthy()) {
			$self->{vector_transmission}->{denovo}->{$chr->id}->{$child->id} = $chr->getNewVector();
			return $self->{vector_transmission}->{denovo}->{$chr->id}->{$child->id};
		}

			$vector -= $self->getMotherVector($chr);
			$vector -= $self->getFatherVector($chr);
			 $self->{vector_transmission_child}->{denovo}->{$chr->id}->{$child->id} = $vector;
		return $vector
}

sub getVectorStrictDenovoTransmission {
	my ($self,$chr,$child) = @_;
	return $self->{vector_transmission_child}->{strict_denovo}->{$chr->id}->{$child->name()} if exists $self->{vector_transmission_child}->{strict_denovo}->{$chr->id}->{$child->name()};
	my $var_model_null = $chr->getNewVector();
	return $var_model_null unless (-d $chr->dir_sqlite_strict_denovo);
	return $var_model_null unless (-e $chr->dir_sqlite_strict_denovo.'/'.$chr->name().'.lite');
	my $nosql = $chr->sqlite_strict_denovo();
	my $h;
	if ($nosql->exists_db($chr->name())) {
		($h) = $nosql->get_bulk($chr->name());
	}
	$nosql->close();
	return $var_model_null unless ($h);
	#die();
	return $var_model_null if $h->{$child->name()}->is_empty();
	$self->{vector_transmission_child}->{strict_denovo}->{$chr->id} = $h;
	return $self->{vector_transmission_child}->{strict_denovo}->{$chr->id}->{$child->name()};
}

sub getModelVector {
	my ($self,$chr,$model) = @_;
	my $chrid = $chr->id;
	return  $self->{vector}->{$chrid}->{$model} if (exists $self->{vector}->{$chrid}->{$model});
	if ($model eq "recessive"){
			$self->{vector}->{$chrid}->{$model} = $self->checkModel_fam_recessif($chr) 
	}
	elsif ($model eq "denovo"){
		$self->{vector}->{$chrid}->{$model} = $self->checkModel_fam_denovo($chr) 
	}
	else {
		confess("mdodel : $model no code ");
	}

	return $self->{vector}->{$chrid}->{$model};
}

sub getCategoryVariantsVector {
	my ($self, $chr, $cat) = @_;
	my $vector = $chr->getNewVector();
	foreach my $patient (@{$self->getPatients()}) {
		$vector += $patient->getCategoryVariantsVector($chr, $cat);
	}
	return $vector;
}

has hash_regionRec_vector => (
	is 		=> 'rw',
	lazy	=> 1,
	default => sub { {} },
);

sub getVectorRegionRec {
	my ($self, $chr, $nbVar, $filter_regionho_sub_only) = @_;
	confess() unless ($nbVar);
	return $self->hash_regionRec_vector->{$chr->id()}->{$nbVar} if (exists $self->hash_regionRec_vector->{$chr->id()}->{$nbVar});
	#warn ' -> TODO';
	my $hregions;				# table de hash pour stocker les regions Ho : patientname -> chrname = spanHo 

	#  creation de la liste de patients a prendre en compte
	my $patients;
	my $familles;
	#my @attic = split(/,/,$liste_familles_attic); #TODO:?
	my @attic = ();
	$familles= $self->project->getFamilies();
	
	if ( scalar(@attic)==0)		# on cherche les regions Ho de tous les patients de toutes les familles
	{
		$patients = $self->getPatients();
	}
	else	# on ne considere que les familles non in the attic
	{
		#next if ( $liste_familles_attic =~ m/$fam_name/ );
		

		foreach my $pat (@{$self->getPatientsIll()}) 
		{
				push(@{$patients},$pat);
		}
		foreach my $pat (@{$self->getPatientsHealthy()}) 
		{
				push(@{$patients},$pat);
		}
	}
	
	# recherche des regions Ho de chaque patients des familles selectionnees
	foreach my $pat (@{$patients})		
	{
		next if $chr->not_used();
		my $chrname = $chr->name();
		my $pname = $pat->name();
			
    	my $spanHo = Set::IntSpan::Fast::XS->new();			# creation du span  : pour stocker la liste des regions Ho du patient sur le chr
	    
		# creation de nouveaux vecteurs
		my $vecSubHe = Bit::Vector->new($chr->size_vector());
		my $vecSubHo = Bit::Vector->new($chr->size_vector());
		my $vecReg = Bit::Vector->new($chr->size_vector());

		# pour ne tenir compte que des substitutions
		if ($filter_regionho_sub_only) {
			$vecSubHe ->Intersection($pat->getHe($chr), $chr->getVectorSubstitutions());
			$vecSubHo ->Intersection($pat->getHo($chr), $chr->getVectorSubstitutions());
		}
		
		# recupération des régions Ho  étendues aux zones sans variants	
		$vecReg->Complement($vecSubHe);
			
        # pour ne conserver que les régions dont le nombre de variants Ho est > a filtre 
        # et redelimiter la region Ho (au premier / dernier variant Ho rencontre dans la region)
            
        my @H_index = split(/,/,$vecReg->to_Enum());
		
		$spanHo = $pat->getRegionHo_intspan($chr, $nbVar, $filter_regionho_sub_only);
		
		$hregions->{$pname}->{$chrname}=$spanHo;
	} # boucle sur les patients
	
	
	#----------------------------------------------------------------------------------------------
	# regions recessives = Ho chez enfants malades He chez parents et enfants sains
	#----------------------------------------------------------------------------------------------
	
	# gestion des familles : recuperation de liste des enfants et liste des parents
	my $hfam;

	my $fam_name=$self->name();
	
	my $children_ill = $self->getChildrenIll(); 								# enfants malades 
	my $children_healthy = $self->getChildrenHealthy(); 			# enfants sains
	my $parents = $self->getParentsHealthy(); 							# parents sains
	my $others = $self->getPatientsHealthy(); 							# pour les patients malades isoles (ni enfant ni parent)
	
	my $spanIntersChild = Set::IntSpan::Fast::XS->new();							# creation du span  : pour stocker intersection des regionsHo des enfants malades
	my $spanIntersParentsEnfantsSains = Set::IntSpan::Fast::XS->new();						# creation du span  : pour stocker intersection des regionsHe des parents et enfants sains
	my $spanRegionsFam = Set::IntSpan::Fast::XS->new();						# creation du span  : pour stocker les regions Ho pour une famille 
	
	my $spanOthers = Set::IntSpan::Fast::XS->new();						# creation du span  : pour stocker union des regions Ho des patients malades isoles
 
	# boucle sur les chromosomes du projet
	next if $chr->not_used();
	my $chrname=$chr->name();
	
	# familles avec enfants malades et 1 ou 2 parents sains
	if ( (scalar(@{$children_ill}) > 0 )  && ( (scalar(@{$parents}) > 0) || (scalar(@{$children_healthy}) > 0 ) )) 
	{
		$spanIntersChild = $hregions->{@{$children_ill}[0]->name()}->{$chrname};
		for( my $i=1; $i<scalar(@{$children_ill}); $i++ )
		{
			my $pname=@{$children_ill}[$i]->name();
			$spanIntersChild = $spanIntersChild->intersection($hregions->{$pname}->{$chrname});
		}
			
		$spanIntersParentsEnfantsSains = $hregions->{@{$parents}[0]->name()}->{$chrname}->complement();
		$spanIntersParentsEnfantsSains = $spanIntersParentsEnfantsSains->intersection($hregions->{@{$parents}[1]->name()}->{$chrname}->complement()) if(scalar(@{$parents}) == 2);
		for( my $i=1; $i<scalar(@{$children_healthy}); $i++ )
		{
			my $pname=@{$children_healthy}[$i]->name();
			$spanIntersParentsEnfantsSains = $spanIntersParentsEnfantsSains->intersection($hregions->{@{$children_healthy}[$i]->name()}->{$chrname}->complement());
		}
		$spanRegionsFam = $spanIntersChild->intersection($spanIntersParentsEnfantsSains);

	}
	else # individus malades isolés
	{
		if ( (scalar(@{$others}) > 0 ) ) 
		{
			$spanOthers = $hregions->{@{$others}[0]->name()}->{$chrname};
			for( my $i=1; $i<scalar(@{$others}); $i++ )
			{
				my $pname=@{$others}[$i]->name();
				$spanOthers = $spanOthers->union($hregions->{$pname}->{$chrname});
			}
		}
	}
	my @lRegionRec;
	foreach my $coord_vec (split(',',$spanRegionsFam->as_string())) {
		next unless ($coord_vec =~ /-/);
		
		my $v_tmp = Bit::Vector->new_Enum($chr->size_vector(), $coord_vec);
		$v_tmp->Intersection($v_tmp, $self->getHo($chr));
		if ($filter_regionho_sub_only) { $v_tmp->Intersection($v_tmp, $chr->getVectorSubstitutions()); }
		my $nb_var = $self->countThisVariants($v_tmp);
		
		if ($nb_var >= $nbVar){
			push(@lRegionRec, $coord_vec);
		}
	}
	#die;
	my $v_rec;
	if (@lRegionRec) { $v_rec = Bit::Vector->new_Enum($chr->size_vector(), join(',', @lRegionRec)); }
	else { $v_rec = $chr->getNewVector(); }
	$self->{hash_regionRec_vector}->{$chr->id()}->{$nbVar} = $v_rec;
	return $v_rec;
}

has hash_stats_regionRec_vector => (
	is 		=> 'rw',
	lazy	=> 1,
	default => sub { {} },
);


1;