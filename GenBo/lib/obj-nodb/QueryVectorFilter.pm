package QueryVectorFilter;
use strict;
use Vcf;
use Moo;

use Data::Dumper;
use Config::Std;
use Carp;


has verbose_debug => (
	is => 'rw',
); 

sub filter_vector_ratio {
	my ($self, $chr, $limit_ratio, $filter_type_ratio) = @_;
	
	
	return unless ($limit_ratio);
	return if ($limit_ratio eq 'all');
	my $vector_ok = $chr->getNewVector();
	foreach my $patient (@{$chr->getPatients()}) {
		next if ($patient->in_the_attic());
		my $vector_ratio_name = "ratio_".$limit_ratio;
		if ($filter_type_ratio eq 'max') { $vector_ratio_name = 'lower_'.$vector_ratio_name; }
		my $vquality = $patient->_getRocksVector($chr, $vector_ratio_name);;
		$vector_ok += $vquality;
	}
	$vector_ok->Intersection($vector_ok, $chr->getVariantsVector());
	$chr->setVariantsVector($vector_ok);
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_ratio_min - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
}

sub filter_vector_ncboost {
	my ($self, $chr, $cat) = @_;
	return unless ($cat);
	if ($chr->id() eq 'MT') {
		$chr->setVariantsVector( $chr->getNewVector() );
		return;
	}
	$chr->project->print_dot(1);
	my $vector_ncboost = $chr->getVectorNcboost($cat);
	$vector_ncboost->Intersection($vector_ncboost, $chr->getVariantsVector());
	$chr->setVariantsVector($vector_ncboost);
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_ncboost - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	return;
}

sub filter_vector_hgmd_dm {
	my ($self, $chr) = @_;
	confess() unless ($chr->project->isUpdate());
	my $vector_hgmd_dm = $chr->getVectorLmdbDm();
	$vector_hgmd_dm->Intersection($vector_hgmd_dm, $chr->getVariantsVector());
	$chr->setVariantsVector($vector_hgmd_dm);
}

sub filter_vector_hgmd_dm_new_version {
	my ($self, $chr) = @_;
	confess unless ($chr->project->isUpdate());
	my $vector_hgmd_dm_new = $chr->getVectorLmdbDm_newVersion();
	$vector_hgmd_dm_new->Intersection($vector_hgmd_dm_new, $chr->getVariantsVector());
	$chr->setVariantsVector($vector_hgmd_dm_new);
}

sub filter_vector_region_ho {
	my ($self, $chr, $filter_nbvar_regionho, $filter_regionho_sub_only, $indiv_or_fam) = @_;
	return unless ($filter_nbvar_regionho);
	return filter_vector_region_ho_individual($chr, $filter_nbvar_regionho, $filter_regionho_sub_only) if ($indiv_or_fam eq 'individual');
	return filter_vector_region_ho_familial($chr, $filter_nbvar_regionho, $filter_regionho_sub_only) if ($indiv_or_fam eq 'familial');
}

sub filter_vector_region_ho_individual {
	my ($self, $chr, $filter_nbvar_regionho, $filter_regionho_sub_only) = @_;
	my $vector_all = $chr->getNewVector();
	my $vectorRegionHo_global = $chr->getNewVector();
	foreach my $patient (@{$chr->getPatients()}) {
		$chr->getProject->print_dot(1);
		$vectorRegionHo_global += $patient->getRegionHo($chr, $filter_nbvar_regionho, $filter_regionho_sub_only);
	}
	foreach my $patient (@{$chr->getPatients()}) {
		$chr->getProject->print_dot(1);
		$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $vectorRegionHo_global);
	}
	$vectorRegionHo_global->Intersection($vectorRegionHo_global, $chr->getVariantsVector());
	$chr->setVariantsVector($vectorRegionHo_global);
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_region_ho - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	return;
}

sub filter_vector_region_ho_familial {
	my ($self, $chr, $filter_nbvar_regionho, $filter_regionho_sub_only) = @_;
	my $vector_all = $chr->getNewVector();
	my $vectorRegionRec_global = $chr->getNewVector();
	foreach my $family (@{$chr->getFamilies()}) {
		$chr->getProject->print_dot(1);
		if ($family->getVectorRegionRec($chr, $filter_nbvar_regionho, $filter_regionho_sub_only)->is_empty()) {
			foreach my $patient (@{$family->getPatients()}){
				$patient->getVariantsVector($chr)->Empty();
			}
		}
		$vectorRegionRec_global += $family->getVectorRegionRec($chr, $filter_nbvar_regionho, $filter_regionho_sub_only);
	}
	foreach my $patient (@{$chr->getPatients()}) {
		$chr->getProject->print_dot(1);
		$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $vectorRegionRec_global);
	}
	$chr->setVariantsVector($vectorRegionRec_global);
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_region_ho - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	return;
}

sub filter_vector_type_variants {
	my ($self, $chr, $hToFilter) = @_;
	return unless ($hToFilter);
	my $vector = $chr->getVariantsVector->Clone();
	$vector->AndNot($vector, $chr->get_vector_type_variants($hToFilter));
	$chr->setVariantsVector($vector);
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_type_variants - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	return;
}

sub filter_vector_predicted_splice_region_global {
	my ($self, $chr, $hToFilter) = @_;
	return unless ($hToFilter);
	return unless (exists $hToFilter->{predicted_splice_site});
	my $vector = $chr->getVariantsVector->Clone();
	$vector->AndNot($vector, $chr->get_vector_predicted_splice_region($hToFilter));
	$chr->setVariantsVector($vector);
	return;
}

sub filter_vector_gnomad_ac {
	my ($self, $chr, $filter_name) = @_;
	return if not $filter_name;
	return if not $filter_name =~ /gnomad/;
	my $no = $chr->flush_rocks_vector();
	my $vector_ok = $chr->getVariantsVector() & $no->get_vector_chromosome($filter_name);
	$chr->setVariantsVector($vector_ok);
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_gnomad_ac - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
}

sub filter_vector_freq_variants {
	my ($self, $chr, $hToFilter) = @_;
	return unless ($hToFilter);
	my $vector = $chr->getVariantsVector->Clone();
	$vector->Intersection($vector, $chr->get_vector_freq_variants($hToFilter));
	$chr->setVariantsVector($vector);
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_freq_variants - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	return;
}

sub filter_vector_gnomad_ho_ac_variants {
	my ($self, $chr, $hToFilter) = @_;
	return unless ($hToFilter);
	my $vector = $chr->getVariantsVector->Clone();
	if (exists $hToFilter->{gnomad_ac_all})    { $vector->Intersection($vector, $chr->get_vector_gnomad_ac_variants($hToFilter)); }
	if (exists $hToFilter->{gnomad_ho_ac_all}) { $vector->Intersection($vector, $chr->get_vector_gnomad_ho_ac_variants($hToFilter)); }
	$chr->setVariantsVector($vector);
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_gnomad_ho_ac_variants - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	return;
}

sub filter_vector_confidence_variants {
	my ($self, $chr, $hToFilter) = @_;
	return unless ($hToFilter);
	my $vector = $chr->getVariantsVector->Clone();
	$vector->AndNot($vector, $chr->get_vector_confidence_variants($hToFilter));
	$chr->setVariantsVector($vector);
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_confidence_variants - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	return;
}

sub filter_vector_cadd_variants {
	my ($self, $chr, $hToFilter, $keep_indels_cadd) = @_;
	return unless ($hToFilter);
	my $can_do;
	foreach my $filter_name (keys %{$chr->project->hash_cadd_filters()}) {
		$can_do = 1 if (exists $hToFilter->{$filter_name});
	}
	return unless ($can_do);
	my $vector_init = $chr->getVariantsVector->Clone();
	my $vector = $chr->getVariantsVector->Clone();
	$vector->AndNot($vector, $chr->get_vector_cadd_variants($hToFilter));
	if ($keep_indels_cadd) {
		my $vector_to_keep = $chr->getNewVector();
		$vector_to_keep += $chr->getVectorInsertions() unless (exists $hToFilter->{insertion});
		$vector_to_keep += $chr->getVectorDeletions() unless (exists $hToFilter->{deletion});
		$vector_to_keep += $chr->getVectorLargeDeletions() unless (exists $hToFilter->{large_deletion});
		$vector_to_keep += $chr->getVectorLargeDuplications() unless (exists $hToFilter->{large_duplication});
		$vector_to_keep->Intersection($vector_to_keep, $vector_init);
		$vector += $vector_to_keep;
	}
	$chr->setVariantsVector($vector);
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_cadd_variants - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	return;
}

sub filter_vector_polyscore {
	my ($self, $chr, $polyscore_value) = @_;
	my $min_value;
	$min_value = 10 if ($polyscore_value eq 'high');
	$min_value = 7 if ($polyscore_value eq 'moderate');
	$min_value = 5 if ($polyscore_value eq 'low');
	my $vector = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		next if ($gene->is_intergenic());
		foreach my $var (@{$gene->getVariants()}) {
			foreach my $family (@{$chr->getProject->getFamilies()}) {
				foreach my $patient (@{$family->getChildrenIll()}) {
					next unless ($patient->getVariantsVector($chr)->contains($var->vector_id()));
					eval {
						my $this_score = $var->scaledScoreVariant($gene, $patient);
						$vector->Bit_On($var->vector_id()) if ($this_score >= $min_value);
					};
					if ($@) {
						next;
					}
				}
			}
		}
	}
	$chr->setVariantsVector($vector);
}

sub filter_vector_dejavu {
	my ($self, $chr, $nb_dejavu, $dejavu_ho, $test) = @_;
	return unless ($nb_dejavu);
	my $vector = $chr->getVariantsVector->Clone();
	my $no = $chr->flush_rocks_vector();
	if ($dejavu_ho and $nb_dejavu eq 'uniq'){
		my $v1 = $no->get_vector_chromosome("shodv_0");
		$vector &= $v1 ;
		foreach my $id (@{$chr->getListVarVectorIds($vector)}) {
			$chr->project->print_dot(200);
			my $v = $chr->getProject->returnVariants($chr->name."!".$id);
			$vector->Bit_Off($id) if ($v->other_patients_ho() > 0);
		}
	}
	elsif ($dejavu_ho) {
		my $v1 = $no->get_vector_chromosome("shodv_".$nb_dejavu);
		$vector &= $v1 ;
		
	}
	elsif ($nb_dejavu eq 'uniq') {
		my $v0 = $no->get_vector_chromosome("sdv_0");
		$vector &= $v0 ;
		foreach my $id (@{$chr->getListVarVectorIds($vector)}) {
			$chr->project->print_dot(200);
			my $v = $chr->getProject->returnVariants($chr->name."!".$id);
			$vector->Bit_Off($id) if ($v->other_patients() > 0);
		}

	}
	elsif (not $nb_dejavu eq 'all') {
		my $v1 = $no->get_vector_chromosome("sdv_".$nb_dejavu);
		$vector &= $v1 ;
	}
	$chr->setVariantsVector($vector);
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_dejavu - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	return;
}

# AT_LEAST par patients en mode VAR et IND 
sub atLeastFilter_var_ind {
	my ($self, $chr, $atLeast) = @_;
	my (@lOk, $hVectPat);
	my @lPat;
	foreach my $patient (@{$chr->getPatients()}) {
		next if ($patient->in_the_attic());
		push(@lPat, $patient);
	}
	if ($atLeast <= scalar(@lPat)) {
		foreach my $patient (@lPat) {
			$hVectPat->{$patient->name()} = $patient->getVariantsVector($chr);
		}
		foreach my $index (@{$chr->getIdsBitOn($chr->getVariantsVector())}) {
			$chr->project->print_dot(50);
			my $ok = 0;
			foreach my $patient (@lPat) {
				if ($hVectPat->{$patient->name()}->contains($index)) {
					$ok++;
					if ($ok == $atLeast) {
						push(@lOk, $index);
						last;
					}
				}
			}
		}
	}
	my $vector_ok = Bit::Vector->new_Enum($chr->getVariantsVector->Size(), join(',', @lOk));
	$chr->setVariantsVector($vector_ok);
}

# AT_LEAST par familles en mode VAR et FAM
sub atLeastFilter_var_fam {
	my ($self, $chr, $atLeast) = @_;
	my ( @lOk, $hVectFam);
	my @lFam;
	foreach my $family (@{$chr->getFamilies()}) {
		next if ($family->in_the_attic());
		push(@lFam, $family);
	}
	if ($atLeast <= scalar(@lFam)) {
		foreach my $family (@lFam) {
			$hVectFam->{$family->name()} = $chr->getNewVector();
			foreach my $patient (@{$family->getPatients()}) {
				next if ($patient->in_the_attic());
				$hVectFam->{$family->name()} += $patient->getVariantsVector($chr);
			}
		}
		foreach my $index (@{$chr->getIdsBitOn($chr->getVariantsVector())}) {
			$chr->project->print_dot(50);
			my $ok = 0;
			foreach my $fam_name (keys %$hVectFam) {
				if ($hVectFam->{$fam_name}->contains($index)) {
					$ok++;
					if ($ok == $atLeast) {
						push(@lOk, $index);
						last;
					}
				}
			}
		}
	}
	my $vector_ok = Bit::Vector->new_Enum($chr->getVariantsVector->Size(), join(',', @lOk));
	$chr->setVariantsVector($vector_ok);
}

sub atLeastFilter_var_som {
	my ($self, $chr, $atLeast) = @_;
	my (@lOk, $hVectFam);
	my @lGroups;
	foreach my $group (@{$chr->getSomaticGroups()}) {
		next if ($group->in_the_attic());
		push(@lGroups, $group);
	}
	if ($atLeast <= scalar(@lGroups)) {
		foreach my $group (@lGroups) {
			$hVectFam->{$group->name()} = $chr->getNewVector();
			foreach my $patient (@{$group->getPatients()}) {
				next if ($patient->in_the_attic());
				$hVectFam->{$group->name()} += $patient->getVariantsVector($chr);
			}
		}
		foreach my $index (@{$chr->getIdsBitOn($chr->getVariantsVector())}) {
			$chr->project->print_dot(50);
			my $ok = 0;
			foreach my $fam_name (keys %$hVectFam) {
				if ($hVectFam->{$fam_name}->contains($index)) {
					$ok++;
					if ($ok == $atLeast) {
						push(@lOk, $index);
						last;
					}
				}
			}
		}
	}
	my $vector_ok = Bit::Vector->new_Enum($chr->getVariantsVector->Size(), join(',', @lOk));
	$chr->setVariantsVector($vector_ok);
}

# AT_LEAST par patients en mode GENES et IND
sub atLeastFilter_genes_ind {
	my ($self, $chr, $atLeast) = @_;
	my $variants_genes  = $chr->getNewVector();
	my $var_tmp_atleast = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		my $nb_ok = 0;
		foreach my $patient (@{$chr->getPatients()}) {
			$var_tmp_atleast->Intersection($gene->getVariantsVector(), $patient->getVariantsVector($chr));
			unless ($var_tmp_atleast->is_empty()) {
				$nb_ok++;
				last if ($nb_ok == $atLeast);
			}
		}
		if ($nb_ok < $atLeast) { $chr->supressGene($gene);  }
		else { $variants_genes += $gene->getVariantsVector(); }
	}
	$chr->getVariantsVector->Intersection($chr->getVariantsVector(), $variants_genes);
	foreach my $patient (@{$chr->getPatients()}) {
		$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $chr->getVariantsVector());
	}
}

# AT_LEAST par familles en mode GENES et FAM
sub atLeastFilter_genes_fam {
	my ($self, $chr, $atLeast) = @_;
	my $variants_genes  = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		my $nb_ok = 0;
		my $var_tmp_atleast = $chr->getNewVector();
		foreach my $family (@{$chr->getFamilies}) {
			$var_tmp_atleast->Intersection($gene->getVariantsVector(), $family->getVariantsVector($chr));
			unless ($var_tmp_atleast->is_empty()) {
				$nb_ok++;
				last if ($nb_ok == $atLeast);
			}
		}
		if ($nb_ok < $atLeast) { $chr->supressGene($gene); }
		else { $variants_genes += $gene->getVariantsVector(); }
	}
	$chr->getVariantsVector->Intersection($chr->getVariantsVector(), $variants_genes);
	foreach my $patient (@{$chr->getPatients()}) {
		$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $chr->getVariantsVector());
	}
}

sub atLeastFilter_genes_som{
	my ($self, $chr, $atLeast) = @_;
	my $variants_genes  = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		my $nb_ok = 0;
		my $var_tmp_atleast = $chr->getNewVector();
		foreach my $group (@{$chr->getSomaticGroups}) {
			my $vector_gr = $chr->getNewVector();
			foreach my $patient (@{$group->getPatients()}) { $vector_gr += $patient->getVariantsVector($chr); }
			$var_tmp_atleast->Intersection($gene->getVariantsVector(), $vector_gr);
			unless ($var_tmp_atleast->is_empty()) {
				$nb_ok++;
				last if ($nb_ok == $atLeast);
			}
		}
		if ($nb_ok < $atLeast) { $chr->supressGene($gene); }
		else { $variants_genes += $gene->getVariantsVector(); }
	}
	$chr->getVariantsVector->Intersection($chr->getVariantsVector(), $variants_genes);
	foreach my $patient (@{$chr->getPatients()}) {
		$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $chr->getVariantsVector());
	}
}


######### FILTRES MODELS CHR ##########

our $hModelsMethodsNames_familial = {
   'denovo' => 'isDenovoTransmission',
   'strict-denovo' => 'isStrictDenovoTransmission',
   'dominant' => 'isDominantTransmission',
   'mosaic' => 'isMosaicTransmission',
   'recessif' => 'isRecessiveTransmission',
#   'incomplete_penetrance' => 'getModelVector_fam_incomplete_penetrance',
   'uniparental_disomy' => 'isUniparentalDisomyTransmission',
   
   'only_tissues_somatic' => 'getModelVector_som_only_tissues_somatic',
   
   
   'recessif_vector' => 'getVector_individual_recessive',
   'dominant_vector' => 'getVector_individual_dominant',
   'denovo_vector' => 'getVector_individual_denovo',
   'uniparental_disomy_vector' => 'getVector_individual_uniparental_disomy',
};

sub filter_models_familial_union {
	my ($self, $chr, $list_models, $h_arguments) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_models_familial_union need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	confess("\n\nERROR: QueryVectorFilter::filter_models_familial_union need a list of models. Die.\n\n") unless ($list_models);
	my $vector_models = $self->get_models_familial_union($chr, $list_models, $h_arguments);
#	foreach my $patient (@{$chr->getPatients()}) {
#		my $vector_patient = $chr->getNewVector();
#		foreach my $model (@$list_models) {
#			confess("\n\nERROR: $model model not used for patient ".$patient->name()." for CHR".$chr->id()."... Die.\n\n") unless (exists $patient->hash_models_genetics_used->{$model}->{$chr->id()});
#			$vector_patient += $patient->hash_models_genetics_used->{$model}->{$chr->id()};
#		}
#		$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $vector_patient);
#	}
	$chr->setVariantsVector($vector_models);
}

sub get_models_familial_union {
	my ($self, $chr, $list_models, $h_arguments) = @_;
	confess("\n\nERROR: QueryVectorFilter::get_models_familial_union need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	confess("\n\nERROR: QueryVectorFilter::get_models_familial_union need a list of models. Die.\n\n") unless ($list_models);
	my $vector_models = $chr->getNewVector();
	foreach my $model (@$list_models) {
		#isMosaicTransmission
		if ($model eq 'compound') {
			$vector_models += $self->get_model_familial_compound($chr, $h_arguments);
		}
		else {
			my $model_method_name = $hModelsMethodsNames_familial->{$model};
			my $model_method_vector_name;
			$model_method_vector_name = $hModelsMethodsNames_familial->{$model.'_vector'} if exists $hModelsMethodsNames_familial->{$model.'_vector'};
			foreach my $family (@{$chr->getProject->getFamilies()}) {
				next if $family->in_the_attic();
				if ($model eq 'dominant') { $family->{isDominant} = 1; }
				foreach my $child (@{$family->getChildrenIll()}) {
					next if $child->in_the_attic();
					my $var_pat = $child->getVectorOrigin($chr);
					$var_pat->Intersection($var_pat, $chr->getVariantsVector());
					if ($model_method_vector_name) {
						$var_pat->Intersection($var_pat, $family->$model_method_vector_name($chr, $child));
						$vector_models += $var_pat;
					}
					else {
						foreach my $var (@{$chr->getListVarObjects($var_pat)}) {
							$vector_models->Bit_On($var->vector_id()) if $var->$model_method_name($family, $child);
						}
					}
				}
			}
		}
	}
	$chr->setVariantsVector($vector_models);
	return $vector_models;
}

sub filter_model_familial_denovo {
	my ($self, $chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_denovo need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return filter_model_familial_common($chr, 'denovo');
}

sub get_model_familial_denovo {
	my ($self, $chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::get_model_familial_denovo need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return get_model_familial_common($chr, 'denovo');
}

sub filter_model_familial_strict_denovo {
	my ($self, $chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_strict_denovo need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return filter_model_familial_common($chr, 'strict-denovo');
}
		
sub get_model_familial_strict_denovo {
	my ($self, $chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::get_model_familial_strict_denovo need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return get_model_familial_common($chr, 'strict-denovo');
}

sub filter_model_familial_dominant {
	my ($self, $chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_dominant need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return filter_model_familial_common($chr, 'dominant');
}
		
sub get_model_familial_dominant {
	my ($self, $chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::get_model_familial_dominant need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return get_model_familial_common($chr, 'dominant');
}

sub filter_model_familial_mosaic {
	my ($self, $chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_mosaic need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return filter_model_familial_common($chr, 'mosaic');
}
		
sub get_model_familial_mosaic {
	my ($self, $chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::get_model_familial_mosaic need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return get_model_familial_common($chr, 'mosaic');
}

sub filter_model_familial_recessif {
	my ($self, $chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_recessif need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return filter_model_familial_common($chr, 'recessif');
}
		
sub get_model_familial_recessif {
	my ($self, $chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::get_model_familial_recessif need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return get_model_familial_common($chr, 'recessif');
}

sub filter_model_familial_uniparental_disomy {
	my ($self, $chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_uniparental_disomy need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return filter_model_familial_common($chr, 'uniparental_disomy');
}
		
sub get_model_familial_uniparental_disomy {
	my ($self, $chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::get_model_familial_uniparental_disomy need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return get_model_familial_common($chr, 'uniparental_disomy');
}

sub filter_model_familial_compound {
	my ($self, $chr, $h_arguments) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_compound need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_compound need a Hash of Arguments. Die.\n\n") unless ($h_arguments);
	return filter_model_familial_common($chr, 'compound');
}

sub get_model_familial_compound {
	my ($self, $chr, $h_arguments) = @_;
	my $variants_genes  = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		my $vector_gene = $self->filter_gene_model_familial_compound($gene, $h_arguments);
		$variants_genes += $vector_gene;
	}
	$variants_genes->Intersection($chr->getVariantsVector(), $variants_genes);
	my $vector_chr = $chr->getNewVector();
	foreach my $family (@{$chr->getFamilies()}) {
		my $vector_fam = $family->get_vector_keep_var_compound($chr)->Clone();
		$vector_fam->Intersection($vector_fam, $variants_genes);
		foreach my $patient (@{$family->getPatients}) {
			my $vector_patient = $patient->getVariantsVector($chr)->Clone();
			$vector_patient->Intersection($vector_patient, $vector_fam);
			$vector_chr += $vector_patient;
			$patient->{hash_models_genetics_used}->{'compound'}->{$chr->id()} = $vector_patient->Clone();
		}
		$family->{hash_models_genetics_used}->{'compound'}->{$chr->id()} = $vector_fam->Clone();
	}
	return $vector_chr;
}

sub get_gene_model_familial_compound_prepare_multi_annot {
	my ($self, $gene, $h_arguments) = @_;
	my $hFiltersChr = $h_arguments->{'filters_1'};
	my $hFiltersChr_var2 = $h_arguments->{'filters_2'};
	my $vector_filtered = $h_arguments->{'vector_filters_1'};
	my $vector_filtered_2 = $h_arguments->{'vector_filters_2'};
	my $var_gene_annot = $gene->getChromosome->getNewVector();
	my $var_gene_annot_1 = $self->get_vector_filter_gene_annotations($gene, $hFiltersChr);
	$var_gene_annot_1->Intersection($var_gene_annot_1, $vector_filtered) if ($vector_filtered);
	if ($var_gene_annot_1->is_empty()) { return; }
	$var_gene_annot += $var_gene_annot_1;
	my $var_gene_annot_2 = $self->get_vector_filter_gene_annotations($gene, $hFiltersChr_var2);
	$var_gene_annot_2->Intersection($var_gene_annot_2, $vector_filtered_2) if ($vector_filtered_2);
	if ($var_gene_annot_2->is_empty()) { return; }
	$var_gene_annot += $var_gene_annot_2;
	return ($var_gene_annot_1, $var_gene_annot_2);
}

sub filter_gene_model_familial_compound {
	my ($self, $gene,$h_arguments) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_gene_model_familial_compound need a GenBoGeneCache. Die.\n\n") unless ($gene);
	my $vector_compound = $gene->getChromosome->getNewVector;
	return $vector_compound if not $gene->getPatients();
	foreach my $fam (@{$gene->getFamilies()}) {
		my $vector_fam_gene_compound = $fam->getModelVector_fam_compound($gene);
#		if ($h_arguments->{'filters_2'}) {
#			my ($self, $var_gene_annot_1, $var_gene_annot_2) = $self->get_gene_model_familial_compound_prepare_multi_annot($gene, $h_arguments);
#			if ($var_gene_annot_1 and $var_gene_annot_2) {
#				$var_gene_annot_1->Intersection($var_gene_annot_1, $vector_fam_gene_compound);
#				$vector_fam_gene_compound->Empty() if ($var_gene_annot_1->is_empty());
#				$var_gene_annot_2->Intersection($var_gene_annot_2, $vector_fam_gene_compound);
#				$vector_fam_gene_compound->Empty() if ($var_gene_annot_2->is_empty());
#				$vector_fam_gene_compound->Empty() if ($gene->getChromosome->countThisVariants($vector_fam_gene_compound) < 2);
#			}
#			else {
#				$vector_fam_gene_compound->Empty();
#			}
#		}
		$vector_compound += $vector_fam_gene_compound;
		
	}
	return $vector_compound;
	
#	
#	
#	warn ' -> '.$gene->getVectorOrigin->Norm();
#	my $vector_gene = $gene->getVectorOrigin() & $gene->getChromosome->getVariantsVector();
#	warn ' -> '.$vector_gene->Norm();
#	
#	
#	
#	
#	if ($vector_gene->Norm() < 2) {
#		$vector_gene->Empty();
#		return $vector_gene;
#	}
#	
#	warn  'before comp -> '.$vector_gene->Norm();
#	$vector_gene = $vector_gene & $gene->getModelGeneVector_fam_compound();
#	warn ' after comp -> '.$vector_gene->Norm();
#	
#	my $vector_ok = $gene->getChromosome->getNewVector();
#	foreach my $family (@{$gene->getChromosome->getFamilies()}) {
#		my $vector_gene_tmp = $vector_gene & $family->get_vector_keep_var_compound($gene->getChromosome());
#		$vector_ok += $vector_gene_tmp;
#	}
#	$vector_gene->Intersection($vector_gene, $vector_ok);
#	
#	
#	warn ' after keep fam -> '.$vector_gene->Norm();
#	
#	die if $vector_gene > 2;
#	
#	if ($vector_gene->is_empty()) { return $vector_gene; }
#	elsif ($hFiltersChr_var2) {
#		my ($self, $var_gene_annot_1, $var_gene_annot_2) = $self->get_gene_model_familial_compound_prepare_multi_annot($gene, $hFiltersChr, $hFiltersChr_var2, $vector_filtered, $vector_filtered_2);
#		if ($var_gene_annot_1 and $var_gene_annot_2) {
#			
#			warn $var_gene_annot_1->Norm();
#			warn $var_gene_annot_2->Norm();
#			die;
#			
#			$var_gene_annot_1->Intersection($var_gene_annot_1, $vector_gene);
#			if ($var_gene_annot_1->is_empty()) {
#				$vector_gene->Empty();
#			}
#			$var_gene_annot_2->Intersection($var_gene_annot_2, $vector_gene);
#			if ($var_gene_annot_2->is_empty()) {
#				$vector_gene->Empty();
#			}
#			if ($gene->getChromosome->countThisVariants($vector_gene) < 2) {
#				$vector_gene->Empty();
#			}
#		}
#		else {
#			$vector_gene->Empty();
#		}
#	}
#	return $vector_gene;
}

sub filter_model_familial_common {
	my ($self, $chr, $model) = @_;
	return unless ($model);
	my $vector_model = get_model_familial_common($chr, $model);
	my $vector_chr = $chr->getNewVector();
	foreach my $patient (@{$chr->getPatients()}) {
		confess("\n\nERROR: $model model not used for patient ".$patient->name()." for CHR".$chr->id()."... Die.\n\n") unless (exists $patient->hash_models_genetics_used->{$model}->{$chr->id()});
		my $vector_patient_model = $patient->hash_models_genetics_used->{$model}->{$chr->id()};
		$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $vector_patient_model);
		$vector_chr += $patient->getVariantsVector($chr);
	}
	$chr->setVariantsVector($vector_chr);
}

sub get_model_familial_common {
	my ($self, $chr, $model, $h_arguments) = @_;
	return unless ($model);
	my $model_method_name = $hModelsMethodsNames_familial->{$model};
	return unless ($model_method_name);
	my $vector_update = $chr->getNewVector();
	foreach my $family (@{$chr->getFamilies()}) {
		my ($self, $has_patients, $has_parents); 
		$family->{hash_models_genetics_used}->{$model}->{$chr->id()} = $family->$model_method_name($chr, $h_arguments);
		foreach my $patient (@{$family->getPatients()}) {
			my $vector_patient = $patient->getVariantsVector($chr)->Clone();
			$vector_patient->Intersection( $vector_patient, $family->{hash_models_genetics_used}->{$model}->{$chr->id()});
			$vector_patient->Intersection($vector_patient, $chr->getVariantsVector());
			$vector_update += $vector_patient;
			$patient->{hash_models_genetics_used}->{$model}->{$chr->id()} = $vector_patient->Clone();
			$has_patients = 1;
		}
	}
	return $vector_update;
}

sub filter_model_somatic_loh {
	my ($self, $chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_model_somatic_loh need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	$chr->setVariantsVector($chr->getModelVector_som_loh());
	return;
	
}

sub filter_model_somatic_only_tissues_somatic {
	my ($self, $chr) = @_;
	return $self->filter_model_familial_common($chr, 'only_tissues_somatic');
}

sub filter_model_individual_recessif {
	my ($self, $chr) = @_;
	my $variants_genes  = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		next if $gene->getCurrentVector->is_empty();
		$self->filter_gene_model_individual_recessif($gene);
		next if ($gene->getCurrentVector->is_empty());
		$variants_genes += $gene->getCurrentVector();
	}
	$chr->setVariantsVector($variants_genes);
	foreach my $patient (@{$chr->getPatients()}) {
		$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $patient->{model_indiv_recessif}->{$chr->id}) if exists ($patient->{model_indiv_recessif}->{$chr->id});
	}
}

sub filter_model_individual_compound {
	my ($self, $chr) = @_;
	my $variants_genes  = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		next if $gene->getCurrentVector->is_empty();
		$self->filter_gene_model_individual_compound($gene);
		
		next if ($gene->getCurrentVector->is_empty());
		
#		if ($gene->external_name eq 'UBIAD1') {
#			warn "\n\n";
#			warn $gene->external_name();
#			warn 'All: '.$gene->getCurrentVector->to_Enum();
#			warn 'NB: '.$chr->countThisVariants($gene->getCurrentVector());
#			die;
#		}
#		
#		
#		my $nb = $chr->countThisVariants($gene->getCurrentVector());
#		if ($nb < 2) {
#			$gene->getCurrentVector->Empty();
#			next;
#		}
#		
		$variants_genes += $gene->getCurrentVector();
	}
	$chr->setVariantsVector($variants_genes);
	#$chr->getVariantsVector->Intersection($chr->getVariantsVector(), $variants_genes);
	foreach my $patient (@{$chr->getPatients()}) {
		$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $patient->{model_indiv_compound}->{$chr->id}) if exists ($patient->{model_indiv_compound}->{$chr->id});
	}
	return;
}


sub filter_model_somatic_dbl_evt {
	my ($self, $chr) = @_;
	my $variants_genes  = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		filter_gene_model_somatic_dbl_evt($gene);
		if ($gene->getVariantsVector->is_empty()) {
			$chr->supressGene($gene);
			next;
		}
		$variants_genes += $gene->getVariantsVector();
	}
	$chr->getVariantsVector->Intersection($chr->getVariantsVector(), $variants_genes);
	my $var = $chr->getNewVector();
	foreach my $patient (@{$chr->getPatients()}) {
		if ($patient->model_vector_var_ok->{$chr->id()}) {
			$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $patient->model_vector_var_ok->{$chr->id()});
			$var += $patient->model_vector_var_ok->{$chr->id()};
			$chr->patients_categories->{$patient->name()} = $patient->model_vector_var_ok->{$chr->id()};
		}
		else {
			$patient->getVariantsVector($chr)->Empty();
		}
	}
	$chr->setVariantsVector($var);
	return;
}



########## FILTRES MODELS GENES ##########



sub filter_genes_from_ids {
	my ($self, $chr, $hGeneIds, $can_use_hgmd) = @_;
	my @lOk;
	foreach my $hash (@{$chr->values_lmdb_genes()}) {
		next unless (exists $hGeneIds->{$hash->{name}});
		push(@lOk, $hash);
	}
	$chr->{values_lmdb_genes} = \@lOk;
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_genes_from_ids - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	return;
}


sub filter_gene_model_individual_recessif {
	my ($self, $gene) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_gene_model_individual_recessif need a GenBoGeneCache. Die.\n\n") unless ($gene);
	my $v = $gene->getCurrentVector;
	$v->Intersection($gene->getCurrentVector(), $gene->getChromosome->getModelGeneVector_indiv_recessif($gene));
	$gene->setCurrentVector($v);
	return;
}

sub filter_gene_model_individual_compound {
	my ($self, $gene) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_gene_model_individual_compound need a GenBoGeneCache. Die.\n\n") unless ($gene);
	my $v = $gene->getCurrentVector;
	$v->Intersection($gene->getCurrentVector(), $gene->getChromosome->getModelGeneVector_indiv_compound($gene));
	$gene->setCurrentVector($v);
	return;
}

sub filter_gene_model_somatic_dbl_evt {
	my ($self, $gene) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_gene_model_somatic_dbl_evt need a GenBoGeneCache. Die.\n\n") unless ($gene);
	my $v = $gene->getCurrentVector;
	$v->Intersection($gene->getCurrentVector(), $gene->getChromosome->getModelGeneVector_som_dbl_evt($gene));
	$gene->setCurrentVector($v);
	return;
}



######### FILTRES PATIENTS / FAMILIES ##########

my $has_intersected_or_excuded_patients;

# method in the attic pour un patient 
sub setInTheAtticPatients {
	my ($self, $chr, $patients) = @_;
	return if (scalar(@$patients) == 0);
	foreach my $patient (@$patients) {
		$patient->in_the_attic(1);
	}
	return;
}

sub setIntersectPatients {
	my ($self, $chr, $patients) = @_;
	return if (scalar(@$patients) == 0);
	foreach my $patient (@$patients) {
		$patient->intersected(1);
	}
	$has_intersected_or_excuded_patients++;
	return;
}

sub setExcludePatients {
	my ($self, $chr, $patients, $status) = @_;
	return if (scalar(@$patients) == 0);
	return if ($status ne 'all' and $status ne 'he' and $status ne 'ho');
	foreach my $patient (@$patients) {
		$patient->excluded($status);
	}
	$has_intersected_or_excuded_patients++;
	return;
}

sub setIntersectFamilies {
	my ($self, $chr, $families) = @_;
	return if (scalar(@$families) == 0);
	foreach my $family (@$families) {
		$family->intersected(1);
	}
	$has_intersected_or_excuded_patients++;
	return;
}

sub setExcludeFamilies {
	my ($self, $chr, $families) = @_;
	return if (scalar(@$families) == 0);
	foreach my $family (@$families) {
		$family->excluded('all');
	}
	$has_intersected_or_excuded_patients++;
	return;
}


sub getVector_project {
	my ($self,  $chr, $vector_chr_gene_init ) = @_;
	return $self->getVector_project_or_fam( $chr, $vector_chr_gene_init, $chr->getProject() );
}

sub getVector_fam {
	my ($self,  $chr, $vector_chr_gene_init, $family ) = @_;
	return $self->getVector_project_or_fam( $chr, $vector_chr_gene_init, $family );
}

sub getVector_project_or_fam {
	my ($self,  $chr, $vector_chr_gene_init, $family_or_project ) = @_;
#	warn 'NB: '.$chr->countThisVariants($vector_chr_gene_init);
	my $var_ok = $chr->getNewVector();
	my $var_excluded = $chr->getNewVector();
	my $var_intersect = $chr->getNewVector();
	return ($vector_chr_gene_init,$var_excluded,$var_intersect) if $vector_chr_gene_init->is_empty();
	my $nb_intersected = 0;
	foreach my $patient (@{$family_or_project->getPatients()}) {
		next if ($patient->in_the_attic());
		my $var_pat = $patient->getVariantsVector($chr)->Clone();
		$var_pat->Intersection($var_pat, $vector_chr_gene_init);
#		warn 'pat init '.$patient->name().': '.$chr->countThisVariants($var_pat);
		next if $var_pat->is_empty();
		if ($patient->excluded()) {
			my $status = $patient->excluded();
#			warn '-> exclude '.$status;
			if ($status eq 'all') { $var_excluded += $patient->getVectorOrigin( $chr ); }
			elsif ($status eq 'he') { $var_excluded += $patient->getVectorOriginHe( $chr ); }
			elsif ($status eq 'ho') { $var_excluded += $patient->getVectorOriginHo( $chr ); }
		}
		if ($patient->intersected() and $nb_intersected == 0) {
#			warn '-> intersected ';
			$nb_intersected++;
			$var_intersect += $var_pat;
			$var_ok += $var_pat;
		}
		elsif ($patient->intersected()) {
#			warn '-> intersected ';
			$var_intersect->Intersection($var_pat, $var_intersect);
			$var_ok += $var_pat;
		}
		$var_ok += $var_pat;
#		warn 'pat '.$patient->name().': '.$chr->countThisVariants($var_pat);
#		warn 'all: '.$chr->countThisVariants($var_ok);
#		warn "\n\n";
	}
	$var_ok -= $var_excluded;
	if ($nb_intersected) {
		$var_ok->Intersection($var_ok, $var_intersect);
	}
#	die;
	
	#renvoie le vector du gene current si level_ind est gene
	return $vector_chr_gene_init if ($chr->getProject->level_ind() eq 'gene');
	
	#renvoie le vector des variants restants si level_ind est ind
	return ($var_ok, $var_excluded, $var_intersect);
}

sub updateVectorPatients {
	my ($self, $lPatients, $chr, $var_local_ok, $var_local_excluded, $var_local_intersect) = @_;
	foreach my $p (@$lPatients) {
		if ($var_local_excluded and not $var_local_excluded->is_empty()) {
			$p->{origin}->{all}->{$chr->name} -= $var_local_excluded;
			$p->{origin}->{he}->{$chr->name} -= $var_local_excluded;
			$p->{origin}->{ho}->{$chr->name} -= $var_local_excluded;
		}
		if ($var_local_intersect) {
			$p->{origin}->{all}->{$chr->name}->Intersection($p->{origin}->{all}->{$chr->name}, $var_local_intersect);
			$p->{origin}->{he}->{$chr->name}->Intersection($p->{origin}->{he}->{$chr->name}, $var_local_intersect);
			$p->{origin}->{ho}->{$chr->name}->Intersection($p->{origin}->{ho}->{$chr->name}, $var_local_intersect);
		}
		$p->{origin}->{all}->{$chr->name}->Intersection($p->{origin}->{all}->{$chr->name}, $var_local_ok);
		$p->{origin}->{he}->{$chr->name}->Intersection($p->{origin}->{he}->{$chr->name}, $var_local_ok);
		$p->{origin}->{ho}->{$chr->name}->Intersection($p->{origin}->{ho}->{$chr->name}, $var_local_ok);
	}
}

sub setIntersectExclude_PAT_FAM {
	my ($self,  $chr ) = @_;
	return if not $has_intersected_or_excuded_patients;
	return setIntersectExclude_PAT_FAM_genes($chr) if ($chr->getProject->level_ind() eq 'gene');
	return setIntersectExclude_PAT_FAM_genes($chr) if ($chr->getProject->level_fam() eq 'gene');
	my $var_ok = $chr->getNewVector();
	my $var_excluded = $chr->getNewVector();
	my $var_intersect = $chr->getNewVector();
	my $nb_intersected = 0;
	if ($chr->getProject->typeFilters() eq 'individual') {
		my ($var_local_ok, $var_local_excluded, $var_local_intersect) = $self->getVector_project($chr, $chr->getVariantsVector());
		$var_ok += $var_local_ok;
		$var_ok -= $var_local_excluded;
		$var_excluded += $var_local_excluded;
		$var_intersect += $var_local_intersect if $var_local_intersect;
		my @lPat = @{$chr->getProject->getPatients()};
		$self->updateVectorPatients(\@lPat, $chr, $var_local_ok, $var_local_excluded, $var_local_intersect);
		
	}
	else {
		foreach my $family (@{$chr->getProject->getFamilies()}) {
			my ($var_local_ok, $var_local_excluded, $var_local_intersect) = $self->getVector_fam($chr, $chr->getVariantsVector(), $family);
			if ($family->in_the_attic()) {
				next;
			}
			elsif ($family->excluded()) {
				$var_excluded += $var_local_ok;
			}
			elsif ($family->intersected() and $nb_intersected == 0) {
				$nb_intersected++;
				$var_intersect += $var_local_ok;
			}
			elsif ($family->intersected()) {
				$var_intersect->Intersection($var_local_ok, $var_intersect);
			}
			$var_ok += $var_local_ok;
			my @lPat = @{$family->getPatients()};
			$self->updateVectorPatients(\@lPat, $chr, $var_local_ok, $var_local_excluded, $var_local_intersect);
		}
	}
	$var_ok -= $var_excluded;
	if ($nb_intersected) {
		$var_ok->Intersection($var_ok, $var_intersect);
	}
	$chr->setVariantsVector($var_ok);
}

sub setIntersectExclude_PAT_FAM_genes {
	my ($self,  $chr ) = @_;
	return if not $has_intersected_or_excuded_patients;
	my $var_ok = $chr->getNewVector();
	my $var_excluded = $chr->getNewVector();
	my $var_intersect = $chr->getNewVector();
	my $nb_intersected = 0;
	
	my $h_genes_fam_intersect;
	
	foreach my $gene (@{$chr->getGenes()}) {
		if ($chr->getProject->typeFilters() eq 'individual') {
			my ($self, $var_local_ok, $var_local_excluded, $var_local_intersect) = $self->getVector_project($chr, $gene->getCurrentVector(), $chr->getProject());
			$var_ok += $var_local_ok;
		}
		else {
			foreach my $family (@{$chr->getProject->getFamilies()}) {
				next if ($family->in_the_attic());
				
				my $var_gene = $gene->getCurrentVector->Clone();
				my ($self, $var_local_ok, $var_local_excluded, $var_local_intersect) = $self->getVector_fam($chr, $var_gene, $family);
				my $var_fam = $var_local_ok;
				my $var_to_add = $chr->getNewVector();
				#ajoute le vecteur du gene level_fam est gene
				
				if ($chr->getProject->level_ind() eq 'gene' or $chr->getProject->level_fam() eq 'gene') { 
					if ($var_local_ok->is_empty) {
						$var_excluded += $var_gene;
					}
					elsif ($family->excluded()) {
						$var_excluded += $var_gene;
					}
					if ($chr->getProject->level_ind() eq 'gene') {
						my $has_pat_intersect;
						my $has_intersect = 1;
						foreach my $patient (@{$family->getPatients()}) {
							next if not $patient->intersected();
							$has_pat_intersect = 1;
							my $var_pat = $patient->getVariantsVector($chr)->Clone();
							$var_pat->Intersection($var_pat, $var_gene);
							$has_intersect = undef if $var_pat->is_empty();
							
						}
						if ($has_intersect) { $var_ok += $var_gene; }
						else { $var_gene->Empty(); }
					}
					if ($chr->getProject->level_fam() eq 'gene') {
						if ($family->intersected() and $nb_intersected == 0) {
							$h_genes_fam_intersect->{$gene->id()}->{$family->name()}++;
							$nb_intersected++;
							$var_intersect += $var_gene;
						}
						elsif ($family->intersected()) {
							$h_genes_fam_intersect->{$gene->id()}->{$family->name()}++;
							$var_intersect->Intersection($var_to_add, $var_gene);
						}
					}
					$var_ok += $var_gene;
				}
				else {
					$var_to_add += $var_fam;
					if ($family->excluded()) {
						$var_excluded += $var_to_add;
					}
					elsif ($family->intersected() and $nb_intersected == 0) {
						$nb_intersected++;
						$var_intersect += $var_to_add;
					}
					elsif ($family->intersected()) {
						$var_intersect->Intersection($var_to_add, $var_intersect);
					}
					$var_ok += $var_to_add;
				}
			}
		}
	}
	
	my @lPat = @{$chr->getProject->getPatients()};
	$var_ok -= $var_excluded;
	if ($nb_intersected) {
		if ($chr->getProject->level_fam() eq 'gene') {
			foreach my $gene (@{$chr->getGenes()}) {
				my $has_fam_intersect;
				my $has_intersect = 1;
				foreach my $family (@{$chr->getProject->getFamilies()}) {
					next if not $family->intersected();
					$has_fam_intersect = 1;
					$has_intersect = undef if not exists $h_genes_fam_intersect->{$gene->id()};
					$has_intersect = undef if not exists $h_genes_fam_intersect->{$gene->id()}->{$family->name()};
				}
				if ($has_fam_intersect and not $has_intersect) {
					$var_excluded += $gene->getCurrentVector();
					$var_ok -= $gene->getCurrentVector();
				}
			}
		}
		else {
			$var_ok->Intersection($var_ok, $var_intersect);
		}
		$self->updateVectorPatients(\@lPat, $chr, $var_ok, $var_excluded, $var_intersect);
	}
	else {
		$self->updateVectorPatients(\@lPat, $chr, $var_ok, $var_excluded);
	}
	$chr->setVariantsVector($var_ok);
}


sub filter_atLeast {
	my ($self, $chr, $atLeast, $typeFilters, $level_fam) = @_;
	if ($atLeast and $atLeast >= 2) {
		if ($typeFilters eq 'familial') {
			if ($level_fam eq 'gene') { atLeastFilter_genes_fam($chr, $atLeast); }
			else { atLeastFilter_var_fam($chr, $atLeast); }
		}
		if ($typeFilters eq 'somatic') {
			if ($level_fam eq 'gene') { atLeastFilter_genes_som($chr, $atLeast); }
			else { atLeastFilter_var_som($chr, $atLeast); }
		}
		else {
			if ($level_fam eq 'gene') { atLeastFilter_genes_ind($chr, $atLeast); }
			else { atLeastFilter_var_ind($chr, $atLeast); }
		}
	}
	return;
}

sub filter_genes_text_search {
	my ($self, $chr, $filter_text) = @_;
	return unless ($filter_text);
	$chr->getProject->filter_text($filter_text);
	my $vector_genes = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		$chr->getProject->print_dot(1);
		unless ($chr->check_filter_text($gene)) {
			next;
		}
		$vector_genes += $gene->getVariantsVector();
	}
	$vector_genes->Intersection($vector_genes, $chr->getVariantsVector());
	$chr->setVariantsVector($vector_genes);
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_genes_text_search - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
}

sub filter_genes_only_genes_names {
	my ($self, $chr, $only_genes) = @_;
	return unless ($only_genes);
	foreach my $name (split(',', $only_genes)) {
		if ($name =~ /capture/) {
			my @lTmp = split(':', $name);
			my $capture_name = $lTmp[-1];
			foreach my $parent_capture_name (keys %{$chr->getProject->buffer->getAllGenesNamesInAllBundle()}) {
				next unless (exists $chr->getProject->buffer->getAllGenesNamesInAllBundle->{$parent_capture_name}->{$capture_name});
				foreach my $gene_name (keys %{$chr->getProject->buffer->getAllGenesNamesInAllBundle->{$parent_capture_name}->{$capture_name}->{genes}}) {
					$chr->getProject->{only_genes}->{$gene_name} = undef;
				}
			}
		}
		else { $chr->getProject->{only_genes}->{$name} = undef; }
	}
	my $vector_genes = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		$chr->getProject->print_dot(1);
		next if ($chr->getProject->only_genes() and not exists $chr->getProject->only_genes->{uc($gene->external_name())});
		$vector_genes += $gene->getVariantsVector();
	}
	$chr->setVariantsVector($vector_genes);
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_genes_only_genes_names - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
}


sub get_vector_filter_gene_annotations {
	my ($self, $gene, $hcat,$rocks) = @_;
	
	my $vector = $gene->getChromosome->getNewVector();
	my $z = $gene->getCurrentVector();
	#warn $z->Min()." ".$z->Max;
	my $st= "";
	my $cs = $gene->getChromosome->rocks_vector('r')->get_vector_gene($gene->id);
 	foreach my $c (values %$hcat){
 		next unless $cs->{$c};
 		$st.= $cs->{$c}.",";
 	}
 	chop($st);
 	my $set = Set::IntSpan::Fast->new($st);
 	$vector->from_Enum($st);
	return $vector;
}
sub return_cat {
	my ($project,@unselect) = @_;
	my %hcat = %{$project->buffer->config->{'functional_annotations'}};

	foreach my $c (@unselect){
		delete $hcat{$c};
	}
	return keys %hcat;
}
sub filter_genes_annotations {
	my ($self, $chr, $hFiltersChr) = @_;
	return unless ($hFiltersChr);
	my $tt =time;
	

	my @all_cat = return_cat($chr->project,keys %$hFiltersChr);
	
	my $variants_genes  = $chr->getNewVector();
	my $vvv = $chr->getVariantsVector();
	#$chr->rocks_vector("r")->get_vector_chromosome("coding") & $chr->getVariantsVector();
	
	my $cc = $chr->getVariantsVector()->Clone;
	my $nb = 0;
	my $nb_genes = 0;
	my $id_genes = {};
	my $t =time;
	#$rocks->prepare([keys %$lh]);
	foreach my $gene ( @{$chr->getGenesFromVector($vvv)}) {
	#foreach my $gene ( @{$chr->getGenesFromVector($chr->getVariantsVector())}) {
	#foreach my $gene ( @{$objs}) {
		$nb_genes ++;
		$chr->getProject->print_dot(1);
		my $debug;
#		if ($gene->external_name eq "RNF115"){
#			$debug =1;
#		}
		my $vsmall = $gene->getCompactVectorOriginCategories(\@all_cat,$debug);
#		if ($debug){
#			warn Dumper @all_cat;
#			warn "\t\t ".$vsmall->Norm;
#			#die();
#		}
		
		my $vchr = $gene->return_compact_vector( $chr->getVariantsVector());
		$vsmall &= $vchr;
		if ($vsmall->is_empty()){
			delete $chr->project->{genes_object}->{$gene->id};
			next;
		}
		$gene->setCurrentVector($vsmall);
		
		$nb ++;
		$id_genes->{$gene->id} ++;
		$variants_genes += $gene->enlarge_compact_vector($vsmall);
	}
	
	$chr->getVariantsVector->Intersection($chr->getVariantsVector(), $variants_genes);
	$chr->{buffer_vector} = $chr->getVariantsVector();
#	warn "---------------- $nb/$nb_genes/".scalar(keys %{$chr->{genes_object}});
#	warn abs(time - $t);
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_genes_annotations - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
}


sub filter_genetics_models {
	my ($self, $chr, $model, $h_args) = @_;
	return if $chr->getVariantsVector->is_empty();
	if ($model and $chr->getProject->typeFilters() eq 'individual') {
		$self->filter_model_individual_recessif($chr) if ($model eq 'recessif');
		$self->filter_model_individual_compound($chr) if ($model eq 'compound');
	}
	elsif ($model and $chr->getProject->typeFilters() eq 'familial') {
		my @lModels;
		if ($model eq 'recessif_compound') {
			push(@lModels, 'recessif');
			push(@lModels, 'compound');
		}
		elsif ($model eq 'uniparental_recessif_compound') {
			push(@lModels, 'recessif');
			push(@lModels, 'compound');
			push(@lModels, 'uniparental_disomy');
		}
		else { push(@lModels, $model); }
		$self->filter_models_familial_union($chr, \@lModels, $h_args);
	}
	elsif ($model and $chr->getProject->typeFilters() eq 'somatic') {
		$self->filter_model_somatic_loh($chr) if ($model eq 'loh');
		$self->filter_model_somatic_dbl_evt($chr) if ($model eq 'dbl_evt');
		$self->filter_model_somatic_only_tissues_somatic($chr) if ($model eq 'only_tissues_somatic');
	}
	foreach my $gene (@{$chr->getGenes()}) {
		next if not $gene->compact_vector();
		my $vv = $chr->getVariantsVector() & $gene->getCurrentVector();
		$gene->setCurrentVector($vv);
	}
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER models - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
}


1;
