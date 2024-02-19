package QueryVectorFilter;
use strict;
use Vcf;
use Moo;

use Data::Dumper;
use Config::Std;
use Carp;



sub filter_vector_ratio {
	my ($chr, $limit_ratio, $filter_type_ratio) = @_;
	return unless ($limit_ratio);
	return if ($limit_ratio eq 'all');
	my $vector_ok = $chr->getNewVector();
	foreach my $patient (@{$chr->getPatients()}) {
		next if ($patient->in_the_attic());
		my $vector_ratio_name = $patient->name()."_ratio_".$limit_ratio;
		my $vquality = $chr->getVectorScore($vector_ratio_name);
		if (lc($filter_type_ratio) eq 'max') {
			$vector_ok->Fill();
			$vector_ok -= $vquality;
		}
		else {
			$vector_ok += $vquality;
		}
	}
	$vector_ok->Intersection($vector_ok, $chr->getVariantsVector());
	$chr->setVariantsVector($vector_ok);
}

sub filter_vector_ncboost {
	my ($chr, $cat) = @_;
	return unless ($cat);
	if ($chr->id() eq 'MT') {
		$chr->setVariantsVector( $chr->getNewVector() );
		return;
	}
	$chr->project->print_dot(1);
	my $vector_ncboost = $chr->getVectorNcboost($cat);
	$vector_ncboost->Intersection($vector_ncboost, $chr->getVariantsVector());
	$chr->setVariantsVector($vector_ncboost);
	return;
}

sub filter_vector_hgmd_dm {
	my ($chr) = @_;
	confess() unless ($chr->project->isUpdate());
	my $vector_hgmd_dm = $chr->getVectorLmdbDm();
	$vector_hgmd_dm->Intersection($vector_hgmd_dm, $chr->getVariantsVector());
	$chr->setVariantsVector($vector_hgmd_dm);
}

sub filter_vector_hgmd_dm_new_version {
	my ($chr) = @_;
	confess unless ($chr->project->isUpdate());
	my $vector_hgmd_dm_new = $chr->getVectorLmdbDm_newVersion();
	$vector_hgmd_dm_new->Intersection($vector_hgmd_dm_new, $chr->getVariantsVector());
	$chr->setVariantsVector($vector_hgmd_dm_new);
}

sub filter_vector_region_ho {
	my ($chr, $filter_nbvar_regionho, $filter_regionho_sub_only, $indiv_or_fam) = @_;
	return unless ($filter_nbvar_regionho);
	return filter_vector_region_ho_individual($chr, $filter_nbvar_regionho, $filter_regionho_sub_only) if ($indiv_or_fam eq 'individual');
	return filter_vector_region_ho_familial($chr, $filter_nbvar_regionho, $filter_regionho_sub_only) if ($indiv_or_fam eq 'familial');
}

sub filter_vector_region_ho_individual {
	my ($chr, $filter_nbvar_regionho, $filter_regionho_sub_only) = @_;
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
	return;
}

sub filter_vector_region_ho_familial {
	my ($chr, $filter_nbvar_regionho, $filter_regionho_sub_only) = @_;
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
	return;
}

sub filter_vector_type_variants {
	my ($chr, $hToFilter) = @_;
	return unless ($hToFilter);
	my $vector = $chr->getVariantsVector->Clone();
	$vector->AndNot($vector, $chr->get_vector_type_variants($hToFilter));
	$chr->setVariantsVector($vector);
	return;
}

sub filter_vector_predicted_splice_region_global {
	my ($chr, $hToFilter) = @_;
	return unless ($hToFilter);
	return unless (exists $hToFilter->{predicted_splice_site});
	my $vector = $chr->getVariantsVector->Clone();
	$vector->AndNot($vector, $chr->get_vector_predicted_splice_region($hToFilter));
	$chr->setVariantsVector($vector);
	return;
}

sub filter_vector_freq_variants {
	my ($chr, $hToFilter) = @_;
	return unless ($hToFilter);
	my $vector = $chr->getVariantsVector->Clone();
	$vector->Intersection($vector, $chr->get_vector_freq_variants($hToFilter));
	$chr->setVariantsVector($vector);
	return;
}

sub filter_vector_gnomad_ho_ac_variants {
	my ($chr, $hToFilter) = @_;
	return unless ($hToFilter);
	my $vector = $chr->getVariantsVector->Clone();
	if (exists $hToFilter->{gnomad_ac_all})    { $vector->Intersection($vector, $chr->get_vector_gnomad_ac_variants($hToFilter)); }
	if (exists $hToFilter->{gnomad_ho_ac_all}) { $vector->Intersection($vector, $chr->get_vector_gnomad_ho_ac_variants($hToFilter)); }
	$chr->setVariantsVector($vector);
	return;
}

sub filter_vector_confidence_variants {
	my ($chr, $hToFilter) = @_;
	return unless ($hToFilter);
	my $vector = $chr->getVariantsVector->Clone();
	$vector->AndNot($vector, $chr->get_vector_confidence_variants($hToFilter));
	$chr->setVariantsVector($vector);
	return;
}

sub filter_vector_cadd_variants {
	my ($chr, $hToFilter, $keep_indels_cadd) = @_;
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
	return;
}

sub filter_vector_polyscore {
	my ($chr, $polyscore_value) = @_;
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
	my ($chr, $nb_dejavu, $dejavu_ho, $test) = @_;
	return unless ($nb_dejavu);
	my $vector = $chr->getVariantsVector->Clone();
	if ($test) {
		$vector->AndNot($vector, $chr->get_vector_dejavu_TESTS_F_ONLY($vector, $nb_dejavu, $dejavu_ho));
	}
	else {
		$vector->AndNot($vector, $chr->get_vector_dejavu($vector, $nb_dejavu, $dejavu_ho));
	}
	$chr->setVariantsVector($vector);
	return;
}

# AT_LEAST par patients en mode VAR et IND 
sub atLeastFilter_var_ind {
	my ($chr, $atLeast) = @_;
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
	my ($chr, $atLeast) = @_;
	my (@lOk, $hVectFam);
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
	my ($chr, $atLeast) = @_;
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
	my ($chr, $atLeast) = @_;
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
	my ($chr, $atLeast) = @_;
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
	my ($chr, $atLeast) = @_;
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

sub filter_vector_gnomad_ac {
	my ($chr, $filter_name) = @_;
	return unless ($filter_name =~ /gnomad_/);
	my $vector_ok = $chr->getVariantsVector->Clone();
	$vector_ok->Intersection($vector_ok, $chr->getVectorScore($filter_name));
	$chr->setVariantsVector($vector_ok);
}


######### FILTRES MODELS CHR ##########

our $hModelsMethodsNames_familial = {
   'denovo' => 'getModelVector_fam_denovo',
   'strict-denovo' => 'getModelVector_fam_strict_denovo',
   'dominant' => 'getModelVector_fam_dominant',
   'mosaic' => 'getModelVector_fam_mosaique',
   'recessif' => 'getModelVector_fam_recessif',
#   'incomplete_penetrance' => 'getModelVector_fam_incomplete_penetrance',
   'uniparental_disomy' => 'getModelVector_fam_uniparental_disomy',
   'only_tissues_somatic' => 'getModelVector_som_only_tissues_somatic',
};

sub filter_models_familial_union {
	my ($chr, $list_models, $h_arguments) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_models_familial_union need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	confess("\n\nERROR: QueryVectorFilter::filter_models_familial_union need a list of models. Die.\n\n") unless ($list_models);
	my $vector_models = get_models_familial_union($chr, $list_models, $h_arguments);
	foreach my $patient (@{$chr->getPatients()}) {
		my $vector_patient = $chr->getNewVector();
		foreach my $model (@$list_models) {
			confess("\n\nERROR: $model model not used for patient ".$patient->name()." for CHR".$chr->id()."... Die.\n\n") unless (exists $patient->hash_models_genetics_used->{$model}->{$chr->id()});
			$vector_patient += $patient->hash_models_genetics_used->{$model}->{$chr->id()};
		}
		$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $vector_patient);
	}
	$chr->setVariantsVector($vector_models);
}

sub get_models_familial_union {
	my ($chr, $list_models, $h_arguments) = @_;
	confess("\n\nERROR: QueryVectorFilter::get_models_familial_union need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	confess("\n\nERROR: QueryVectorFilter::get_models_familial_union need a list of models. Die.\n\n") unless ($list_models);
	my $vector_models = $chr->getNewVector();
	foreach my $model (@$list_models) {
		if ($model eq 'compound') { $vector_models += get_model_familial_compound($chr, $h_arguments); }
		else { $vector_models += get_model_familial_common($chr, $model, $h_arguments); }
	}	
	return $vector_models;
}

sub filter_model_familial_denovo {
	my ($chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_denovo need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return filter_model_familial_common($chr, 'denovo');
}

sub get_model_familial_denovo {
	my ($chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::get_model_familial_denovo need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return get_model_familial_common($chr, 'denovo');
}

sub filter_model_familial_strict_denovo {
	my ($chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_strict_denovo need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return filter_model_familial_common($chr, 'strict-denovo');
}
		
sub get_model_familial_strict_denovo {
	my ($chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::get_model_familial_strict_denovo need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return get_model_familial_common($chr, 'strict-denovo');
}

sub filter_model_familial_dominant {
	my ($chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_dominant need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return filter_model_familial_common($chr, 'dominant');
}
		
sub get_model_familial_dominant {
	my ($chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::get_model_familial_dominant need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return get_model_familial_common($chr, 'dominant');
}

sub filter_model_familial_mosaic {
	my ($chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_mosaic need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return filter_model_familial_common($chr, 'mosaic');
}
		
sub get_model_familial_mosaic {
	my ($chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::get_model_familial_mosaic need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return get_model_familial_common($chr, 'mosaic');
}

sub filter_model_familial_recessif {
	my ($chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_recessif need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return filter_model_familial_common($chr, 'recessif');
}
		
sub get_model_familial_recessif {
	my ($chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::get_model_familial_recessif need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return get_model_familial_common($chr, 'recessif');
}

sub filter_model_familial_uniparental_disomy {
	my ($chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_uniparental_disomy need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return filter_model_familial_common($chr, 'uniparental_disomy');
}
		
sub get_model_familial_uniparental_disomy {
	my ($chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::get_model_familial_uniparental_disomy need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	return get_model_familial_common($chr, 'uniparental_disomy');
}

sub filter_model_familial_compound {
	my ($chr, $h_arguments) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_compound need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_compound need a Hash of Arguments. Die.\n\n") unless ($h_arguments);
	return filter_model_familial_common($chr, 'compound');
}

sub get_model_familial_compound {
	my ($chr, $h_arguments) = @_;
	my $hFiltersChr = $h_arguments->{'filters_1'};
	my $hFiltersChr_var2 = $h_arguments->{'filters_2'};
	my $vector_filtered = $h_arguments->{'vector_filters_1'};
	my $vector_filtered_2 = $h_arguments->{'vector_filters_2'};
	my $variants_genes  = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		my $vector_gene = filter_gene_model_familial_compound($gene, $hFiltersChr, $hFiltersChr_var2, $vector_filtered, $vector_filtered_2);
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
	my ($gene, $hFiltersChr, $hFiltersChr_var2, $vector_filtered, $vector_filtered_2) = @_;
	my $var_gene_annot = $gene->getChromosome->getNewVector();
	my $var_gene_annot_1 = $gene->getCategoriesVariantsVector($hFiltersChr);
	$var_gene_annot_1->Intersection($var_gene_annot_1, $vector_filtered) if ($vector_filtered);
	if ($var_gene_annot_1->is_empty()) { return; }
	$var_gene_annot += $var_gene_annot_1;
	my $var_gene_annot_2 = $gene->getCategoriesVariantsVector($hFiltersChr_var2);
	$var_gene_annot_2->Intersection($var_gene_annot_2, $vector_filtered_2) if ($vector_filtered_2);
	if ($var_gene_annot_2->is_empty()) { return; }
	$var_gene_annot += $var_gene_annot_2;
	return ($var_gene_annot_1, $var_gene_annot_2);
}

sub filter_gene_model_familial_compound {
	my ($gene, $hFiltersChr, $hFiltersChr_var2, $vector_filtered, $vector_filtered_2) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_gene_model_familial_compound need a GenBoGeneCache. Die.\n\n") unless ($gene);
	my $vector_gene = $gene->getVariantsVector->Clone();
	$vector_gene->Intersection($vector_gene, $gene->getChromosome->getVariantsVector());
	$vector_gene->Intersection($vector_gene, $gene->getModelGeneVector_fam_compound());
	my $vector_ok = $gene->getChromosome->getNewVector();
	foreach my $family (@{$gene->getChromosome->getFamilies()}) {
		my $vector_gene_tmp = $vector_gene->Clone();
		$vector_gene_tmp->Intersection($vector_gene_tmp, $family->get_vector_keep_var_compound($gene->getChromosome()));
		$vector_ok += $vector_gene_tmp;
	}
	$vector_gene->Intersection($vector_gene, $vector_ok);
	
	if ($vector_gene->is_empty()) { $vector_gene->Empty(); }
	elsif ($hFiltersChr_var2) {
		my ($var_gene_annot_1, $var_gene_annot_2) = get_gene_model_familial_compound_prepare_multi_annot($gene, $hFiltersChr, $hFiltersChr_var2, $vector_filtered, $vector_filtered_2);
		if ($var_gene_annot_1 and $var_gene_annot_2) {
			$var_gene_annot_1->Intersection($var_gene_annot_1, $vector_gene);
			if ($var_gene_annot_1->is_empty()) {
				$vector_gene->Empty();
			}
			$var_gene_annot_2->Intersection($var_gene_annot_2, $vector_gene);
			if ($var_gene_annot_2->is_empty()) {
				$vector_gene->Empty();
			}
			if ($gene->getChromosome->countThisVariants($vector_gene) < 2) {
				$vector_gene->Empty();
			}
		}
		else {
			$vector_gene->Empty();
		}
	}
	return $vector_gene;
}

sub filter_model_familial_common {
	my ($chr, $model) = @_;
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
	my ($chr, $model, $h_arguments) = @_;
	return unless ($model);
	my $model_method_name = $hModelsMethodsNames_familial->{$model};
	return unless ($model_method_name);
	my $vector_update = $chr->getNewVector();
	foreach my $family (@{$chr->getFamilies()}) {
		my ($has_patients, $has_parents); 
		$family->{hash_models_genetics_used}->{$model}->{$chr->id()} = $family->$model_method_name($chr, $h_arguments);
		foreach my $patient (@{$family->getPatients()}) {
			my $vector_patient = $patient->getVariantsVector($chr)->Clone();
			$vector_patient->Intersection( $vector_patient, $family->{hash_models_genetics_used}->{$model}->{$chr->id()});
			$vector_update += $vector_patient;
			$patient->{hash_models_genetics_used}->{$model}->{$chr->id()} = $vector_patient->Clone();
			$has_patients = 1;
		}
	}
	return $vector_update;
}

sub filter_model_somatic_loh {
	my ($chr) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_model_somatic_loh need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
	$chr->setVariantsVector($chr->getModelVector_som_loh());
	return;
	
}

sub filter_model_somatic_only_tissues_somatic {
	my $chr = shift;
	return filter_model_familial_common($chr, 'only_tissues_somatic');
}

sub filter_model_individual_recessif {
	my $chr = shift;
	my $variants_genes  = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		filter_gene_model_individual_recessif($gene);
		if ($gene->getVariantsVector->is_empty()) {
			$chr->supressGene($gene);
			next;
		}
		$variants_genes += $gene->getVariantsVector();
	}
	$chr->getVariantsVector->Intersection($chr->getVariantsVector(), $variants_genes);
	foreach my $patient (@{$chr->getPatients()}) {
		$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $patient->{model_indiv_recessif}->{$chr->id}) if exists ($patient->{model_indiv_recessif}->{$chr->id});
	}
	return;
}

sub filter_model_individual_compound {
	my $chr = shift;
	my $variants_genes  = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		filter_gene_model_individual_compound($gene);
		if ($gene->getVariantsVector->is_empty()) {
			$chr->supressGene($gene);
			next;
		}
		$variants_genes += $gene->getVariantsVector();
	}
	$chr->getVariantsVector->Intersection($chr->getVariantsVector(), $variants_genes);
	foreach my $patient (@{$chr->getPatients()}) {
		$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $patient->{model_indiv_compound}->{$chr->id}) if exists ($patient->{model_indiv_compound}->{$chr->id});
	}
	return;
}


sub filter_model_somatic_dbl_evt {
	my $chr = shift;
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
	my ($chr, $hGeneIds, $can_use_hgmd) = @_;
	my @lOk;
	foreach my $hash (@{$chr->values_lmdb_genes()}) {
		next unless (exists $hGeneIds->{$hash->{name}});
		push(@lOk, $hash);
	}
	$chr->{values_lmdb_genes} = \@lOk;
	return;
}

sub filter_usefull_genes_ids {
	my ($chr, $hFiltersChr, $can_use_hgmd) = @_;
	# Filter list genes to construct
	my $hKeep_annot;
	foreach my $filter_name (keys %{$chr->getProject->ensembl_annotations()}) {
		next if (exists $hFiltersChr->{$filter_name});
		$hKeep_annot->{$filter_name} = undef;
	}
	my @lOk;
	foreach my $hash (@{$chr->values_lmdb_genes()}) {
		my $is_usefull;
		foreach my $cat (keys %{$hKeep_annot}) {
			if (exists $hash->{bitvector}->{$cat}) {
				$is_usefull = 1;
				push(@lOk, $hash);
				last;
			}
		}
		# TODO: a decocher - KEEP IF HAS HGMD DM VAR
		if ($can_use_hgmd and $chr->getProject->isUpdate()) {
			if (not $is_usefull and not exists $hash->{bitvector}->{intergenic}) {
				if ($chr->variations_hgmd_dm_tree->fetch_window(int($hash->{start}), (int($hash->{end}) + 1)) > 0) {
					push(@lOk, $hash);
				}
			}
		}
	}
	$chr->{values_lmdb_genes} = \@lOk;
	return;
}

sub filter_gene_model_individual_recessif {
	my ($gene) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_gene_model_individual_recessif need a GenBoGeneCache. Die.\n\n") unless ($gene);
	$gene->getVariantsVector->Intersection($gene->getVariantsVector(), $gene->getChromosome->getModelGeneVector_indiv_recessif($gene));
	return;
}

sub filter_gene_model_individual_compound {
	my ($gene) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_gene_model_individual_compound need a GenBoGeneCache. Die.\n\n") unless ($gene);
	$gene->getVariantsVector->Intersection($gene->getVariantsVector(), $gene->getChromosome->getModelGeneVector_indiv_compound($gene));
	return;
}

sub filter_gene_model_familial_recessif_compound {
	my ($gene, $vector_global_fam_recessif) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_gene_model_familial_recessif_compound need a GenBoGeneCache. Die.\n\n") unless ($gene);
	my $vector_fam_recessif = $gene->getChromosome->getNewVector();
	$vector_fam_recessif->Intersection($vector_global_fam_recessif, $gene->getVariantsVector());
	my $vector_fam_compound = $gene->getChromosome->getNewVector();
	$vector_fam_compound->Intersection($gene->getModelGeneVector_fam_compound(), $gene->getVariantsVector());
	my $vector_fam_recessif_compound = $gene->getChromosome->getNewVector();
	$vector_fam_recessif_compound += $vector_fam_recessif;
	$vector_fam_recessif_compound += $vector_fam_compound;
	$gene->getVariantsVector->Intersection($gene->getVariantsVector(), $vector_fam_recessif_compound);
	return;
}

sub filter_gene_model_somatic_dbl_evt {
	my ($gene) = @_;
	confess("\n\nERROR: QueryVectorFilter::filter_gene_model_somatic_dbl_evt need a GenBoGeneCache. Die.\n\n") unless ($gene);
	$gene->getVariantsVector->Intersection($gene->getVariantsVector(), $gene->getChromosome->getModelGeneVector_som_dbl_evt($gene));
	return;
}



######### FILTRES PATIENTS / FAMILIES ##########



# method in the attic pour un patient 
sub setInTheAttic {
	my ($chr, $patients) = @_;
	return if (scalar(@$patients) == 0);
	foreach my $patient (@$patients) { $patient->setInTheAttic($chr, 1); }
	$chr->update_from_patients();
}

sub setExcludePatient {
	my ( $chr, $patients, $status, $force_ind) = @_;
	return if (scalar(@$patients) == 0);
	if ($chr->getProject->typeFilters() eq 'familial' and not $force_ind) { return setExcludePatient_FAM($chr, $patients, $status); }
	return setExcludePatient_IND($chr, $patients, $status);
}

# method exclude (all, ho, he) pour un patient (IND)
sub setExcludePatient_IND {
	my ( $chr, $patients, $status ) = @_;
	my $var = $chr->getNewVector();
	my $var_pat = $chr->getNewVector();
	my $var_excluded = $chr->getNewVector();
	foreach my $pat (@$patients) {
		$pat->excluded($status);
		if ($status eq 'all') { $var_excluded += $pat->getVectorOrigin( $chr ); }
		elsif ($status eq 'he' or $status eq 'ho') { $var_excluded += $pat->getCategoryVariantsVector( $chr, $status ); }
	}
	unless ($chr->variants_excluded()) { $chr->variants_excluded($var_excluded); }
	else { $chr->{variants_excluded} += $var_excluded; }
	foreach my $pat (@{$chr->getPatients()}) {
		next if ($pat->in_the_attic());
		$var_pat->Intersection($chr->patients_categories->{$pat->name()}, $chr->getVariantsVector());
		$var_pat -= $chr->variants_excluded();
		#$chr->patients_categories->{$pat->name()}->Intersection($chr->patients_categories->{$pat->name()}, $var_pat);
		$var += $var_pat;
	}
	$chr->setVariantsVector( $var );
}

# method exclude (all, ho, he) pour un patient (FAM)
sub setExcludePatient_FAM {
	my ( $chr, $patients, $status ) = @_;
	my $hPat;
	foreach my $pat (@$patients) { $hPat->{$pat->name()} = undef; }
	my $var = $chr->getNewVector();
	foreach my $family (@{$chr->getFamilies()}) {
		my @lPat;
		foreach my $patient (@{$family->getPatients()}) {
			push(@lPat, $patient) if (exists $hPat->{$patient->name()});
		}
		my $var_tmp = $chr->getNewVector();
		if (@lPat) {
			$var_tmp = $family->setExcludePatients(\@lPat, $status, $chr);
		}
		else {
			foreach my $patient (@{$family->getPatients()}) {
				$var_tmp += $patient->getVariantsVector($chr);
			}
		}
		$var += $var_tmp;
	}
	$chr->setVariantsVector($var);
}

# method exclude (all, ho, he) pour une famille
sub setExcludeFamily {
	my ($chr, $families) = @_;
	return if (scalar(@$families) == 0);
	my $hFamVar;
	my @lPatients;
	foreach my $family (@$families) {
		$hFamVar->{$family->name()} = $chr->getNewVector();
		foreach my $patient (@{$family->getPatients()}) {
			$hFamVar->{$family->name()} += $patient->getVectorOrigin($chr);
			push(@lPatients, $patient);
		}
		$family->excluded('all');
	}
	foreach my $patient (@{$chr->getPatients()}) {
		foreach my $fam_name (keys %$hFamVar) {
			$patient->getVariantsVector($chr)->AndNot($patient->getVariantsVector($chr), $hFamVar->{$fam_name});
		}
	}
	$chr->update_from_patients();
}

# method intersect exclude ho region patients
sub setExcludePatient_HO_REGIONS {
	my ( $chr, $patients, $filter_nbvar_regionho) = @_;
	return if ($chr->getProject->typeFilters() eq 'familial');
	return if ($filter_nbvar_regionho == 0);
	my $var = $chr->getNewVector();
	my $vector_del_ho_regions = $chr->getVector_ho_regions_after_exclude($patients, $filter_nbvar_regionho);
	foreach my $pat (@{$chr->getPatients()}) {
		my $vector_patient = $pat->getVariantsVector( $chr )->Clone();
		$vector_patient -= $vector_del_ho_regions;
		$var += $vector_patient;
	}
	$chr->setVariantsVector( $var );
}

# method intersect exclude ho region patients
sub setExcludeFamily_HO_REGIONS {
	my ( $chr, $families, $filter_nbvar_regionho) = @_;
	return if ($chr->getProject->typeFilters() eq 'individual');
	return if ($filter_nbvar_regionho == 0);
	my $var = $chr->getNewVector();
	my $vector_del_rec_regions = $chr->getVector_rec_regions_after_exclude($families, $filter_nbvar_regionho);
	foreach my $pat (@{$chr->getPatients()}) {
		my $vector_patient = $pat->getVariantsVector( $chr )->Clone();
		$vector_patient -= $vector_del_rec_regions;
		$var += $vector_patient;
	}
	$chr->setVariantsVector( $var );
}

# method intersect common ho region patients
sub setIntersectPatient_HO_REGIONS {
	my ($chr, $patients, $filter_nbvar_regionho) = @_;
	return if (scalar(@$patients) == 0);
	return if ($chr->getProject->typeFilters() eq 'familial');
	return if ($filter_nbvar_regionho == 0);
	my $var = $chr->getNewVector();
	my $vector_common_ho_regions = $chr->getVector_ho_regions_after_intersect($patients, $filter_nbvar_regionho);
	foreach my $pat (@{$chr->getPatients()}) {
		if ($vector_common_ho_regions->is_empty()) {
			$chr->patients_categories->{$pat->name()}->Empty();
		}
		else {
			$chr->patients_categories->{$pat->name()}->Intersection( $chr->patients_categories->{$pat->name()}, $vector_common_ho_regions );
			$var += $pat->getVariantsVector( $chr );
		}
	}
	$chr->setVariantsVector( $var );
}



# method intersect common recessive region families
sub setIntersectFamily_REC_REGIONS {
	my ($chr, $families, $filter_nbvar_regionho) = @_;
	return if (scalar(@$families) == 0);
	return if ($chr->getProject->typeFilters() eq 'individual');
	return if ($filter_nbvar_regionho == 0);
	my $var = $chr->getNewVector();
	my $vector_common_rec_regions = $chr->getVector_rec_regions_after_intersect($families, $filter_nbvar_regionho);
	foreach my $pat (@{$chr->getPatients()}) {
		$chr->patients_categories->{$pat->name()}->Intersection( $chr->patients_categories->{$pat->name()}, $vector_common_rec_regions );
		$var += $pat->getVariantsVector( $chr );
	}
	$chr->setVariantsVector( $var );
}

# method intersect pour un patient
sub setIntersectPatient {
	my ($chr, $patients) = @_;
	return if ($chr->getProject->level_ind() eq 'gene');
	return if (scalar(@$patients) == 0);
	if ($chr->getProject->typeFilters() eq 'familial') { return setIntersectPatient_FAM($chr, $patients); }
	return setIntersectPatient_IND($chr, $patients);
}

# method intersect pour un patient (IND)
sub setIntersectPatient_IND {
	my ($chr, $patients) = @_;
	return if (scalar(@$patients) == 0);
	my $var = $chr->getNewVector();
	my $var_intersected = $chr->getNewVector();
	my $nb_pat;
	foreach my $pat (@$patients) {
		$nb_pat++;
		$pat->intersected(1);
		if ($nb_pat == 1) { $var_intersected += $pat->getVariantsVector( $chr ); }
		else { $var_intersected->Intersection($var_intersected, $pat->getVariantsVector( $chr )) }
	}
	unless ($chr->variants_intersected()) { $chr->variants_intersected($var_intersected); }
	else { $chr->{variants_intersected} += $var_intersected; }
	foreach my $pat (@{$chr->getPatients()}) {
		$chr->patients_categories->{$pat->name()}->Intersection( $chr->patients_categories->{$pat->name()}, $chr->variants_intersected() );
		$var += $pat->getVariantsVector( $chr );
	}
	$chr->setVariantsVector( $var );
}

# method intersect pour un patient (FAM)
sub setIntersectPatient_FAM {
	my ($chr, $patients) = @_;
	return if (scalar(@$patients) == 0);
	my $hPat;
	foreach my $pat (@$patients) { push(@{$hPat->{$pat->family()}}, $pat); }
	my $var = $chr->getNewVector();
	foreach my $family (@{$chr->getFamilies()}) {
		my $var_tmp = $chr->getNewVector();
		if (exists $hPat->{$family->name()}) {
			$var_tmp = $family->setIntersectPatients($hPat->{$family->name()}, $chr);
		}
		else {
			foreach my $patient (@{$family->getPatients()}) {
				$var_tmp += $patient->getVariantsVector($chr);
			}
		}
		$var += $var_tmp;
	}
	$chr->setVariantsVector($var);
}

# methode intersect par genes (IND)
sub setIntersectPatient_GENES {
	my ( $chr, $patients ) = @_;
	return unless ($chr->getProject->level_ind() eq 'gene');
	return if (scalar(@$patients) == 0);
	if ($chr->getProject->typeFilters() eq 'familial') {  return setIntersectPatient_GENES_FAM($chr, $patients); }
	return setIntersectPatient_GENES_IND($chr, $patients);
}

# methode intersect par genes (IND)
sub setIntersectPatient_GENES_IND {
	my ($chr, $patients) = @_;
	return if (scalar(@$patients) == 0);
	my $var_ok = $chr->getNewVector();
	my $var_not = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		next if ($gene->is_intergenic());
		if ($gene->hasVariantsForAllPatients($patients)) { $var_ok += $gene->getVariantsVector(); }
		else {
			$var_not += $gene->getVariantsVector();
			$chr->supressGene($gene);
		}
	}
	$var_not -= $var_ok;
	$chr->getVariantsVector->AndNot($chr->getVariantsVector(), $var_not);
	foreach my $patient (@$patients) { $patient->intersected(1); }
	foreach my $patient (@{$chr->getPatients()}) { $patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $chr->getVariantsVector()); }
}

# methode intersect par genes (IND each FAM)
sub setIntersectPatient_GENES_FAM {
	my ( $chr, $patients ) = @_;
	return if (scalar(@$patients) == 0);
	my $hPatNames;
	my $hVarToSuppress;
	my $var_tmp = $chr->getNewVector();
	my $var_gene = $chr->getNewVector();
	foreach my $patient (@$patients) {
		$hPatNames->{$patient->getFamily->name()}->{$patient->name()} = undef;
		$hVarToSuppress->{$patient->getFamily->name()} = $chr->getNewVector();
	}
	my @lGenesImpacted;
	foreach my $gene (@{$chr->getGenes()}) {
		foreach my $family (@{$chr->getFamilies()}) {
			next unless (exists $hPatNames->{$family->name()});
			my $to_supress;
			foreach my $patient (@{$family->getPatients()}) {
				next unless (exists $hPatNames->{$family->name()}->{$patient->name()});
				$patient->intersected(1);
				$var_tmp->Intersection($gene->getVariantsVector(), $patient->getVariantsVector($chr));
				if ($var_tmp->is_empty()) {
					$hVarToSuppress->{$family->name()} += $gene->getVariantsVector();
					push(@lGenesImpacted, $gene);
					$to_supress = 1;
					last;
				}
			}
		}
	}
	foreach my $family (@{$chr->getFamilies()}) {
		next unless (exists $hVarToSuppress->{$family->name()});
		foreach my $patient (@{$family->getPatients()}) {
			$patient->getVariantsVector($chr)->AndNot($patient->getVariantsVector($chr), $hVarToSuppress->{$family->name()});
		}
	}
	$chr->update_from_patients();
	foreach my $gene (@lGenesImpacted) {
		$gene->getVariantsVector()->Intersection($gene->getVariantsVector(), $chr->getVariantsVector());
		if ($gene->getVariantsVector->is_empty()) { $chr->supressGene($gene); }
	}
	return;
}

# method intersect pour une famille
sub setIntersectFamily {
	my ( $chr, $families ) = @_;
	return if ($chr->getProject->level_fam() eq 'gene');
	return if (scalar(@$families) == 0);
	my $hFamVar;
	foreach my $family (@$families) {
		my $fam_name = $family->name();
		$hFamVar->{$fam_name} = $chr->getNewVector();
		my $family = $chr->getFamily($fam_name);
		foreach my $patient (@{$family->getPatients()}) {
			$hFamVar->{$fam_name} += $patient->getVariantsVector($chr);
		}
		$family->intersected(1);
	}
	my $var_inter = $chr->getNewVector();
	my $i = 0;
	foreach my $fam_name (keys %$hFamVar) {
		if ($i == 0) { $var_inter += $hFamVar->{$fam_name}; }
		else { $var_inter->Intersection( $var_inter, $hFamVar->{$fam_name} ) }
		$i++;
	}
	foreach my $patient (@{$chr->getPatients()}) {
		$patient->getVariantsVector($chr)->Intersection( $var_inter, $patient->getVariantsVector($chr) );
	}
	$chr->update_from_patients();
}

# methode intersect par genes (FAM)
sub setIntersectFamily_GENES {
	my ($chr, $families) = @_;
	return unless ($chr->getProject->level_fam() eq 'gene');
	return if (scalar(@$families) == 0);
	my $var_ok  = $chr->getNewVector();
	my $var_tmp = $chr->getNewVector();
	my ($hVarFam, $hGenesids);
	foreach my $family (@$families) {
		$hVarFam->{$family->name()} = $chr->getNewVector();
		foreach my $patient (@{$family->getPatients()}) {
			$hVarFam->{$family->name()} += $patient->getVariantsVector($chr);
		}
		$family->intersected(1);
	}
	foreach my $gene (@{$chr->getGenes()}) {
		my $is_del;
		foreach my $family (@$families) {
			$var_tmp->Intersection($gene->getVariantsVector(), $hVarFam->{$family->name()});
			if ($var_tmp->is_empty()) {
				$is_del = 1;
				last;
			}
		}
		if ($is_del) { $chr->supressGene($gene); }
		else { $var_ok += $gene->getVariantsVector(); }
	}
	foreach my $patient (@{$chr->getPatients()}) { $patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $var_ok); }
	$chr->update_from_patients();
}

sub filter_atLeast {
	my ($chr, $atLeast, $typeFilters, $level_fam) = @_;
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
	my ($chr, $filter_text) = @_;
	return unless ($filter_text);
	$chr->getProject->filter_text($filter_text);
	my $vector_genes = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		$chr->getProject->print_dot(1);
		unless ($chr->check_filter_text($gene)) {
			$chr->supressGene($gene);
			next;
		}
		$vector_genes += $gene->getVariantsVector();
	}
	$vector_genes->Intersection($vector_genes, $chr->getVariantsVector());
	$chr->setVariantsVector($vector_genes);
}

sub filter_genes_only_genes_names {
	my ($chr, $only_genes) = @_;
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
		if ($chr->getProject->only_genes() and not exists $chr->getProject->only_genes->{uc($gene->external_name())}) {
			$chr->supressGene($gene);
			next;
		}
		else {
			$vector_genes += $gene->getVariantsVector();
		}
	}
#	$chr->setVariantsVector($vector_genes);
}

sub filter_genes_annotations {
	my ($chr, $hFiltersChr) = @_;
	return unless ($hFiltersChr);
	my $variants_genes  = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		$chr->getProject->print_dot(1);
		my $var_gene_annot = $gene->getCategoriesVariantsVector($hFiltersChr);
		$gene->getVariantsVector->Intersection($var_gene_annot, $gene->getVariantsVector());
		if ($gene->getVariantsVector->is_empty()) {
			$chr->supressGene($gene);
			next;
		}
		$variants_genes += $gene->getVariantsVector();
	}
	$chr->getVariantsVector->Intersection($chr->getVariantsVector(), $variants_genes);
}


1;
