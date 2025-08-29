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


has is_filter_gene_id => (
	is => 'rw',
	default	=> sub {
		return undef;
	},
); 

has hash_filter_gene => (
	is => 'rw',
	default	=> sub {
		return {};
	},
); 

sub getGenes {
	my ($self,$chr) = @_;
	confess() unless $chr;
	my @t = values %{$self->{genes}->{$chr->name}};
	return \@t;
}

sub setGenes {
	my ($self,$chr) = @_;
	unless (exists $self->{genes}->{$chr->name}){
			my $vvv = $chr->getVariantsVector();
		foreach my $g (@{$chr->getGenesFromVector($vvv)}){
			$self->{genes}->{$chr->name}->{$g->id} = $g;
		}
	 }
	
}
sub refreshGenes {
	my ($self,$chr,$vector) = @_;
	$vector = $chr->getVariantsVector unless $vector; 
	foreach my $g (@{$self->getGenes($chr)}){
		my $v1 = $g->getCurrentVector();
		$v1 &= $vector;
		if ($v1->is_empty){
			$self->deleteGene($chr,$g);
		}
		else {
			$g->getCurrentVector($v1);
		}
	}
}

sub deleteGene {
	my ($self,$chr,$gene) = @_;
	delete $self->{genes}->{$chr->name}->{$gene->id} ;
	return 1;
}

sub existsGene {
	my ($self,$chr,$gene) = @_;
	return exists  $self->{genes}->{$chr->name}->{$gene->id} ;
}

sub filter_vector_ratio {
	my ($self, $chr, $limit_ratio, $filter_type_ratio) = @_;
	
	
	return unless ($limit_ratio);
	return if ($limit_ratio eq 'all');
	my $vector_ok = $chr->getNewVector();
	foreach my $patient (@{$chr->project->getPatients()}) {
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
	my @lPat = @{$chr->getPatients};
	if ($atLeast > scalar(@lPat)) {
		$chr->setVariantsVector($chr->getNewVector());
		return;
	}
	my $vector_ok = $chr->getNewVector();
	foreach my $patient (@lPat) {
		$hVectPat->{$patient->name()} = $patient->getVariantsVector($chr);
		$hVectPat->{$patient->name()} &= $chr->getVariantsVector();
	}
	foreach my $index (@{$chr->getIdsBitOn($chr->getVariantsVector())}) {
		$chr->project->print_dot(50);
		my $ok = 0;
		foreach my $patient (@lPat) {
			if ($hVectPat->{$patient->name()}->contains($index)) {
				$ok++;
				if ($ok == $atLeast) {
					$vector_ok->Bit_On($index);
					last;
				}
			}
		}
	}
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
	if ($atLeast > scalar(@lFam)) {
		$chr->setVariantsVector($chr->getNewVector());
		return;
	}
	my $vector_ok = $chr->getNewVector();
	my $vo = $chr->getVariantsVector();
	foreach my $family (@lFam) {
		$hVectFam->{$family->name()} = $family->getCurrentVariantsVector($chr) & $vo; 
	}
	foreach my $index (@{$chr->getIdsBitOn($chr->getVariantsVector())}) {
		$chr->project->print_dot(50);
		my $ok = 0;
		foreach my $fam_name (keys %$hVectFam) {
			if ($hVectFam->{$fam_name}->contains($index)) {
				$ok++;
				if ($ok == $atLeast) {
					$vector_ok->Bit_On($index);
					last;
				}
			}
		}
	}
	foreach my $fam (@{$chr->getFamilies()}) {
		my $v_fam = $fam->getCurrentVariantsVector($chr) & $vector_ok;
		$fam->setCurrentVariantsVector($chr, $vector_ok);
	}
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
foreach my $gene (@{$self->getGenes($chr)}) {
		my $nb_ok = 0;
		foreach my $patient ($self->getPatients) {
			$var_tmp_atleast->Intersection($gene->getVariantsVector(), $patient->getVariantsVector($chr));
			unless ($var_tmp_atleast->is_empty()) {
				$nb_ok++;
				last if ($nb_ok == $atLeast);
				
			}
		}
		if ($nb_ok < $atLeast) { 
			$self->deleteGene($chr,$gene);
			delete $chr->{genes_object}->{$gene->id}; 
		}
		else { $variants_genes += $gene->getVariantsVector();}
		}
	my $v = $chr->getVariantsVector() & $variants_genes;
	$chr->setVariantsVector($v);
		
	foreach my $patient (@{$chr->getPatients()}) {
		$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $chr->getVariantsVector());
	}
}

# AT_LEAST par familles en mode GENES et FAM
sub atLeastFilter_genes_fam {
	my ($self, $chr, $atLeast) = @_;
	my $variants_genes  = $chr->getNewVector();
	foreach my $gene (@{$self->getGenes($chr)}) {
		my $nb_ok = 0;
		my $var_tmp_atleast = $chr->getNewVector();
		foreach my $family (@{$chr->getFamilies}) {
			$var_tmp_atleast->Intersection($gene->getVariantsVector(), $family->getVariantsVector($chr));
			unless ($var_tmp_atleast->is_empty()) {
				$nb_ok++;
				last if ($nb_ok == $atLeast);
			}
		}
		if ($nb_ok < $atLeast) { 
			$self->deleteGene($chr,$gene);
			$chr->supressGene($gene); 
			}
		else { $variants_genes += $gene->getVariantsVector(); }
	}
	my $v = $chr->getVariantsVector() & $variants_genes;
	$chr->setVariantsVector($v);
	
}

sub atLeastFilter_genes_som{
	my ($self, $chr, $atLeast) = @_;
	my $variants_genes  = $chr->getNewVector();
	foreach my $gene (@{$self->getGenes($chr)}) {
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
		if ($nb_ok < $atLeast){ 
			$self->deleteGene($chr,$gene);
			$chr->supressGene($gene); 
			}
		else { $variants_genes += $gene->getVariantsVector(); }
	}
	my $v = $chr->getVariantsVector() & $variants_genes;
	$chr->setVariantsVector($v);
	
	foreach my $patient (@{$chr->getPatients()}) {
		$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $chr->getVariantsVector());
	}
}


######### FILTRES MODELS CHR ##########

#our $hModelsMethodsNames_familial = {
#   'denovo' => 'isDenovoTransmission',
#   'strict-denovo' => 'isStrictDenovoTransmission',
#   'dominant' => 'isDominantTransmission',
#   'mosaic' => 'isMosaicTransmission',
#   'recessif' => 'isRecessiveTransmission',
##   'incomplete_penetrance' => 'getModelVector_fam_incomplete_penetrance',
#   'uniparental_disomy' => 'isUniparentalDisomyTransmission',
#   
#   'only_tissues_somatic' => 'getModelVector_som_only_tissues_somatic',
#   
#   
#   'recessif_vector' => 'getVector_individual_recessive',
#   'dominant_vector' => 'getVector_individual_dominant',
#   'denovo_vector' => 'getVector_individual_denovo',
#   'uniparental_disomy_vector' => 'getVector_individual_uniparental_disomy',
#};
#
#sub filter_models_familial_union {
#	my ($self, $chr, $list_models, $h_arguments) = @_;
#	confess("\n\nERROR: QueryVectorFilter::filter_models_familial_union need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
#	confess("\n\nERROR: QueryVectorFilter::filter_models_familial_union need a list of models. Die.\n\n") unless ($list_models);
#	my $vector_models = $self->get_models_familial_union($chr, $list_models, $h_arguments);
##	foreach my $patient (@{$chr->getPatients()}) {
##		my $vector_patient = $chr->getNewVector();
##		foreach my $model (@$list_models) {
##			confess("\n\nERROR: $model model not used for patient ".$patient->name()." for CHR".$chr->id()."... Die.\n\n") unless (exists $patient->hash_models_genetics_used->{$model}->{$chr->id()});
##			$vector_patient += $patient->hash_models_genetics_used->{$model}->{$chr->id()};
##		}
##		$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $vector_patient);
##	}
#	$chr->setVariantsVector($vector_models);
#}
#
#sub get_models_familial_union {
#	my ($self, $chr, $list_models, $h_arguments) = @_;
#	confess("\n\nERROR: QueryVectorFilter::get_models_familial_union need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
#	confess("\n\nERROR: QueryVectorFilter::get_models_familial_union need a list of models. Die.\n\n") unless ($list_models);
#	my $vector_models = $chr->getNewVector();
#	foreach my $model (@$list_models) {
#		#isMosaicTransmission
#		if ($model eq 'compound') {
#			$vector_models += $self->get_model_familial_compound($chr, $h_arguments);
#		}
#		else {
#			my $model_method_name = $hModelsMethodsNames_familial->{$model};
#			my $model_method_vector_name;
#			$model_method_vector_name = $hModelsMethodsNames_familial->{$model.'_vector'} if exists $hModelsMethodsNames_familial->{$model.'_vector'};
#			foreach my $family (@{$chr->getProject->getFamilies()}) {
#				next if $family->in_the_attic();
#				if ($model eq 'dominant') { $family->{isDominant} = 1; }
#				foreach my $child (@{$family->getChildrenIll()}) {
#					next if $child->in_the_attic();
#					my $var_pat = $child->getVectorOrigin($chr);
#					$var_pat->Intersection($var_pat, $chr->getVariantsVector());
#					if ($model_method_vector_name) {
#						$var_pat->Intersection($var_pat, $family->$model_method_vector_name($chr, $child));
#						$vector_models += $var_pat;
#					}
#					else {
#						foreach my $var (@{$chr->getListVarObjects($var_pat)}) {
#							$vector_models->Bit_On($var->vector_id()) if $var->$model_method_name($family, $child);
#						}
#					}
#				}
#			}
#		}
#	}
#	$chr->setVariantsVector($vector_models);
#	return $vector_models;
#}
#
#sub filter_model_familial_denovo {
#	my ($self, $chr) = @_;
#	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_denovo need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
#	return filter_model_familial_common($chr, 'denovo');
#}
#
#sub get_model_familial_denovo {
#	my ($self, $chr) = @_;
#	confess("\n\nERROR: QueryVectorFilter::get_model_familial_denovo need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
#	return get_model_familial_common($chr, 'denovo');
#}
#
#sub filter_model_familial_strict_denovo {
#	my ($self, $chr) = @_;
#	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_strict_denovo need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
#	return filter_model_familial_common($chr, 'strict-denovo');
#}
#		
#sub get_model_familial_strict_denovo {
#	my ($self, $chr) = @_;
#	confess("\n\nERROR: QueryVectorFilter::get_model_familial_strict_denovo need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
#	return get_model_familial_common($chr, 'strict-denovo');
#}
#
#sub filter_model_familial_dominant {
#	my ($self, $chr) = @_;
#	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_dominant need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
#	return filter_model_familial_common($chr, 'dominant');
#}
#		
#sub get_model_familial_dominant {
#	my ($self, $chr) = @_;
#	confess("\n\nERROR: QueryVectorFilter::get_model_familial_dominant need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
#	return get_model_familial_common($chr, 'dominant');
#}
#
#sub filter_model_familial_mosaic {
#	my ($self, $chr) = @_;
#	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_mosaic need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
#	return filter_model_familial_common($chr, 'mosaic');
#}
#		
#sub get_model_familial_mosaic {
#	my ($self, $chr) = @_;
#	confess("\n\nERROR: QueryVectorFilter::get_model_familial_mosaic need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
#	return get_model_familial_common($chr, 'mosaic');
#}
#
#sub filter_model_familial_recessif {
#	my ($self, $chr) = @_;
#	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_recessif need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
#	return filter_model_familial_common($chr, 'recessif');
#}
#		
#sub get_model_familial_recessif {
#	my ($self, $chr) = @_;
#	confess("\n\nERROR: QueryVectorFilter::get_model_familial_recessif need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
#	return get_model_familial_common($chr, 'recessif');
#}
#
#sub filter_model_familial_uniparental_disomy {
#	my ($self, $chr) = @_;
#	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_uniparental_disomy need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
#	return filter_model_familial_common($chr, 'uniparental_disomy');
#}
#		
#sub get_model_familial_uniparental_disomy {
#	my ($self, $chr) = @_;
#	confess("\n\nERROR: QueryVectorFilter::get_model_familial_uniparental_disomy need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
#	return get_model_familial_common($chr, 'uniparental_disomy');
#}
#
#sub filter_model_familial_compound {
#	my ($self, $chr, $h_arguments) = @_;
#	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_compound need a GenBoChromosomeCache. Die.\n\n") unless ($chr);
#	confess("\n\nERROR: QueryVectorFilter::filter_model_familial_compound need a Hash of Arguments. Die.\n\n") unless ($h_arguments);
#	return filter_model_familial_common($chr, 'compound');
#}
#
#sub get_model_familial_compound {
#	my ($self, $chr, $h_arguments) = @_;
#	my $variants_genes  = $chr->getNewVector();
#	foreach my $gene (@{$chr->getGenes()}) {
#		my $vector_gene = $self->filter_gene_model_familial_compound($gene, $h_arguments);
#		$variants_genes += $vector_gene;
#	}
#	$variants_genes->Intersection($chr->getVariantsVector(), $variants_genes);
#	my $vector_chr = $chr->getNewVector();
#	foreach my $family (@{$chr->getFamilies()}) {
#		my $vector_fam = $family->get_vector_keep_var_compound($chr)->Clone();
#		$vector_fam->Intersection($vector_fam, $variants_genes);
#		foreach my $patient (@{$family->getPatients}) {
#			my $vector_patient = $patient->getVariantsVector($chr)->Clone();
#			$vector_patient->Intersection($vector_patient, $vector_fam);
#			$vector_chr += $vector_patient;
#			$patient->{hash_models_genetics_used}->{'compound'}->{$chr->id()} = $vector_patient->Clone();
#		}
#		$family->{hash_models_genetics_used}->{'compound'}->{$chr->id()} = $vector_fam->Clone();
#	}
#	return $vector_chr;
#}

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

#sub filter_model_familial_common {
#	my ($self, $chr, $model) = @_;
#	return unless ($model);
#	my $vector_model = get_model_familial_common($chr, $model);
#	my $vector_chr = $chr->getNewVector();
#	foreach my $patient (@{$chr->getPatients()}) {
#		confess("\n\nERROR: $model model not used for patient ".$patient->name()." for CHR".$chr->id()."... Die.\n\n") unless (exists $patient->hash_models_genetics_used->{$model}->{$chr->id()});
#		my $vector_patient_model = $patient->hash_models_genetics_used->{$model}->{$chr->id()};
#		$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $vector_patient_model);
#		$vector_chr += $patient->getVariantsVector($chr);
#	}
#	$chr->setVariantsVector($vector_chr);
#}
#
#sub get_model_familial_common {
#	my ($self, $chr, $model, $h_arguments) = @_;
#	return unless ($model);
#	my $model_method_name = $hModelsMethodsNames_familial->{$model};
#	return unless ($model_method_name);
#	my $vector_update = $chr->getNewVector();
#	foreach my $family (@{$chr->getFamilies()}) {
#		my ($self, $has_patients, $has_parents); 
#		$family->{hash_models_genetics_used}->{$model}->{$chr->id()} = $family->$model_method_name($chr, $h_arguments);
#		foreach my $patient (@{$family->getPatients()}) {
#			my $vector_patient = $patient->getVariantsVector($chr)->Clone();
#			$vector_patient->Intersection( $vector_patient, $family->{hash_models_genetics_used}->{$model}->{$chr->id()});
#			$vector_patient->Intersection($vector_patient, $chr->getVariantsVector());
#			$vector_update += $vector_patient;
#			$patient->{hash_models_genetics_used}->{$model}->{$chr->id()} = $vector_patient->Clone();
#			$has_patients = 1;
#		}
#	}
#	return $vector_update;
#}




##########################################


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






######### FILTRES PATIENTS / FAMILIES ##########

my $has_intersected_or_excuded_patients;

sub setPatients {
	my ($self,$project) = @_;
	foreach my $p (@{$project->getPatients}){
		$self->{patients}->{$p->id} = $p;
	}
}

sub patients {
	my ($self,$project) = @_;
	confess unless $self->{patients};
	return values %{$self->{patients}};
}
sub delPatient {
	my ($self,$patient) = @_;
	delete $self->{patients}->{$patient->id};
	
	return values %{$self->{patients}};
}

# method in the attic pour un patient 

sub setIntheAttic {
	my ($self, $patients) = @_;
	return if (scalar(@$patients) == 0);
	foreach my $patient (@$patients) {
		$self->delPatient($patient);
		$patient->in_the_attic(1);
		push(@{$self->{in_the_attic}->{patients}},$patient);
	}
}

sub FilterInTheAtticPatients {
	my ($self, $chr) = @_;
	return if exists $self->{in_the_attic}->{patients};
	
	my $v1 = $chr->getVariantsVector();
	my $vfinal = $chr->getNewVector(); ;
	foreach my $patient ($self->patients) {
		my $v = $patient->getVectorOrigin($chr);
		$vfinal += $v;
	}
	$chr->setVariantsVector($vfinal);
	
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
			warn $patient->name;
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

#sub setIntersectExclude_PAT_FAM {
#	my ($self,  $chr ) = @_;
#	return if not $has_intersected_or_excuded_patients;
#	return $self->setIntersectExclude_PAT_FAM_genes($chr) if ($chr->getProject->level_ind() eq 'gene');
#	return $self->setIntersectExclude_PAT_FAM_genes($chr) if ($chr->getProject->level_fam() eq 'gene');
#	my $var_ok = $chr->getNewVector();
#	my $var_excluded = $chr->getNewVector();
#	my $var_intersect = $chr->getNewVector();
#	my $nb_intersected = 0;
#	if ($chr->getProject->typeFilters() eq 'individual') {
#		my ($var_local_ok, $var_local_excluded, $var_local_intersect) = $self->getVector_project($chr, $chr->getVariantsVector());
#		$var_ok += $var_local_ok;
#		$var_ok -= $var_local_excluded;
#		$var_excluded += $var_local_excluded;
#		$var_intersect += $var_local_intersect if $var_local_intersect;
#		my @lPat = @{$chr->getProject->getPatients()};
#		$self->updateVectorPatients(\@lPat, $chr, $var_local_ok, $var_local_excluded, $var_local_intersect);
#		
#	}
#	else {
#		foreach my $family (@{$chr->getProject->getFamilies()}) {
#			my ($var_local_ok, $var_local_excluded, $var_local_intersect) = $self->getVector_fam($chr, $chr->getVariantsVector(), $family);
#			if ($family->in_the_attic()) {
#				next;
#			}
#			elsif ($family->excluded()) {
#				$var_excluded += $var_local_ok;
#			}
#			elsif ($family->intersected() and $nb_intersected == 0) {
#				$nb_intersected++;
#				$var_intersect += $var_local_ok;
#			}
#			elsif ($family->intersected()) {
#				$var_intersect->Intersection($var_local_ok, $var_intersect);
#			}
#			$var_ok += $var_local_ok;
#			my @lPat = @{$family->getPatients()};
#			$self->updateVectorPatients(\@lPat, $chr, $var_local_ok, $var_local_excluded, $var_local_intersect);
#		}
#	}
#	$var_ok -= $var_excluded;
#	if ($nb_intersected) {
#		$var_ok->Intersection($var_ok, $var_intersect);
#	}
#	$chr->setVariantsVector($var_ok);
#}

#setPatients 
# intersect patient 


sub setIntersectedPatients {
	my ($self, $patients, $typeFilters) = @_;
	return if (scalar(@$patients) == 0);
	foreach my $patient (@$patients) {
		$patient->intersected(1);
		if ($patient->project->level_ind() eq 'gene'){
			if ($typeFilters eq 'familial') {
				push(@{$self->{genes}->{intersected}->{patients}->{$patient->getFamily->name()}},$patient);
			}
			else {
				push(@{$self->{genes}->{intersected}->{patients}},$patient);
			}
		}
		else {
			if ($typeFilters eq 'familial') {
				push(@{$self->{variants}->{intersected}->{patients}->{$patient->getFamily->name()}},$patient);
			}
			else {
				push(@{$self->{variants}->{intersected}->{patients}},$patient);
			}
		}
		
	}
	return;
}


# excluded  patient 

sub setExcludedPatients {
	my ($self, $patients, $typeFilters) = @_;
	return if (scalar(@$patients) == 0);
	
	foreach my $patient (@$patients) {
		$patient->excluded('all');
		if ($patient->getProject->level_ind() eq 'gene'){
			#push(@{$self->{genes}->{excluded}->{patients}},$patient);
			if ($typeFilters eq 'familial') {
				push(@{$self->{genes}->{excluded}->{patients}->{$patient->getFamily->name()}},$patient);
			}
			else {
				push(@{$self->{genes}->{excluded}->{patients}},$patient);
			}
		}
		else {
			#push(@{$self->{variants}->{excluded}->{patients}},$patient);
			if ($typeFilters eq 'familial') {
				push(@{$self->{variants}->{excluded}->{patients}->{$patient->getFamily->name()}},$patient);
			}
			else {
				push(@{$self->{variants}->{excluded}->{patients}},$patient);
			}
		}
		
	}
	return;
}
#he 

sub setFilteredHoPatients {
	my ($self, $patients, $typeFilters) = @_;
	return if (scalar(@$patients) == 0);
	
	foreach my $patient (@$patients) {
		$patient->excluded('ho');
		#push(@{$self->{variants}->{filtered_ho}->{patients}},$patient);
		if ($typeFilters eq 'familial') {
			push(@{$self->{variants}->{filtered_ho}->{patients}->{$patient->getFamily->name()}},$patient);
		}
		else {
			push(@{$self->{variants}->{filtered_ho}->{patients}},$patient);
		}
		
	}
	return;
}

#ho 
sub setFilteredHePatients {
	
	my ($self,  $patients, $typeFilters) = @_;
	return if (scalar(@$patients) == 0);
	
	foreach my $patient (@$patients) {
		$patient->excluded('he');
		#push(@{$self->{variants}->{filtered_he}->{patients}},$patient);
		if ($typeFilters eq 'familial') {
			push(@{$self->{variants}->{filtered_he}->{patients}->{$patient->getFamily->name()}},$patient);
		}
		else {
			push(@{$self->{variants}->{filtered_he}->{patients}},$patient);
		}
		
	}
	return;
}

# intersect families 
sub setIntersectedFamilies {
	my ($self, $families) = @_;
	return if (scalar(@$families) == 0);
	
	foreach my $fam (@$families) {
		
		$fam->intersected(1);
		if ($fam->getProject->level_fam() eq 'gene'){
				push(@{$self->{genes}->{intersected}->{families}},$fam);
		}
		elsif ($fam->getProject->level_fam() eq 'variation'){
				push(@{$self->{variants}->{intersected}->{families}},$fam);
		}
	}
}


sub setExcludedFamilies {
	my ($self, $families) = @_;
	return if (scalar(@$families) == 0);
	
	foreach my $fam (@$families) {
		$fam->excluded(1);
		if ($fam->getProject->level_fam() eq 'gene'){
				push(@{$self->{genes}->{excluded}->{families}},$fam);
		}
		elsif ($fam->getProject->level_fam() eq 'variation'){
				push(@{$self->{variants}->{excluded}->{families}},$fam);
		}
	}
}




#########
# END SET
###########

#filter by variants 
#    intersect Patient 
sub intersectVariantsPatients {
	my ($self, $chr, $typeFilters) = @_;
	return unless $self->{variants}->{intersected}->{patients};
	if ($typeFilters eq 'familial') {
		foreach my $fam_name (keys %{$self->{variants}->{intersected}->{patients}}) {
			my $fam = $chr->getFamily($fam_name);
			my $var_ok = $fam->getVariantsVector($chr);
			foreach my $patient (@{$self->{variants}->{intersected}->{patients}->{$fam_name}}) {
				$var_ok &= $patient->getVectorOrigin($chr);
			}
			$fam->setCurrentVariantsVector($chr, $var_ok);
		}
	}
	else {
		my $var_ok = $chr->getVariantsVector();
		foreach my $patient (@{$self->{variants}->{intersected}->{patients}}) {
			$var_ok &= $patient->getVectorOrigin($chr);
		}
		$chr->setVariantsVector($var_ok);
	}
}



#    intersect Fam 
sub intersectVariantsFamilies {
	my ($self, $chr ) = @_;
	return unless $self->{variants}->{intersected}->{families};
	my $var_ok = $chr->getVariantsVector();
	foreach my $fam (@{$self->{variants}->{intersected}->{families}}) {
		$var_ok &= $fam->getCurrentVariantsVector($chr);
	}
	foreach my $fam (@{$self->{variants}->{intersected}->{families}}) {
		my $v_fam = $fam->getCurrentVariantsVector($chr);
		$v_fam &= $var_ok;
		$fam->setCurrentVariantsVector($chr, $v_fam);
	}
	$chr->setVariantsVector($var_ok);
}




#    exclude Patient 
sub excludeVariantsPatients {
	my ($self, $chr, $typeFilters) = @_;
	return unless $self->{variants}->{excluded}->{patients};
	if ($typeFilters eq 'familial') {
		my $var_del = $chr->getNewVector();
		foreach my $fam_name (keys %{$self->{variants}->{excluded}->{patients}}) {
			my $fam = $chr->getFamily($fam_name);
			my $var_ok = $fam->getVariantsVector($chr);
			foreach my $patient (@{$self->{variants}->{excluded}->{patients}->{$fam_name}}) {
				$var_ok -= $patient->getVectorOrigin($chr);
			}
			$fam->setCurrentVariantsVector($chr, $var_ok);
		}
	}
	else {
		my $var_ok = $chr->getVariantsVector();
		foreach my $patient (@{$self->{variants}->{excluded}->{patients}}) {
			$var_ok -= $patient->getVectorOrigin($chr);
		}
		$chr->setVariantsVector($var_ok);
	}
}

#ho

sub filterHoVariantsPatients {
	my ($self,  $chr, $typeFilters ) = @_;
	return unless $self->{variants}->{filtered_ho}->{patients};
	if ($typeFilters eq 'familial') {
		foreach my $fam_name (keys %{$self->{variants}->{filtered_ho}->{patients}}) {
			my $fam = $chr->getFamily($fam_name);
			my $v_fam = $fam->getCurrentVariantsVector($chr);
			foreach my $pat (@{$self->{variants}->{filtered_ho}->{patients}->{$fam_name}}) {
				$v_fam -= $pat->getHo($chr);
			}
			$fam->setCurrentVariantsVector($chr, $v_fam);
		}
	}
	else {
		my $var_ok = $chr->getVariantsVector();
		foreach my $patient (@{$self->{variants}->{filtered_ho}->{patients}}) {
			$var_ok -= $patient->getVectorOriginHo($chr);
		}
		$chr->setVariantsVector($var_ok);
	}
}

sub filterHeVariantsPatients {
	my ($self,  $chr, $typeFilters ) = @_;
	return unless $self->{variants}->{filtered_he}->{patients};
	if ($typeFilters eq 'familial') {
		foreach my $fam_name (keys %{$self->{variants}->{filtered_he}->{patients}}) {
			my $fam = $chr->getFamily($fam_name);
			my $v_fam = $fam->getCurrentVariantsVector($chr);
			foreach my $pat (@{$self->{variants}->{filtered_he}->{patients}->{$fam_name}}) {
				$v_fam -= $pat->getHe($chr);
			}
			$fam->setCurrentVariantsVector($chr, $v_fam);
		}
	}
	else {
		my $var_ok = $chr->getVariantsVector();
		foreach my $patient (@{$self->{variants}->{filtered_he}->{patients}}) {
			$var_ok -= $patient->getVectorOriginHe($chr);
		}
		$chr->setVariantsVector($var_ok);
	}
}

#    exclude Families 
sub excludeVariantsFamilies {
	my ($self,  $chr ) = @_;
	return unless $self->{variants}->{excluded}->{families};
	my $var_del = $chr->getNewVector();
	foreach my $fam (@{$self->{variants}->{excluded}->{families}}) {
		foreach my $patient (@{$fam->getPatients()}) {
			$var_del += $patient->getVectorOrigin($chr);
		}
		$fam->setCurrentVariantsVector($chr, $chr->getNewVector());
	}
	foreach my $fam (@{$chr->getFamilies()}) {
		my $v_fam = $fam->getCurrentVariantsVector($chr);
		$v_fam -= $var_del;
		$fam->setCurrentVariantsVector($chr, $v_fam);
	}
}


#filter by Genes 
#	intersect patient

sub intersectGenesPatients {
	my ($self,  $chr ) = @_;
	return unless $self->{genes}->{intersected}->{patients};
	my $var_ok = $chr->getNewVector();
	my (@intersect) = @{$self->{genes}->{intersected}->{patients}};
	
	foreach my $gene (@{$self->getGenes($chr)}) {
		my $ok;
		foreach my $patient (@intersect) {
				my $v = $gene->getCurrentVector() & $patient->getVectorOrigin($chr);
				if ($v->is_empty){
					$ok =0;
					last;
				}
				$ok ++;
		}
		unless ($ok){
			delete $chr->{genes_object}->{$gene->id};
			$self->deleteGene($chr,$gene);
			next;
		}
		$var_ok += $gene->getCurrentVector();
		
	}
	$chr->setVariantsVector($var_ok);
}
#	exclude patient
sub excludeGenesPatients {
	my ($self,  $chr ) = @_;
	return unless $self->{genes}->{excluded}->{patients};
	my $var_ok = $chr->getNewVector();
	my (@intersect) = @{$self->{genes}->{excluded}->{patients}};
	foreach my $gene (@{$self->getGenes($chr)}) {
		my $ok;
		foreach my $patient (@intersect) {
				my $v = $gene->getCurrentVector() & $patient->getVectorOrigin($chr);
				unless ($v->is_empty){
					$ok =0;
					last;
				}
				$ok ++;
		}
		unless ($ok){
			delete $chr->{genes_object}->{$gene->id};
			$self->deleteGene($chr,$gene);
			next;
		}
		$var_ok += $gene->getCurrentVector();
		
	}
	$chr->setVariantsVector($var_ok);
}


#	intersect Fam

sub intersectGenesFamilies {
	my ($self,  $chr ) = @_;
	return unless $self->{genes}->{intersected}->{families};
	my $var_ok = $chr->getNewVector();
	my (@intersect) = @{$self->{genes}->{intersected}->{families}};
	foreach my $gene (@{$self->getGenes($chr)}) {
		my $ok;
		foreach my $fam (@intersect) {
				my $v = $gene->getCurrentVector() & $fam->getCurrentVariantsVector($chr);
				if ($v->is_empty){
					$ok =0;
					last;
				}
				$ok ++;
		}
		unless ($ok){
			delete $chr->{genes_object}->{$gene->id};
			$self->deleteGene($chr,$gene);
			next;
		}
		$var_ok += $gene->getCurrentVector();
		
	}
	$chr->setVariantsVector($var_ok);
}

#	exclud  Fam


sub excludeGenesFamilies {
	my ($self,  $chr ) = @_;
	return unless $self->{genes}->{excluded}->{families};
	my $var_ok = $chr->getNewVector();
	my (@intersect) = @{$self->{genes}->{excluded}->{families}};
	
	foreach my $gene (@{$self->getGenes($chr)}) {
		my $ok;
		foreach my $fam (@intersect) {
				my $v = $gene->getCurrentVector() & $fam->getVectorOrigin($chr);
				unless ($v->is_empty){
					$ok =0;
					last;
				}
				$ok ++;
		}
		unless ($ok){
			delete $chr->{genes_object}->{$gene->id};
			$self->deleteGene($chr,$gene);
			next;
		}
		$var_ok += $gene->getCurrentVector();
	}
	$chr->setVariantsVector($var_ok);
}



sub filter_atLeast {
	my ($self, $chr, $atLeast, $typeFilters, $level_fam) = @_;
	if ($atLeast and $atLeast >= 2) {
		if ($typeFilters eq 'familial') {
			if ($level_fam eq 'gene') { $self->atLeastFilter_genes_fam($chr, $atLeast); }
			else { $self->atLeastFilter_var_fam($chr, $atLeast); }
		}
		if ($typeFilters eq 'somatic') {
			if ($level_fam eq 'gene') { $self->atLeastFilter_genes_som($chr, $atLeast); }
			else { $self->atLeastFilter_var_som($chr, $atLeast); }
		}
		else {
			if ($level_fam eq 'gene') { $self->atLeastFilter_genes_ind($chr, $atLeast); }
			else { $self->atLeastFilter_var_ind($chr, $atLeast); }
		}
	}
	return;
}

sub filter_genes_text_search {
	my ($self, $chr, $filter_text) = @_;
	return unless ($filter_text);
	$chr->getProject->filter_text($filter_text);
	my $vector_genes = $chr->getNewVector();
	foreach my $gene (@{$self->getGenes($chr)}) {
		$chr->getProject->print_dot(1);
		next if not $chr->check_filter_text($gene);
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
	my $vector_genes = $chr->getVariantsVector();
	$self->{genes}->{$chr->name} ={};
	foreach my $g (keys %{$chr->getProject->only_genes}){
		$chr->getProject->print_dot(1);
		
		my $gene= $chr->project->newGene($g);

		next if $gene->getChromosome->name ne $chr->name;
#		warn $g."=====";
		my $vgene =  $gene->getVariantsVector();
		$vgene &=  $vector_genes;
		next if $vgene->is_empty();
		$vector_genes &= $gene->getVariantsVector();
		$chr->{genes_object}->{$gene->id} ++;
		$self->{genes}->{$chr->name}->{$gene->id} = $gene;
	
	}
	$chr->setVariantsVector($vector_genes);
	return;
	
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
	my %hcat = %{$project->buffer->config->{'genbo_annotations_names'}};

	foreach my $c (@unselect){
		delete $hcat{$c};
	}
	return keys %hcat;
}
sub filter_genes_annotations {
	my ($self, $chr, $hFiltersChr) = @_;
	
	unless ($hFiltersChr) {
		$chr->getVariantsVector();
		$self->setGenes($chr);
		return;
	}

	my $tt =time;
	

	my @all_cat = return_cat($chr->project,keys %$hFiltersChr);
	
#	warn Dumper @all_cat;
#	die;
	
	my $variants_genes  = $chr->getNewVector();
	my $vvv = $chr->getVariantsVector();
	#$chr->rocks_vector("r")->get_vector_chromosome("coding") & $chr->getVariantsVector();
	
	my $cc = $chr->getVariantsVector()->Clone;
	my $nb = 0;
	my $nb_genes = 0;
	my $id_genes = {};
	my $t =time;
	#$rocks->prepare([keys %$lh]);
	$chr->{genes_object} = {};
	my @genes ;
	
	$self->setGenes($chr);
	 	
	
	foreach my $gene (@{$self->getGenes($chr)}) {
		$nb_genes ++;
		$chr->getProject->print_dot(1);
		my $debug;
		my $vsmall = $gene->getCompactVectorOriginCategories(\@all_cat,$debug);
		
		my $vchr = $gene->return_compact_vector( $chr->getVariantsVector());
		$vsmall &= $vchr;
		
		
#		warn $gene->external_name().' -> '.$vsmall->Norm();
#		warn Dumper @all_cat if $gene->external_name() eq 'AEN';		
#		die if $gene->external_name() eq 'AEN';
		#die if $gene->external_name() eq 'ACYP2';
		
		if ($vsmall->is_empty()){
			delete $chr->project->{genes_object}->{$gene->id};
			$self->deleteGene($chr,$gene);
			next;
		}
		$gene->setCurrentVector($vsmall);
		$chr->{genes_object}->{$gene->id} ++;
		$nb ++;
		$id_genes->{$gene->id} ++;
		$variants_genes += $gene->enlarge_compact_vector($vsmall);
	}
	$chr->getVariantsVector->Intersection($chr->getVariantsVector(), $variants_genes);
	$chr->{buffer_vector} = $chr->getVariantsVector();
	
#	die;
	
#	warn "---------------- $nb/$nb_genes/".scalar(keys %{$chr->{genes_object}});
#	warn abs(time - $t);
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_genes_annotations - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
}


sub filter_model_individual_recessif {
	my ($self, $chr) = @_;
	my $variants_genes  = $chr->getNewVector();
	foreach my $gene (@{$self->getGenes($chr)}) {
		next if $gene->getCurrentVector->is_empty();
		$self->filter_gene_model_individual_recessif($gene);
		
		if ($gene->getCurrentVector->is_empty()){
			$self->deleteGene($chr,$gene);
			next;
		}
		$variants_genes += $gene->getCurrentVector();
	}
	$chr->setVariantsVector($variants_genes);
	foreach my $patient (@{$chr->getPatients()}) {
		$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $patient->{model_indiv_recessif}->{$chr->id}) if exists ($patient->{model_indiv_recessif}->{$chr->id});
	}
}




sub filter_genetics_models {
	my ($self, $chr, $level_ind, $level_fam, $model, $h_args) = @_;
	return if $chr->getVariantsVector->is_empty();
	return if not $model;
	my $vector_models = $chr->getNewVector();
	if ($model and $chr->getProject->typeFilters() eq 'individual') {
		$vector_models += $self->filter_genetics_models_individual($chr, $level_ind, $level_fam, $model, $h_args);
	}
	elsif ($model and $chr->getProject->typeFilters() eq 'familial') {
		$vector_models += $self->filter_genetics_models_familial($chr, $level_ind, $level_fam, $model, $h_args);
	}
	elsif ($model and $chr->getProject->typeFilters() eq 'somatic') {
		$vector_models += $self->filter_genetics_models_somatic($chr, $level_ind, $level_fam, $model, $h_args);
	}
	
	if ($level_ind eq 'variation' and $level_fam eq 'variation') {
		$chr->setVariantsVector($vector_models);
	}
	else {
		#TODO: voir avec en genes level_ind ou fam
		$chr->setVariantsVector($vector_models);
	}
	
	my $vc = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		next if not $gene->compact_vector();
		my $vv = $chr->getVariantsVector() & $gene->getCurrentVector();
		$gene->setCurrentVector($vv);
#		if ($model eq 'compound' and $vv->Norm() == 1) { $gene->setCurrentVector($chr->getNewVector()); }
		$vc += $vv;
	}
	$chr->setVariantsVector($vc);
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER models - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
}

sub filter_genetics_models_individual {
	my ($self, $chr, $level_ind, $level_fam, $model, $h_args) = @_;
	my $vector_models = $chr->getNewVector();
	if ($model eq 'recessif') {
		foreach my $pat (@{$chr->getPatients()}) {
			next if $pat->in_the_attic();
			my $v_pat = $chr->getNewVector();
			my $v_all = $chr->getVariantsVector() & $pat->getVariantsVector($chr);
			my $v_ho = $chr->getVariantsVector() & $pat->getHo($chr);
			my $v_he = $chr->getVariantsVector() & $pat->getHe($chr);
			foreach my $gene (@{$chr->getGenes()}) {
				my $v_gene = $gene->getVariantsVector() & $v_all;
				my $v_gene_ho = $gene->getVariantsVector() & $v_ho;
				my $v_gene_he = $gene->getVariantsVector() & $v_he;
				if ($v_gene_ho->Norm() >= 1 or $v_gene_he->Norm() >= 2) {
					$v_pat += $gene->getVariantsVector();
				}
			}
			$vector_models += $v_pat;
		}
	}
	if ($model eq 'compound') {
		foreach my $pat (@{$chr->getPatients()}) {
			next if $pat->in_the_attic();
			my $v_pat = $chr->getNewVector();
			my $v_all = $chr->getVariantsVector() & $pat->getVariantsVector($chr);
			my $v_he = $chr->getVariantsVector() & $pat->getHe($chr);
			foreach my $gene (@{$chr->getGenes()}) {
				my $v_gene = $gene->getVariantsVector() & $v_all;
				my $v_gene_he = $gene->getVariantsVector() & $v_he;
				if ($v_gene_he->Norm() >= 2) {
					$v_pat += $gene->getVariantsVector();
				}
			}
			$vector_models += $v_pat;
		}
	}
	return $vector_models;
}

sub filter_genetics_models_familial {
	my ($self, $chr, $level_ind, $level_fam, $model, $h_args) = @_;
	my $vector_models = $chr->getNewVector();
	if ($model eq 'compound' or $model eq 'recessif_compound' or $model eq 'uniparental_recessif_compound') {
		foreach my $fam (@{$chr->getProject->getFamilies()}) {
			$fam->{keep_var_compound}->{$chr->id()} = $chr->getNewVector();
			$fam->{hash_models_genetics_used}->{'compound'}->{$chr->id()} = $chr->getNewVector();
			if (not $fam->getMother() or not $fam->getFather()) {
				$fam->setCurrentVariantsVector($chr, $chr->getNewVector());
				next;
			}
			
			my $v_fam = $fam->getVariantsVector($chr);
			my $v_children_ill = $chr->getNewVector();
			foreach my $child (@{$fam->getChildrenIll()}) {
				$v_children_ill += $child->getHe($chr);
			}
			foreach my $child (@{$fam->getChildrenHealthy()}) {
				$v_children_ill -= $child->getVariantsVector($chr);
			}
			my $v_common_parents = $chr->getNewVector();
			$v_common_parents += $fam->getFather->getVariantsVector($chr);
			$v_common_parents &= $fam->getMother->getVariantsVector($chr);
			$v_children_ill -= $v_common_parents;
			
			my $v_mother = $fam->getVectorMotherTransmission($chr);
			$v_mother &= $fam->getMother->getHe($chr);
			$v_mother &= $v_children_ill;
			
			my $v_father = $fam->getVectorFatherTransmission($chr);
			$v_father &= $fam->getFather->getHe($chr);
			$v_father &= $v_children_ill;
			
			foreach my $gene (@{$chr->getGenes()}) {
				my $v_gene_m = $gene->getVariantsVector()->Clone();
				$v_gene_m &= $chr->getVariantsVector();
				$v_gene_m &= $v_mother;
				
				my $v_gene_f = $gene->getVariantsVector()->Clone();
				$v_gene_f &= $chr->getVariantsVector();
				$v_gene_f &= $v_father;
				
				if ($v_gene_m->Norm() >= 1 and $v_gene_f->Norm() >= 1) {
					$fam->{keep_var_compound}->{$chr->id()} += $v_gene_m;
					$fam->{keep_var_compound}->{$chr->id()} += $v_gene_f;
					#warn '--- '.$gene->external_name().' --- '.($v_gene_f->Norm() + $v_gene_m->Norm()).' in '.$fam->name;
				}
			}
#			warn '--- FAM '.$fam->name().' found '.$fam->keep_var_compound->{$chr->id()}->Norm().' ---';
			$fam->setCurrentVariantsVector($chr, $fam->{keep_var_compound}->{$chr->id()});
			$vector_models += $fam->getCurrentVariantsVector($chr);
		}
		foreach my $gene (@{$chr->getGenes()}) {
			my $v_gene = $gene->getVariantsVector()->Clone();
			$v_gene &= $vector_models;
			$gene->{variants} = $chr->getNewVector() if (not $v_gene->Norm >= 2);
		}
	}
	if ($model eq 'recessif' or $model eq 'recessif_compound' or $model eq 'uniparental_recessif_compound') {
		foreach my $fam (@{$chr->getProject->getFamilies()}) {
			my $v_fam = $chr->getNewVector();
			foreach my $child (@{$fam->getChildrenIll()}) {
				my $v_child = $fam->getVector_individual_recessive($chr, $child) & $child->getVariantsVector($chr);
				$v_fam += $v_child;
			}
			$v_fam &= $chr->getVariantsVector();
			if ($model eq 'recessif') { $fam->setCurrentVariantsVector($chr, $v_fam); }
			else { $fam->addCurrentVariantsVector($chr, $v_fam); }
			$vector_models += $v_fam;
		}
	}
	if ($model eq 'dominant') {
		foreach my $fam (@{$chr->getProject->getFamilies()}) {
			my $v_fam = $fam->getVariantsVector($chr);
			$v_fam &= $fam->getModelVector_fam_dominant($chr);
			$v_fam &= $chr->getVariantsVector();
			$fam->setCurrentVariantsVector($chr, $v_fam);
			$vector_models += $v_fam;
		}
	}
	if ($model eq 'denovo') {
		foreach my $fam (@{$chr->getProject->getFamilies()}) {
			my $v_fam = $fam->getCurrentVariantsVector($chr);
			$v_fam &= $fam->getModelVector_fam_denovo($chr);
			$v_fam &= $chr->getVariantsVector();
			$fam->setCurrentVariantsVector($chr, $v_fam);
			$vector_models += $v_fam;
		}
	}
	if ($model eq 'strict-denovo') {
		foreach my $fam (@{$chr->getProject->getFamilies()}) {
			my $v_fam = $fam->getVariantsVector($chr);
			$v_fam &= $fam->getModelVector_fam_denovo($chr);
			my $v_strict = $chr->getNewVector();
			foreach my $child (@{$fam->getChildren()}) {
				$v_strict += $fam->getVector_individual_strict_denovo($chr, $child);
			}
			$v_fam &= $v_strict;
			$v_fam &= $chr->getVariantsVector();
			$fam->setCurrentVariantsVector($chr, $v_fam);
			$vector_models += $v_fam;
		}
	}
	
	if ($model eq 'uniparental_disomy' or $model eq 'uniparental_recessif_compound') {
		foreach my $fam (@{$chr->getProject->getFamilies()}) {
			my $v_fam = $chr->getNewVector();
			foreach my $child (@{$fam->getChildren()}) {
				my $v_child = $fam->getVector_individual_uniparental_disomy($chr, $child) & $child->getVariantsVector($chr);
				$v_fam += $v_child;
			}
			$v_fam &= $chr->getVariantsVector();
			if ($model eq 'uniparental_disomy') { $fam->setCurrentVariantsVector($chr, $v_fam); }
			else { $fam->addCurrentVariantsVector($chr, $v_fam); }
			$vector_models += $v_fam;
		}
	}
	
	#TODO: a terminer	
	if ($model eq 'mosaic') {
		foreach my $fam (@{$chr->getProject->getFamilies()}) {
			my $v_mosaic = $fam->getModelVector_fam_mosaique($chr);
			$v_mosaic &= $fam->getVariantsVector($chr);
			$fam->setCurrentVariantsVector($chr, $v_mosaic);
			$vector_models += $v_mosaic;
		}
	}
	return $vector_models;
}

sub filter_genetics_models_somatic {
	my ($self, $chr, $level_ind, $level_fam, $model, $h_args) = @_;
	my $vector_models = $chr->getNewVector();
	if ($model eq 'loh') {
		$vector_models += $chr->getModelVector_som_loh();
	}
	if ($model eq 'dbl_evt') {
		my $variants_genes  = $chr->getNewVector();
		foreach my $gene (@{$self->getGenes($chr)}) {
			my $v = $gene->getCurrentVector;
			$v &= $gene->getChromosome->getModelGeneVector_som_dbl_evt($gene);
			$gene->setCurrentVector($v);
			next if ($gene->getVariantsVector->is_empty());
			$variants_genes += $gene->getVariantsVector();
		}
		$chr->getVariantsVector->Intersection($chr->getVariantsVector(), $variants_genes);
		foreach my $patient (@{$chr->getPatients()}) {
			if ($patient->model_vector_var_ok->{$chr->id()}) {
				$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $patient->model_vector_var_ok->{$chr->id()});
				$vector_models += $patient->model_vector_var_ok->{$chr->id()};
				$chr->patients_categories->{$patient->name()} = $patient->model_vector_var_ok->{$chr->id()};
			}
			else {
				$patient->getVariantsVector($chr)->Empty();
			}
		}
		
	}
	if ($model eq 'only_tissues_somatic') {
		foreach my $fam (@{$chr->getProject->getFamilies()}) {
			my $v_fam = $fam->getModelVector_som_only_tissues_somatic($chr);
			$v_fam &= $fam->getVariantsVector($chr);
			$fam->setCurrentVariantsVector($chr, $v_fam);
			$vector_models += $v_fam;
		}
	}
	return $vector_models;
}

1;
