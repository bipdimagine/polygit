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
		$vvv &= $self->{regions_filtred}->{$chr->id()} if exists $self->{regions_filtred}->{$chr->id()};
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

sub filter_regions {
	my ($self, $chr, $filter_region) = @_;
	return unless ($filter_region);
	my @lFilters;
	my $isIntersect;
	my $hFilters;
	
	my $use_filter;
	my $var_del = $chr->getNewVector();
	my $var_inter = $chr->getNewVector();
	foreach my $this_filter (split(" ", $filter_region)) {
		print $chr->getProject->print_dot(1000);
		next if (exists $hFilters->{$this_filter});
		$hFilters->{$this_filter} = undef;
		my ($chrId, $start, $end, $include) = split(":", $this_filter);
		next unless ($chrId eq $chr->id());
		print 'b_r:'.$chr->getVariantsVector()->Norm();
		#$chr->check_each_var_filter_region($this_filter, $first_launch);
		my $var_filter = $chr->getFilterRegionVector($this_filter);
		if ($include eq '0' or $include eq '99') {
			$isIntersect = 1;
			$var_inter += $var_filter;
		}
		elsif ($include eq '-1') {
			$var_del += $var_filter;
		}
		$use_filter = 1;
	}
	return if not $use_filter;
	
	my $var_ok = $chr->getNewVector();
#	foreach my $fam (@{$chr->getProject->getFamilies()}) {
#		my $v_fam = $fam->getVariantsVector($chr);
#		$v_fam &= $var_inter;
#		$v_fam -= $var_del;
#		$var_ok += $v_fam;
#		$fam->setCurrentVariantsVector($chr, $v_fam);
#	}
	$var_ok += $chr->getVariantsVector();
	$var_ok &= $var_inter;
	$var_ok -= $var_del;
	$chr->setVariantsVector($var_ok);
	print 'a_r:'.$chr->getVariantsVector()->Norm();
	
	foreach my $fam (@{$chr->getFamilies()}) {
		my $v_fam = $fam->getCurrentVariantsVector($chr);
		$v_fam &= $var_ok;
		$fam->setCurrentVariantsVector($chr, $v_fam);
	}
	
	$self->{regions_filtred}->{$chr->id()} = $var_ok->Clone();
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
	my $var_ok = $chr->getNewVector();
	my $var_excluded = $chr->getNewVector();
	my $var_intersect = $chr->getNewVector();
	return ($vector_chr_gene_init,$var_excluded,$var_intersect) if $vector_chr_gene_init->is_empty();
	my $nb_intersected = 0;
	foreach my $patient (@{$family_or_project->getPatients()}) {
		next if ($patient->in_the_attic());
		my $var_pat = $patient->getVariantsVector($chr)->Clone();
		$var_pat->Intersection($var_pat, $vector_chr_gene_init);
		next if $var_pat->is_empty();
		if ($patient->excluded()) {
			my $status = $patient->excluded();
			if ($status eq 'all') { $var_excluded += $patient->getVectorOrigin( $chr ); }
			elsif ($status eq 'he') { $var_excluded += $patient->getVectorOriginHe( $chr ); }
			elsif ($status eq 'ho') { $var_excluded += $patient->getVectorOriginHo( $chr ); }
		}
		if ($patient->intersected() and $nb_intersected == 0) {
			warn $patient->name;
			$nb_intersected++;
			$var_intersect += $var_pat;
			$var_ok += $var_pat;
		}
		elsif ($patient->intersected()) {
			$var_intersect->Intersection($var_pat, $var_intersect);
			$var_ok += $var_pat;
		}
		$var_ok += $var_pat;
	}
	$var_ok -= $var_excluded;
	if ($nb_intersected) {
		$var_ok->Intersection($var_ok, $var_intersect);
	}
	
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
	return;
}


sub get_vector_filter_gene_annotations {
	my ($self, $chr, $hFiltersChr, $not_set_vector_gene) = @_;
	my $nb = 0;
	my $id_genes = {};
	my @all_cat = return_cat($chr->project,keys %$hFiltersChr);
	my $variants_genes  = $chr->getNewVector();
	foreach my $gene (@{$self->getGenes($chr)}) {
		if (exists $chr->getProject->{only_genes}) {
			my $found;
			$found = 1 if exists $chr->getProject->{only_genes}->{$gene->id()};
			$found = 1 if exists $chr->getProject->{only_genes}->{$gene->external_name()};
			if (not $found) {
				delete $chr->project->{genes_object}->{$gene->id};
				$self->deleteGene($chr,$gene);
				next;
			}
		}
		
		$chr->getProject->print_dot(1);
		my $debug;
		my $vsmall = $gene->getCompactVectorOriginCategories(\@all_cat,$debug);
		
		my $vchr = $gene->return_compact_vector( $chr->getVariantsVector());
		$vsmall &= $vchr;
		
		if ($vsmall->is_empty()){
			delete $chr->project->{genes_object}->{$gene->id};
			$self->deleteGene($chr,$gene);
			next;
		}
		$gene->setCurrentVector($vsmall) if not $not_set_vector_gene;
		$chr->{genes_object}->{$gene->id} ++;
		$nb ++;
		$id_genes->{$gene->id} ++;
		$variants_genes += $gene->enlarge_compact_vector($vsmall);
	}
	return $variants_genes;
}

sub return_cat {
	my ($project,@unselect) = @_;
	my %hcat = %{$project->buffer->config->{'genbo_annotations_names'}};
	foreach my $c (@unselect){
		delete $hcat{$c};
		if ($c eq 'predicted_splice_site') {
			delete $hcat{'predicted_promoter_ai'} if exists $hcat{'predicted_promoter_ai'};
		}
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
	my $vvv = $chr->getVariantsVector();
	#$chr->rocks_vector("r")->get_vector_chromosome("coding") & $chr->getVariantsVector();
	
	my $cc = $chr->getVariantsVector()->Clone;
	my $t =time;
	#$rocks->prepare([keys %$lh]);
	$chr->{genes_object} = {};
	my @genes ;
	
	$self->setGenes($chr);
	 
	#here
	my $variants_genes = $self->get_vector_filter_gene_annotations($chr, $hFiltersChr);
	
	$chr->getVariantsVector->Intersection($chr->getVariantsVector(), $variants_genes);
	$chr->{buffer_vector} = $chr->getVariantsVector();
	
	
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
	my ($self, $chr, $model, $h_args_compound_second_variant) = @_;
	return if $chr->getVariantsVector->is_empty();
	return if not $model;
	my $vector_models = $chr->getNewVector();
	if ($model and $chr->getProject->typeFilters() eq 'individual') {
		$vector_models += $self->filter_genetics_models_individual($chr, $model);
	}
	elsif ($model and $chr->getProject->typeFilters() eq 'familial') {
		$vector_models += $self->filter_genetics_models_familial($chr, $model, $h_args_compound_second_variant);
	}
	elsif ($model and $chr->getProject->typeFilters() eq 'somatic') {
		$vector_models += $self->filter_genetics_models_somatic($chr, $model);
	}
	$chr->setVariantsVector($vector_models);
	
	my $vc = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		next if not $gene->compact_vector();
			my $vv = $chr->getVariantsVector() & $gene->getCurrentVector();
			$gene->setCurrentVector($vv);
			$vc += $vv;
	}
	$chr->setVariantsVector($vc);
	
	if ($self->verbose_debug) { warn "\nCHR ".$chr->id()." -> AFTER models - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
}

sub filter_genetics_models_individual {
	my ($self, $chr, $model, $h_args) = @_;
	my $vector_models = $chr->getNewVector();
	if ($model eq 'recessif') {
		foreach my $pat (@{$chr->getPatients()}) {
			next if $pat->in_the_attic();
			my $v_pat = $chr->getNewVector();
			my $v_all = $chr->getVariantsVector() & $pat->getVariantsVector($chr);
			my $v_ho = $chr->getVariantsVector() & $pat->getHo($chr);
			my $v_he = $chr->getVariantsVector() & $pat->getHe($chr);
			foreach my $gene (@{$chr->getGenes()}) {
				my $v_gene = $gene->getCurrentVector() & $v_all;
				my $v_gene_ho = $v_gene & $v_ho;
				my $v_gene_he = $v_gene & $v_he;
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
				my $v_gene = $gene->getCurrentVector() & $v_all;
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
	my ($self, $chr, $model, $h_args_compound_second_variant) = @_;
	my $vector_models = $chr->getNewVector();
				
	if ($model eq 'compound' or $model eq 'recessif_compound' or $model eq 'uniparental_recessif_compound') {
		my $no = $chr->flush_rocks_vector();
		my $gnomad_ac_2 = $h_args_compound_second_variant->{gnomad_ac_2};
		my $dejavu_2 = $h_args_compound_second_variant->{dejavu_2};
		
		my $not_update_gene_vector = 1;
		my $variants_genes_2 = $self->get_vector_filter_gene_annotations($chr, $h_args_compound_second_variant, $not_update_gene_vector);
		$variants_genes_2 &= $no->get_vector_chromosome("sdv_".$dejavu_2) if not $dejavu_2 eq 'all';
		$variants_genes_2 &= $no->get_vector_chromosome($gnomad_ac_2) if $gnomad_ac_2;
		$no->close();
		
		my $h_fam_vectors;
		foreach my $fam (@{$chr->getProject->getFamilies()}) {
			next if $fam->in_the_attic();
			next if (not $fam->getMother() or not $fam->getFather());
			my $v_fam = $fam->getVariantsVector($chr);
			my $v_fam_2 = $v_fam & $variants_genes_2;
			
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
			
			$h_fam_vectors->{$chr->id()}->{$fam->name()}->{mother_1} = $fam->getVectorMotherTransmission($chr);
			$h_fam_vectors->{$chr->id()}->{$fam->name()}->{mother_1} &= $fam->getMother->getHe($chr);
			$h_fam_vectors->{$chr->id()}->{$fam->name()}->{mother_1} &= $v_children_ill;
			$h_fam_vectors->{$chr->id()}->{$fam->name()}->{father_1} = $fam->getVectorFatherTransmission($chr);
			$h_fam_vectors->{$chr->id()}->{$fam->name()}->{father_1} &= $fam->getFather->getHe($chr);
			$h_fam_vectors->{$chr->id()}->{$fam->name()}->{father_1} &= $v_children_ill;
			$h_fam_vectors->{$chr->id()}->{$fam->name()}->{mother_2} = $h_fam_vectors->{$chr->id()}->{$fam->name()}->{mother_1} & $variants_genes_2;
			$h_fam_vectors->{$chr->id()}->{$fam->name()}->{father_2} = $h_fam_vectors->{$chr->id()}->{$fam->name()}->{father_1} & $variants_genes_2;
			$h_fam_vectors->{$chr->id()}->{$fam->name()}->{compound} = $chr->getNewVector();
		}
			
		foreach my $gene (@{$chr->getGenes()}) {
			my $v_gene_after = $chr->getNewVector();
			if ($chr->{genes_object}->{$gene->id()} == 1) {
				delete $chr->project->{genes_object}->{$gene->id};
				$self->deleteGene($chr,$gene);
				next;
			}
			my $v_gene_1 = $gene->getVectorOrigin() & $chr->getVariantsVector();
			my $v_gene_2 = $gene->getVectorOrigin() & $variants_genes_2;
			
			foreach my $fam (@{$chr->getFamilies()}) {
				next if not exists $h_fam_vectors->{$chr->id()}->{$fam->name()};
				
				my $v_gene_m = $v_gene_1 & $h_fam_vectors->{$chr->id()}->{$fam->name()}->{mother_1};
				my $v_gene_f = $v_gene_1 & $h_fam_vectors->{$chr->id()}->{$fam->name()}->{father_1};
				my $v_gene_m_2 = $v_gene_2 & $h_fam_vectors->{$chr->id()}->{$fam->name()}->{mother_2};
				my $v_gene_f_2 = $v_gene_2 & $h_fam_vectors->{$chr->id()}->{$fam->name()}->{father_2};
					
				if ($v_gene_m->Norm() >= 1 and $v_gene_f_2->Norm() >= 1) {
					$h_fam_vectors->{$chr->id()}->{$fam->name()}->{compound} += $v_gene_m;
					$h_fam_vectors->{$chr->id()}->{$fam->name()}->{compound} += $v_gene_f_2;
				}
				if ($v_gene_m_2->Norm() >= 1 and $v_gene_f->Norm() >= 1) {
					$h_fam_vectors->{$chr->id()}->{$fam->name()}->{compound} += $v_gene_m_2;
					$h_fam_vectors->{$chr->id()}->{$fam->name()}->{compound} += $v_gene_f;
				}
			}
		}
			
		foreach my $fam (@{$chr->getFamilies()}) {
			if (exists $h_fam_vectors->{$chr->id()}->{$fam->name()} and exists $h_fam_vectors->{$chr->id()}->{$fam->name()}->{compound}) {
				$vector_models += $h_fam_vectors->{$chr->id()}->{$fam->name()}->{compound};
			}
			$fam->setCurrentVariantsVector($chr, $vector_models);
		}
	}
	if ($model eq 'recessif' or $model eq 'recessif_compound' or $model eq 'uniparental_recessif_compound') {
		foreach my $fam (@{$chr->getProject->getFamilies()}) {
			next if $fam->in_the_attic();
			my $v_fam = $fam->getVariantsVector($chr);
			foreach my $pat (@{$fam->getPatients()}) {
				if ($pat->isHealthy()) { $v_fam -= $pat->getVectorOriginHo($chr); }
				else { $v_fam -= $pat->getVectorOriginHe($chr); }
			}
			my $v_fam_ok = $chr->getNewVector();
			foreach my $child (@{$fam->getPatientsIll()}) {
				my $v_child = $v_fam & $fam->getVector_individual_recessive($chr, $child);
				$v_fam_ok += $v_child;
				$v_fam_ok -= $child->getVectorOriginHe($chr);
			}
			$v_fam_ok &= $chr->getVariantsVector();
			if ($model eq 'recessif') { $fam->setCurrentVariantsVector($chr, $v_fam_ok); }
			else { $fam->addCurrentVariantsVector($chr, $v_fam_ok); }
			$vector_models += $v_fam_ok;
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
	my ($self, $chr, $model) = @_;
	my $vector_models = $chr->getNewVector();
	if ($model eq 'loh') {
		foreach my $group (@{$chr->getProject->getSomaticGroups}) {
			foreach my $pat (@{$group->getPatients()}) {
				$vector_models += $group->getVector_individual_loh($chr, $pat);
			}
		}
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
