package hvariant;
use strict;
use Data::Dumper;

sub hash_variant_2 {
	my ($patient,$vid) = @_;
	my $project = $patient->getProject();
	my $v = $project->_newVariant($vid);
	my $global_genes;
	my $chr = $v->getChromosome();
	my $fam = $patient->getFamily();
	my $buffer = $v->buffer();
	my $h =  update_variant_editor::construct_hash_variant ( $project, $v,undef,$patient,1);
	$h->{vector_id} = $v->vector_id;
	$h->{global_vector_id} = $v->global_vector_id;
	update_variant_editor::vvarsome($h,$patient);
	foreach my $g (@{$v->getGenes}){
		my $score_variant = $v->scaledScoreVariant($g,$patient);
		$h->{scaled_score}->{$g->id} = $score_variant;
	}
	return $h;
}

sub hash_variant {
	my ($patient,$vid) = @_;
	my $project = $patient->getProject();
	my $v = $project->returnVariants($vid);
	my $global_genes;
	my $chr = $v->getChromosome();
	my $fam = $patient->getFamily();
	my $buffer = $v->buffer();
#	my $h =  update_variant_editor::construct_hash_variant ( $project, $v,undef,$patient,1);
#	$h->{vector_id} = $v->vector_id;
#	$h->{global_vector_id} = $v->global_vector_id;
#	update_variant_editor::vvarsome($h,$patient);
		
	my $mother = $patient->getFamily->getMother();
	my $father = $patient->getFamily->getFather();
	my $agenes =[];
	my $score_genes;
	my $javascript_id = time;
	foreach my $g (@{$v->getGenes}){
		my $hgenes;
		$javascript_id ++;
		my $score_variant = $v->scaledScoreVariant($g,$patient);
		my $parent_score = 0;
		$parent_score = $score_variant  if $score_variant >0 &&  $v->isConsequence($g,"medium");
		my $nid = "s:".$v->{vector_id}.":".$g->id."@".$patient->name;
#		$h->{scaled_score}->{$g->id} = $score_variant;
		unless (exists $global_genes->{$patient->name}->{$g->id}){
					$global_genes->{$patient->name}->{$g->id}->{id} = $g->id;
					$global_genes->{$patient->name}->{$g->id}->{name} = $g->external_name;
					$global_genes->{$patient->name}->{$g->id}->{score} = $g->score;
					$global_genes->{$patient->name}->{$g->id}->{omim_inheritance} = $g->omim_inheritance;
					$global_genes->{$patient->name}->{$g->id}->{omim_inheritance} = $g->omim_inheritance;
					$global_genes->{$patient->name}->{$g->id}->{external_name} = $g->external_name;
					$global_genes->{$patient->name}->{$g->id}->{pLI} = $g->pLI;
					$global_genes->{$patient->name}->{$g->id}->{omim_id} = $g->omim_id;
					$global_genes->{$patient->name}->{$g->id}->{panels} = $buffer->queryPanel()->getPanelsForGeneName($g->external_name);
					$global_genes->{$patient->name}->{$g->id}->{js_id} = $javascript_id."_".$g->id;
					$global_genes->{$patient->name}->{$g->id}->{father}->{score} = -100;
					$global_genes->{$patient->name}->{$g->id}->{mother}->{score} = -100;
					$global_genes->{$patient->name}->{$g->id}->{description} = $g->description;
					$global_genes->{$patient->name}->{$g->id}->{phenotypes} = $g->phenotypes;
					$global_genes->{$patient->name}->{$g->id}->{vector} = $g->getCurrentVector();
					if ($patient->isMale && $chr->name eq "X" && $chr->isAutosomal($g->start,$g->end) ){
							$global_genes->{$patient->name}->{$g->id}->{no_compound} ++;
					}
					unless($mother){
						$global_genes->{$patient->name}->{$g->id}->{no_compound} ++;
					}
					unless($father){
							$global_genes->{$patient->name}->{$g->id}->{no_compound} ++;
					}
				
				}
				unless (exists $hgenes->{$g->id}){
					$hgenes->{score_mother} = -100;
					$hgenes->{score_father} = -100;
					$hgenes->{score_biallelic} = -10;
					$hgenes->{id} = $g->id;
			}
			$hgenes->{score_variant}->{$v->id} = $score_variant;
			die() if exists $hgenes->{variant}->{$v->id};
			
			$hgenes->{variant}->{$v->id}->{score} = $score_variant;
			$hgenes->{pathogenic} ++  if  $v->isDM or $v->isClinvarPathogenic;
			$hgenes->{clinvar_hgmd} ++ if  $v->hgmd or $v->clinvar;
			
			$hgenes->{cnv_del} ++ if  $v->isLargeDeletion();
			$hgenes->{cnv_dup} ++ if  $v->isLargeDuplication();

			$hgenes->{denovo_rare} ++ if $v->getGnomadAC< 10 && $v->isDenovoTransmission($fam,$patient) ;#&& $v->effectImpact($g) eq "moderate";
			$global_genes->{$patient->name}->{$g->id}->{denovo_rare} ++ if $v->getGnomadAC< 10 && $v->isDenovoTransmission($fam,$patient);
			$global_genes->{$patient->name}->{$g->id}->{cnv_del} ++ if  $v->isLargeDeletion();;
			$global_genes->{$patient->name}->{$g->id}->{cnv_dup} ++ if  $v->isLargeDuplication();
			
			
			my $val = $v->score_validations($g);
			#my $val = 0;
			if ($val){
				$hgenes->{validations} = $val->{validation} ;#if $val->{validation} > $hgenes->{$g->id}->{validations};
			}
			 
			push (@{$hgenes->{variants}}, $v->id);
			my $debug;
				
			if ($patient->isChild && $patient->getFamily->isTrio()){
				
				if  ($v->isUniparentalDisomyTransmission($fam,$patient)){
					$hgenes->{score_biallelic} = $score_variant;
					$hgenes->{biallelic} ++;
					$hgenes->{variant}->{$v->id}->{biallelic} ++;
				}
				
				if  ($v->isMosaicTransmission($fam,$patient)){
					$hgenes->{score_biallelic} = $score_variant;
					$hgenes->{biallelic} ++;
					$hgenes->{variant}->{$v->id}->{biallelic} ++;
				}
				elsif  ($v->isDominantTransmission($fam,$patient)){
					$hgenes->{score_biallelic} = $score_variant;
					$hgenes->{biallelic} ++;
					$hgenes->{variant}->{$v->id}->{biallelic} ++;
				}
				 elsif  ($v->isMotherTransmission($fam,$patient)){
				 	$hgenes->{score_mother} = $score_variant;
				 	$hgenes->{mother} ++;
				 	$hgenes->{variant}->{$v->id}->{mother} ++;
				 	if ($global_genes->{$patient->name}->{$g->id}->{mother}->{score} < $score_variant  ){
				 		$global_genes->{$patient->name}->{$g->id}->{mother}->{id} = $v->id;
				 		$global_genes->{$patient->name}->{$g->id}->{mother}->{score} = $score_variant;
				 		$global_genes->{$patient->name}->{$g->id}->{mother}->{vector_id} .= $v->{vector_id};
				 		$global_genes->{$patient->name}->{$g->id}->{mother}->{global_vector_id} .= $v->{global_vector_id};
				 	}
				 	
				 }
				 elsif  ($v->isFatherTransmission($fam,$patient)){
				 	$hgenes->{score_father} = $score_variant;#>=$hgenes->{score_father}?$score_variant:$hgenes->{score_father});
				 	$hgenes->{father} ++;
				 	$hgenes->{variant}->{$v->id}->{father} ++;
				 	if ($global_genes->{$patient->name}->{$g->id}->{father}->{score} < $score_variant){
				 		$global_genes->{$patient->name}->{$g->id}->{father}->{id} = $v->id;
				 		$global_genes->{$patient->name}->{$g->id}->{father}->{score} = $score_variant;
				 		$global_genes->{$patient->name}->{$g->id}->{father}->{vector_id} .= $v->{vector_id};
				 		$global_genes->{$patient->name}->{$g->id}->{father}->{global_vector_id} .= $v->{global_vector_id};
				 	}
				 }
				 
				 else {
				 	warn "else" if $debug;
				 	$hgenes->{variant}->{$v->id}->{biallelic} ++;
				 	$hgenes->{score_biallelic} = $score_variant;
				 }
				
			}
			else {
				$hgenes->{biallelic} ++;	
				$hgenes->{score_biallelic} = $score_variant;
				$hgenes->{variant}->{$v->id}->{biallelic} ++;
			}
			my $mask = $v->annotation()->{$g->id}->{mask};
			$hgenes->{mask} = $mask;
			push(@$agenes,$hgenes);
		}#end Genes 
	#	warn Dumper $h;
	#	die("--------");
		return {array=>$agenes};
	#	 	$no->put($v->{global_vector_id}."-hashvariants@".$patient->name,$h);
	#	 	$no->put($v->{global_vector_id}."-mask@".$patient->name,{value=>$vmask,patient=>$patient->name,id1=>$v->vector_id,type=>"mask"}) if $vmask;
	#	 	$no->put($v->{global_vector_id}."-genes@".$patient->name,{array=>$agenes,patient=>$patient->name,id1=>$v->{global_vector_id},type=>"genes"}) if $agenes;
	#	 	$no->put($v->{global_vector_id}."@".$patient->name,{id=>$v->id,id1=>$v->{global_vector_id},patient=>$patient->name,type=>"variants",});
	#	 	push(@{$zh->{$patient->name}},$v->{global_vector_id});
		

	
	
	
}





1;