#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;

package get_variations;
use strict;
use Data::Dumper;
#use KyotoCabinet;
use Storable qw/freeze thaw nfreeze nstore_fd nstore retrieve/;
use JSON;
use GenBoNoSql;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/utility.pm";

require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";


sub getIds_byCache_onlive_polyviewer {
	my ($buffer, $project, $chr, $hIds, $user) = @_;
	
	my $can_use_hgmd = $buffer->hasHgmdAccess($user);
	my @headers_validations = ("varsome","igv","alamut","var_name","trio","gnomad","deja_vu","table_validation","table_transcript");
	my @header_transcripts = ("consequence","enst","nm","ccds","appris","exon","nomenclature","codons","codons_AA", "polyphen","sift","ncboost","cadd","revel",);
	my $fsize = "font-size:10px";
	my $class;
	$class->{rowspan} -= 1;
	$class->{rowspan} = 1 if $class->{rowspan} <=0;
	$class->{style} = "min-width:10%;padding:1px";
	my $value_red = 10;
	my $value_coral = 9;
	my $value_orange = 8;
	my $value_yellow = 5;
	my $color_red = '#CE0000';
	my $color_coral = 'coral';
	my $color_orange = '#EFB73E';
	my $color_yellow = 'yellow';
	
	my ($h_scaled_score, $h_scaled_score_gene);
	my $hash_res;
	my $proj_name = $project->name();
	
	my $cgi = new CGI();
	my ($h_html_genes, $hgenes, $h_var, $h_var_trio);
	
	foreach my $var_id (keys %$hIds) {
		my ($h_local_scaled_score_done);
		my $v;
		eval { $v = $project->getVariant($var_id); };
		if ($@) { last; }
		foreach my $g (@{$v->getGenes()}) {
			next if ($g->getVariantsVector->is_empty());
			$h_html_genes->{$g->id()}->{obj} = $g;
			push(@{$h_html_genes->{$g->id()}->{var_ids}}, $var_id);
			my $h_families_done;
			foreach my $family (@{$project->getFamilies()}) {
				$v = $project->getVariant($var_id) unless ($v);
				next unless ($family->getVariantsVector($v->getChromosome())->contains($v->vector_id()));
				my @lTrios;
				foreach my $patient (@{$family->getChildrenIll()}) {
					next unless ($patient);
					#delete $v->{scale_score} if (exists $v->{scale_score});
					my $hvariation;
					eval { $hvariation = update_variant_editor::construct_hash_variant( $project, $v, undef, $patient); };
					if ($@) { next; }
					# NO CSS ADDED from others scripts, like view.pl (HGMD).
					$hvariation->{html}->{no_css_polydiag} = 1;
					$hvariation->{obj} = $v;
					my $max_score = -99;
					my $max_tr;
					foreach my $tr (@{$g->getTranscripts()}) {
						my $this_score = $v->scaledScoreVariant($tr, $patient, undef);
						if ($this_score > $max_score) {
							$max_score = $this_score;
							$max_tr = $tr;
						}
					}
					my $this_score = $max_score;
					my $tr = $max_tr;
					$hvariation->{scaled_score}->{$g->id()} = $this_score;
					$hvariation->{scaled_score_gene_used} = $g->external_name();
					$hvariation->{scaled_score_transcript_used} = $tr->id();
					unless (exists $hgenes->{$g->id}){
						$hgenes->{$g->id}->{project_name} = $proj_name;
						$hgenes->{$g->id}->{name} = $g->external_name;
						$hgenes->{$g->id}->{description} = $g->description;
						$hgenes->{$g->id}->{phenotypes} = $g->phenotypes;
						$hgenes->{$g->id}->{score_mother} = 0;
						$hgenes->{$g->id}->{score_father} = 0;
						$hgenes->{$g->id}->{score_biallelic} = -10;
						$hgenes->{$g->id}->{score} = $g->score;
						$hgenes->{$g->id}->{id} = $g->id;
						$hgenes->{$g->id}->{omim_inheritance} = $g->omim_inheritance;
						$hgenes->{$g->id}->{external_name} = $g->external_name;
						$hgenes->{$g->id}->{pLI} = $g->pLI;
						$hgenes->{$g->id}->{omim_id} = $g->omim_id;
						$hgenes->{$g->id}->{panels} = $buffer->queryPanel()->getPanelsForGeneName($g->external_name);
						$hgenes->{$g->id}->{js_id} = $proj_name."_".$g->id;
					}
					$hgenes->{$g->id}->{score_variant}->{$v->id} = $this_score;
					$hgenes->{$g->id}->{pathogenic} ++  if  $v->isDM or $v->isClinvarPathogenic;
					$hgenes->{$g->id}->{clinvar_hgmd} ++ if  $v->hgmd or $v->clinvar;
					#$hgenes->{$g->id}->{denovo_rare} ++ if $v->getGnomadAC< 10 &&  $v->isDenovoTransmission($patient->getFamily,$patient);
					my $val = $v->score_validations($g);
					if ($val){
						$hgenes->{$g->id}->{validations} = $val->{validation};
					}
					push (@{$hgenes->{$g->id}->{variants}}, $v->id);
					if ($patient->isChild && $patient->getFamily->isTrio()){
						 if  ($v->isMotherTransmission($family, $patient)){
						 	$hgenes->{$g->id}->{score_mother} = $this_score;
						 }
						 elsif  ($v->isFatherTransmission($family, $patient)){
						 	$hgenes->{$g->id}->{score_father} = $this_score;
						 }
						 else {
						 	$hgenes->{$g->id}->{score_biallelic} = $this_score;
						 }
					}
					else {
						$hgenes->{$g->id}->{score_biallelic} = $this_score;
					}
					my $this_alert = 0;
					if ($hgenes->{$g->id}->{validations} > 2) { $this_alert = 3; }
					if ($hgenes->{$g->id}->{validations} > 4) { $this_alert = 4; }
					update_variant_editor::vhgmd($v,$hvariation);
					update_variant_editor::vclinvar($v,$hvariation);
					if ($can_use_hgmd and $v->hgmd() and $v->hgmd_phenotype()) {
						$hvariation->{hgmd_phenotype} = $v->hgmd_phenotype();
					}
					$hvariation->{html}->{hgmd_no_access} = qq{<span class="glyphicon glyphicon-ban-circle" aria-hidden="true" style='font-size:12px;color:black;'></span>};
					update_variant_editor::table_validation($patient,$hvariation,$g);
					$h_local_scaled_score_done->{$hvariation->{scaled_score}} = undef;
					unless (exists $h_families_done->{$family->name()}) {
						update_variant_editor::trio($v, $hvariation, $patient);
						push(@{$h_var_trio->{$var_id}->{$g->id()}}, "<b><span style='font-size:10px;'>Fam ".$patient->getFamily->name().'</b></span><br>'.$hvariation->{html}->{trio});
					}
					$hvariation->{html}->{$g->id()}->{'table_transcript'} = update_variant_editor::table_transcripts($hvariation->{genes}->{$g->id()}, \@header_transcripts);
					$h_families_done->{$family->name()} = undef;
					delete $hvariation->{obj};
					$h_scaled_score->{$g->id()}->{$this_score}->{$var_id} = $hvariation;
					$v = undef;
				}
			}
		}
		$v = undef;
	}
	
	
	my $h_genes_max_score;
	foreach my $gene_id (keys %$h_html_genes) {
		my $max_score = -999;
		foreach my $var_id (keys %{$hgenes->{$gene_id}->{score_variant}}) {
			my $this_score = $hgenes->{$gene_id}->{score_variant}->{$var_id};
			$max_score = $this_score if ($this_score > $max_score);
		}
		push(@{$h_genes_max_score->{$max_score}}, $gene_id);
	}
	my @l_max_score = reverse sort {$a <=> $b} keys %$h_genes_max_score;
	
	my $out2 = '';
	foreach my $max_score (@l_max_score) {
		foreach my $gene_id (@{$h_genes_max_score->{$max_score}}) {
			my $g = $h_html_genes->{$gene_id}->{obj};
			#next unless ($hgenes->{$gene_id}->{variants});
			# PANEL Gene de polydiag avec modif local de css
			
			my $collapse_id = 'tr_'.$gene_id;
			$out2 .= "<tr class='collapse' id='$collapse_id'>";
			my $bgcolor = '#60798B';
			$out2 .=  $cgi->start_div({class=>"panel-heading panel-face panel-grey",style=>"background-color:$bgcolor;min-height:13px;border:0px"});
			my $panel_id = "panel_".$hgenes->{$gene_id}->{js_id};
			
			$hgenes->{$gene_id}->{max_score} = $max_score;
			my $this_html = update_variant_editor::panel_gene($hgenes->{$gene_id}, $panel_id);
			$this_html =~ s/text\-shadow\:1px 1px 2px black/color:white/g;
			$out2 .= $this_html;
			$out2 .= $cgi->end_div();
			$out2 .= "</tr>";
			
			$out2 .=  $cgi->start_div({style=>"background-color:$bgcolor;"});
			$out2 .= "<table class='table' style='font-size:13px;'>";
			$out2 .= "<thead>";
			$out2 .= $cgi->start_Tr({style=>"background-color:#E9DEFF;$fsize;border: solid 1px black;"});
			foreach my $h (@headers_validations){
				if ($h eq 'scaled_score') { $out2 .= $cgi->th('Old_diag_score'); }
				#elsif ($h eq 'trio') { $out2 .= $cgi->th('Families / Patients'); }
				else { $out2 .=  $cgi->th(ucfirst($h)); }
			}
			$out2 .= $cgi->end_Tr();
			$out2 .= "</thead>";
			$out2 .= "<tbody>";
			
			my $h_var_done;
			foreach my $this_score (reverse sort {$a <=> $b} keys %{$h_scaled_score->{$g->id()}}) {
				my @lVarIds = keys %{$h_scaled_score->{$g->id()}->{$this_score}};
				foreach my $var_id (@lVarIds) {
					next if (exists $h_var_done->{$var_id});
					my $hvariation = $h_scaled_score->{$g->id()}->{$this_score}->{$var_id};
					$hvariation->{html}->{trio} = qq{<div style='max-height:180px;overflow-x:hidden;overflow-y:auto;'>};
					$hvariation->{html}->{trio} .= join('<br>', @{$h_var_trio->{$var_id}->{$g->id()}});
					$hvariation->{html}->{trio} .= '</div>';
					unless ($can_use_hgmd) {
						$hvariation->{html}->{hgmd} = qq{<span class="glyphicon glyphicon-ban-circle" aria-hidden="true" style='font-size:12px;color:black;'></span>};
					}
					my $class_tr;
					$class_tr->{style} = "background:#D3D6D3;border: solid 1px grey;";
					$out2 .= $cgi->start_Tr($class_tr);
					foreach my $h (@headers_validations){
						if ($h eq "trio" or "table_transcript"){
							$class->{style} = "min-width:10%;max-height:200px;vertical-align:middle;padding:5px;border-top: solid 1px grey;";
						}
						elsif ($h eq "igv" or "alamut"){
							$class->{style} = "max-width:50px;vertical-align:middle;padding:5px;border-top: solid 1px grey;";
						}
						else {
							$class->{style} = "min-width:5%;vertical-align:middle;padding:5px;border-top: solid 1px grey;";
						}
						if ($h eq 'scaled_score') {
							my $color = 'black';
							if ($hvariation->{$h}->{$gene_id} >= $value_red) { $color = $color_red; }
							elsif ($hvariation->{$h}->{$gene_id} >= $value_coral) { $color = $color_coral; }
							elsif ($hvariation->{$h}->{$gene_id} >= $value_orange) { $color = $color_orange; }
							elsif ($hvariation->{$h}->{$gene_id} >= $value_yellow) { $color = $color_yellow; }
							my $score = sprintf("%.1f", $hvariation->{$h}->{$gene_id});
							my $b = "<button type='button' class='btn ' style='background-color:white;font-size:10px;'>";
							$b .= "<span style='color:$color;'>".$score."<span>";
							$b .= "</button>";
							$out2.= $cgi->td($class,$b);
						}
						elsif ($h eq 'table_validation') {
							if ($can_use_hgmd) { $out2.= $cgi->td($class,$hvariation->{html}->{$h}); }
							else { $out2.= $cgi->td($class,$hvariation->{html}->{table_validation_hgmd_no_access}); }
							
						}
						elsif ($h eq 'table_transcript') {
							$out2 .= $cgi->td($class, $hvariation->{html}->{$gene_id}->{'table_transcript'});
						}
						else {
							$out2.= $cgi->td($class,$hvariation->{html}->{$h});
						}
					}
					$h_var_done->{$var_id} = undef;
				}
				$out2 .= $cgi->end_Tr();
				$out2 .= $cgi->end_Tr();
			}
			
			
			$out2 .= "</tbody>";
			$out2 .= "</table>";
			$out2 .= "</div>";
			$out2 .= "<br>";
		}
	}
	$hash_res->{html} = $out2;
	
	#die;
	
	my @lRes;
	push(@lRes, $hash_res);
	return \@lRes;
}

sub getIds_byCache_onlive {
	my ($buffer, $project, $chr, $hIds, $user) = @_;
	my @lRes;
	my $hAnnotImpacts;
	foreach my $impact (keys %{$project->impacts_ensembl_annotations()}) {
		foreach my $annot (keys %{$project->impacts_ensembl_annotations->{$impact}}) {
			$hAnnotImpacts->{$annot} = $impact;
		}
	}
	foreach my $id (keys %$hIds) {
		my $v_id = $hIds->{$id}->{'v_id'};
		my $key_id = $id;
		my @lTmp = split('_', $id);
		unless ($lTmp[0] =~ /[0-9XYMT]+/) {
			shift(@lTmp);
			$id = join('_', @lTmp);
		}
		my $var = $project->getVariant($id);
		my @lGenesObj = @{$var->getGenes()};
		my $hash;
		$hash->{'name'} = $var->name();
#		$hash->{'name'} = $var->rs_name() if ($var->rs_name());
		$hash->{'id'} = $id;
		$hash->{'v_id'} = $v_id;
		$hash->{'position'} = $chr->name().':'.$var->start().'-'.$var->end();
		$hash->{'chromosome'} = $chr->name();
		my $pos = $var->position($chr);
		$hash->{'start'} = $pos->start();
		$hash->{'end'} = $pos->end();
		$hash->{'text'} = $var->alleles();
		$hash->{'alamut'} = $var->alamut_id();
		$hash->{'ho_region'} = $hIds->{$id}->{ho_region};
		my @lPatNames;
		foreach my $patient (@{$var->getPatients()}) { push(@lPatNames, $patient->name()); }
		$hash->{'patients'} = join(';', @lPatNames);
#		if ($var->isVariation()) {
#			$hash->{'structural_type'} = 'snp';
#			$hash->{'type'} = 'variations';
#		}
#		elsif ($var->isInsertion()) {
#			$hash->{'structural_type'} = 'ins';
#			$hash->{'type'} = 'insertions';
#		}
#		elsif ($var->isDeletion()) {
#			$hash->{'structural_type'} = 'del';
#			$hash->{'type'} = 'deletions';
#		}
		
		#$hash->{'structural_type'} = $var->structural_type();
		$hash->{'type'} = $var->type();
		
		$hash->{'gnomad_ac'} = '-';
		$hash->{'gnomad_ho'} = '-';
		$hash->{'gnomad_an'} = '-';
		if ($var->getGnomadAC()) {
			my $vn = $var->id();
			$vn =~ s/_/-/g;
			$vn =~ s/chr//;
			my $link = "https://gnomad.broadinstitute.org/variant/$vn";
			$hash->{'gnomad_ac'} = $var->getGnomadAC().';'.$link;
			$hash->{'gnomad_ho'} = $var->getGnomadHO().';'.$link;
			$hash->{'gnomad_an'} = $var->getGnomadAN().';'.$link;
		}
		
		$hash->{'isClinical'} = 'No';
		$hash->{'isClinical'} = 'Yes' if ($var->isClinical());
		$hash->{'polyphen_status'} = $var->polyphenStatus();
		$hash->{'sift_status'} = $var->siftStatus();
		$hash->{'cadd_score'} = $var->cadd_score();
		$hash->{'cadd_score'} = '.' if $hash->{'cadd_score'} == -1;
		$hash->{'ncboost'} = $var->ncboost_category().';'.$var->ncboost_score();
		my @l_dbscsnv;
		push(@l_dbscsnv, 'rf:'.$var->dbscsnv->{'rf'}) if exists $var->dbscsnv->{'rf'};
		push(@l_dbscsnv, 'ada:'.$var->dbscsnv->{'ada'}) if exists $var->dbscsnv->{'ada'};
		$hash->{'dbscsnv'} = join(' ', @l_dbscsnv);
		
		my $release = 'hg38';
		$release = 'hg19' if ($var->getProject->annotation_genome_version() =~ /HG19/);
		my $varsome_url = qq{https://varsome.com/variant/$release/}.$var->gnomad_id();
		$hash->{'varsome'} = $varsome_url;
		
		$hash->{'promoterAI'} = '-';
		$hash->{'spliceAI'} = '-';
		
		my (@lGenesIds, @GenesNames, @lGenesCons, @lGenesOmim, @lGenesPolywebScore);
		foreach my $gene (@lGenesObj) {
			if ($gene->omim_id()) { push(@lGenesOmim, $gene->omim_id()); }
			else { push(@lGenesOmim, undef); }
			my ($gene_id, $id_chr) = split('_', $gene->id());
			push (@lGenesIds, $gene_id);
			push (@GenesNames, $gene->external_name()); 
			push (@lGenesCons, $var->variationTypeInterface($gene));

			my @l_score_spliceAI;
			my $h_score_spliceAI = $var->spliceAI_score($gene);
			if ($h_score_spliceAI) {
				foreach my $cat (sort keys %$h_score_spliceAI) {
					push(@l_score_spliceAI, $cat.'='.$h_score_spliceAI->{$cat});
				}
			}
			else {
				push(@l_score_spliceAI, '-');
			}
			$hash->{'spliceAI'} = join(' ', @l_score_spliceAI);
			if ($var->promoterAI) {
				foreach my $tr (@{$gene->getTranscripts()}) {
					my $score_promoterAI = $var->promoterAI_score($tr);
					$hash->{'promoterAI'} = $score_promoterAI if $score_promoterAI;
					last if $hash->{'promoterAI'} ne '-';
				}
			}
		}
		unless (@lGenesIds) {
			push (@lGenesCons, 'Intergenic');
			push (@GenesNames, 'intergenic_'.$hash->{'chromosome'}.'_'.$var->start().'_'.$var->end());
			push (@lGenesIds, 'intergenic_'.$hash->{'chromosome'}.'_'.$var->start().'_'.$var->end());
			push (@lGenesPolywebScore, '.');
		}
		$hash->{genes_name} = join(";", @lGenesIds);
		$hash->{external_names} = join(";", @GenesNames);
		$hash->{'gene_in_omim'} = join(';', @lGenesOmim);
		#TODO: PolyScore
#		$hash->{'polyweb_score'} = join(';', @lGenesPolywebScore); 
		my $var_cons = join(";", @lGenesCons);
#		my @lTmp2;
#		foreach my $this_cons (split(',', $var_cons)) {
#			my ($name, $name2) = split(';', $buffer->config->{ensembl_annotations}->{$this_cons});
#			#$name .= ';'.$hAnnotImpacts->{$name} if (exists $hAnnotImpacts->{$name});
#			push(@lTmp2, $name)
#		}
#		$var_cons = join(', ', @lTmp2);
		if ($project->hasHgmdAccess($user)) {
			$hash->{'hgmd_class'} = '';
			if ($var->hgmd()) {
				my $h;
				$h->{obj} = $var;
				my ($class, $cmd) = update::hgmd($project,$h,'return_cmd');
				$hash->{'hgmd_class'} = $class.';'.$cmd;
				if ($var->isNewHgmd()) {
					$hash->{'hgmd_class'} .= ";New!";
				}
			}
		}
		$hash->{'consequence!all'} = $var_cons;
		$hash->{'homo_hetero'} += 1 if ($hIds->{$key_id}->{ho});
		$hash->{'homo_hetero'} += 2 if ($hIds->{$key_id}->{he});
		$hash->{'nb_patients'} = $hIds->{$key_id}->{nb_patients};
		my $freq = $var->frequency();
		if ($freq == -1) { $hash->{'freq'} = undef; }
		else { $hash->{'freq'} = sprintf("%.3f", ($freq*100)).' %'; }
		my ($is_genome, $type, $db_gnomad);
		$is_genome = 1 if ($project->isGenome());
		if ($hash->{'type'} eq 'substitution') {
			$type = 'snps';
			$hash->{'structural_type'} = 'snp';
		}
		elsif ($hash->{'type'} eq 'insertion') {
			$type = 'insertions';
			$hash->{'structural_type'} = 'ins';
		}
		elsif ($hash->{'type'} eq 'deletion') {
			$type = 'deletions';
			$hash->{'structural_type'} = 'del';
		}
		elsif ($hash->{'type'} eq 'large_duplication') {
			$type = 'insertions';
			$hash->{'structural_type'} = 'cnv:dup';
		}
		elsif ($hash->{'type'} eq 'large_deletion') {
			$type = 'deletions';
			$hash->{'structural_type'} = 'cnv:del';
		}
		if ($is_genome) { $db_gnomad = 'gnomad-genome'; }
		else { $db_gnomad = 'gnomad-exome'; }
		foreach my $fam_name (keys %{$hIds->{$id}->{json_ped}}) {
			$hash->{'json_ped_'.$fam_name} = $hIds->{$id}->{json_ped}->{$fam_name};
			my $fam = $project->getFamily($fam_name);
			foreach my $patient (@{$fam->getPatients()}) {
				$hash->{'json_ped_'.$patient->name()} = $hash->{'json_ped_'.$fam_name};
				if ($patient->isChild() and $patient->getVariantsVector($chr)->contains($v_id)) {
					my $transmission = $var->getTransmissionModel( $fam, $patient);
					push(@{$hash->{'patient_transmission'}}, $transmission);
				}
				
			}
		}
		push(@lRes, $hash);
	}
	return \@lRes;
}

sub getIds_onlive {
	my ($buffer, $project, $chr, $list_select_ids) = @_;
	my @lRes;
	my $hRegionsHo;
	foreach my $v_id (@$list_select_ids) {
		my $var = $chr->getVarObject($v_id);
		my $var_id = $var->id();
		$var->{project} = $project;
		$var->{chromosome_object} = $chr;
		my @lGenesObj = @{$var->getGenes()};
		my @lTransObj = @{$var->getTranscripts()};
		my $hash;
		$hash->{'name'} = $var_id;
		$hash->{'name'} = $var->rs_name() if ($var->rs_name());
		$hash->{'id'} = $var_id;
		$hash->{'chromosome'} = $chr->name();
		my $pos = $var->position($chr);
		$hash->{'start'} = $pos->start();
		$hash->{'end'} = $pos->end();
		$hash->{'text'} = $chr->sequence($pos->start(), $pos->end())."/".$var->sequence();
		$hash->{'heterozygote'} = 0;
		$hash->{'homozygote'} = 0;
		my $hVarPatInfos = $var->{annex};
		
		my @lVar_same_pos;
		my ($check_before_done, $check_after_done);
		my $i = $v_id;
		my $j = $v_id;
		$check_before_done = 1 if $i == 0;
		$check_after_done = 1 if $j+1 == $chr->size_vector();
		while ($check_before_done == undef) {
			$i--;
			my $var_before = $chr->getVarObject($i);
			if ($var_before->start() == $var->start()) { push(@lVar_same_pos, $var_before); }
			else { $check_before_done = 1; }
		}
		while ($check_after_done == undef) {
			$j++;
			my $var_after = $chr->getVarObject($j);
			if ($var_after->start() == $var->start()) { push(@lVar_same_pos, $var_after); }
			else { $check_after_done = 1; }
		}
		
		my $hPat;
		foreach my $patient (@{$var->getPatients()}) {
			$hPat->{$patient->id()} = $patient->name();
		}
		foreach my $patient (@{$project->getPatients()}) {
			my $pat_id = $patient->id();
			my $vector_regionho = $patient->getRegionHo($chr, 25)->Clone();
			#if ($project->filter_nbvar_regionho()) { $vector_regionho = $patient->getRegionHo($chr)->Clone(); } 
			#else { $vector_regionho = $patient->getRegionHo($chr, 25)->Clone(); }
			if ($vector_regionho->contains($v_id)) {
				push(@{$hash->{'patient_inregionho'} }, 'true');
				my @lReg = split(',', $vector_regionho->to_Enum());
				my $hThisRegion;
				foreach my $region (@lReg) {
					my $vector_region = Bit::Vector->new_Enum($chr->getVariantsVector->Size(), $region);
					$vector_region->contains($v_id);
					next if ($vector_region->is_empty());
					my ($v_id_start, $v_id_end) = split('-', $region);
					my $start = $chr->getVarObject($v_id_start)->start();
					my $end   = $chr->getVarObject($v_id_end)->end();
					$vector_region->Intersection($vector_region, $patient->getHo($chr));
					my $nb_var = $chr->countThisVariants($vector_region);
					my $region_id = $patient->name().';chr'.$chr->id().';'.$start.';'.$end.';'.$nb_var;
					unless (exists $hRegionsHo->{$region_id}) {
						$hRegionsHo->{$region_id} = 1;
						my $length = $end - $start + 1;
						$hThisRegion->{$nb_var}->{id} = $region_id;
						$hThisRegion->{$nb_var}->{patient_name} = $patient->name();
						$hThisRegion->{$nb_var}->{start} = $start;
						$hThisRegion->{$nb_var}->{end} = $end;
						$hThisRegion->{$nb_var}->{length} = $length;
						$hThisRegion->{$nb_var}->{nb_var} = $nb_var;
					}
				}
				if ($hThisRegion) {
					my @lNbVar = sort {$a <=> $b} keys %$hThisRegion;
					@lNbVar = reverse @lNbVar;
					foreach my $nb_var (@lNbVar) {
						if ($var->start() >= $hThisRegion->{$nb_var}->{start} and $var->start() <= $hThisRegion->{$nb_var}->{end}) {
							push(@{$hash->{'regions_ho'}}, $hThisRegion->{$nb_var});
							last;
						}
					}
				}
			}
			else { push(@{$hash->{'patient_inregionho'} }, 'false'); }
			
			my $h_infos_var = $var->sequencing_infos->{$patient->id()};
			if ($h_infos_var) {
				push(@{$hash->{'patient_id'}}, $pat_id);
				push(@{$hash->{'patient_name'}}, $hPat->{$pat_id});
				#push(@{$hash->{'patient_unifiedgenotyper'}  }, $hVarPatInfos->{$pat_id}->{'score'});
				push(@{$hash->{'patient_unifiedgenotyper'}  }, '.');
				my $nb_all_ref = $var->getNbAlleleRef($patient);
				my $nb_all_mut = $var->getNbAlleleAlt($patient);
				push(@{$hash->{'patient_unifiedgenotyper4'} }, $var->getMeanDP($patient));
				push(@{$hash->{'patient_unifiedgenotyper5'} }, 'uni');
				my @lBases;
				if ($var->isHeterozygote($patient)) {
					$hash->{'heterozygote'}++;
					my $is_he_ho = 'he';
					my $base = $var->getGenotype($patient);
					push(@lBases, $base);
					
					
					foreach my $v_other (@lVar_same_pos) {
						my $h_infos_var_other = $v_other->sequencing_infos->{$patient->id()};
						next unless ($h_infos_var_other);
						my $base_alt = $v_other->getGenotype($patient);
						if ($v_other->isHeterozygote($patient)) {
							$is_he_ho .= '  he';
						}
						if ($v_other->isHomozygote($patient)) {
							$is_he_ho .= '  ho';
						}
						push(@lBases, $base_alt);
					}
					push(@{$hash->{'patient_unifiedgenotyperheho'}}, $is_he_ho);
					$base = join(' / ', @lBases) if (scalar @lBases > 1);
					push(@{$hash->{'patient_base'}}, $base);
				}
				if ($var->isHomozygote($patient)) {
					$hash->{'homozygote'}++;
					my $is_he_ho = 'ho';
					my $base = $var->getGenotype($patient);
					push(@lBases, $base);
					foreach my $v_other (@lVar_same_pos) {
						my $h_infos_var_other = $v_other->sequencing_infos->{$patient->id()};
						next unless ($h_infos_var_other);
						if ($v_other->isHeterozygote($patient)) { $is_he_ho .= '  he'; };
						if ($v_other->isHomozygote($patient)) { $is_he_ho .= '  ho'; };
						my $base_alt = $v_other->getGenotype($patient);
						push(@lBases, $base_alt);
					}
					push(@{$hash->{'patient_unifiedgenotyperheho'}}, $is_he_ho);
					$base = join(' / ', @lBases) if (scalar @lBases > 1);
					push(@{$hash->{'patient_base'}}, $base);
				}
				push(@{$hash->{'patient_unifiedgenotyper2'} }, $nb_all_ref);
				push(@{$hash->{'patient_unifiedgenotyper3'} }, $nb_all_mut);
				$hash->{'score'} = $hVarPatInfos->{$pat_id}->{'score'} if ($hVarPatInfos->{$pat_id}->{'score'} > $hash->{'score'});
				push(@{$hash->{'patient_found'}}, 'yes');
			}
			else {
				push(@{$hash->{'patient_id'}}, $pat_id);
				push(@{$hash->{'patient_name'}}, $patient->name());
				my $nb_all_ref;
				my $nb_all_mut;
				my $first_v_alt = 1;
				my (@lBases, $is_he_ho, $base);
				foreach my $v_other (@lVar_same_pos) {
					my $h_infos_var_other = $v_other->sequencing_infos->{$patient->id()};
					next unless ($h_infos_var_other);
					if ($first_v_alt == 1) {
						$is_he_ho = '';
						$nb_all_ref = '';
						$nb_all_mut = '';
						$first_v_alt = undef;
					}
					my $base_alt = $v_other->getGenotype($patient);
					if ($v_other->isHeterozygote($patient)) {
						$is_he_ho .= '  he';
					};
					if ($v_other->isHomozygote($patient)) {
						$is_he_ho .= '  ho';
					};
					if ($nb_all_ref and $nb_all_mut) {
						$nb_all_ref .= ';';
						$nb_all_mut .= ';';
					}
#					if (exists $v_other->sequencing_infos->{$patient->id()}->{'nb_all_ref'} and exists $v_other->sequencing_infos->{$patient->id()}->{'nb_all_other_mut'}) {
#						$nb_all_ref .= $v_other->sequencing_infos->{$patient->id()}->{'nb_all_ref'} + $v_other->sequencing_infos->{$patient->id()}->{'nb_all_other_mut'};
#					}
#					else {
#						$nb_all_ref .= $v_other->sequencing_infos->{$patient->id()}->{'nb_all_ref'};
#					}
#					$nb_all_mut .= $v_other->sequencing_infos->{$patient->id()}->{'nb_all_mut'};
					
					push(@lBases, $base_alt);
				}
				my $nb_bases = scalar(@lBases);
				if ($nb_bases == 0) {
					$base = $var->getGenotype($patient);
					push(@lBases, $base);
					$is_he_ho = '.';
				}
				elsif (scalar @lBases == 1) {
					$base = $lBases[0];
				}
				else {
					my @lBases2;
					foreach my $this_base (@lBases) {
						my @t = split('/', $this_base);
						push(@lBases2, $t[1]);
					}
					$base = join('/', @lBases2);
				}
				
				push(@{$hash->{'patient_base'}}, $base);
				push(@{$hash->{'patient_unifiedgenotyper2'} }, $nb_all_ref);
				push(@{$hash->{'patient_unifiedgenotyper3'} }, $nb_all_mut);
				
				push(@{$hash->{'patient_unifiedgenotyper'}  }, '.');
				my $thisCov = $patient->meanDepth($var->getChromosome->name, $var->start, $var->end);
				#	push(@{$hash->{'patient_unifiedgenotyper4'} }, );
				if ($thisCov) { push(@{$hash->{'patient_unifiedgenotyper4'} }, int($thisCov)); }
				else { push(@{$hash->{'patient_unifiedgenotyper4'} }, '.'); }
				push(@{$hash->{'patient_unifiedgenotyper5'} }, '.');
				$hash->{'not'}++;
				push(@{$hash->{'patient_unifiedgenotyperheho'}}, $is_he_ho);
				$hash->{'score'} = '.';
				push(@{$hash->{'patient_found'}}, 'no');
			}
			if ($patient->isChild() and $patient->getVariantsVector($chr)->contains($v_id)) {
				my $transmission = $var->getTransmissionModelType( $patient->getFamily(), $patient);
				push(@{$hash->{'patient_transmission'}}, $transmission);
				my @l_scale_score;
				my $max_score = -999;
#				foreach my $g (@lGenesObj) {
#					#my $score = $g->external_name().':'.$var->scaledScoreVariant($g, $patient);
#					#push(@l_scale_score, $score);
#					my $score = $var->scaledScoreVariant($g, $patient);
#					$max_score = $score if ($score > $max_score)
#				}
				#push(@{$hash->{'patient_scale_score'}}, join(' ', @l_scale_score));
				$max_score = ' ' if ($max_score == -999);
				push(@{$hash->{'patient_scale_score'}}, $max_score);
			}
			else {
				push(@{$hash->{'patient_transmission'}}, '');
				push(@{$hash->{'patient_scale_score'}}, '');
			}
			my @lTmpBam = split('ngs', $patient->getBamFiles->[0]); 
			push(@{$hash->{'patient_bam'}}, '/NGS'.$lTmpBam[1]);
		}
		$hash->{'allpatients'} = join(';', @{$hash->{patient_name}});
		$hash->{'nb_patients'} = scalar(@{$hash->{patient_name}});
		$hash->{'homo_hetero'} += 1 if ($hash->{'homozygote'} > 0);
		$hash->{'homo_hetero'} += 2 if ($hash->{'heterozygote'} > 0);
		$hash->{'polyphen_status'} = $var->polyphenStatus();
		$hash->{'sift_status'} = $var->siftStatus();
		$hash->{'cadd_score'} = $var->cadd_score();
		my @l_dbscsnv;
		push(@l_dbscsnv, 'rf:'.$var->dbscsnv->{'rf'}) if exists $var->dbscsnv->{'rf'};
		push(@l_dbscsnv, 'ada:'.$var->dbscsnv->{'ada'}) if exists $var->dbscsnv->{'ada'};
		$hash->{'dbscsnv'} = join(' ', @l_dbscsnv);
		
#		if ($var->isStop()) {
#			$hash->{'polyphen_status'} = 5;
#			$hash->{'sift_status'} = 5;
#		}
#		if ($var->isPhase()) {
#			$hash->{'polyphen_status'} = 4;
#			$hash->{'sift_status'} = 4;
#		}
		if ($var->isVariation()) {
			$hash->{'structural_type'} = 'snp';
			$hash->{'type'} = 'variations';
		}
		elsif ($var->isInsertion()) {
			$hash->{'structural_type'} = 'ins';
			$hash->{'type'} = 'insertions';
		}
		elsif ($var->isDeletion()) {
			$hash->{'structural_type'} = 'del';
			$hash->{'type'} = 'deletions';
		}
		elsif ($var->isCnv()) {
			$hash->{'structural_type'} = 'cnv';
			$hash->{'type'} = 'cnv';
		}
		foreach my $gene (@lGenesObj) {
			my ($gene_id, $id_chr) = split('_', $gene->id());
			push(@{$hash->{'genes'}}, $gene_id);
		}
		#$hash->{reference} = $var->getReference->name();
		my @lGenesOmim;
		foreach my $gene (@lGenesObj) {
			if ($gene->omim_id()) { push(@lGenesOmim, $gene->omim_id()); }
			else { push(@lGenesOmim, undef); }
			my ($gene_id, $this_chr) = split('_', $gene->id());
			my $polyphen = $var->polyphenStatus($gene);
			my $sift = $var->siftStatus($gene);
			eval { $hash->{$gene_id.'_consequence'} = $var->variationTypeInterface($gene); };
			if ($@) { $hash->{$gene_id.'_consequence'} = 'Error...'; }
		}
		$hash->{'gene_in_omim'} = join(';', @lGenesOmim);
		foreach my $transcript (@lTransObj) {
			eval { $hash->{'consequence!'.$transcript->id()} = $var->variationTypeInterface($transcript); };
			if ($@) {
				$hash->{'consequence!'.$transcript->id()} = 'Error...';
				next;
			}
			$hash->{'nomenclature!'.$transcript->id()} = $var->getNomenclature($transcript);
			my $gene = $transcript->getGene();
			my $prot = $transcript->getProtein();
			my $hTrans;
			$hTrans->{'gene'} = $gene->name()." (".$gene->external_name().")";
			$hTrans->{'transcript'} = $transcript->name()."+".$transcript->{'external_name'};
			$hTrans->{'cdna_position'} = $transcript->translate_position($var->start());
			$hTrans->{'exon'} = $transcript->findExonNumber($var->start());
			$hTrans->{'exon'} = $transcript->findNearestExon($var->start()) if $hTrans->{'exon'} == -1;
			$hTrans->{'name'} = $var->name();
			$hTrans->{'description'} = $transcript->getGene()->description();
			$hTrans->{'protein'} = $prot->name()."+".$transcript->{'external_protein_name'} if ($prot);
			if ($var->isCoding($transcript)) {
				$hTrans->{'nomenclature'} = $var->getNomenclature($transcript);
				my ($polyphen, $sift);
				if ($prot) {
					$hTrans->{'cds_position'} = $var->getOrfPosition($prot);
					$hTrans->{'AA_variation'} = $var->changeAA($prot);
					$hTrans->{'AA_protein'} = $var->getProteinAA($prot);
					$hTrans->{'protein_position'} = $var->getProteinPosition($prot);
					$hTrans->{'AA_position'}  = $hTrans->{'protein_position'};
					$hTrans->{'consequence'} = $hTrans->{'AA_protein'}."/".$hTrans->{'AA_variation'}." (".$var->getCodons($transcript).")";
					$hTrans->{'consequence'} = "ins:".$var->sequence() if ($var->isInsertion());
					$hTrans->{'consequence'} = "del:".$var->sequence() if ($var->isDeletion());
					$polyphen = $var->polyphenStatus($prot);
					$sift = $var->siftStatus($prot);
				}
				if ($var->isStop($transcript)) {
					$polyphen = 5;
					$sift = 5;
				}
				if ($var->isPhase($transcript)) {
					$polyphen = 4;
					$sift = 4;
				}
				$hTrans->{'polyphen_status'} = $polyphen."+".$var->polyphenScore($prot);
				$hTrans->{'sift_status'} = $sift."+".$var->siftScore($prot);
			}
			else { $hash->{'consequence'} = $var->variationType($transcript); }
			$hTrans->{'polyphen_html'} = '-';
			$hTrans->{'revel_score'} = $var->revel_score($transcript);
			$hTrans->{'revel_score'} = '-' if $hTrans->{'revel_score'} eq '-99';
			$hTrans->{'amissense'} = $var->alphamissense($transcript);
			
			push(@{$hash->{'tab_consequences'}}, $hTrans);
			
			
		}
		#$hash->{'public'} = $var->origin_database();
		$hash->{'filter'} = $var->getNGSScore();
		eval { $hash->{'consequence!all'} = $var->variationType(); };
		if ($@) { $hash->{'consequence!all'} = 'Error...'; }
		push(@lRes, $hash);
#		warn Dumper $hash; die;
	}
	return \@lRes;
}

sub getIds {
	my ($buffer,$project,$chr,$temp_type_name,$select_ids) = @_;
	my $file_cache_variation = $project->getCacheVariationsFile($chr->name);
	my $data;
	my $flag;
	
	($data,$flag) = getData($project,$chr,$temp_type_name);
	
	if ($select_ids ){
	my @res;
	my %hids;
	my @keep;
	
	foreach my $id (@$select_ids){
	if (exists $data->{$id}) {
		my $hh; 
		if ($flag ==1){
			$hh = thaw $data->{$id};
		}
		else {
			$hh =  $data->{$id};
		}
	 push(@res,$hh);
	}
	 # if exists $data->{$id};
	 #next  if exists $data->{$id};
	 push(@keep,$id);
	}
	@res = sort{$a->{start} <=> $b->{start}} @res;
	
	@$select_ids = @keep;
	return (\@res);
	}
	else {
		
		return $data;
	}
	
}

sub getData{
	my ($project,$chr,$type_name) = @_;
	my $file_cache_variation = $project->getCacheVariationsFile($chr->name);
	my $kyoto_cache_variation = $project->getCacheVariationsKyotoFile($chr->name);
	my $data;
	
	if (-e $kyoto_cache_variation){
		my $hvars;
		
		#tie(%{$hvars}, 'KyotoCabinet::DB', $kyoto_cache_variation , KyotoCabinet::DB::ONOLOCK || KyotoCabinet::DB::OREADER) || die(" not file ".$kyoto_cache_variation);
		return ($hvars,1);
		
	}
	elsif (-e $file_cache_variation){
	 	 return ($project->buffer->getStore($file_cache_variation),0);
	 	warn "file cache !!!!!!!!!!!!!!";

	 }
	 else {
		my $id_storable = GenBoStorable::getStoreId($project->buffer->dbh, $project->id, $chr->id,$type_name );
		unless ($id_storable) {
			exit(0);
		}
		return (GenBoStorable::getStore( $project->buffer->dbh, $id_storable ),0);
 	
	 }
}



sub getIdsByPosition {
	my ($buffer,$project,$chr,$start,$end,$type_name) = @_;
	my $data = getData($project,$chr,$type_name);
 	my $set = Set::IntSpan::Fast->new($start."-".$end);
	my @res;
	foreach my $id (keys %$data){
		next unless $set->contains($data->{$id}->{start});
		 push(@res,$data->{$id});
		
	}
	@res = sort{$a->{start} <=> $b->{start}} @res;
	return \@res;
}

1;
