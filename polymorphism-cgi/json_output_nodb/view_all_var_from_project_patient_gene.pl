#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;

use strict;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../cache_nodb/scripts/";

use GBuffer;
use export_data;
use JSON;
use VcfMerge;
use GenBoNoSql;
use Set::IntervalTree;
use Spreadsheet::WriteExcel;
use Compress::Snappy;
use Storable qw(store retrieve freeze dclone thaw);
use POSIX qw(strftime);
use CGI::Session;
use html; 
use Carp;
use Cache_Commons;
use QueryVectorFilter;
use IO::Handle;

require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";

my $io = IO::Handle->new();
$io->autoflush(1);

my $cgi = new CGI();
my $project_name = $cgi->param('project');
my $patient_name = $cgi->param('patient');
my $gene_name = $cgi->param('gene');
my $gene_transcript = $cgi->param('gene_transcript');
my $only_model = $cgi->param('only_model');
my $only_type = $cgi->param('only_type');
my $keep_var_id = $cgi->param('keep');
my $skip_header = $cgi->param('skip_header');
my $hide_filtred_model = $cgi->param('hide_filtred_model');

confess("\n\nERROR: -project option mandatory. Die.\n\n") unless ($project_name);
confess("\n\nERROR: -patient option mandatory. Die.\n\n") unless ($patient_name);
my $transcript_name;
if ($gene_transcript) {
	my ($this_gene_id, $this_tr_id) = split(':', $gene_transcript);
	$transcript_name = $this_tr_id;
	$gene_name = $this_gene_id;
}

confess("\n\nERROR: -gene option mandatory. Die.\n\n") unless ($gene_name);
#confess("\n\nERROR: -keep option mandatory. Die.\n\n") unless ($keep_var_id);

my $class;
my $fsize = "font-size:10px";

my $value_red = 10;
my $value_coral = 9;
my $value_orange = 8;
my $value_yellow = 5;

my @headers_validations = ("igv","alamut","var_name","trio","gnomad","deja_vu","table_validation","table_transcript");
my @header_transcripts = ("consequence","enst","nm","ccds","appris","exon","nomenclature","codons","codons_AA", "polyphen","sift","ncboost","cadd","revel","dbscsnv","spliceAI");
my @header_transcripts_cnv = ("consequence", "enst", "nm", "ccds", "appris", "start", "end");

my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => $project_name );
my $patient = $project->getPatient($patient_name);

my $h_new = $buffer->queryHgmd->get_hash_last_released_DM();
my $hAllVarIds;
foreach my $acc_num (keys %{$h_new}) {
	$project->print_dot(5);
	next unless $h_new->{$acc_num}->{tag} eq 'DM';
	my $chr_id = $h_new->{$acc_num}->{chrom};
	my $pos = $h_new->{$acc_num}->{pos};
	my $ref = $h_new->{$acc_num}->{ref};
	my $alt = $h_new->{$acc_num}->{alt};
	my $var_id = $chr_id.'_'.$pos.'_'.$ref.'_'.$alt;
	$hAllVarIds->{$var_id}->{hgmd} = 1;
}
foreach my $var_id (@{$buffer->queryClinvarPathogenic->getAllVarIds_onlyLastRelease()}) {
	$hAllVarIds->{$var_id}->{clinvar} = 1;
}

my @lGenes;
if ($keep_var_id) {
	my $v = $project->getVariant($keep_var_id);
	@lGenes = (@{$v->getGenes()});
}
else {
	push(@lGenes, $project->newGene($gene_name));
} 
eval {
	foreach my $gene (@lGenes) {
		if (($gene->id() eq $gene_name) or ($gene->external_name() eq $gene_name)) {
			my $chr = $gene->getChromosome();
			my $vector_chr = $chr->getVariantsVector();
			my $vector_gene = $gene->getVariantsVector();
			my @lVar;
			if ($patient eq 'all') {
				@lVar = @{$chr->getListVarObjects($vector_gene)};
			}
			else {
				my $vector_patient = $patient->getVariantsVector($chr);
				$vector_patient->Intersection($vector_patient, $vector_gene);
				if ($only_type eq 'cnv') {
					my $v_cnv = $chr->getNewVector();
					$v_cnv += $chr->global_categories->{'large_deletion'} if (exists $chr->global_categories->{'large_deletion'});
					$v_cnv += $chr->global_categories->{'large_duplication'} if (exists $chr->global_categories->{'large_duplication'});
					$vector_patient->Intersection($vector_patient, $v_cnv);
				}
				elsif ($only_type eq 'nocnv') {
					$vector_patient -= $chr->global_categories->{'large_deletion'} if (exists $chr->global_categories->{'large_deletion'});
					$vector_patient -= $chr->global_categories->{'large_duplication'} if (exists $chr->global_categories->{'large_duplication'});
				}
				@lVar = @{$chr->getListVarObjects($vector_patient)};
			}
			view_html_variants($patient, $gene, \@lVar);
			exit(0);
		}
	}
};
if($@) {
	
	warn Dumper $@;
	
	my $hashRes;
	$hashRes->{html_table} = "Not Available -> Analyse too old for project $project_name (PolyViewer not available for $project_name / patient $patient_name)";
	print $cgi->header('text/json-comment-filtered');
	my $json_encode = encode_json $hashRes;
	print $json_encode;
	exit(0);
}

sub view_html_variants {
	my ($patient, $g, $lVar) = @_;
	my $tr_id = 'tr_'.$patient->project->name();
	my $color = "#9BC8A5";
	my $out2;
	my ($h_scaled_score, $hgenes, $hvariation_selected, $hvariation_by_model);
	
	my $hTrInGene;
	foreach my $tr (@{$g->getTranscripts()}) {
		$hTrInGene->{$tr->id()} = undef;
	}
	
	# CAS recherche avec variant Down/Up stream
	if ($keep_var_id) {
		my $found_keep_var_id;
		foreach my $v (@$lVar) {
			if ($v->id() eq $keep_var_id) {
				$found_keep_var_id = 1;
				last;
			}
		}
		unless ($found_keep_var_id) {
			my $v = $project->getVariant($keep_var_id);
			push(@$lVar, $v);
		}
	}
	
	foreach my $v (@$lVar) {
		my $var_id = $v->id();
		my @lTrios;
		$v->{isDM} = 1 if (exists $hAllVarIds->{$var_id}->{hgmd});
		$v->{isClinvarPathogenic} = 1 if (exists $hAllVarIds->{$var_id}->{clinvar});
		my $hvariation = update_variant_editor::construct_hash_variant_global($project, $v);
		update_variant_editor::construct_hash_variant_patient( $project, $v, $patient, $hvariation);
		$hvariation->{obj} = $v;
		$hvariation->{html}->{no_css_polydiag} = 1;
		
		my $this_score;
		if ($transcript_name and $var_id eq $keep_var_id) {
			foreach my $tr (@{$v->getTranscripts()}) {
				if ($tr->id() eq $transcript_name) {
					eval { $this_score = $v->scaledScoreVariant($tr, $patient, undef); };
					if ($@) { die; }
					last;
				}
			}
		}
		else {
			my $max_score = -99;
			my $max_tr;
			foreach my $tr (@{$v->getTranscripts()}) {
				my $this_score = $v->scaledScoreVariant($tr, $patient, undef);
				if ($this_score > $max_score) {
					$max_score = $this_score;
					$max_tr = $tr;
				}
			}
			$this_score = $max_score;
		}
		
		$hvariation->{scaled_score}->{$g->id()} = $this_score;
		$hgenes->{project_name} = $patient->project->name();
		$hgenes->{name} = $g->external_name;
		$hgenes->{description} = $g->description;
		$hgenes->{phenotypes} = $g->phenotypes;
		$hgenes->{score_mother} = 0;
		$hgenes->{score_father} = 0;
		$hgenes->{score_biallelic} = -10;
		$hgenes->{score} = $g->score;
		$hgenes->{id} = $g->id;
		$hgenes->{omim_inheritance} = $g->omim_inheritance;
		$hgenes->{external_name} = $g->external_name;
		$hgenes->{pLI} = $g->pLI;
		$hgenes->{omim_id} = $g->omim_id;
		$hgenes->{panels} = $buffer->queryPanel()->getPanelsForGeneName($g->external_name);
		$hgenes->{js_id} = $project_name."_".$g->id;
		$hgenes->{score_variant}->{$v->id} = $this_score;
		$hgenes->{pathogenic} ++  if  $v->isDM or $v->isClinvarPathogenic;
		$hgenes->{clinvar_hgmd} ++ if  $v->hgmd or $v->clinvar;
		#$hgenes->{denovo_rare} ++ if $v->getGnomadAC< 10 &&  $v->isDenovoTransmission($patient->getFamily,$patient);
		my $val = $v->score_validations($g);
		if ($val){
			$hgenes->{validations} = $val->{validation};
		}
		push (@{$hgenes->{variants}}, $v->id);
		if ($patient->isChild && $patient->getFamily->isTrio()){
			 if  ($v->isMotherTransmission($patient->getFamily(), $patient)){
			 	$hgenes->{score_mother} = $this_score;
			 }
			 elsif  ($v->isFatherTransmission($patient->getFamily(), $patient)){
			 	$hgenes->{score_father} = $this_score;
			 }
			 else {
			 	$hgenes->{score_biallelic} = $this_score;
			 }
		}
		else {
			$hgenes->{score_biallelic} = $this_score;
		}
		
		delete $v->{hgmd};
		delete $v->{isNewHgmd};
		delete $v->{hgmd_details};
		delete $v->{hgmd_id};
		delete $v->{hgmd_class};
#		delete $v->{isDM};
		delete $v->{hgmd_inheritance};
		delete $v->{hgmd_hgvs};
		delete $v->{hgmd_phenotype};
		delete $v->{hgmd_disease};
		delete $v->{hgmd_releases};
		delete $v->{hgmd_pubmed};
		delete $v->{clinvar};
		delete $v->{text_clinvar};
		delete $v->{score_clinvar};
#		delete $v->{isClinvarPathogenic};
		delete $v->{clinical};
		delete $hvariation->{html}->{hgmd};
		delete $hvariation->{html}->{clinvar};
		update_variant_editor::vhgmd($v,$hvariation);
		update_variant_editor::vclinvar($v,$hvariation);
		if (exists $hAllVarIds->{$var_id}->{hgmd}) { $hvariation->{html}->{hgmd} .= "<img src='images/polyicons/new.jpeg'>"; }
		if (exists $hAllVarIds->{$var_id}->{clinvar}) { $hvariation->{html}->{clinvar} .= "<img src='images/polyicons/new.jpeg'>"; }
		delete $hvariation->{html}->{table_validation};
		update_variant_editor::table_validation($patient,$hvariation,$g);
		#update_variant_editor::trio($v, $hvariation, $patient);
		
		my $model = lc($v->getTransmissionModel($patient->getFamily(), $patient));
		if ($keep_var_id and $keep_var_id eq $v->id()) {
			$this_score += 100000;
			$h_scaled_score->{$this_score}->{$var_id} = $hvariation;
			$hvariation_selected = $hvariation;
		}
		elsif ($only_model eq 'denovo' and $model =~ /denovo/) {
			$this_score += 900000;
			$hvariation_by_model->{denovo}->{$var_id} = undef;
		}
		elsif ($only_model eq 'dominant' and $model eq 'dominant') {
			$this_score += 900000;
			$hvariation_by_model->{dominant}->{$var_id} = undef;
		}
		elsif ($only_model eq 'recessive' and $model eq 'recessive') {
			$this_score += 900000;
			$hvariation_by_model->{recessive}->{$var_id} = undef;
		}
		elsif ($only_model eq 'mosaic' and lc($model) eq 'mosaic') {
			$this_score += 900000;
			$hvariation_by_model->{mosaic}->{$var_id} = undef;
		}
		elsif ($only_model eq 'uniparental' and lc($model) eq 'uniparental disomy') {
			$this_score += 900000;
			$hvariation_by_model->{uniparental}->{$var_id} = undef;
		}
		elsif ($only_model eq 'mother' and ($model eq 'mother' or ($model eq 'compound' and $v->isMotherTransmission($patient->getFamily(),$patient)))) {
			$this_score += 50000;
			$hvariation_by_model->{mother}->{$var_id} = undef;
		}
		elsif ($only_model eq 'father' and ($model eq 'father' or ($model eq 'compound' and $v->isFatherTransmission($patient->getFamily(),$patient)))) {
			$this_score += 50000;
			$hvariation_by_model->{father}->{$var_id} = undef;
		}
		elsif (lc($model) eq 'both') {
			$this_score += 20000;
			$hvariation_by_model->{both}->{$var_id} = undef;
		}
		else {
			$hvariation_by_model->{others}->{$var_id} = undef;
		}
		$h_scaled_score->{$this_score}->{$var_id} = $hvariation;
	}
	
	my $nb_red = 0;
	my $nb_orange = 0;
	my $nb_grey = 0;
	foreach my $scale_score (keys %{$h_scaled_score}) {
		if ($scale_score >= $value_red) { $nb_red += scalar keys(%{$h_scaled_score->{$scale_score}}); }
		elsif ($scale_score >= $value_orange) { $nb_orange += scalar keys(%{$h_scaled_score->{$scale_score}}); }
		else { $nb_grey += scalar keys(%{$h_scaled_score->{$scale_score}}); }
	}
	
	my $panel_gene_done;
	my @lScores = sort {$a <=> $b} keys %{$h_scaled_score};
	foreach my $hash_scale_score (reverse @lScores) {
		my ($out3, $out3_header);
		foreach my $var_id (sort keys %{$h_scaled_score->{$hash_scale_score}}) {
			my $scale_score = $hash_scale_score;
			unless ($out2) {
				$out2 = "<center>";
				if ($skip_header) {
					if ($skip_header == 1) {
						$out2 .= "<div style='width:97%;background-color:#F3F3F3;border: 1px solid black;box-shadow: 1px 1px 6px black;'>";
						my $my_style_ped = "";
						if ($patient->getFamily->isTrio && $patient->isChild()){
							my $t_tab = update_variant_editor::print_table_nav_trio($patient, $my_style_ped);
							$t_tab =~ s/position:relative;left:10%;//;
							$out2 .= $t_tab;
							$out2 .= "<br></div><br>";
						}
						else {
							$out2 .= "<br>";
							my $t_tab = update_variant_editor::print_table_nav_solo($patient, $my_style_ped);
							$t_tab =~ s/position:relative;left:10%;//;
							$out2 .= $t_tab;
							$out2 .= "</div><br>";
						}
					}
				}
				else {
					$out2 .= "<div style='width:97%;background-color:#F3F3F3;border: 1px solid black;box-shadow: 1px 1px 6px black;'>";
					my $my_style_ped = "";
					if ($patient->getFamily->isTrio && $patient->isChild()){
						$out2 .= update_variant_editor::print_table_nav_trio($patient, $my_style_ped);
					}
					else {
						$out2 .= update_variant_editor::print_table_nav_solo($patient, $my_style_ped);
					}
					$out2 .= "<br>";
					
					my $sex_status;
					if ($patient->isChild()) {
						if ($patient->isIll()) {
							if ($patient->sex() eq '1') { $sex_status = "<img src='/icons/Polyicons/baby-boy-d.png'>"; }
							else { $sex_status = "<img src='/icons/Polyicons/baby-girl-d.png'>"; }
						}
						else {
							if ($patient->sex() eq '1') { $sex_status = "<img src='/icons/Polyicons/baby-boy-s.png'>"; }
							else { $sex_status = "<img src='/icons/Polyicons/baby-girl-s.png'>"; }
						}
					}
					elsif ($patient->isFather()) {
						if ($patient->isIll()) { $sex_status = "<img src='/icons/Polyicons/male-d.png'>"; }
						else { $sex_status = "<img src='/icons/Polyicons/male-s.png'>"; }
					}
					elsif ($patient->isMother()) {
						if ($patient->isIll()) { $sex_status = "<img src='/icons/Polyicons/female-d.png'>"; }
						else { $sex_status = "<img src='/icons/Polyicons/female-s.png'>"; }
					}
					$out2 .= "<span style='padding-left:50px'><b>Project:</b> ".$project->name()."</span>, <span style:'padding-left:30px;'><b>Description:</b> ".$project->description()."</span>";
					$out2 .= "<br><br>";
					$out2 .= "<span style='padding-left:50px'><b>Family:</b> ".$patient->getFamily->name()."</span> , <span style:'padding-left:30px;'><b>Name:</b> ".$patient->name()."</span> , <span style:'padding-left:30px;'><b>Status:</b> ".$sex_status."</span>";
					$out2 .= "<br><br>";
					$out2 .= "<span style='padding-left:50px'><b>Gene Name:</b> ".$g->external_name()."</span> , <span style:'padding-left:30px;'><b>Gene ID:</b> ".$g->ensg()."</span> , <span style:'padding-left:30px;'><b>Transmission:</b> ".$g->omim_inheritance()."</span> ";
					$out2 .= "<br>";
					my $phenotypes = $g->short_phenotypes();
					$out2 .= "<span style='padding-left:50px;padding-bottom:20px;'><b>Phenotype(s):</b> ".$phenotypes."<br></span>";
					$out2 .= "<br>";
					$out2 .= "</div>";
					$out2 .= "</center>";
					$out2 .= "<br>";
					$out2 .= "</tr>";
					$out2 .= "</center>";
					
					$out2 .= "<tr>";
					$out2 .= "<td>";
					# PANEL Gene de polydiag avec modif local de css
					unless ($panel_gene_done) {
						my $bgcolor = '#60798B';
						$out2 .= $cgi->start_div({class=>"panel-heading panel-face panel-grey",style=>"background-color:$bgcolor;min-height:40px;border:0px;width:100%;"});
						my $panel_id = "panel_".$hgenes->{js_id};
						$hgenes->{max_score} = $scale_score;
						if ($hgenes->{max_score} >= 90000) { $hgenes->{max_score} -= 100000; }
						elsif ($hgenes->{max_score} >= 40000) { $hgenes->{max_score} -= 50000; }
						$hgenes->{max_score} = sprintf("%.1f", $hgenes->{max_score});
						my $this_html = update_variant_editor::panel_gene($hgenes, $panel_id, $project_name, $patient);
						#$this_html =~ s/607[Dd]8[Bb]/5A6268/g;
						$this_html =~ s/text\-shadow\:1px 1px 2px black/color:white/g;
						$out2 .= $this_html;
						$out2 .= $cgi->end_div();
						$out2 .= "</td>";
						$out2 .= "</tr>";
						$out2 .= "<tr>";
						$out2 .= "<td>";
						$panel_gene_done = 1;
					}
				}
				$out2 .= "<div>";
				
				$out2 .= "<center>";
				my $color_div_1 = 'black';
				$color_div_1 = '#779ECB' if ($only_model eq 'mother');
				$color_div_1 = '#F7C9C9' if ($only_model eq 'father');
				$out2 .= qq{<div style='font-size:13px;width:97%;background-color:white;max-height:300px;overflow-y:scroll;border: 1px solid $color_div_1;box-shadow: 1px 1px 6px $color_div_1;'>};
				$out2 .= qq{<table class='table' style='padding:5px;font-size:13px;'>};
				$out2 .= get_header_table();
			}
			$out2 .= "</div>";
			$out2 .= "</div>";
			$out2 .= "</td>";
			$out2 .= "</tr>";
			
			my $id_tr = 'tr_'.$var_id;
			my $class_tr = 'tr_variant';
			if (exists $hvariation_by_model->{both}->{$var_id}) { $class_tr .= ' tr_both'; }
			if (exists $hvariation_by_model->{mother}->{$var_id}) { $class_tr .= ' tr_mother'; }
			if (exists $hvariation_by_model->{father}->{$var_id}) { $class_tr .= ' tr_father'; }
			if (exists $hvariation_by_model->{others}->{$var_id}) { $class_tr .= ' tr_other'; }
			if (exists $hvariation_by_model->{denovo}->{$var_id}) { $class_tr .= ' tr_denovo'; }
			if (exists $hvariation_by_model->{dominant}->{$var_id}) { $class_tr .= ' tr_dominant'; }
			if (exists $hvariation_by_model->{recessive}->{$var_id}) { $class_tr .= ' tr_recessive'; }
			if (exists $hvariation_by_model->{mosaic}->{$var_id}) { $class_tr .= ' tr_mosaic'; }
			if (exists $hvariation_by_model->{uniparental}->{$var_id}) { $class_tr .= ' tr_uniparental'; }
			
			if ($hide_filtred_model and not exists $hvariation_by_model->{$only_model}->{$var_id}) { $class_tr .= ' hidden'; }
			
			my $stlye_tr = "background-color:white;";
			$out2 .= qq{<tr id="$id_tr" class="$class_tr" style="$stlye_tr">};
			my $hvariation = $h_scaled_score->{$hash_scale_score}->{$var_id};
			
			if ($scale_score >= 90000) { $scale_score -= 100000; }
			elsif ($scale_score >= 40000) { $scale_score -= 50000; }
			foreach my $h (@headers_validations){
				if (lc($h) eq "trio" or lc($h) eq "table_transcript"){
					$class->{style} = "background-color:white;min-width:10%;vertical-align:middle;padding:5px;";
				}
				elsif (lc($h) eq "igv" or lc($h) eq "alamut"){
					$class->{style} = "background-color:white;max-width:50px;vertical-align:middle;padding:5px;";
				}
				else {
					$class->{style} = "background-color:white;min-width:5%;vertical-align:middle;padding:5px;";
				}
				if ($h eq 'scaled_score') {
					my $color = 'black';
					if ($scale_score >= $value_red) { $color = '#CE0000'; }
					elsif ($scale_score >= $value_coral) { $color = 'coral'; }
					elsif ($scale_score >= $value_orange) { $color = '#EFB73E'; }
					my $score = sprintf("%.1f", $scale_score);
					my $b = "<button type='button' class='btn ' style='background-color:white;font-size:10px;'>";
					$b .= "<span style='color:$color;'>".$score."<span>";
					$b .= "</button>";
					$out2.= $cgi->td($class,$b);
				}
				elsif ($h eq 'table_transcript') {
					if ($hvariation->{value}->{type} eq 'large_deletion' or $hvariation->{value}->{type} eq 'large_duplication') {
						$out2 .= $cgi->td($class, update_variant_editor::table_transcripts_cnv($hvariation, $hvariation->{genes}->{ $g->{id} }, $g->{id}, \@header_transcripts_cnv));
					}
					else {
						$out2 .= $cgi->td($class, update_variant_editor::table_transcripts($hvariation->{genes}->{ $g->{id} }, \@header_transcripts));
					}
				}
				else {
					$out2.= $cgi->td($class,$hvariation->{html}->{$h});
				}
			}
			$out2 .= $cgi->end_Tr();
			
			if ($keep_var_id eq $var_id) {
				$out2 .= "</table>";
				$out2 .= "</div>";
				$out2 .= "</center>";
				
				$out2 .= "<div style='float:left;'>";
				$out2 .= "<span style='text-align:left;padding-left:25px;'><span class='glyphicon glyphicon-menu-up' aria-hidden='true'></span> <b>Selected Variant</b></span>";
				$out2 .= "</div>";
				
				$out2 .= "<br>";
				
				$out2 .= "<div style='float:right;'>";
				my $color_div_2 = 'black';
				my $text_nb_others;
				if ($only_model) {
					my $nb_others = scalar (keys %{$hvariation_by_model->{others}});
					if ($nb_others > 0) {
						$text_nb_others = qq{$nb_others VAR Other(s) Transmissions};
					}
				}
				if ($only_model eq 'mother') {
					my $nb_others = scalar (keys %{$hvariation_by_model->{mother}});
					$color_div_2 = '#F7C9C9';
					$out2 .= "<span style='padding-right:25px'><b>Found: </b><button onClick=\"select_tr('table_variants', 'tr_variant', 'tr_mother')\"><font color='red'>$nb_others VAR Mother Transmission</font></button> AND <button onClick=\"select_tr('table_variants', 'tr_variant', 'tr_other')\">$text_nb_others</button> <span class='glyphicon glyphicon-menu-down' aria-hidden='true'></span></span>";
				}
				elsif ($only_model eq 'father') {
					my $nb_others = scalar (keys %{$hvariation_by_model->{father}});
					$color_div_2 = '#779ECB';
					$out2 .= "<span style='padding-right:25px'><b>Found: </b><button onClick=\"select_tr('table_variants', 'tr_variant', 'tr_father')\"><font color='blue'>$nb_others VAR Father Transmission</font></button> AND <button onClick=\"select_tr('table_variants', 'tr_variant', 'tr_other')\">$text_nb_others</button> <span class='glyphicon glyphicon-menu-down' aria-hidden='true'></span></span>";
				}
				else {
					my $nb_others = scalar (keys %{$hvariation_by_model->{others}});
					$out2 .= "<span style='padding-right:25px'><b>Found: $nb_others Variant(s)</b> <span class='glyphicon glyphicon-menu-down' aria-hidden='true'></span></span>";
				}
				$out2 .= "</div>";
				$out2 .= "<center>";
				$out2 .= qq{<div style='font-size:13px;width:97%;background-color:white;max-height:300px;overflow-y:scroll;border: 1px solid $color_div_2;box-shadow: 1px 1px 6px $color_div_2;'>};
				$out2 .= "<table id='table_variants' class='table' style='padding:5px;font-size:13px;'>";
				$out2 .= get_header_table();
			}
		}
	}
	$out2 .= "</table>";
	$out2 .= "</div>";
	$out2 .= "</center>";
	$out2 .= "</div>";
	$out2 .= "</td>";
	$out2 .= "</tr>";

	my  $CSS = qq{
	<style type="text/css"> 
	.bs-callout {
	    padding: 20px;
	    margin: 20px 0;
	    border: 1px solid #eee;
	    border-left-width: 5px;
	    border-radius: 3px;
	}
	.bs-callout h4 {
	    margin-top: 0;
	    margin-bottom: 5px;
	}
	.bs-callout p:last-child {
	    margin-bottom: 0;
	}
	.bs-callout code {
	    border-radius: 3px;
	}
	.bs-callout+.bs-callout {
	    margin-top: -5px;
	}
	.bs-callout-default {
	    border-left-color: #777;
	}
	.bs-callout-default h4 {
	    color: #777;
	}
	.bs-callout-primary {
	    border-left-color: #428bca;
	}
	.bs-callout-primary h4 {
	    color: #428bca;
	}
	.bs-callout-success {
	    border-left-color: #5cb85c;
	}
	.bs-callout-success h4 {
	    color: #5cb85c;
	}
	.bs-callout-danger {
	    border-left-color: #d9534f;
	}
	.bs-callout-danger h3 {
	    color: #d9534f;
	}
	.bs-callout-danger h4 {
	    color: black;
	}
	.bs-callout-warning {
	    border-left-color: #f0ad4e;
	}
	.bs-callout-warning h4 {
	    color: #f0ad4e;
	}
	.bs-callout-info {
	    border-left-color: #5bc0de;
	}
	.bs-callout-info h4 {
	    color: #5bc0de;
	}
	</style>
	};
	
	my $hashRes;
	$hashRes->{is_trio} = $out2;
	$hashRes->{html_table} = $out2;
	print $cgi->header('text/json-comment-filtered');
	my $json_encode = encode_json $hashRes;
	print $json_encode;
	exit(0);
	
	#print $out;
	#exit(0);
}

sub get_header_table {
	my $header = "<thead>";
	$header .= $cgi->start_Tr({style=>"background-color:#E9DEFF;$fsize"});
	foreach my $h (@headers_validations){
		if (lc($h) eq 'scaled_score') { $header .= $cgi->th('Diag_score'); }
 		else { $header .=  $cgi->th(ucfirst($h)); }
	}
	$header .= $cgi->end_Tr();
	$header .= "</thead>";
	return $header;
}