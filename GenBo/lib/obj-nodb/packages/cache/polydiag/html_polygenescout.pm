$| = 1;
use CGI qw/:standard :html3/;

package html_polygenescout;
use strict;
use FindBin qw($Bin);
use Data::Dumper;
use lib "$Bin/";
use lib "$Bin/../../../";
use GBuffer;
use GenBoProject;
use GenBoProjectCache;


my $cgi = new CGI();



sub print_hgmd_total_new {
	my ($total_new, $external_name) = @_;
	return qq{<td><a class="btn btn-xs" onClick="get_variants_from_hgmd_gene_filters('$external_name', '1');" style="min-width:30px"><button style="background-color:green;color:white;font-size:12px;padding-left:5px;padding-right:5px;">$total_new</button></a></td>} if ($total_new);
	return $cgi->td({rowspan => 1,style => "max-height:60px;overflow-y: auto;width:90px;"}, '<i class="fa fa-minus"></i>');
}

sub print_hgmd_total_mut {
	my ($total_mut, $external_name) = @_;
	return qq{<td><a class="btn btn-xs" onClick="get_variants_from_hgmd_gene_filters('$external_name');" style="min-width:30px"><button style="background-color:green;color:white;font-size:12px;padding-left:5px;padding-right:5px;">$total_mut</button></a></td>} if ($total_mut);
	return $cgi->td({rowspan => 1,style => "max-height:60px;overflow-y: auto;width:90px;"}, '<i class="fa fa-minus"></i>');
}

sub print_hgmd_concepts {
	my ($nb_concepts, $concept_list, $used_concept) = @_;
	if ($nb_concepts) {
		if ($used_concept and $concept_list =~ /$used_concept/) {
			$concept_list =~ s/$used_concept/$used_concept <b><i>Found !<\/b><\/i>/;
			return qq{<td><a class="btn btn-xs" onClick="document.getElementById('span_list_panels').innerHTML='$concept_list';dijit.byId('dialog_list_panels').show();" style="min-width:30px"><button style="background-color:red;color:white;font-size:12px;padding-left:5px;padding-right:5px;">$nb_concepts</button></a></td>}
		}
		else {
			return qq{<td><a class="btn btn-xs" onClick="document.getElementById('span_list_panels').innerHTML='$concept_list';dijit.byId('dialog_list_panels').show();" style="min-width:30px"><button style="background-color:green;color:white;font-size:12px;padding-left:5px;padding-right:5px;">$nb_concepts</button></a></td>}
		}
	}
	return $cgi->td({rowspan => 1,style => "max-height:60px;overflow-y: auto;width:90px;"}, '<i class="fa fa-minus"></i>');
}

sub print_tab_variants_project {
	my ($chr_id, $external_name, $trio_project, $trio_patient, $h_genes,$project_trio) = @_;
	my $buffer_trio;
	if ($project_trio){
		#$buffer_trio = $project_trio;
		
	}
	else {
		my $buffer_trio = GBuffer->new();
		$project_trio = $buffer_trio->newProjectCache( -name => $trio_project );
	}
	# $buffer_trio = GBuffer->new() unless $buffer_trio;
	# die();
	#my $project_trio = $buffer_trio->newProjectCache( -name => $trio_project );
	my $patient_trio = $project_trio->getPatient($trio_patient);
	
	my $chr_trio = $project_trio->getChromosome($chr_id);
	my $gene_trio = $project_trio->newGene($h_genes->{$external_name});
	my $vector_gene = $gene_trio->getVariantsVector();
	my $vector_patient = $patient_trio->getVectorOrigin($chr_trio)->Clone();
	$vector_patient->Intersection($vector_patient, $vector_gene);
	my @lVar = @{$chr_trio->getListVarObjects($vector_patient)};
	my $total = 0;
	my $denovo = 0;
	my $dominant = 0;
	my $recessif = 0;
	my $mosaic = 0;
	my $uniparental = 0;
	my $mother = 0;
	my $father = 0;
	my $both = 0;
	
	foreach my $v (@lVar) {
		$total++;
		my $v_model = lc($v->getTransmissionModelType($patient_trio->getFamily(), $patient_trio));
		if ($v_model =~ /denovo/) { $denovo++; }
		if ($v_model eq 'dominant') { $dominant++; }
		if ($v_model =~ /recessi/) { $recessif++; }
		if ($v_model eq 'mosaic') { $mosaic++; }
		if ($v_model =~ /uniparental/ ) { $uniparental++; }
		if ($v_model =~ /mother/) { $mother++; }
		if ($v_model =~ /father/) { $father++; }
		if ($v_model eq 'both') { $both++; }
	}
	
	my $color_father = '#779ecb';
	my $color_mother = '#f7c9c9';
	my $color_denovo = '#e74c3c';
	my $color_dominant = '#e74c3c';
	my $color_recessif = '#ee82ee';
	my $color_mosaic = '#f9885c';
	my $color_uniparental = '#45b8ac';
	my $color_both = 'grey';
	
	my (@l_td_header, @l_td_values);
	foreach my $cat ('Denovo', 'Dominant', 'Recessive', 'Mosaic', 'Uniparental', 'Mother', 'Father', 'Both') {
		my $use_patient = $trio_patient;
		if (lc($cat) eq 'father') {
			next unless $patient_trio->getFamily->getFather();
			$use_patient = $patient_trio->getFamily->getFather->name();
		}
		if (lc($cat) eq 'mother') {
			next unless $patient_trio->getFamily->getMother();
			$use_patient = $patient_trio->getFamily->getMother->name();
		}
		
		my ($value, $badge_color);
		$badge_color = 'grey';
		if (lc($cat) eq 'total') { $value = $total; }
		if (lc($cat) eq 'denovo') {
			$value = $denovo;
			$badge_color = $color_denovo;
		}
		if (lc($cat) eq 'dominant') {
			$value = $dominant;
			$badge_color = $color_dominant;
		}
		if (lc($cat) eq 'recessive') {
			$value = $recessif;
			$badge_color = $color_recessif;
		}
		if (lc($cat) eq 'mosaic') {
			$value = $mosaic;
			$badge_color = $color_mosaic;
		}
		if (lc($cat) eq 'uniparental') {
			$value = $uniparental;
			$badge_color = $color_uniparental;
		}
		if (lc($cat) eq 'mother') {
			$value = $mother;
			$badge_color = $color_mother;
		}
		if (lc($cat) eq 'father') {
			$value = $father;
			$badge_color = $color_father;
		}
		if (lc($cat) eq 'both') { $value = $both; }
		
		next if (lc($cat) ne 'total' and $value == 0);
		
		my ($td_header, $td_value);
		if ($value > 0) {
			my $only_model = lc($cat);
			my $ensg = $h_genes->{$external_name};
			my $cmd_var_base = qq{view_var_from_proj_gene_pat('$trio_project', '$ensg', '$use_patient', '', '$only_model')};
			$td_header = qq{<td style="color:black;"><b>$cat</b></td>};
			$td_value = qq{<td onClick="$cmd_var_base"><span class="badge badge-success badge-xs" style="border-color:$badge_color;background-color:$badge_color;color:white;font-size:9px;">$value</span></td>};
			
			
			
		}
		else {
			$td_header = qq{<td style="color:black;opacity:0.5;"><b>$cat</b></td>};
			$td_value = qq{<td style="color:black;opacity:0.5;">$value</td>};
		}
		push(@l_td_header, $td_header);
		push(@l_td_values, $td_value);
	}
	
	my $color_text = 'black';
#	$color_text = 'red' if ($nb_dm_text > 0);
	my $nb_var_text = $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"text-align:center;max-width:300px;box-shadow: 1px 1px 6px $color_text;font-size: 8px;font-family:  Verdana;margin-bottom:0px"});
	$nb_var_text .= $cgi->start_Tr();
	foreach my $td (@l_td_header) { $nb_var_text .= $td; }
	$nb_var_text .= $cgi->end_Tr();
	
	$nb_var_text .= $cgi->start_Tr();
	foreach my $td (@l_td_values) { $nb_var_text .= $td; }
	
	$nb_var_text .= $cgi->end_Tr();
	$nb_var_text .= $cgi->end_table();
	return qq{<td><a class="btn btn-primary btn-xs" style="min-width:40px;border-color:transparent;background-color:transparent;color:black;"><center><u>$nb_var_text</u></center></a></td>};
}

sub print_link_dejavu {
	my ($cmd_link) = @_;
	return qq{<td><a class="btn btn-primary btn-xs" onClick="$cmd_link" style="min-width:40px;border-color:transparent;background-color:transparent;color:black;"><center><button>DejaVu</button></center></a></td>};
}

sub print_link_dejavu_with_locus {
	my ($cmd_link, $hash_locus_intervals) = @_;
	my $html = qq{<center><td>};
	$html .= qq{<table style="width:100%;">};
	$html .= qq{<tr><td>};
	$html .= qq{<a class="btn btn-primary btn-xs" onClick="$cmd_link" style="min-width:40px;border-color:transparent;background-color:transparent;color:black;"><center><button>DejaVu</button></center></a>};
	$html .= qq{</td></tr>};
	foreach my $locus (sort keys %$hash_locus_intervals) {
		my $cmd_link_locus = $hash_locus_intervals->{$locus};
		$html .= qq{<tr><td>};
		$html .= qq{<a class="btn btn-primary btn-xs" onClick="$cmd_link_locus" style="min-width:40px;border-color:transparent;background-color:transparent;color:black;"><center><u>Only $locus</u></center></a>};
		$html .= qq{</td></tr>};
	}
	$html .= qq{</table>};
	$html .= qq{</td></center>};
	return $html;
}

sub get_project_phenotype {
	my ($trio_project) = @_;
	my $buffer_trio = GBuffer->new();
	my $project_trio = $buffer_trio->newProject( -name => $trio_project );
	return $project_trio->phenotypes->[0];
}

sub print_gene_score {
	my ($gene_score) = @_;
	my $bcolor = "grey";
	$bcolor = "yellow" if $gene_score >= 1;
	$bcolor = "orange" if $gene_score >= 2;
	$bcolor = "coral" if $gene_score >= 3;
	$bcolor = "red" if $gene_score >= 4;
	my $color = '#FFFFFF';
	$color = 'grey' if $bcolor eq 'yellow';
 	return qq{<td style="border-color:$bcolor;background-color:$bcolor;border-color:black;"><span class="badge badge-success badge-xs" style="background-color:$bcolor;color:$color;font-size:25px;">$gene_score</span></td>};
}

sub print_locus {
	my ($gene) = @_;
	my $chr_locus = $gene->getChromosome->id();
	my $start_locus = $gene->start();
	my $end_locus = $gene->end();
	my $gene_locus = $chr_locus.':'.$start_locus.'-'.$end_locus;
	return $cgi->td({ rowspan => 1, style => "max-height:60px;overflow-y:auto;font-size:12px;" }, $gene_locus );
}

sub print_gene_basic_tables_infos {
	my ($gene_id, $use_phenotype_score) = @_;
	return print_table_base_line_gene($gene_id, $use_phenotype_score);
}

sub print_table_base_line_gene {
	my ($hgene, $used_phenotype) = @_;
	my $project_name;
	my $cgi = new CGI();
	my $out;
	my $gene_id = $hgene->{id};
	my $max_score = $hgene->{max_score};
	
	my $glyph = "";
	$glyph = qq{<span class="glyphicon glyphicon-star-empty text-default" aria-hidden="true"></span>} if ($hgene->{nb_clinvar} and $hgene->{nb_clinvar} > 0);
	$glyph = qq{<span class="glyphicon  glyphicon-alert text-alert" aria-hidden="true" style="color:red"></span>} if (($hgene->{nb_clinvar_alert} and $hgene->{nb_clinvar_alert} > 0) or $hgene->{dm});
	
	my $gene_name = $hgene->{external_name};
	my $in = $hgene->{omim_inheritance};
	$in ="" if $in eq "-";
	$in = "X-linked " if $in =~/X-linked/;
	
	my $pli = 0;
	$pli = $hgene->{pLI}*1.0 if ($hgene->{pLI} and not $hgene->{pLI} eq '-');
				
	my $bcolor = "grey";
	$bcolor = "green" if $max_score >= 0;
	$bcolor = "yellow" if $max_score >= 5;
	$bcolor = "orange" if $max_score >= 8;
	$bcolor = "coral" if $max_score >= 12;
	$bcolor = "red" if $max_score >= 14;
			
	my $cnv_status = "cnv_none";
	if ($hgene->{level_dude} and $hgene->{level_dude} ne '-1') {
		my ($l,$t) = split(":",$hgene->{level_dude});
		if ($t eq "del"){
			$cnv_status = "cnv_del_medium" if $l eq "medium";
			$cnv_status = "cnv_del_high" if $l eq "high";
			$cnv_status = "cnv_del_low" if $l eq "low";
		}
		if ($t eq "dup"){ 
			$cnv_status = "cnv_dup_medium" if $l eq "medium";
			$cnv_status = "cnv_dup_high" if $l eq "high";
			$cnv_status = "cnv_dup_low" if $l eq "low";
		}
	}
	my $this_b_cmd;
	if (exists $hgene->{specific_cmd}) {
		$this_b_cmd = $hgene->{specific_cmd};
	}
	if ($this_b_cmd) {
		$out .= qq{<td><button class="btn btn-primary btn-xs" style="border-color:transparent;background-color:transparent;color:black;font-size:20px;" onClick='$this_b_cmd'><u>$gene_name</u></button></td>};
	}
	else {
		$out .= qq{<td><a class="btn btn-primary btn-xs" style="border-color:transparent;background-color:transparent;color:black;font-size:20px;"  href="https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=$gene_name" target="_blank"><u>$gene_name</u></a></td>};
	}
	
	my $pheno = $hgene->{phenotypes}->{pheno};
	my $nb_other_terms = $hgene->{phenotypes}->{nb_other_terms};
	my $color ;
	$color = qq{ style = "color:#E74C3C"} if $pheno =~/intellectual/ or $pheno =~/mental/ or $pheno =~/retar/;
	my $jtid = 'zz'.time.rand(500);
   	my $div_pheno = qq{<a aria-disabled="true" class="btn btn-primary btn-xs" href="#" role="button" style="text-align:left;font-family: proxima-nova, sans-serif;font-style:normal;border-color:transparent;background-color:transparent;color:black;">};
   	if ($pheno) {
		if (length($pheno) > 70) { $pheno = substr($pheno, 0, 70).'...'; }
   		if ($nb_other_terms > 0) {
	   		if ($project_name) {
	   			$pheno .= qq{</u><span style='color:blue' onclick="update_grid_gene_phenotypes(\'$gene_id\', \'$project_name\')";">+ $nb_other_terms terms</span>};
	   		}
	   		else {
	   			$pheno .= qq{</u><span style='color:blue' onclick="update_grid_gene_phenotypes(\'$gene_id\')";">+ $nb_other_terms terms</span>};
	   		}
   		}
   		if ($project_name) {
   			$div_pheno .= qq{<span onclick="update_grid_gene_phenotypes(\'$gene_id\', \'$project_name\')";">$pheno</span>};
   		}
   		else {
   			$div_pheno .= qq{<span onclick="update_grid_gene_phenotypes(\'$gene_id\')";">$pheno</span>};
   		}
   	}
   	$div_pheno .= qq{</a>};
	$out .= '<td><u>'.$div_pheno.'</u></td>';

	my $nbv = $hgene->{nb};
	my $omim = $hgene->{omim_id};
	
	$out .= qq{<td><a class="btn btn-primary btn-xs" href="http://www.omim.org/entry/$omim" role="button" target="_blank" style="min-width:40px;border-color:transparent;background-color:transparent;color:black;"><u>Omim</u></a></td>} if ($omim and $omim ne "");
	$out .= qq{<td><i class="fa fa-minus"></i></td>} if (not $omim or $omim eq "");
				
	my ($gid,$t) = split("_",$hgene->{id});
	$out .=qq{<td><a class="btn btn-primary btn-xs" href="https://gtexportal.org/home/gene/$gid" role="button" target="_blank" style="min-width:40px;border-color:transparent;background-color:transparent;color:black;"><u>Gtex</u></a></td>};
	
	my $oid = $hgene->{external_name};
	my $type ="green";
	$type = "orange" if $pli >= 0.75;
	$type = "red" if $pli >= 0.9;
	my $m = $hgene->{max_score};
	$out .=qq{<td><center><a class="btn btn-xs" href="https://gnomad.broadinstitute.org/gene/$oid" target="_blank" style="min-width:30px"><button style="background-color:$type;color:white;font-size:11px;padding-left:5px;padding-right:5px;">$pli</button></a></center></td>};
		
	my $panel_name1 = join("-",keys %{$hgene->{panels}});
	my $hPanels;
	foreach my $panel_name (keys %{$hgene->{panels}}) {
		my $pheno_name = $hgene->{panels}->{$panel_name}->{phenotype};
		$hPanels->{$pheno_name}->{$panel_name} = undef;
	}
	my $found_phenotype;
	my $panel_list;
	foreach my $pheno_name (sort keys %$hPanels) {
		$found_phenotype = 1 if ($used_phenotype and $pheno_name eq $used_phenotype);
		$panel_list .= "<br><center><b><u>Phenotype: ".$pheno_name."</b></u></center><br>";
		foreach my $panel_name (sort keys %{$hPanels->{$pheno_name}}) {
			$panel_list .= '<center>'.$panel_name."</center>";
		}
	}
	$panel_list .= "<br>";
	$panel_name1 = scalar (keys %{$hgene->{panels}});
	if ($panel_name1 and $found_phenotype and $found_phenotype == 1) {
		$out .=qq{<td><center><a class="btn btn-primary btn-xs" href="#" role="button" style="border-color:transparent;background-color:transparent;color:black;font-size:11px;" onclick="document.getElementById('span_list_panels').innerHTML='$panel_list';dijit.byId('dialog_list_panels').show();"><button style="background-color:red;color:white;font-size:11px;padding-left:5px;padding-right:5px;">$panel_name1</button></a></center></td>};
	}
	elsif ($panel_name1) {
		$out .=qq{<td><center><a class="btn btn-primary btn-xs" href="#" role="button" style="border-color:transparent;background-color:transparent;color:black;font-size:11px;" onclick="document.getElementById('span_list_panels').innerHTML='$panel_list';dijit.byId('dialog_list_panels').show();"><button style="background-color:green;color:white;font-size:11px;padding-left:5px;padding-right:5px;">$panel_name1</button></a></center></td>};
	}
	else { $out .= qq{<td><center><i class="fa fa-minus"></i></center></td>}; }

#	$out .= '</tr></table>';
	return $out;	   		
}

sub print_alamut_variant_button {
	my ($alamut_id) = @_;
	return qq{	<button class="alamutView3" onClick ="displayInAlamutVector('$alamut_id');"></button>};
}


1;