#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;

use strict;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../cache_nodb/scripts/";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use GBuffer;
use export_data;
use JSON;
use polyweb_dude;
use VcfMerge;
use GenBoNoSql;
use Set::IntervalTree;
use Spreadsheet::WriteExcel;
use Bit::Vector;
use Bit::Vector::Overload;
use Compress::Snappy;
use Storable qw(store retrieve freeze dclone thaw);
use POSIX qw(strftime);
use List::MoreUtils qw(natatime);
use CGI::Session;
use html; 
use Carp;
use Cache_Commons;
use QueryVectorFilter;
use IO::Handle;
use xls_export;

require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";

my $io = IO::Handle->new();
$io->autoflush(1);

my $fsize = "font-size:10px";
my $value_red = 14;
my $value_coral = 12;
my $value_orange = 8;
my $value_yellow = 5;
my $color_red = '#CE0000';
my $color_coral = 'coral';
my $color_orange = '#EFB73E';
my $color_yellow = 'yellow';

my ($class, $class_default, $class_tr);
$class->{rowspan} -= 1;
$class->{rowspan} = 1 if $class->{rowspan} <=0;
$class->{style} = "min-width:10%;padding:1px";
$class_default->{style} = "max-width:350px;overflow-x:auto;vertical-align:middle;padding:5px;";
			


my $cgi = new CGI();
my $hgmd_version = $cgi->param('hgmd_version');
my $gene_name = $cgi->param('gene_name');
my $only_new = $cgi->param('only_new');
my $only_dm = $cgi->param('only_dm');
my $max_dejavu = $cgi->param('dejavu');
my $max_dejavu_ho = $cgi->param('dejavu_ho');
my $max_gnomad = $cgi->param('gnomad');
my $max_gnomad_ho = $cgi->param('gnomad_ho');
my $filters_cons = $cgi->param('filters_cons');

my $h_filters_cons;
if ($filters_cons) {
	foreach my $cons (split(',', $filters_cons)) {
		$h_filters_cons->{lc($cons)} = undef;
	}
}
my @headers_validations = ("var_name","locus","gnomad","deja_vu","table_validation","table_transcript");
my @header_transcripts = ("consequence","enst","nm","ccds","appris","exon","nomenclature","codons","codons_AA", "polyphen","sift","ncboost","cadd","revel","dbscsnv");
my $buffer = new GBuffer;

my $project_name = $buffer->get_random_project_name_with_this_annotations_and_genecode();
my $project = $buffer->newProject( -name => $project_name);
my $gene = $project->newGene($gene_name);

my @list_lines;
my $nb_var = 0;
my $h_v_hgmd_ids = $buffer->queryHgmd->search_variant_for_gene($gene_name);
my $h_new_dm = $buffer->queryHgmd->get_hash_last_released_DM();
my $h_no_coord = $buffer->queryHgmd->search_variant_for_gene_without_coord($gene_name);


my $h_by_coord;
my $hVarErrors;
my @lVar;
foreach my $v_hgmd_id (keys %{$h_v_hgmd_ids}) {
	next if (exists $h_no_coord->{$v_hgmd_id});
	next if ($only_new and not exists $h_new_dm->{$v_hgmd_id});
	
	$nb_var++;
	my $tag = $h_v_hgmd_ids->{$v_hgmd_id}->{tag};
	next if ($only_dm and $tag ne 'DM');
	
	my $chromosome = $h_v_hgmd_ids->{$v_hgmd_id}->{chromosome};
	my $start = $h_v_hgmd_ids->{$v_hgmd_id}->{coordSTART};
	my $end = $h_v_hgmd_ids->{$v_hgmd_id}->{coordEND};
	my $locus = 'chr'.$chromosome.':'.$start.'-'.$end;
	my $hgvs = $h_v_hgmd_ids->{$v_hgmd_id}->{hgvs};
	my $strand = $h_v_hgmd_ids->{$v_hgmd_id}->{strand};
	my $rs_name = $h_v_hgmd_ids->{$v_hgmd_id}->{dbsnp};
	my $ref_allele = $h_v_hgmd_ids->{$v_hgmd_id}->{'ref'};
	my $mut_allele = $h_v_hgmd_ids->{$v_hgmd_id}->{'alt'};
#	my ($ref_allele, $mut_allele) = get_ref_mut_alleles_from_hgvs($hgvs);
	
	my $polyweb_id = $chromosome.'_'.$start.'_'.$ref_allele.'_'.$mut_allele;
	my $gnomad_id = $chromosome.'-'.$start.'-'.$ref_allele.'-'.$mut_allele;
	$gnomad_id =~ s/chr//;

	my $var_id = $polyweb_id;
	my $v = $project->_newVariant($var_id);
	
	my $not_ok;
	my $var_gnomad = $v->getGnomadAC();
	$not_ok++ if ($var_gnomad and $max_gnomad and $var_gnomad > $max_gnomad);
	next if $not_ok;
	
	my $var_gnomad_ho = $v->getGnomadHO();
	$not_ok++ if ($var_gnomad_ho and $max_gnomad_ho and $var_gnomad_ho > $max_gnomad_ho);
	next if $not_ok;
	
	my $is_ok_annot;
	eval { $v->annotation(); };
	if ($@) {
		$hVarErrors->{$v->id()} = undef;
		next;
	}
	unless (exists $v->annotation->{$gene->id()}) {
		$hVarErrors->{$v->id()} = undef;
		next;
	}
	my $var_annot = $v->variationTypeInterface($gene);
	foreach my $this_annot (split(',', $var_annot)) {
		$this_annot =~ s/ /_/g;
		$is_ok_annot ++ if (exists $h_filters_cons->{lc($this_annot)});
	}
	next unless ($is_ok_annot);
	
	my $not_ok;
	my $var_dejavu = $v->nb_dejavu();
	$not_ok++ if ($max_dejavu and $var_dejavu > $max_dejavu);
	next if $not_ok;
	
	my $var_dejavu_ho = $v->nb_dejavu();
	$not_ok++ if ($max_dejavu_ho and $var_dejavu_ho > $max_dejavu_ho);
	next if $not_ok;
	
	
	push(@lVar, $v);

	$h_v_hgmd_ids->{$v_hgmd_id}->{class} = $tag;
	$h_v_hgmd_ids->{$v_hgmd_id}->{fields}->{class} = $tag;
	$h_v_hgmd_ids->{$v_hgmd_id}->{hgmd_id} = $v_hgmd_id;
	$h_v_hgmd_ids->{$v_hgmd_id}->{fields}->{hgmd_id} = $v_hgmd_id;
	$v->{hgmd} = $h_v_hgmd_ids->{$v_hgmd_id};

	my $hvariation;
	$hvariation->{value}->{id} =  $v->id;
	$hvariation->{html}->{id} =  $v->id;
	$hvariation->{value}->{type} = $v->type;
	$hvariation->{html}->{type} = $v->type;
	$hvariation->{value}->{locus} = $locus;
	$hvariation->{html}->{locus} = qq{<span style="font-size:8px;">$locus</span>};
	$hvariation->{value}->{gnomad_id} = $gnomad_id;
	$hvariation->{html}->{gnomad_id} = $gnomad_id;	
	$hvariation->{value}->{is_cnv} = 0;	
	update_variant_editor::vgnomad($v,$hvariation);
	update_variant_editor::vname($v,$hvariation);
	update_variant_editor::vspliceAI($v,$hvariation);
	update_variant_editor::vdivers($v,$hvariation);
	update_variant_editor::vclinvar($v,$hvariation);
	update_variant_editor::vhgmd($v,$hvariation);
	update_variant_editor::vdejavu($v,$hvariation);
	eval {
		$hvariation->{genes}->{$gene->id} = update_variant_editor::construct_hash_transcript($v, $cgi, \@header_transcripts, 2, $gene);
		$hvariation->{html}->{table_transcript} = update_variant_editor::table_transcripts($hvariation->{genes}->{$gene->id}, \@header_transcripts, 1);
	};
	if ($@) {
		$hvariation->{html}->{table_transcript} = "Problem...";
		$hvariation->{html}->{table_transcript} .= ' -> HGMD DB infos:  '.$h_v_hgmd_ids->{$v_hgmd_id}->{mutype};
		
	}
	if (exists $h_new_dm->{$v_hgmd_id}) {
		$hvariation->{html}->{hgmd} .= qq{ <b><i><font color='red'>New!</font></b></i>};
	}
	my $disease_hgmd = $v->hgmd_details->{disease};
	if ($disease_hgmd) {
		$hvariation->{html}->{hgmd} .= qq{<br><i><font color='green'>$disease_hgmd</font></i>};
	}
	update_variant_editor::table_validation_without_local($project, $hvariation, $gene);
	
	my $nb_other_project = $hvariation->{value}->{other_project};
	
	my $out;
	if ($nb_other_project and $nb_other_project > 0) { $out = $cgi->start_Tr(); }
	else { $out = $cgi->start_Tr({class=>"tr_hgmd_not_found"}); }
	foreach my $h (@headers_validations){
		if ($h eq "trio" or "table_transcript"){
			$class->{style} = "min-width:200px;max-width:450px;max-height:200px;overflow-x:auto;vertical-align:middle;padding:5px;white-space: nowrap;";
		}
		elsif ($h eq 'var_name') {
			$class->{style} = "max-width:130px;overflow-x:auto;vertical-align:middle;padding:5px;";
		}
		else {
			$class->{style} = "vertical-align:middle;padding:5px;";
		}
		
		my $this_html = $hvariation->{html}->{$h};
		
		if ($h eq 'table_transcript') {
			$out .= qq{<div style="max-height:200px;overflow-y:auto;max-width:220px;overflow-x:auto;">};
			$out .= $cgi->td($class, $this_html);
			$out .= qq{</div>};
		}
		elsif ($h eq '#') {
			$out.= $cgi->td($class_default, $nb_var);
		}
		elsif ($h eq 'var_name') {
			$out .= $cgi->td($class_default, "<center>".$this_html."</center>");
		}
		elsif ($h eq 'varsome') {
			$out .= $cgi->td($class_default, "<center>".$this_html."</center>");
		}
		elsif ($h eq 'locus') {
			$out .= $cgi->td($class_default, "<center>".$this_html."</center>");
		}
		else {
			$out.= $cgi->td($class_default,$this_html);
		}
	}
	$out .= $cgi->end_Tr();
	$out =~ s/zoomHgmd/zoomHgmdWithoutCss/;
	$h_by_coord->{$v->start()} = $out;
}

foreach my $start (sort {$a <=> $b} keys %$h_by_coord) {
	my $out = $h_by_coord->{$start};
	push(@list_lines, $out);
}

my ($hResGene, $h_panels_found);
my $gene_id = $gene->id();
my $external_name = $gene->external_name();
$hResGene->{$gene_id}->{id} = $gene_id;
$hResGene->{$gene_id}->{external_name} = $external_name;
$hResGene->{$gene_id}->{pLI} = $gene->pLI();
$hResGene->{$gene_id}->{omim_id} = $gene->omim_id();
$hResGene->{$gene_id}->{omim_inheritance} = $gene->omim_inheritance();
$hResGene->{$gene_id}->{variants} = \@lVar;
my ($pheno,$nb_other_terms) = $gene->polyviewer_phentotypes();
$hResGene->{$gene_id}->{phenotypes}->{pheno} = $pheno;
$hResGene->{$gene_id}->{phenotypes}->{nb_other_terms} = $nb_other_terms;
$hResGene->{$gene_id}->{specific_cmd} = '';
$hResGene->{$gene_id}->{max_score} = $gene->score();
my $description_gene = $gene->description();
$class = 'hgmd_gene';
eval {
	foreach my $panel (@{$gene->getPanels()}) {
		my $name = lc($panel->name());
		$name =~ s/ /_/g;
		unless (exists $hResGene->{$gene_id}->{panels}->{$name}) {
			$hResGene->{$gene_id}->{panels}->{$name}->{phenotype} = $panel->getPhenotypes()->[0]->name();
			$h_panels_found->{$name}++;
			my $value = 'panel_'.$name;
			$class .= ' '.$value;
		}
	}
};
if ($@) { $hResGene->{$gene_id}->{panels} = undef; }
my $panel_id = 'p_'.$gene_id.'_hgmd';
$hResGene->{$gene_id}->{uid} = $panel_id;
my $html_gene = update_variant_editor::panel_gene($hResGene->{$gene_id},$panel_id);
$html_gene =~ s/float:right;/float:right;display:none;/;
$html_gene =~ s/glyphicon-triangle-right//;
my $b_hgmd_hide = qq{</span></span></div><div style="float:right;" class="form-check"><input class="form-check-input" type="checkbox" value="" onClick="show_hide_hgmd_not_in_polyweb();" id="b_found_hgmd_in_polyweb"><label style="color:white;padding:5px;font-size:12px;" class="form-check-label" for="b_found_hgmd_in_polyweb"> only variant(s) found in PolyWeb</label></div>};
$html_gene =~ s/<\/span><\/span><\/div>/$b_hgmd_hide/;
#warn Dumper $html_gene; die;


### EXPORT HTML


my $out2 =  $cgi->start_div();
$out2 .= qq{<table data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-v-align="both" data-pagination-pre-text="Previous" data-pagination-next-text="Next"data-pagination="true" data-page-size="100" data-page-list="[25, 50, 100, 200, 300]" data-resizable='true' id='table_variants' class='table table-striped' style='font-size:13px;'>};
$out2 .= "<thead>";
$out2 .= $cgi->start_Tr({style=>"background-color:#E9DEFF;$fsize"});
foreach my $h (@headers_validations) {
	my $cat = ucfirst($h);
	if (lc($h) eq '#' or lc($h) eq 'gnomad' or lc($h) eq 'deja_vu') {
		$out2 .= qq{<th data-field="$h">$cat</th>};
	}
	else {
		if (lc($cat) eq '#' or lc($cat) eq 'varsome' or lc($cat) eq 'gnomad' or lc($cat) eq 'deja_vu') {
			$out2 .= qq{<th data-field="$h">$cat</th>};
		}
		elsif (lc($cat) eq 'table_validation') {
			$out2 .= qq{<th data-field="$h" data-filter-control="input" data-filter-control-placeholder="DM / Pathogenic">HGMD / Clinvar / Local</th>};
		}
		elsif (lc($cat) eq 'table_transcript') {
			$out2 .= qq{<th data-field="$h" data-filter-control="input" data-filter-control-placeholder="Stop / c.73A>C / exon2 / ENST00000145855">Annotations</th>};
		}
		elsif (lc($cat) eq 'var_name') {
			$out2 .= qq{<th data-field="$h" data-filter-control="input" data-sortable="true" data-filter-control-placeholder="7-15456-A-G">Variant ID</th>};
		}
		elsif (lc($cat) eq 'locus') {
			$out2 .= qq{<th data-field="$h" data-filter-control="input" data-sortable="true" data-filter-control-placeholder="Position">Locus</th>};
		}
		else {
			$out2 .= qq{<th data-field="$h" data-filter-control="input">$cat</th>};
		}
	}
}
$out2 .= $cgi->end_Tr();
$out2 .= "</thead>";
$out2 .= "<tbody>";
foreach my $line (@list_lines) {
	$out2 .= $line;
}

my $hgmd_version_used = $buffer->queryHgmd->database();
my $nb_no_coord = scalar keys %{$h_no_coord};
if ($nb_no_coord > 0) {
	my @l_infos;
	foreach my $v_hgmd_id (sort keys %$h_no_coord) {
		my $h_infos = $buffer->queryHgmd->getDataHGMDPro($v_hgmd_id);
		my $line = "<center>HGMD ID: <b>$v_hgmd_id</b>";
		if (exists $h_new_dm->{$v_hgmd_id}) {
			$line .= " - class found: <b>".$h_infos->{'tag'}.' <i><u>New!</i></u></b>';
		}
		else {
			$line .= " - class found: <b>".$h_infos->{'tag'}.'</b>';
		}
		$line .= " - HGMD annotation: <b>".$h_infos->{'mutype'}.'</b></center>';
		push(@l_infos, $line);
	}
	my $infos_list = "<b><u><center>Variants without coordinates in HGMD ($hgmd_version_used) for gene $gene_name</center></b></u><br>";
	$infos_list .= join('<br>', @l_infos);
	$out2 .= $cgi->start_Tr({style=>"background-color:#f5a260;font-size:14px;"});
	my $nb_rowspan = scalar(@headers_validations);
	$out2 .= $cgi->td({style=>"vertical-align:middle;padding:1px;text-align: center;vertical-align:middle", colspan=>$nb_rowspan},qq{<b>+ <a class="btn btn-primary btn-xs" role="button" onclick="document.getElementById('span_list_panels').innerHTML='$infos_list';dijit.byId('dialog_list_panels').show();" style="min-width:40px;background-color:white;text-shadow:0.5px 0.5px 0.5px green;color:red;font-size:9px;">$nb_no_coord</a> variants found without coordinates in HGMD Database...</b>} );
	$out2 .= $cgi->end_Tr();
}

$out2 .= "</tbody>";
$out2 .= "</table>";
$out2 .= "</div>";

my $hRes;
$hRes->{html_variants} = $out2;
$hRes->{html_gene} = $html_gene;
$hRes->{html_page_title} = 'HGMD '.$gene->external_name();
$hRes->{html_source} = "HGMD DataBase $hgmd_version_used";
$hRes->{html_title} = "Gene ".$gene->external_name();

my $session_id = save_export_xls($gene, \@lVar);
save_html($session_id, $hRes);

my $hRes2;
$hRes2->{session_id} = $session_id;
my $json_encode = encode_json $hRes2;
print $cgi->header('text/json-comment-filtered');
print $json_encode;
exit(0);


sub save_export_xls {
	my ($gene, $list_var) = @_;
	my $buffer = new GBuffer;
	my $proj_name = $buffer->get_random_project_name_with_this_annotations_and_genecode();
	my $project = $buffer->newProject( -name => $proj_name );
	$project->cgi_object(1);
	my @lVarObj;
	my $xls_export = new xls_export();
	$xls_export->title_page('HGMD_'.$gene->external_name().'.xls');
	$xls_export->store_variants_infos(\@$list_var, $project);
	
	my ($h_patients, $h_row_span);
	foreach my $chr_id (keys %{$xls_export->{hash_variants_global}}) {
		foreach my $var_id (keys %{$xls_export->{hash_variants_global}->{$chr_id}}) {
			my $project_name = 'HGMD';
			my $pat_name = 'HGMD';
			$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'variation'} = $var_id;
			$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'project'} = $project_name;
			$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'description'} = '';
			$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'phenotypes'} = '';
			$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'fam'} = '';
			$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'name'} = '';
			$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'sex'} = '';
			$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'status'} = '';
			$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'parent_child'} = '';
			$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'sex_status_icon'} = '';
			$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'perc'} = '';
			$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'model'} = '';
		}
	}
	
	$xls_export->store_specific_infos('projects_patients_infos', $h_patients);
	my $session_id = $xls_export->save();
	return $session_id;
}

sub save_html {
	my ($session_id, $hRes) = @_;
	my $session = new session_export();
	$session->load_session( $session_id );
	$session->save( 'html_page_title', $hRes->{'html_page_title'} );
	$session->save( 'html_source', $hRes->{'html_source'} );
	$session->save( 'html_title', $hRes->{'html_title'} );
	$session->save( 'html_gene', $hRes->{'html_gene'} );
	$session->save( 'html_variants', $hRes->{'html_variants'} );
}

sub get_ref_mut_alleles_from_hgvs {
	my ($hgvs, $strand) = shift;
	$hgvs =~ s/[-+]+//g;
	$hgvs =~ s/[0-9]+//g;
	my ($ref, $mut) = split('>', $hgvs);
	return ($ref, $mut) if $strand eq '+';
	my $ref2 = convert_allele_if_strand_negatif($ref);
	my $mut2 = convert_allele_if_strand_negatif($mut);
	return ($ref2, $mut2);
}

sub convert_allele_if_strand_negatif {
	my $allele = shift;
	return 'A' if $allele eq 'T';
	return 'T' if $allele eq 'A';
	return 'G' if $allele eq 'C';
	return 'C' if $allele eq 'G';
	confess("\n\nERROR: pb in method convert_allele_if_strand_negatif with allele $allele. Die.\n\n");
}




