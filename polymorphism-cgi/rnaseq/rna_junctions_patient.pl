#!/usr/bin/perl
$|=1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo";
use lib "$Bin/../../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";

use connect;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use JSON;
use xls_export;
use session_export;

my $cgi    = new CGI;
my $project_name = $cgi->param('project');
my $use_patient = $cgi->param('patient');
my $view_all_junctions = $cgi->param('view_not_filtred_junctions');
my $isChrome = $cgi->param('is_chrome');
my $max_dejavu = $cgi->param('dejavu');
my $use_percent_dejavu = $cgi->param('dejavu_percent');
my $min_score = $cgi->param('min_score');
my $only_gene_name = $cgi->param('only_gene');
my $only_positions = $cgi->param('only_positions');
my $only_dejavu_ratio_10 = $cgi->param('only_dejavu_ratio_10');


$min_score = 20 unless (defined($min_score));

my $nb_group_junctions_colors = 0;

my $buffer = GBuffer->new;
my $project = $buffer->newProjectCache( -name => $project_name );
$project->getChromosomes();
my $patient = $project->getPatient($use_patient);
$patient->use_not_filtred_junction_files(1);
my $only_gene;
$only_gene = $project->newGene($only_gene_name) if ($only_gene_name);

print $cgi->header('text/json-comment-filtered');
print "{\"progress\":\".";

my $patient_name = $patient->name();
my $table_id = 'table_'.$patient_name;
my $html;
my $release = $project->getVersion();
my $gencode = '-';
$gencode = $project->gencode_version() if ($release =~ /HG19/);

my ($default_filter_dejavu_project, $default_filter_score);
my ($nb_max_patient, $nb_limit);
my $h_desc = $project->get_hash_patients_description_rna_seq_junction_analyse();
if ($h_desc) {
	foreach my $pat_name (keys %$h_desc) { $nb_max_patient++ if (exists $h_desc->{$pat_name}->{pat}); }
}
else { $nb_max_patient = scalar(@{$project->getPatients()}); }
if ($nb_max_patient == 1) { $nb_limit = 0; }
elsif ($nb_max_patient <= 3) { $nb_limit = 1; }
else { $nb_limit = int($nb_max_patient / 4); }
$default_filter_dejavu_project = "data-filter-default='<=$nb_limit'" if $nb_limit > 0;
$default_filter_score = "data-filter-default='>=1'";

my $score_slider = 0;
$score_slider = $min_score / 10 if ($min_score and $min_score > 0);
my (@lJunctions, $h_varids, $h_var_linked_ids);
foreach my $chr (@{$project->getChromosomes()}) {
	#next if $chr->id ne '4';
	my $vector_patient = $patient->getJunctionsVector($chr);
	if ($only_positions) {
		my @lTmp = split(':', $only_positions);
		my $chr_id_limit = lc($lTmp[0]);
		$chr_id_limit =~ s/chr//;
		next if lc($chr->id) ne $chr_id_limit;
		my ($from, $to) = split('-', $lTmp[1]);
		my $vector_postions = $chr->getVectorByPosition($from, $to);
		$vector_patient->Intersection($vector_patient, $vector_postions);
	}
	foreach my $junction (@{$chr->getListVarObjects($vector_patient)}) {
		$h_varids->{$junction->id()} = undef;
		if ($junction->is_junctions_linked($patient)) {
			if (exists $h_var_linked_ids->{$junction->id()}) {
				$h_var_linked_ids->{$junction->id()}->{vector_id} = $chr->id().'-'.$junction->vector_id();
			}
			else {
				$h_var_linked_ids->{$junction->id()}->{vector_id} = $chr->id().'-'.$junction->vector_id();
				my $hdone;
				my @lOthers = (keys %{$junction->get_hash_junctions_linked_to_me->{$patient->name()}});
				my $nb_others = scalar(@lOthers);
				while ($nb_others > 0) {
					my $other_id = $lOthers[0];
					#if (not exists $hdone->{$other_id}) {
						foreach my $id (keys %{$junction->get_hash_junctions_linked_to_me->{$patient->name()}}) {
							$h_var_linked_ids->{$junction->id()}->{linked_to}->{$id} = undef;
							$h_var_linked_ids->{$id}->{linked_to}->{$junction->id()} = undef;
							push(@lOthers, $id) if not (exists $hdone->{$other_id});
						}
						$hdone->{$other_id} = undef;
					#}
					my $supress = shift(@lOthers);
					$nb_others = scalar(@lOthers);
				}
				
				my $hallids = $h_var_linked_ids->{$junction->id()}->{linked_to};
				$hallids->{$junction->id()} = undef;
				foreach my $jid (keys %$hallids) {
					next unless exists $h_var_linked_ids->{$jid}->{linked_to};
					$h_var_linked_ids->{$jid}->{linked_to} = $hallids;
				}
			}
			
		}
		push(@lJunctions, $junction);
	}
}

my $percent_dejavu = 0;
if (defined($use_percent_dejavu)) {
	$percent_dejavu = $use_percent_dejavu - 90;
}
my $nb_percent_dejavu_value = 90 + $percent_dejavu;

my $checked_only_dejavu_ratio_10;
$checked_only_dejavu_ratio_10 = qq{checked="checked"} if $only_dejavu_ratio_10;
my $max_dejavu_value = 51;
$max_dejavu_value = $max_dejavu if defined($max_dejavu);
my $html_dejavu = qq{
	<table style="width:100%;">
		<tr>
			<td style="padding-top:5px;">
				<center>
					<label for="slider_dejavu_percent" class="form-label" style="font-size:10px;font-weight:300;"><i>if <span id="nb_percent_dejavu" style="color:blue;">$nb_percent_dejavu_value</span> % similar coord</i></label>
					<input style="width:180px;" type="range" class="form-range" min="0" max="10" value="$percent_dejavu" step="1" id="slider_dejavu_percent" onchange="update_dejavu_percent_span_value()">
				</center>
			</td>
			<td style="padding-top:5px;">
				<center>
					<label for="slider_dejavu" class="form-label" style="font-size:10px;font-weight:300;"><i>in max <span id="nb_max_dejavu_patients" style="color:blue;">$max_dejavu_value</span> Patient(s)</i></label>
					<input style="width:180px;" type="range" class="form-range" min="0" max="51" value="$max_dejavu_value" step="1" id="slider_dejavu" onchange="update_dejavu_span_value()">
				</center>
			</td>
			<td style="padding-top:5px;">
				<center>
					<div class="form-check"><input $checked_only_dejavu_ratio_10 class="form-check-input" type="checkbox" value="" id="b_dejavu_min_ratio_10"><label class="form-check-label" for="b_dejavu_min_ratio_10" style="padding-left:10px;font-size:11px;">only ratio >10%</label></div>
				</center>		
			</td>
					
		</tr>
		<tr>
			<td style="padding-top:5px;font-size:11px;"><center><b>DejaVu Min Similarity</center></b></td>
			<td style="padding-top:5px;font-size:11px;"><center><b>DejaVu Max patients</center></b></td>
			<td style="padding-top:5px;font-size:11px;"><center><b>DejaVu Ratio</center></b></td>
		</tr>
	</table>};
	

my $html_gene_select = qq{<input type="text" style="font-size:11px;width:180px;" placeholder="COL4A5, 1:100000-1500000" class="form-control" id="input_gene_positions">};
if ($only_gene_name) {
	$html_gene_select = qq{<input type="text" style="font-size:11px;width:180px;" value="$only_gene_name" class="form-control" id="input_gene_positions">};
}
elsif ($only_positions) {
	$html_gene_select = qq{<input type="text" style="font-size:11px;width:180px;" value="$only_positions" class="form-control" id="input_gene_positions">};
}

	
my $html_filters = qq{
	<table style="width:100%;">
		<tr>
			<td style="padding-top:5px;">
				<center>
					<label for="slider_score" class="form-label" style="font-size:10px;font-weight:300;"><i>Ratio >= <span id="nb_score" style="color:blue;">$min_score%</span></i></label>
					<input style="width:150px;" type="range" class="form-range" min="0" max="10" value="$score_slider" step="1" id="slider_score" onchange="update_score_span_value()">
				</center>
			</td>
			<td style="padding-top:5px;">
				<center>
					$html_gene_select
				</center>
			</td>
					
		</tr>
		<tr>
			<td style="padding-top:5px;font-size:11px;"><center><b>Filter Min Ratio</center></b></td>
			<td style="padding-top:5px;font-size:11px;"><center><b>Filter Only Gene / Positions</center></b></td>
		</tr>
	</table>};

	
my $html_refresh = qq{
	<table style="width:100%;">
		<tr>
			<td style="padding-top:12px;">
				<center>
					<button style="font-size:25px;" type="button" class="btn btn-sm btn-primary" id="b_update_dejavu" onclick="launch_dejavu_span_value()"><span class="glyphicon glyphicon-refresh" aria-hidden="true"></span></button>
				</center>
			</td>
		<tr>
		</tr>
			<td style="padding-top:5px;">
				<center>
					<b><span style="color:#363945">REFRESH</span></b>
				</center>
			</td>
		</tr>
	</table>};
	
my $html_patient = qq{
	<table style="width:100%;">
		<tr>
			<td style="padding-top:8px;"><center><b>Name</b></center></td>
			<td style="padding-top:8px;"><center><span style='color:red;font-size:13px'>$patient_name</span></center></td>
		</tr>
		<tr>
			<td><center><b>Release</b></center></td>
			<td><center>$release</center></td>
		</tr>
		<tr>
			<td><center><b>Gencode</b></center></td>
			<td><center>$gencode</center></td>
		</tr>
	</table>};

$html .= qq{<br>};
$html .= qq{<div class="container" style="width:100%;height:86px;vertical-align:middle;"><div class="row">};
$html .= qq{<div class="col-sm-2"><div style="height:82px;border:3px #363945 double;">$html_patient</div></div>};
$html .= qq{<div class="col-sm-5"><div style="height:82px;border:1px solid black;text-align:center;"><center>$html_dejavu</center></div></div>};
$html .= qq{<div class="col-sm-4"><div style="height:82px;border:1px solid black;text-align:center;"><center>$html_filters</center></div></div>};
$html .= qq{<div class="col-sm-1"><div style="height:82px;border:1px solid black;text-align:center;"><center>$html_refresh</center></div></div>};
$html .= qq{</div></div>};
$html .= qq{<br>};

my $_tr_lines_by_genes;
my $h_junctions_color;
my $h_dejavu_cnv;
my $n = 0;


my $h_junctions_scores;
foreach my $junction (@lJunctions) {
	$n++;
	print '.' if ($n % 100);
	
	my $is_junction_linked_filtred;
	next if ($junction->isCanonique($patient));
	next if ($junction->get_ratio_new_count($patient) == 1);
	
	my $gene_name = $junction->annex->{$patient->name()}->{ensid};
	my $gene_name2 = $junction->annex->{$patient->name()}->{gene};
	if ($only_gene) {
		my $keep;
		$keep = 1 if ($only_gene->id() eq $gene_name);
		$keep = 1 if ($only_gene->external_name() eq $gene_name);
		$keep = 1 if ($only_gene->id() eq $gene_name2);
		$keep = 1 if ($only_gene->external_name() eq $gene_name2);
		next unless $keep;
	}
	
	if ($junction->get_percent_new_count($patient) < $min_score) {
		next;
		if (exists $h_var_linked_ids->{$junction->id()}) { $is_junction_linked_filtred = 1; }
		else { next; }
	}
	
	$junction->dejavu_percent_coordinate_similar($nb_percent_dejavu_value);
	my $nb_dejavu_pat = $junction->dejavu_nb_others_patients();
	$nb_dejavu_pat = $junction->dejavu_nb_other_patients_min_ratio_10($patient) if ($only_dejavu_ratio_10);
	if ($nb_dejavu_pat > $max_dejavu_value) {
		next;
		if (exists $h_var_linked_ids->{$junction->id()}) { $is_junction_linked_filtred = 1; }
		else { next; }
	}
	my $html_sashimi = get_sashimi_plot($junction, $patient);
	my $html_igv = get_igv_button($junction, $patient);
	my $html_id = get_html_id($junction);
	my $html_patients = get_html_patients($junction, $patient);
	my $html_dv = get_html_dejavu($junction, $patient);
	my $html_validation = '';
	my ($html_trans, $has_linked_junctions) = get_html_transcripts($junction, $patient);
	my $html_to_validate = '';
	
	my $score = $junction->junction_score($patient);
	#$score += 1 if ($junction->isSE($patient));
	#$score += 0.1 if ($junction->isRI($patient) and $has_linked_junctions);
	my $score_details_text = get_html_score_details($junction, $patient);
	
	$h_junctions_scores->{all}->{$gene_name}->{$score} = $junction->id();

	my $html_tr;
	if ($is_junction_linked_filtred) { $html_tr .= qq{<tr style="text-align:center;font-size:11px;opacity:0.55;">}; }
	else { $html_tr .= qq{<tr style="text-align:center;font-size:11px;">}; }
	$html_tr .= qq{<td style="width:230px;">$html_sashimi</td>};
	$html_tr .= qq{<td>$html_igv</td>};
	$html_tr .= qq{<td>$html_id</td>};
	$html_tr .= qq{<td>$html_patients</td>};
#	$html_tr .= qq{<td></td>};
	$html_tr .= qq{<td>$html_dv</td>};
#	$html_tr .= qq{<td>$html_validation</td>};
	$html_tr .= qq{<td>$html_trans</td>};
#	$html_tr .= qq{<td>$html_to_validate</td>};
	$html_tr .= qq{<td>$score</td>};
	$html_tr .= qq{<td>$score_details_text</td>};
	$html_tr .= qq{</tr>};
	push(@{$_tr_lines_by_genes->{$gene_name}}, $html_tr);
}
foreach my $gene_name (keys %{$h_junctions_scores->{all}}) {
	my @lscores = sort {$a <=> $b} keys %{$h_junctions_scores->{all}->{$gene_name}};
	$h_junctions_scores->{max}->{$lscores[-1]}->{$gene_name} = undef;
}
my @lGenesNames;
foreach my $score (sort {$b <=> $a} keys %{$h_junctions_scores->{max}}) {
	foreach my $gene_name (sort keys %{$h_junctions_scores->{max}->{$score}}) {
		push(@lGenesNames, $gene_name);
	}
}

if (scalar(@lGenesNames) == 0) {
	$html .= qq{<tr><td><center><b><i>No Result Found...</b></i></center></td></tr>};
	$html .= qq{</table>};
	my $hash;
	$hash->{html} = $html;
	$hash->{table_id} = '';
	if ($release =~ /HG19/) {
		$hash->{used_dejavu} = $max_dejavu_value;
		$hash->{used_dejavu} = $max_dejavu if $max_dejavu;
	}
	else {
		$hash->{used_dejavu} = 'not';
	}
	printJson($hash);
	exit(0);
}

my @lTablesIds;
$html .= "<table style='width:100%;>";
foreach my $gene_name (@lGenesNames) {
	my $class;
	my @lScores = sort {$a <=> $b} keys %{$h_junctions_scores->{all}->{$gene_name}};
	my $max_score = $lScores[-1];
	
	my $this_panel_gene_id = 'panel_'.$table_id.'_'.$gene_name;
	my $this_table_id = $table_id.'_'.$gene_name;
	push(@lTablesIds, $this_table_id);
	my $g = $project->newGene($gene_name);
	my $hgene;
	$hgene->{description} = $g->description();
	$hgene->{external_name} = $g->external_name();
	$hgene->{max_score} = $max_score;
	$hgene->{js_id} = $g->id();
	$hgene->{id} = $g->id();
	$hgene->{omim_id} = $g->omim_id();
	$hgene->{pLI} = $g->pLI();
	$hgene->{omim_inheritance} = $g->omim_inheritance();
	$hgene->{variants} = $_tr_lines_by_genes->{$gene_name};
	
	foreach my $panel (@{$g->getPanels()}) {
		$hgene->{panels}->{$panel->name()} = undef;
	}
	
	my $div_id = 'div_'.$this_table_id;
	$hgene->{collapse_with_id} = $div_id;
	
	my $class_tr_gene->{style} = "background-color:#607D8B;border:1px black solid;padding:0px;white-space: nowrap;height:50px;";
	
	$html .= qq{<tr>};
	$html .= "<table style='width:100%;background-color:#F3F3F3;'>";
	my $html_gene_panel = update_variant_editor::panel_gene($hgene, $this_panel_gene_id,$project->name(), $patient);
	$html_gene_panel =~ s/<b>pLI<\/b> Score//;
	$html_gene_panel =~ s/bottom:5px;/bottom:-5px;/;
	$html .= $cgi->td($class_tr_gene, $html_gene_panel);
	$html .= "</table>";
	$html .= $cgi->end_Tr();
	
	$html .= qq{<tr>};
	$html .= qq{<div class="collapse" id="$div_id">};
	$html .= qq{<table id='$this_table_id' data-sort-name='locus' data-sort-order='desc' data-filter-control='true' data-toggle="table" data-pagination-v-align="both" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="10" data-page-list="[10]" data-resizable='true' class='table table-striped' style='font-size:11px;'>};
	$html .= qq{<thead style="text-align:center;">};
	$html .= qq{<th data-field="plot"><b><center>Sashimi Plot</center></b></th>};
	$html .= qq{<th data-field="igv"><b><center>IGV</center></b></th>};
	$html .= qq{<th data-sortable='true' data-field="var_name"><b><center>Var name</center></b></th>};
	$html .= qq{<th data-field="trio"><b><center>Trio</center></b></th>};
#	$html .= qq{<th data-field="gnomad"><b><center>Gnomad</center></b></th>};
	$html .= qq{<th data-field="deja_vu"><b><center>DejaVu (Nb Samples)</center></b></th>};
#	$html .= qq{<th data-field="validations"><b><center>Validations</center></b></th>};
	$html .= qq{<th data-field="transcripts"><b><center>Transcripts</center></b></th>};
#	$html .= qq{<th data-field="to_validation"><b><center></center></b></th>};
	$html .= qq{<th data-sortable='true' data-field="jscore"><b><center>Score</center></b></th>};
	$html .= qq{<th data-field="noise"><b><center>Details Score</center></b></th>};
	$html .= qq{</thead>};
	$html .= qq{<tbody>};
	foreach my $l (@{$_tr_lines_by_genes->{$gene_name}}) {
		$html .= $l;
	}
	$html .= qq{</tbody>};
	$html .= qq{</table>};
	$html .= qq{</div>};
	$html .= $cgi->end_Tr();
	$html .= "<br>";
}
$html .= qq{</table>};


#$no_cnv->close();
$project->dejavuJunctions->close() if ($release =~ /HG19/);


my $hash;
$hash->{html} = $html;
$hash->{tables_id} = join(',', @lTablesIds);
#$hash->{table_id} = $table_id;
if ($release =~ /HG19/) {
	$hash->{used_dejavu} = $max_dejavu_value;
	$hash->{used_dejavu} = $max_dejavu if $max_dejavu;
}
else {
	$hash->{used_dejavu} = 'not';
}
printJson($hash);
exit(0);



sub transform_results_file {
	my ($type_analyse_name, $file_name, $path_analysis) = @_;
	my ($h_header, $list_res) = parse_results_file($file_name);
	#my ($list_sort) = sort_results($list_res);
	my ($table_id, $html) = add_table_results($type_analyse_name, $h_header, $list_res, $path_analysis);
	return ($table_id, $html);
}

sub get_sashimi_plot_list_files {
	my ($path_analysis, $ensg, $patient, $locus, $score) = @_;
	my $sashimi_plot_file = get_sashimi_plot_file($path_analysis, $ensg, $patient, $locus, $score);
	my @lFiles;
	if (-e $sashimi_plot_file) {
		push(@lFiles, $sashimi_plot_file);
		my $locus_text = $locus;
		$locus_text =~ s/chr//;
		$locus_text =~ s/:/-/;
		my ($chr_id, $start, $end) = split('-', $locus_text);
		my $i = 1;
		while ($i < 5) {
			$start -= (1000*$i);
			$end += (1000*$i);
			my $locus_extended = $chr_id.':'.$start.'-'.$end;
			my $sashimi_plot_file = get_sashimi_plot_file($path_analysis, $ensg, $patient, $locus_extended, $score);
			push(@lFiles, $sashimi_plot_file);
			$i++;
		}
		return \@lFiles;
	}
	return;
}

sub get_sashimi_plot_file {
	my ($path_analysis, $ensg, $patient, $locus, $score) = @_;
	my $locus_text = $locus;
	$locus_text =~ s/chr//;
	$locus_text =~ s/:/-/;
	my $patient_name = $patient->name();
	my $path = $patient->getProject->getProjectPath.'/align/sashimi_plots/';
	my $outfile = $path.'/sashimi_'.$patient_name.'.'.$ensg.'.'.$locus_text.'.pdf';
	return $outfile if (-e $outfile);
}


sub get_html_dejavu {
	my ($junction, $patient) = @_;
	my $color = 'black';
	my $my_ratio = $junction->get_percent_new_count($patient);
	my $project_name = $patient->getProject->name();
	my $patient_name = $patient->name();
	my $vector_id = $junction->getChromosome->id().'-'.$junction->vector_id();
	my $cmd_all = qq{view_deja_vu_rna_junction(\"$project_name\",\"$patient_name\",\"$vector_id\")};
	my $cmd_inthisrun = qq{view_dejavu_nb_int_this_run_patients(\"$project_name\",\"$patient_name\",\"$vector_id\")};
	my $html = $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
	$html.= $cgi->start_Tr();
	$html.= $cgi->th("");
	$html.= $cgi->th("<center><b>Ratio All</b></center>");
#	if ($my_ratio >= 70) {
#		$html.= $cgi->th("<center><b>Ratio >70%</b></center>");
#		$html.= $cgi->th("<center><b>Ratio >90%</b></center>");
#	}
#	else {
		$html.= $cgi->th("<center><b>Ratio >10%</b></center>");
		$html.= $cgi->th("<center><b>Ratio >20%</b></center>");
#	}
	$html.= $cgi->end_Tr();
	
	$html.= $cgi->start_Tr();
	$html.= $cgi->td("<center><b>DejaVu</b></center>");
	$html.= $cgi->td(obutton($cmd_all, $junction->dejavu_nb_other_patients($patient)));
#	if ($my_ratio >= 70) {
#		$html.= $cgi->td(obutton($cmd_all, $junction->dejavu_nb_other_patients_min_ratio_70($patient)));
#		$html.= $cgi->td(obutton($cmd_all, $junction->dejavu_nb_other_patients_min_ratio_90($patient)));
#	}
#	else {
		$html.= $cgi->td(obutton($cmd_all, $junction->dejavu_nb_other_patients_min_ratio_10($patient)));
		$html.= $cgi->td(obutton($cmd_all, $junction->dejavu_nb_other_patients_min_ratio_20($patient)));
#	}
	$html.= $cgi->end_Tr();
	
	$html.= $cgi->start_Tr();
	$html.= $cgi->td("<center><b>InThisRun</b></center>");
	$html.= $cgi->td(obutton($cmd_inthisrun,$junction->dejavu_nb_int_this_run_patients($patient)));
#	if ($my_ratio >= 70) {
#		$html.= $cgi->td(obutton($cmd_inthisrun,$junction->dejavu_nb_int_this_run_patients($patient,70)));
#		$html.= $cgi->td(obutton($cmd_inthisrun,$junction->dejavu_nb_int_this_run_patients($patient,90)));
#	}
#	else {
		$html.= $cgi->td(obutton($cmd_inthisrun,$junction->dejavu_nb_int_this_run_patients($patient,10)));
		$html.= $cgi->td(obutton($cmd_inthisrun,$junction->dejavu_nb_int_this_run_patients($patient,20)));
#	}
	$html.= $cgi->end_Tr();
	$html.=$cgi->end_table();
	return $html;
}

sub get_html_transcripts {
	my ($junction, $patient) = @_;
	my $has_linked_junctions;
	my $color = 'black';
	my $project_name = $patient->getProject->name();
	my $patient_name = $patient->name();
	my $is_junction_linked = $junction->is_junctions_linked($patient);
	my (@junctions_ids_linked, $bcolor);
	@junctions_ids_linked = keys %{$junction->get_hash_junctions_linked_to_me->{$patient->name()}} if ($is_junction_linked);
	my @l_group_junctions_colors = ('#5D3EFF', '#FF4571', '#8FFF49', '#FF495F');
	my $h_exons_introns = $junction->get_hash_exons_introns();
	
	my $html_tr = $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
	$html_tr.= "<tr style='background-color:#FFA81E;'>";
	$html_tr.= $cgi->th("<center><b>enst</b></center>");
	$html_tr.= $cgi->th("<center><b>nm</b></center>");
	$html_tr.= $cgi->th("<center><b>ccds</b></center>");
	$html_tr.= $cgi->th("<center><b>appris</b></center>");
	$html_tr.= $cgi->th("<center><b>start</b></center>");
	$html_tr.= $cgi->th("<center><b>end</b></center>");
	$html_tr.= $cgi->end_Tr();
	
	foreach my $tid (sort keys %{$h_exons_introns}) {
		my ($h_junctions_linked, $h_junctions_exons_introns);
		if (exists $h_junctions_color->{$junction->id()}) {
			$bcolor = $h_junctions_color->{$junction->id()};
		}
		else {
			$bcolor = $l_group_junctions_colors[$nb_group_junctions_colors];
			$nb_group_junctions_colors ++;
			$nb_group_junctions_colors = 0 if ($nb_group_junctions_colors+1 == scalar(@l_group_junctions_colors));
		}
		
		my $t = $patient->getProject->newTranscript($tid);
		$html_tr .= $cgi->start_Tr();
		$html_tr .= $cgi->td("<center>$tid</center>");
		$html_tr .= $cgi->td("<center>".$t->external_name()."</center>");
		$html_tr .= $cgi->td("<center>".$t->ccds_name()."</center>");
		$html_tr .= $cgi->td("<center>".$t->appris_type()."</center>");
		my @lPos = (sort keys %{$h_exons_introns->{$tid}->{by_pos}});
		my $first_exon_intron = $h_exons_introns->{$tid}->{by_pos}->{$lPos[0]};
		my $last_exon_intron = $h_exons_introns->{$tid}->{by_pos}->{$lPos[-1]};
		my ($first_style_color, $last_style_color);
		if ($is_junction_linked) {
			foreach my $other_j_id (@junctions_ids_linked) {
				$first_style_color = "style='background-color:$bcolor;'" if (exists $junction->get_hash_junctions_linked_to_me->{$patient->name()}->{$other_j_id}->{$tid}->{$first_exon_intron});
				$last_style_color = "style='background-color:$bcolor;'" if (exists $junction->get_hash_junctions_linked_to_me->{$patient->name()}->{$other_j_id}->{$tid}->{$last_exon_intron});
				if ($first_style_color or $last_style_color) {
					
					#NEW
					$h_junctions_linked->{$junction->id()} = $h_var_linked_ids->{$junction->id()}->{vector_id};
					foreach my $other_id (keys %{$h_var_linked_ids->{$junction->id()}->{linked_to}}) {
						$h_junctions_linked->{$other_j_id} = $h_var_linked_ids->{$other_j_id}->{vector_id};
					}
					$h_junctions_exons_introns->{$tid}->{$first_exon_intron} = undef;
					$h_junctions_exons_introns->{$tid}->{$last_style_color} = undef;
					
					my $this_chr = $junction->getChromosome();
					my $i = $junction->vector_id();
					my $min = $junction->vector_id() - 15;
					$min = 0 if ($min < 0);
					while ($i >= $min) {
						my $junction2 = $this_chr->getVarObject($i);
						$i--;
						next if (not $junction2->get_hash_junctions_linked_to_me());
						
						my @lPos = (sort keys %{$junction2->get_hash_exons_introns->{$tid}->{by_pos}});
						my $first_exon_intron_2 = $junction2->get_hash_exons_introns->{$tid}->{by_pos}->{$lPos[0]};
						my $last_exon_intron_2 = $junction2->get_hash_exons_introns->{$tid}->{by_pos}->{$lPos[-1]};
						next if (not exists $h_junctions_exons_introns->{$tid}->{$first_exon_intron_2} and not exists $h_junctions_exons_introns->{$tid}->{$last_exon_intron_2});
						
						foreach my $other_jid (keys %$h_junctions_linked) {
							if (exists $junction2->get_hash_junctions_linked_to_me->{$patient->name()}->{$other_jid}->{$tid}) {
								$h_junctions_linked->{$junction2->id()} = $this_chr->id().'-'.$junction2->vector_id();
								$h_junctions_exons_introns->{$tid}->{$first_exon_intron_2} = undef;
								$h_junctions_exons_introns->{$tid}->{$last_exon_intron_2} = undef;
								$min -= 5;
								$min = 0 if ($min < 0);
							}
						}
					}
					$i = $junction->vector_id();
					my $max = $junction->vector_id() + 15;
					$max = $this_chr->size_vector() if ($max >= $this_chr->size_vector());
					while ($i < $max) {
						my $junction2 = $this_chr->getVarObject($i);
						$i++;
						next if (not $junction2->get_hash_junctions_linked_to_me());
						
						my @lPos = (sort keys %{$junction2->get_hash_exons_introns->{$tid}->{by_pos}});
						my $first_exon_intron_2 = $junction2->get_hash_exons_introns->{$tid}->{by_pos}->{$lPos[0]};
						my $last_exon_intron_2 = $junction2->get_hash_exons_introns->{$tid}->{by_pos}->{$lPos[-1]};
						
						next if (not exists $h_junctions_exons_introns->{$tid}->{$first_exon_intron_2} and not exists $h_junctions_exons_introns->{$tid}->{$last_exon_intron_2});
						
						foreach my $other_jid (keys %$h_junctions_linked) {
							if (exists $junction2->get_hash_junctions_linked_to_me->{$patient->name()}->{$other_jid}->{$tid}) {
								next if (scalar keys %{$junction2->get_hash_junctions_linked_to_me->{$patient->name()}->{$other_jid}->{$tid}} ==0 );
								$h_junctions_linked->{$junction2->id()} = $this_chr->id().'-'.$junction2->vector_id();
								$h_junctions_exons_introns->{$tid}->{$first_exon_intron_2} = undef;
								$h_junctions_exons_introns->{$tid}->{$last_exon_intron_2} = undef;
								$max += 5;
								$max = $this_chr->size_vector() if ($max >= $this_chr->size_vector());
							}
						}
					}
				}
			}
		}
		foreach my $jid (keys %$h_junctions_linked) {
			$h_junctions_color->{$jid} = $bcolor;
		}
		my $j_linked = join(',', sort values %$h_junctions_linked);
		$h_junctions_linked = undef;
		my $my_junction_id = $junction->id();
		my $cmd_linked = qq{view_linked_junctions(\"$patient_name\",\"$tid\",\"$j_linked\",\"$my_junction_id\",\"$min_score\")};
		if (scalar(@lPos) == 1) {
			if ($first_style_color) { $html_tr .= "<td colspan='2' $first_style_color>".obutton($cmd_linked,$first_exon_intron)."</td>"; }
			else { $html_tr .= "<td colspan='2' $first_style_color>$first_exon_intron</td>"; }
		}
		else {
			if ($first_style_color) { $html_tr .= "<td $first_style_color>".obutton($cmd_linked, $first_exon_intron)."</td>"; }
			else { $html_tr .= "<td $first_style_color>$first_exon_intron</td>"; }
			if ($last_style_color) { $html_tr .= "<td $last_style_color>".obutton($cmd_linked, $last_exon_intron)."</td>"; }
			else { $html_tr .= "<td $last_style_color>$last_exon_intron</td>"; }
		}
		$has_linked_junctions = 1 if ($first_style_color or $last_style_color);
	}
	$html_tr.= $cgi->end_Tr();
	$html_tr.= qq{</table>};
	return ($html_tr, $has_linked_junctions);
}

sub obutton {
	my ($onclick,$value) = @_;
	return qq{<a class="btn btn-xs btn-primary" onclick=\'$onclick\' target="_blank" style="background-color: #D0D0D0;font-size: 7px;font-family:  Verdana;color:black" role="button">}.$value."</a>";
}

sub get_html_patients {
	my ($junction, $patient) = @_;
	my $h_by_pat;
	foreach my $pat (@{$patient->getFamily->getPatients()}) {
		next if (not $junction->get_dp_count($pat));
		next if (not $junction->get_nb_new_count($pat));
		my $fam_name = $pat->getFamily->name();
		if ($pat->isFather()) {
			if ($pat->isIll()) { $h_by_pat->{$fam_name}->{$pat->name()}->{status} = "<center><img src='/icons/Polyicons/male-d.png'></center>"; }
			else { $h_by_pat->{$fam_name}->{$pat->name()}->{status} = "<center><img src='/icons/Polyicons/male-s.png'></center>"; }
		}
		if ($pat->isMother()) {
			if ($pat->isIll()) { $h_by_pat->{$fam_name}->{$pat->name()}->{status} = "<center><img src='/icons/Polyicons/female-d.png'></center>"; }
			else { $h_by_pat->{$fam_name}->{$pat->name()}->{status} = "<center><img src='/icons/Polyicons/female-s.png'></center>"; }
		}
		if ($pat->isChild()) {
			if ($pat->sex() eq '1') { 
				if ($pat->isIll()) { $h_by_pat->{$fam_name}->{$pat->name()}->{status} = "<center><img src='/icons/Polyicons/baby-boy-d.png'></center>"; }
				else { $h_by_pat->{$fam_name}->{$pat->name()}->{status} = "<center><img src='/icons/Polyicons/baby-boy-s.png'></center>"; }
			}
			else {
				if ($pat->isIll()) { $h_by_pat->{$fam_name}->{$pat->name()}->{status} = "<center><img src='/icons/Polyicons/baby-girl-d.png'></center>"; }
				else { $h_by_pat->{$fam_name}->{$pat->name()}->{status} = "<center><img src='/icons/Polyicons/baby-girl-s.png'></center>"; }
			}
		}
		$h_by_pat->{$fam_name}->{$pat->name()}->{dp} = $junction->get_dp_count($pat);
		$h_by_pat->{$fam_name}->{$pat->name()}->{nb_new} = $junction->get_nb_new_count($pat);
		$h_by_pat->{$fam_name}->{$pat->name()}->{nb_normal} = $junction->get_nb_normal_count($pat);
		$h_by_pat->{$fam_name}->{$pat->name()}->{percent} = sprintf("%.3f", $junction->get_percent_new_count($pat)).'%';
	}
	my $color = 'black';
	my $html_patients= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
	$html_patients .= qq{<thead style="text-align:center;">};
	#$html_patients .= qq{<th data-field="famname"><b><center>Family</center></b></th>};
	$html_patients .= qq{<th data-field="patname"><b><center>Patient</center></b></th>};
	$html_patients .= qq{<th data-field="status"><b><center>Status</center></b></th>};
	$html_patients .= qq{<th data-field="percent"><b><center>Ratio (%)</center></b></th>};
	$html_patients .= qq{<th data-field="nb_new"><b><center>Nb New</center></b></th>};
	$html_patients .= qq{<th data-field="nb_normal"><b><center>Nb Normal</center></b></th>};
	$html_patients .= qq{<th data-field="dp"><b><center>DP</center></b></th>};
	$html_patients .= qq{<th data-field=""><b><center></center></b></th>};
	$html_patients .= qq{</thead>};
	$html_patients .= qq{<tbody>};
	foreach my $fam_name (sort keys %{$h_by_pat}) {
		foreach my $pat_name (sort keys %{$h_by_pat->{$fam_name}}) {
			$html_patients .= qq{<tr>};
			#$html_patients .= qq{<td>}.$fam_name.qq{</td>};
			$html_patients .= qq{<td>}.$pat_name.qq{</td>};
			$html_patients .= qq{<td>}.$h_by_pat->{$fam_name}->{$pat_name}->{status}.qq{</td>};
			$html_patients .= qq{<td>}.$h_by_pat->{$fam_name}->{$pat_name}->{percent}.qq{</td>};
			$html_patients .= qq{<td>}.$h_by_pat->{$fam_name}->{$pat_name}->{nb_new}.qq{</td>};
			$html_patients .= qq{<td>}.$h_by_pat->{$fam_name}->{$pat_name}->{nb_normal}.qq{</td>};
			$html_patients .= qq{<td>}.$h_by_pat->{$fam_name}->{$pat_name}->{dp}.qq{</td>};
			$html_patients .= qq{</tr>};
		}
	}
	$html_patients .= qq{</tbody>};
	$html_patients .= qq{</table>};
	return $html_patients;
}

sub get_igv_button {
	my ($junction, $patient) = @_;
	my $color = 'lightgrey';
	my $bam_file = "https://www.polyweb.fr/".$patient->bamUrl();
	my $list_patients_ctrl = $patient->getPatients_used_control_rna_seq_junctions_analyse();
	if ($list_patients_ctrl) {
		my $nb_control;
		foreach my $other_pat (@$list_patients_ctrl) {
			$bam_file .= ',https://www.polyweb.fr/'.$other_pat->bamUrl();
			$nb_control++;
			last if $nb_control == 3;
		}
	}
	else {
		my $np = 0;
		foreach my $other_pat (@{$project->getPatients()}) {
			next if ($other_pat->name() eq $patient->name());
			$bam_file .= ',https://www.polyweb.fr/'.$other_pat->bamUrl();
			$np++;
			last if $np == 3;
		}
	}
	my $gtf = $patient->getProject->get_gtf_genes_annotations_igv();
	$gtf =~ s/\/data-isilon//;
	$gtf = "https://www.polyweb.fr/".$gtf;
	my $locus = $junction->getChromosome->id().':'.($junction->start()-100).'-'.($junction->end()+100);
	my $igv_link = qq{<button class='igvIcon2' onclick='launch_igv_tool_rna("", "$bam_file,$gtf","$locus")' style="color:black"></button>};
	return $igv_link;
}

sub get_sashimi_plot {
	my ($junction, $patient) = @_;
	my $sashimi_button;
	my $list_sashimi_plot_files = $junction->getListSashimiPlotsPathFiles($patient);
	if ($list_sashimi_plot_files and -e $list_sashimi_plot_files->[0]) {
		$sashimi_button .= qq{<center>};
		my @lFiles;
		foreach my $sashimi_plot_file (@$list_sashimi_plot_files) {
			$sashimi_plot_file =~ s/\/\//\//g;
			$sashimi_plot_file =~ s/\/data-isilon\/sequencing\/ngs/\/NGS/;
			#$sashimi_plot_file = "https://www.polyweb.fr/".$sashimi_plot_file;
			push(@lFiles, $sashimi_plot_file);
		}
		my $files = join(';', @lFiles);
		my $pdf = $lFiles[0].'#toolbar=0&embedded=true';
		$sashimi_button .= qq{<button type="button" class="btn btn-default" style="border:2px black double;overflow:hidden;text-align:center;background-color:white;padding-right:20px;padding-left:4px;padding-top:4px;padding-bottom:4px;" onClick="view_pdf_list_files('$files')"><table><td>};
		$sashimi_button .= qq{<image style="position:relative;width:200px;" loading="lazy" src="$pdf"></image>};
		#$sashimi_button .= qq{<image style="position:relative;z-index:2;width:200px;object-position: -50% -50%;transform: scale(1.7) translate(-27px, 20px);" loading="lazy" src="$pdf"></image>};
		$sashimi_button .= qq{</td><td style="padding-left:1px;"><span style="writing-mode:vertical-lr !important; font: 12px Verdana, sans-serif;letter-spacing: 1px;">Zoom</span></td></table> </button>};
		$sashimi_button .= qq{</center></};
	}
	else {
		$sashimi_button .= qq{<center>N.A.</center>};
	}
	return $sashimi_button;
}

sub get_html_id {
	my ($junction) = @_;
	my $chr_id = $junction->getChromosome->id();
	my $start = $junction->start();
	my $end = $junction->end();
	my $junction_locus = $chr_id.':'.$start.'-'.$end;
	my $length = $junction->length();
	my @lTypes;
	push(@lTypes, 'RI') if $junction->isRI($patient);
	push(@lTypes, 'SE') if $junction->isSE($patient);
	my $type_junction = join('+', @lTypes);
	my $type_junction_description = $junction->getTypeDescription($patient);
	my $html_id = "<center><table>";
	$html_id .= "<tr><td><center><b>$junction_locus</b></center></td></tr>";
	$html_id .= "<tr><td><center>$length nt</td></tr>";
	$html_id .= "<tr><td><center>";
	$html_id .= "$type_junction";
	$html_id .= " - $type_junction_description" if ($type_junction_description and $type_junction_description ne '---');
	$html_id .= "</center></td></tr>";
	$html_id .= "</table></center>";
	return $html_id;
}

sub get_html_score_details {
	my ($junction, $patient) = @_;
	my $noise = $junction->get_noise_score($patient);
	my $score_pen_ratio = $junction->junction_score_penality_ratio($patient);
	my $score_pen_dp = $junction->junction_score_penality_dp($patient);
	my $score_pen_new = $junction->junction_score_penality_new_junction($patient);
	my $score_pen_noise = $junction->junction_score_penality_noise($patient);
	my $score_pen_dvrun = $junction->junction_score_penality_dejavu_inthisrun($patient);
	my $score_pen_dv = $junction->junction_score_penality_dejavu($patient);
	my $score_details_text = $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px black;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
	if ($score_pen_ratio > 0) {
		$score_details_text.= $cgi->start_Tr();
		$score_details_text.= $cgi->td("<center><b>Ratio</b></center>");
		$score_details_text.= $cgi->td("<center>- $score_pen_ratio</center>");
		$score_details_text.= $cgi->end_Tr();
	}
	if ($score_pen_dp > 0) {
		$score_details_text.= $cgi->start_Tr();
		$score_details_text.= $cgi->td("<center><b>DP</b></center>");
		$score_details_text.= $cgi->td("<center>- $score_pen_dp</center>");
		$score_details_text.= $cgi->end_Tr();
	}
	if ($score_pen_new > 0) {
		$score_details_text.= $cgi->start_Tr();
		$score_details_text.= $cgi->td("<center><b>New Junc</b></center>");
		$score_details_text.= $cgi->td("<center>- $score_pen_new</center>");
		$score_details_text.= $cgi->end_Tr();
	}
	if ($score_pen_dvrun > 0) {
		$score_details_text.= $cgi->start_Tr();
		$score_details_text.= $cgi->td("<center><b>Inthisrun</b></center>");
		$score_details_text.= $cgi->td("<center>- $score_pen_dvrun</center>");
		$score_details_text.= $cgi->end_Tr();
	}
	if ($score_pen_dv > 0) {
		$score_details_text.= $cgi->start_Tr();
		$score_details_text.= $cgi->td("<center><b>DejaVu</b></center>");
		$score_details_text.= $cgi->td("<center>- $score_pen_dv</center>");
		$score_details_text.= $cgi->end_Tr();
	}
	if ($score_pen_noise > 0) {
		$score_details_text.= $cgi->start_Tr();
		$score_details_text.= $cgi->td("<center><b>Noise<br></b>$noise found</center>");
		$score_details_text.= $cgi->td("<center>- $score_pen_noise</center>");
		$score_details_text.= $cgi->end_Tr();
	}
	$score_details_text .= "</table>";
	return $score_details_text;
}

sub printJson {
	my ($hashRes) = @_;
	my $json_encode = encode_json $hashRes;
	print ".\",";
	$json_encode =~ s/{//;
	print $json_encode;
	exit(0);
}
