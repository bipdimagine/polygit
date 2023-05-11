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
use List::MoreUtils qw{ natatime };
use Parallel::ForkManager;

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
my $only_html_cache = $cgi->param('only_html_cache');

#$only_positions = '10:17000000-17080000';

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
if ($release =~ /HG19/) { $gencode = $project->gencode_version(); }
elsif ($release eq 'MM38') { $gencode = 'M25'; }
elsif ($release eq 'MM39') { $gencode = 'M32'; }


my $max_dejavu_value = 51;
$max_dejavu_value = $max_dejavu if defined($max_dejavu);

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

my $n = 0;
my $score_slider = 0;
$score_slider = $min_score / 10 if ($min_score and $min_score > 0);
my (@lJunctions, $h_var_linked_ids);

my $j_total = 0;
my $j_selected = 0;
my $h_chr_vectors;
my $h_chr_vectors_counts;
foreach my $chr (@{$project->getChromosomes()}) {
	#next if $chr->id ne '1' and $patient_name eq 'KUCerc';
	my $vector_patient = $patient->getJunctionsVector($chr);
	$j_total += $chr->countThisVariants($vector_patient);
	if ($only_positions) {
		my @lTmp = split(':', $only_positions);
		my $chr_id_limit = lc($lTmp[0]);
		$chr_id_limit =~ s/chr//;
		next if lc($chr->id) ne $chr_id_limit;
		my ($from, $to) = split('-', $lTmp[1]);
		my $vector_postions = $chr->getVectorByPosition($from, $to);
		$vector_patient->Intersection($vector_patient, $vector_postions);
	}
	if ($only_gene) {
		next if ($chr->id() ne $only_gene->getChromosome->id());
		my $vector_postions = $chr->getVectorByPosition(($only_gene->start - 1000), ($only_gene->end + 1000));
		$vector_patient->Intersection($vector_patient, $vector_postions);
	}
	
	$h_chr_vectors->{$chr->id} = $vector_patient->Clone();
	$h_chr_vectors_counts->{$chr->id} = $chr->countThisVariants($h_chr_vectors->{$chr->id});
	$j_selected += $h_chr_vectors_counts->{$chr->id};
}


my $cache_id = 'splices_linked_'.$patient->name().'_'.$j_total;
#warn $cache_id;
my $no_cache = $patient->get_lmdb_cache("r");
my $h_res = $no_cache->get_cache($cache_id);
$no_cache->close();
if ($h_res) {
	$h_var_linked_ids = $h_res;
}
else {
	add_linked_hash_in_cache($cache_id, $patient); 
	my $no_cache = $patient->get_lmdb_cache("r");
	$h_var_linked_ids = $no_cache->get_cache($cache_id);
	$no_cache->close();
}

exit(0) if ($only_html_cache);

if ($j_selected >= 10000) {
	my $no_cache = $patient->get_lmdb_cache("r");
	my $cache_vectors_enum_id = 'splices_linked_'.$patient->name().'_'.$j_total.'_chr_vectors_enum';
	my $h_res_v_enum = $no_cache->get_cache($cache_vectors_enum_id);
	my $use_cat = countMinCatToUse($h_res_v_enum); 
	print ".$use_cat.";
	foreach my $chr_id (keys %$h_res_v_enum) {
		my $v_filters = $project->getChromosome($chr_id)->getNewVector();
		$v_filters->from_Enum($h_res_v_enum->{$chr_id}->{$use_cat});
		$h_chr_vectors->{$chr_id}->Intersection($h_chr_vectors->{$chr_id}, $v_filters);
	}
}

my $percent_dejavu = 0;
if (defined($use_percent_dejavu)) {
	$percent_dejavu = $use_percent_dejavu - 90;
}
my $nb_percent_dejavu_value = 90 + $percent_dejavu;

my $checked_only_dejavu_ratio_10;
$checked_only_dejavu_ratio_10 = qq{checked="checked"} if $only_dejavu_ratio_10;
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
	

my $html_gene_select = qq{<input type="text" style="font-size:11px;width:180px;" placeholder="COL4A5, 1:100000-150000" class="form-control" id="input_gene_positions">};
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
$n = 0;

my $h_junctions_scores;

my $fork = 4;
my $pm = new Parallel::ForkManager($fork);
my $nbErrors = 0;
$pm->run_on_finish(
	sub {
		my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $hres ) = @_;
		unless ( defined($hres) or $exit_code > 0 ) {
			$nbErrors++;
			print qq|No message received from child process $exit_code $pid!\n|;
			warn Dumper $hres;
			return;
		}
		foreach my $gene_name (keys %{$hres->{genes}}) {
			foreach my $html_tr (@{$hres->{genes}->{$gene_name}}) {
				push(@{$_tr_lines_by_genes->{$gene_name}}, $html_tr);
			}
		}
		
		foreach my $gene_name (keys %{$hres->{score}->{all}}) {
			foreach my $score (keys %{$hres->{score}->{all}->{$gene_name}}) {
				$h_junctions_scores->{all}->{$gene_name}->{$score} = undef;
			}
		}
	}
);
		
		
#my $nb_elems = int(scalar(@lJunctions) / $fork);
#$nb_elems += 20;
#
#my $iter = natatime $nb_elems, @lJunctions;
#while( my @tmp = $iter->() ) {
	
	
foreach my $chr_id (sort keys %{$h_chr_vectors}) {
	my $chr = $patient->getProject->getChromosome($chr_id);
	next if $chr->countThisVariants($h_chr_vectors->{$chr_id}) == 0;
	
	$patient->getProject->buffer->dbh_deconnect();
	$pm->start and next;
	$patient->getProject->buffer->dbh_reconnect();
	
	my @lJunctionsChr = @{$chr->getListVarObjects($h_chr_vectors->{$chr_id})};
	my $hres;
	foreach my $junction (@lJunctionsChr) {
		$n++;
		print '.' if ($n % 1000);
		
		my $is_junction_linked_filtred;
		next if ($junction->isCanonique($patient));
		
		next if ($junction->junction_score_without_dejavu_global($patient) < 0);
		
		#next if ($junction->get_ratio_new_count($patient) == 1);
		
		my $gene_name = $junction->annex->{$patient->name()}->{ensid};
		my $gene_name2 = $junction->annex->{$patient->name()}->{gene};
		my @lGenesNames;
		
		if ($only_gene) {
			my $keep;
			$keep = 1 if ($only_gene->id() eq $gene_name);
			$keep = 1 if ($only_gene->external_name() eq $gene_name);
			$keep = 1 if ($only_gene->id() eq $gene_name2);
			$keep = 1 if ($only_gene->external_name() eq $gene_name2);
			next unless $keep;
		}
		
		if ($gene_name) { push (@lGenesNames, $gene_name); }
		else {
			foreach my $gene (@{$junction->getGenes()}) {
				push (@lGenesNames, $gene->id());
			}
		}
		next unless @lGenesNames;
		
		if (not $only_gene) {
			if ($junction->get_percent_new_count($patient) < $min_score) {
				next;
				if (exists $h_var_linked_ids->{$junction->id()}) { $is_junction_linked_filtred = 1; }
				else { next; }
			}
			
			eval {
				$junction->dejavu_percent_coordinate_similar($nb_percent_dejavu_value);
				my $nb_dejavu_pat = $junction->dejavu_nb_others_patients();
				$nb_dejavu_pat = $junction->dejavu_nb_other_patients_min_ratio_10($patient) if ($only_dejavu_ratio_10);
				next if ($nb_dejavu_pat > $max_dejavu_value);
			};
			if ($@) { next; }
		}
		
	#	my $nb_dejavu_pat = $junction->dejavu_nb_others_patients();
	#	$nb_dejavu_pat = $junction->dejavu_nb_other_patients_min_ratio_10($patient) if ($only_dejavu_ratio_10);
	#	if ($nb_dejavu_pat > $max_dejavu_value) {
	#		next;
	#		if (exists $h_var_linked_ids->{$junction->id()}) { $is_junction_linked_filtred = 1; }
	#		else { next; }
	#	}
		my $html_sashimi = get_sashimi_plot($junction, $patient);
		my $html_igv = get_igv_button($junction, $patient);
		my $html_id = get_html_id($junction);
		my $html_patients = get_html_patients($junction, $patient);
		my $html_dv = get_html_dejavu($junction, $patient);
		
		
		my $html_validation = '';
		my $html_to_validate = '';
		
		foreach my $gene_name (@lGenesNames) {
			
			my $ht = $junction->get_hash_exons_introns();
			next unless $ht;
#			next if scalar keys %$ht == 0;
			
			my ($html_trans, $has_linked_junctions) = get_html_transcripts($junction, $patient);
			my $score = $junction->junction_score($patient);
			my $score_details_text = get_html_score_details($junction, $patient);
			
			$hres->{score}->{all}->{$gene_name}->{$score} = $junction->id();
		
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
			push(@{$hres->{genes}->{$gene_name}}, $html_tr);
		}
	}
	$hres->{done} = 1;
	$pm->finish( 0, $hres );
}
$pm->wait_all_children();

die if $nbErrors > 0;


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
	
	if ($only_gene) {
		my $b_igv = get_igv_button($only_gene, $patient);
		$html .= qq{<br><tr><td><center><b><i>...but you can view gene <span style="color:orange;">$only_gene_name</span> in IGV :)</b></i></center></td></tr>};
		$html .= qq{<br><tr><td><center><b><i><span class="glyphicon glyphicon-arrow-down"></span></b></i></center></td></tr><br><tr><td><center><b><i>$b_igv</b></i></center></td></tr>};
	}
	else {
		$html .= qq{<br><tr><td><center><b><i><span class="glyphicon glyphicon-exclamation-sign" style="color:red;font-size:18px;"></span> your gene <span style="color:orange;">$only_gene_name</span> is unknown in my database...</b></i></center></td></tr><br><br><br>};
	}
	
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
elsif (scalar(@lGenesNames) > 0 and $only_gene_name and not $only_gene) {
	$html .= qq{<br><tr><td><b><i><span class="glyphicon glyphicon-exclamation-sign" style="color:red;font-size:18px;"></span> your gene <span style="color:orange;">$only_gene_name</span> is unknown in my database...</b></i></td></tr><br><br><br>};
}

my @lTablesIds;

my $filters = qq{data-filter-control="input" data-filter-control-placeholder="Gene Name / Description / Panel"};
my $numbers_per_page = "100";

my $html_table_1 = qq{<table id='table_major' data-filter-control='true' data-toggle="table" data-pagination-v-align="bottom" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="$numbers_per_page" data-page-list="[$numbers_per_page]" data-resizable='true' class='table' style='font-size:11px;'>};
$html_table_1 .= qq{<thead style="text-align:center;">};
$html_table_1 .= qq{<th data-field="gene" $filters><b><center></center></b></th>};
$html_table_1 .= qq{</thead>};
$html_table_1 .= qq{<tbody>};

my $html_table_2;

#$html .= qq{<table id='table_major' data-filter-control='true' data-toggle="table" data-pagination-v-align="bottom" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="$numbers_per_page" data-page-list="[$numbers_per_page]" data-resizable='true' class='table' style='font-size:11px;'>};
#$html .= qq{<thead style="text-align:center;">};
#$html .= qq{<th data-field="gene" $filters><b><center></center></b></th>};
#$html .= qq{</thead>};
#$html .= qq{<tbody>};
push(@lTablesIds, 'table_major');

my (@l_html, @l_tables_ids);
my $nb_g = 0;
foreach my $gene_name (@lGenesNames) {
	my $class;
	my @lScores = sort {$a <=> $b} keys %{$h_junctions_scores->{all}->{$gene_name}};
	my $max_score = $lScores[-1];
	
	my $this_panel_gene_id = 'panel_'.$table_id.'_'.$gene_name;
	my $this_table_id = $table_id.'_'.$gene_name;
#	push(@lTablesIds, $this_table_id);
	my $g = $project->newGene($gene_name);
	my $hgene;
	$hgene->{description} = $g->description();
	$hgene->{external_name} = $g->external_name();
	$hgene->{max_score} = $max_score;
	$hgene->{js_id} = $g->id();
	$hgene->{id} = $g->id();
	
	$hgene->{pLI} = $g->pLI();
	$hgene->{omim_id} = $g->omim_id();
	$hgene->{omim_inheritance} = $g->omim_inheritance();
	foreach my $panel (@{$g->getPanels()}) {
		$hgene->{panels}->{$panel->name()} = undef;
	}
	
	$hgene->{variants} = $_tr_lines_by_genes->{$gene_name};
	
	
	my $div_id = 'div_'.$this_table_id;
	$hgene->{collapse_with_id} = $div_id;
	
	my $class_tr_gene->{style} = "background-color:#607D8B;border:1px black solid;padding:0px;white-space: nowrap;height:50px;";
	
	$html_table_2 .= qq{<tr><td>};
	$html_table_2 .= "<table style='width:100%;background-color:#F3F3F3;'>";
	my $html_gene_panel = update_variant_editor::panel_gene($hgene, $this_panel_gene_id,$project->name(), $patient);
	$html_gene_panel =~ s/<b>pLI<\/b> Score//;
	$html_gene_panel =~ s/bottom:5px;/bottom:-5px;/;
	$html_table_2 .= $cgi->td($class_tr_gene, $html_gene_panel);
	$html_table_2 .= "</table>";
	
	$html_table_2 .= qq{<div class="collapse" id="$div_id">};
	$html_table_2 .= qq{<table id='$this_table_id' data-sort-name='locus' data-sort-order='desc' data-filter-control='true' data-toggle="table" data-pagination-v-align="both" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="10" data-page-list="[10]" data-resizable='true' class='table table-striped' style='font-size:11px;'>};
	$html_table_2 .= qq{<thead style="text-align:center;">};
	$html_table_2 .= qq{<th data-field="plot"><b><center>Sashimi Plot</center></b></th>};
	$html_table_2 .= qq{<th data-field="igv"><b><center>IGV</center></b></th>};
	$html_table_2 .= qq{<th data-sortable='true' data-field="var_name"><b><center>Var name</center></b></th>};
	$html_table_2 .= qq{<th data-field="trio"><b><center>Trio</center></b></th>};
	$html_table_2 .= qq{<th data-field="deja_vu"><b><center>DejaVu (Nb Samples)</center></b></th>};
	$html_table_2 .= qq{<th data-field="transcripts"><b><center>Transcripts</center></b></th>};
	$html_table_2 .= qq{<th data-sortable='true' data-field="jscore"><b><center>Score</center></b></th>};
	$html_table_2 .= qq{<th data-field="noise"><b><center>Details Score</center></b></th>};
	$html_table_2 .= qq{</thead>};
	$html_table_2 .= qq{<tbody>};
	foreach my $l (@{$_tr_lines_by_genes->{$gene_name}}) {
		$html_table_2 .= $l;
	}
	$html_table_2 .= qq{</tbody>};
	$html_table_2 .= qq{</table>};
	$html_table_2 .= qq{</div>};
	$html_table_2 .= qq{</td></tr>};
	$nb_g++;
}
if ($html_table_2) {
	my $this_html = $html;
	$this_html .= $html_table_1;
	$this_html .= $html_table_2;
	$this_html .= qq{</tbody>};
	$this_html .= qq{</table>};
	push(@l_html, $this_html);
	push(@l_tables_ids, join(',', @lTablesIds));
}


#$no_cnv->close();
$project->dejavuJunctions->close() if ($release =~ /HG19/);


my $hash;
$hash->{html} = $l_html[0];
$hash->{tables_id} = $l_tables_ids[0];
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
	
	my $dv_other_pat = $junction->dejavu_nb_other_patients($patient);
	my $dv_other_pat_ratio_10 = 0;
	my $dv_other_pat_ratio_20 = 0;
	if ($dv_other_pat > 0) {
		$dv_other_pat_ratio_10 = $junction->dejavu_nb_other_patients_min_ratio_10($patient);
		$dv_other_pat_ratio_20 = $junction->dejavu_nb_other_patients_min_ratio_20($patient);
	}
	
	$html.= $cgi->start_Tr();
	$html.= $cgi->td("<center><b>DejaVu</b></center>");
	$html.= $cgi->td(obutton($cmd_all, $dv_other_pat));
#	if ($my_ratio >= 70) {
#		$html.= $cgi->td(obutton($cmd_all, $junction->dejavu_nb_other_patients_min_ratio_70($patient)));
#		$html.= $cgi->td(obutton($cmd_all, $junction->dejavu_nb_other_patients_min_ratio_90($patient)));
#	}
#	else {
		$html.= $cgi->td(obutton($cmd_all, $dv_other_pat_ratio_10));
		$html.= $cgi->td(obutton($cmd_all, $dv_other_pat_ratio_20));
#	}
	$html.= $cgi->end_Tr();
	
	my $dv_run_other_pat = $junction->dejavu_nb_int_this_run_patients($patient);
	my $dv_run_other_pat_ratio_10 = 0;
	my $dv_run_other_pat_ratio_20 = 0;
	if ($dv_run_other_pat > 0) {
		$dv_run_other_pat_ratio_10 = $junction->dejavu_nb_int_this_run_patients($patient, 10);
		$dv_run_other_pat_ratio_20 = $junction->dejavu_nb_int_this_run_patients($patient, 20);
	}
	
#	foreach my $pheno_name (@{$patient->getProject->phenotypes()}) {
#		$html.= $cgi->start_Tr();
#		$html.= $cgi->td("<center><b>$pheno_name</b></center>");
#		my $dv_other_pat_pheno = $junction->dejavu_nb_other_patients_phenotype($patient, $pheno_name);
#		my $dv_other_pat_pheno_ratio_10 = 0;
#		my $dv_other_pat_pheno_ratio_20 = 0;
#		if ($dv_other_pat_pheno > 0) {
#			$dv_other_pat_pheno_ratio_10 = $junction->dejavu_nb_other_patients_phenotype_min_ratio_10($patient, $pheno_name);
#			$dv_other_pat_pheno_ratio_20 = $junction->dejavu_nb_other_patients_phenotype_min_ratio_20($patient, $pheno_name);
#		}
#		$html.= $cgi->td(obutton($cmd_all, $dv_other_pat_pheno));
#		$html.= $cgi->td(obutton($cmd_all, $dv_other_pat_pheno_ratio_10));
#		$html.= $cgi->td(obutton($cmd_all, $dv_other_pat_pheno_ratio_20));
#		$html.= $cgi->end_Tr();
#	}
	
	$html.= $cgi->start_Tr();
	$html.= $cgi->td("<center><b>InThisProject</b></center>");
	$html.= $cgi->td(obutton($cmd_inthisrun, $dv_run_other_pat));
#	if ($my_ratio >= 70) {
#		$html.= $cgi->td(obutton($cmd_inthisrun,$junction->dejavu_nb_int_this_run_patients($patient,70)));
#		$html.= $cgi->td(obutton($cmd_inthisrun,$junction->dejavu_nb_int_this_run_patients($patient,90)));
#	}
#	else {
		$html.= $cgi->td(obutton($cmd_inthisrun, $dv_run_other_pat_ratio_10));
		$html.= $cgi->td(obutton($cmd_inthisrun, $dv_run_other_pat_ratio_20));
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
	
	my $intspan_junction = Set::IntSpan::Fast::XS->new();
	$intspan_junction->add_range($junction->start()-25, $junction->start()+25);
	$intspan_junction->add_range($junction->end()-25, $junction->end()+25);
	
#	my (@l_cov_all_pat, @l_cov_all_pat_gene);
#	foreach my $pat (@{$project->getPatients()}) {
#		foreach my $pcov (@{$pat->depthIntspan($junction->getChromosome->id(), $intspan_junction)}) {
#			push(@l_cov_all_pat, $pcov);
#		}
#		foreach my $pcov (@{$pat->depth($gene->getChromosome->id(), $gene->start(), $gene->end())}) {
#			push(@l_cov_all_pat_gene, $pcov);
#		}
#		
#	}
#	my $cov_mean_all_pat = 0;
#	foreach my $cov (@l_cov_all_pat) { $cov_mean_all_pat += $cov; }
#	$cov_mean_all_pat = $cov_mean_all_pat / scalar(@l_cov_all_pat);
#	
#	my $cov_mean_all_pat_gene = 0;
#	foreach my $cov (@l_cov_all_pat_gene) { $cov_mean_all_pat_gene += $cov; }
#	$cov_mean_all_pat_gene = $cov_mean_all_pat_gene / scalar(@l_cov_all_pat_gene);
    	
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
#		my @l_cov = @{$pat->depthIntspan($junction->getChromosome->id(), $intspan_junction)};
#		my $cov_mean = 0;
#		foreach my $pcov (@l_cov) { $cov_mean += $pcov; }
#		$cov_mean = $cov_mean / scalar(@l_cov);
		
#		my $cov_norm_mean = 0;
#		my @lcov_norm = @{$patient->normalize_depth($gene->getChromosome->id(), $gene->start(), $gene->end())};
#		foreach my $pcov (@lcov_norm) { $cov_norm_mean += $pcov; }
#		$cov_norm_mean = $cov_norm_mean / scalar(@lcov_norm);
		
#		$h_by_pat->{$fam_name}->{$pat->name()}->{dp_pat_gene} = sprintf("%.2f",$cov_mean);
#		$h_by_pat->{$fam_name}->{$pat->name()}->{dp_all_pat_gene} = sprintf("%.2f",$cov_mean_all_pat);
#		$h_by_pat->{$fam_name}->{$pat->name()}->{dp_pat_normalize_gene} = sprintf("%.2f",$cov_norm_mean);
#		$h_by_pat->{$fam_name}->{$pat->name()}->{dp_pat_mean_gene} = sprintf("%.2f",$patient->meanDepth($gene->getChromosome->id(), $gene->start(), $gene->end()));
#		$h_by_pat->{$fam_name}->{$pat->name()}->{dp_all_pat_mean_gene} = sprintf("%.2f",$cov_mean_all_pat_gene);
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
#	$html_patients .= qq{<th data-field="dp_pat_gene"><b><center>DP Pat.</center></b></th>};
#	$html_patients .= qq{<th data-field="dp_pat_gene"><b><center>DP All Pat.</center></b></th>};
#	$html_patients .= qq{<th data-field="dp_mean_gene"><b><center>DP Normalize Pat Gene</center></b></th>};
#	$html_patients .= qq{<th data-field="dp_mean_gene"><b><center>DP Mean Pat Gene</center></b></th>};
#	$html_patients .= qq{<th data-field="dp_mean_gene"><b><center>DP Mean All Pat Gene</center></b></th>};
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
#			$html_patients .= qq{<td>}.$h_by_pat->{$fam_name}->{$pat_name}->{dp_pat_gene}.qq{</td>};
#			$html_patients .= qq{<td>}.$h_by_pat->{$fam_name}->{$pat_name}->{dp_all_pat_gene}.qq{</td>};
#			$html_patients .= qq{<td>}.$h_by_pat->{$fam_name}->{$pat_name}->{dp_pat_normalize_gene}.qq{</td>};
#			$html_patients .= qq{<td>}.$h_by_pat->{$fam_name}->{$pat_name}->{dp_pat_mean_gene}.qq{</td>};
#			$html_patients .= qq{<td>}.$h_by_pat->{$fam_name}->{$pat_name}->{dp_all_pat_mean_gene}.qq{</td>};
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
	
	my $fasta = "";
	if (not $patient->getProject->is_human_genome()) {
		$fasta = $patient->getProject->genomeFasta();
		my $release = $patient->getProject->getVersion();
		my $fasta_named = $fasta;
		$fasta_named =~ s/all/$release/;
		if (-e $fasta_named) { $fasta = $fasta_named; }
		$fasta =~ s/\/data-isilon//;
		$fasta = "https://www.polyweb.fr/".$fasta;
	}
	
	my $igv_link = qq{<button class='igvIcon2' onclick='launch_igv_tool_rna("$fasta", "$bam_file,$gtf","$locus")' style="color:black"></button>};
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
#		my $vid = $junction->vector_id();
#		$sashimi_button .= qq{<center>$vid</center>};
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

sub add_linked_hash_in_cache {
	my ($cache_id, $patient) = @_;
	my $h_vector;
	foreach my $chr (@{$patient->getProject->getChromosomes()}) {
		my $h_vector_chr;
		my $vector_patient = $patient->getJunctionsVector($chr);
		$h_vector_chr->{min0} = $vector_patient->Clone();
		$h_vector_chr->{min2} = $vector_patient->Clone();
		$h_vector_chr->{min4} = $vector_patient->Clone();
		$h_vector_chr->{min6} = $vector_patient->Clone();
		$h_vector_chr->{min8} = $vector_patient->Clone();
		foreach my $junction (@{$chr->getListVarObjects($vector_patient)}) {
			$n++;
			next if ($junction->isCanonique($patient));
			print '.' if ($n % 5000);
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
			my $score_this_j = $junction->junction_score_without_dejavu_global($patient);
			if ($score_this_j < 0) {
				$h_vector_chr->{min0}->Bit_Off($junction->vector_id());
				$h_vector_chr->{min2}->Bit_Off($junction->vector_id());
				$h_vector_chr->{min4}->Bit_Off($junction->vector_id());
				$h_vector_chr->{min6}->Bit_Off($junction->vector_id());
				$h_vector_chr->{min8}->Bit_Off($junction->vector_id());
			}
			if ($score_this_j < 2) {
				$h_vector_chr->{min2}->Bit_Off($junction->vector_id());
				$h_vector_chr->{min4}->Bit_Off($junction->vector_id());
				$h_vector_chr->{min6}->Bit_Off($junction->vector_id());
				$h_vector_chr->{min8}->Bit_Off($junction->vector_id());
			}
			if ($score_this_j < 4) {
				$h_vector_chr->{min4}->Bit_Off($junction->vector_id());
				$h_vector_chr->{min6}->Bit_Off($junction->vector_id());
				$h_vector_chr->{min8}->Bit_Off($junction->vector_id());
			}
			if ($score_this_j < 6) {
				$h_vector_chr->{min6}->Bit_Off($junction->vector_id());
				$h_vector_chr->{min8}->Bit_Off($junction->vector_id());
			}
			if ($score_this_j < 8) {
				$h_vector_chr->{min8}->Bit_Off($junction->vector_id());
			}
		}
		$h_vector->{$chr->id()}->{min0} = $h_vector_chr->{min0}->to_Enum();
		$h_vector->{$chr->id()}->{min2} = $h_vector_chr->{min2}->to_Enum();
		$h_vector->{$chr->id()}->{min4} = $h_vector_chr->{min4}->to_Enum();
		$h_vector->{$chr->id()}->{min6} = $h_vector_chr->{min6}->to_Enum();
		$h_vector->{$chr->id()}->{min8} = $h_vector_chr->{min8}->to_Enum();
		print '.chr'.$chr->id.':'.$chr->countThisVariants($h_vector_chr->{min0}).'|'.$chr->countThisVariants($h_vector_chr->{min2}).'|'.$chr->countThisVariants($h_vector_chr->{min4}).'|'.$chr->countThisVariants($h_vector_chr->{min6}).'|'.$chr->countThisVariants($h_vector_chr->{min8}).'.';
	}
	my $no_cache = $patient->get_lmdb_cache("w");
	$no_cache->put_cache_hash($cache_id, $h_var_linked_ids);
	$no_cache->put_cache_hash($cache_id.'_chr_vectors_enum', $h_vector);
	$no_cache->close();
}

sub countMinCatToUse {
	my ($h_res_v_enum) = @_;
	my $hCount;
	foreach my $chr_id (keys %$h_res_v_enum) {
		foreach my $cat (keys %{$h_res_v_enum->{$chr_id}}) {
			my $v_filters = $project->getChromosome($chr_id)->getNewVector();
			$v_filters->from_Enum($h_res_v_enum->{$chr_id}->{$cat});
			$hCount->{$cat} += $project->getChromosome($chr_id)->countThisVariants($v_filters);
		}
	}
	my @lCat = sort keys %{$hCount};
	foreach my $cat (@lCat) {
		return $cat if $hCount->{$cat} <= 50000;
	}
	return $lCat[-1];
}
	


sub printJson {
	my ($hashRes) = @_;
	my $json_encode = encode_json $hashRes;
	print ".\",";
	$json_encode =~ s/{//;
	print $json_encode;
	exit(0);
}
