#!/usr/bin/perl
$|=1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo";
use lib "$Bin/../../GenBo/lib/obj-nodb";
use lib "$Bin/../../GenBo/lib/obj-nodb/polyviewer/";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";


use connect;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use JSON;
use xls_export;
use session_export;
use Parallel::ForkManager;
use Compress::Snappy;
use Storable qw/thaw freeze/;
use List::MoreUtils qw{ natatime };
use PolyviewerJunction;
use polysplice_html;

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
my $view_polyviewer = $cgi->param('view_polyviewer');

#$only_positions = '10:17000000-17080000';

$min_score = 20 unless (defined($min_score));

my $nb_group_junctions_colors = 0;
my $fork = 5;

my $buffer = GBuffer->new;
my $project = $buffer->newProjectCache( -name => $project_name );
$project->getChromosomes();
my $patient = $project->getPatient($use_patient);
$patient->use_not_filtred_junction_files(1);
my $only_gene;
$only_gene = $project->newGene($only_gene_name) if ($only_gene_name);

if ($view_polyviewer) {
	print $cgi->header;
	print qq{<div hidden>};
}
else {
	print $cgi->header('text/json-comment-filtered');
	print "{\"progress\":\".";
}

my $patient_name = $patient->name();
my $html;
my $release = $project->getVersion();
my $gencode = '-';
$gencode = $project->gencode_version() if ($release =~ /HG19/);

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

my $percent_dejavu = 0;
if (defined($use_percent_dejavu)) {
	$percent_dejavu = $use_percent_dejavu - 90;
}
my $nb_percent_dejavu_value = 90 + $percent_dejavu;

my (@lJunctions, $h_varids, $h_var_linked_ids, $h_genes_by_score);
my $pm1 = new Parallel::ForkManager($fork);
$pm1->run_on_finish(
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $h) = @_;
		my $chr_id_done = $h->{chr_id};
		unless (defined($h) or $exit_code > 0) {
			print "\nERROR for $chr_id_done...";
			print qq|No message received from child process $exit_code $pid!\n|;
			die();
		}
		foreach my $j (@{$h->{list_junctions}}) {
			push(@lJunctions, $j);
			if (exists $h_genes_by_score->{$j->ensid()}) {
				$h_genes_by_score->{$j->ensid()} = $j->score() if ($j->score() > $h_genes_by_score->{$j->ensid()});
			}
			else {
				$h_genes_by_score->{$j->ensid()} = $j->score();
			}
		}
		foreach my $jid (keys %{$h->{hash_varids}}) { $h_varids->{$jid} = $h->{hash_varids}->{$jid}; }
		foreach my $jid (keys %{$h->{hash_var_linked}}) { $h_var_linked_ids->{$jid} = $h->{hash_var_linked}->{$jid}; }
	}
);


$buffer->dbh_deconnect();
foreach my $chr (@{$project->getChromosomes()}) {
	next if ($only_gene and $only_gene->getChromosome->id() ne $chr->id());
	my $pid = $pm1->start and next;
	$buffer->dbh_reconnect();
	#next if $chr->id ne '4';
	my (@lJunctions_tmp, $h_varids_tmp, $h_var_linked_ids_tmp);
	my $res;
	$res->{chr_id} = $chr->id();
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
	if ($only_gene) {
		my $vector_postions = $chr->getVectorByPosition($only_gene->start(), $only_gene->end());
		$vector_patient->Intersection($vector_patient, $vector_postions);
	}
	foreach my $junction (@{$chr->getListVarObjects($vector_patient)}) {
		$n++;
		print '.' if ($n % 50000);
		
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
		eval {
			$junction->dejavu_percent_coordinate_similar($nb_percent_dejavu_value);
			my $nb_dejavu_pat = $junction->dejavu_nb_others_patients();
			$nb_dejavu_pat = $junction->dejavu_nb_other_patients_min_ratio_10($patient) if ($only_dejavu_ratio_10);
			next if ($nb_dejavu_pat > $max_dejavu_value);
		};
		if ($@) {  next; }
		my $vp = PolyviewerJunction->new();
		$vp->setJunction($junction,$patient,$min_score);
		push(@lJunctions_tmp, $vp);
	}
	$res->{list_junctions} = \@lJunctions_tmp;
	$res->{hash_varids} = $h_varids_tmp;
	$res->{hash_var_linked} = $h_var_linked_ids_tmp;
	$pm1->finish(0, $res);
}
$pm1->wait_all_children();
$buffer->dbh_reconnect();

my $print_html = polysplice_html->new( project => $project, patient => $patient );
$print_html->init();

$print_html->set_view_polyviewer($view_polyviewer);
$print_html->set_filter_check_only_dejavu_ratio_10($only_dejavu_ratio_10);
$print_html->set_filter_only_gene_name($only_gene_name);
$print_html->set_filter_only_positions($only_positions);
$print_html->set_filter_percent_similar_dejavu($use_percent_dejavu);
$print_html->set_filter_max_dejavu($max_dejavu);
$print_html->set_filter_min_score($min_score);

if (not $view_polyviewer) {
	$html .= "<br>";
	$html .= $print_html->get_html_header_panel_patient_filters();
}

my $h_junctions_by_genes;
foreach my $junction (@lJunctions) {
	my $ensid = $junction->ensid();
	$h_junctions_by_genes->{$ensid}->{$junction->id()} = $print_html->print_html_line_in_polysplice($junction);
}

my $h_genes_scoring;
foreach my $ensid (keys %$h_genes_by_score) {
	my $score = $h_genes_by_score->{$ensid};
	$h_genes_scoring->{$score}->{$ensid} = undef;
}

my $table_id = 'table_'.$patient_name;

if (not $h_genes_scoring) {
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
elsif ($only_gene_name and not $only_gene) {
	$html .= qq{<br><tr><td><b><i><span class="glyphicon glyphicon-exclamation-sign" style="color:red;font-size:18px;"></span> your gene <span style="color:orange;">$only_gene_name</span> is unknown in my database...</b></i></td></tr><br><br><br>};
}

my @lScores = sort {$b <=> $a} %$h_genes_scoring;
foreach my $score (@lScores) {
	foreach my $ensid (sort keys %{$h_genes_scoring->{$score}}) {
		
		my @list_html_lines_junctions;
		foreach my $jid (sort keys %{$h_junctions_by_genes->{$ensid}}) {
			push(@list_html_lines_junctions, $h_junctions_by_genes->{$ensid}->{$jid});
		}
		my $this_panel_gene_id = 'panel_'.$table_id.'_'.$ensid;
		my $this_table_id = $table_id.'_'.$ensid;
		my $g = $project->newGene($ensid);
		my $hgene;
		$hgene->{description} = $g->description();
		$hgene->{external_name} = $g->external_name();
		$hgene->{max_score} = $score;
		$hgene->{js_id} = $g->id();
		$hgene->{uid} = $g->id();
		$hgene->{id} = $g->id();
		$hgene->{omim_id} = $g->omim_id();
		$hgene->{pLI} = $g->pLI();
		$hgene->{omim_inheritance} = $g->omim_inheritance();
		$hgene->{variants} = \@list_html_lines_junctions;
		foreach my $panel (@{$g->getPanels()}) {
			$hgene->{panels}->{$panel->name()} = undef;
		}
		my $div_id = 'div_'.$this_table_id;
		$hgene->{collapse_with_id} = $div_id;
		my $class_tr_gene->{style} = "background-color:#607D8B;border:1px black solid;padding:0px;white-space: nowrap;height:50px;";
		$html .= qq{<tr>};
		$html .= "<table style='width:100%;background-color:#F3F3F3;'>";
		my $html_gene_panel = $print_html->panel_gene($hgene, $this_panel_gene_id,$project->name(), $patient);
		$html_gene_panel =~ s/<b>pLI<\/b> Score//;
		$html_gene_panel =~ s/bottom:5px;/bottom:-5px;/;
		$html .= $cgi->td($class_tr_gene, $html_gene_panel);
		$html .= "</table>";
		$html .= $cgi->end_Tr();
		$html .= qq{<tr>};
		$html .= qq{<div class="collapse" id="$div_id">};
		$html .= qq{<table id='$this_table_id' data-sort-name='locus' data-sort-order='desc' data-filter-control='true' data-toggle="table" data-pagination-v-align="both" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="10" data-page-list="[10]" data-resizable='true' class='table table-striped' style='font-size:11px;'>};
		$html .= $print_html->print_html_header_in_polysplice();
		$html .= qq{<tbody>};
		$html .= join('', @list_html_lines_junctions);
		$html .= qq{</tbody>};
		$html .= qq{</table>};
		$html .= qq{</div>};
		$html .= $cgi->end_Tr();
		$html .= "<br>";
	}
}
$html .= "</table>";


if ($view_polyviewer) {
	print qq{</div>};
	print $html;
	exit(0);
}

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





sub get_igv_button {
	my ($gene, $patient) = @_;
	my $color = 'lightgrey';
	my $bam_file = $patient->bamUrl();
	my $list_patients_ctrl = $patient->getPatients_used_control_rna_seq_junctions_analyse();
	if ($list_patients_ctrl) {
		my $nb_control;
		foreach my $other_pat (@$list_patients_ctrl) {
			$bam_file .= ','.$other_pat->bamUrl();
			$nb_control++;
			last if $nb_control == 3;
		}
	}
	else {
		my $np = 0;
		foreach my $other_pat (@{$project->getPatients()}) {
			next if ($other_pat->name() eq $patient->name());
			$bam_file .= ','.$other_pat->bamUrl();
			$np++;
			last if $np == 3;
		}
	}
	my $gtf = $patient->getProject->get_gtf_genes_annotations_igv();
	$gtf =~ s/\/data-isilon//;
	$gtf = 'https://'.$ENV{HTTP_HOST}.'/'.$gtf;
	my $locus = $gene->getChromosome->id().':'.($gene->start()-100).'-'.($gene->end()+100);
	my $igv_link = qq{<button class='igvIcon2' onclick='launch_igv_tool_rna("", "$bam_file,$gtf","$locus")' style="color:black"></button>};
	return $igv_link;
}

sub printJson {
	my ($hashRes) = @_;
	my $json_encode = encode_json $hashRes;
	print ".\",";
	$json_encode =~ s/{//;
	print $json_encode;
	exit(0);
}
