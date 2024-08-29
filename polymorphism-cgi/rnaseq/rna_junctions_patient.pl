#!/usr/bin/perl
$| = 1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo";
use lib "$Bin/../../GenBo/lib/obj-nodb";
use lib "$Bin/../../GenBo/lib/obj-nodb/polyviewer/";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
require
  "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";

use connect;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use JSON;
use xls_export;
use session_export;
use List::MoreUtils qw{ natatime };
use Parallel::ForkManager;
use PolyviewerJunction;
use polysplice_html;

my $cgi                  = new CGI;
my $project_name         = $cgi->param('project');
my $use_patient          = $cgi->param('patient');
my $view_all_junctions   = $cgi->param('view_not_filtred_junctions');
my $isChrome             = $cgi->param('is_chrome');
my $max_dejavu           = $cgi->param('dejavu');
my $use_percent_dejavu   = $cgi->param('dejavu_percent');
my $min_score            = $cgi->param('min_score');
my $only_gene_name       = $cgi->param('only_gene');
my $only_positions       = $cgi->param('only_positions');
my $only_dejavu_ratio_10 = $cgi->param('only_dejavu_ratio_10');
my $only_html_cache      = $cgi->param('only_html_cache');
my $view_polyviewer 	 = $cgi->param('view_polyviewer');
my $export_xls			 = $cgi->param('export_xls');
my $session_id			 = $cgi->param('session_id');

my $only_junctions_NDA   = $cgi->param('only_junctions_NDA');
my $only_junctions_DA    = $cgi->param('only_junctions_DA');
my $only_junctions_A     = $cgi->param('only_junctions_A');
my $only_junctions_D     = $cgi->param('only_junctions_D');
my $only_junctions_N     = $cgi->param('only_junctions_N');


if ($export_xls and $session_id) {
	my $xls_export = new xls_export();
	$xls_export->load($session_id);
	my ($list_datas_junctions) = $xls_export->prepare_generic_datas_junctions();
	$xls_export->add_page('Junctions', $xls_export->list_generic_header_junctions(), $list_datas_junctions);
	$xls_export->export();
	exit(0);
}

#$only_positions = '10:17000000-17080000';

$min_score = 20 unless ( defined($min_score) );

my $nb_group_junctions_colors = 0;

my $buffer  = GBuffer->new;
my $project = $buffer->newProjectCache( -name => $project_name );
$project->getChromosomes();
my $patient = $project->getPatient($use_patient);
$patient->use_not_filtred_junction_files(1);
my $only_gene;
$only_gene = $project->newGene($only_gene_name) if ($only_gene_name);

if ( not $only_html_cache ) {
	print $cgi->header('text/json-comment-filtered');
	if ($view_polyviewer) {
		print "<div style='display:none;'>.";
	}
	else {
		print "{\"progress\":\".";
	}
}

my $patient_name = $patient->name();
my $table_id     = 'table_' . $patient_name;
my $html;
my $release = $project->getVersion();
my $gencode = '-';
if    ( $release =~ /HG19/ ) { $gencode = $project->gencode_version(); }
elsif ( $release eq 'MM38' ) { $gencode = 'M25'; }
elsif ( $release eq 'MM39' ) { $gencode = 'M32'; }

my $max_dejavu_value = 51;
$max_dejavu_value = $max_dejavu if defined($max_dejavu);

my ( $default_filter_dejavu_project, $default_filter_score );
my ( $nb_max_patient,                $nb_limit );
my $h_desc = $project->get_hash_patients_description_rna_seq_junction_analyse();
if ($h_desc) {
	foreach my $pat_name ( keys %$h_desc ) {
		$nb_max_patient++ if ( exists $h_desc->{$pat_name}->{pat} );
	}
}
else { $nb_max_patient = scalar( @{ $project->getPatients() } ); }
if    ( $nb_max_patient == 1 ) { $nb_limit = 0; }
elsif ( $nb_max_patient <= 3 ) { $nb_limit = 1; }
else                           { $nb_limit = int( $nb_max_patient / 4 ); }
$default_filter_dejavu_project = "data-filter-default='<=$nb_limit'"
  if $nb_limit > 0;
$default_filter_score = "data-filter-default='>=1'";

my $n            = 0;
my $score_slider = 0;
$score_slider = $min_score / 10 if ( $min_score and $min_score > 0 );
my ( @lJunctions, $h_var_linked_ids );

my $j_total    = 0;
my $j_selected = 0;
my $h_chr_vectors;
my $h_chr_vectors_counts;

my $has_regtools_vectors;
my $exists_vector_ratio_10;
foreach my $chr ( @{ $project->getChromosomes() } ) {
#	warn $chr->id;
#	$chr->global_categories;
	
	#next if $chr->id ne '1' and $patient_name eq 'KUCerc';
	my $vector_patient = $patient->getJunctionsVector($chr);
	$j_total += $chr->countThisVariants($vector_patient);
	if ($only_positions) {
		my @lTmp         = split( ':', $only_positions );
		my $chr_id_limit = lc( $lTmp[0] );
		$chr_id_limit =~ s/chr//;
		next if lc( $chr->id ) ne $chr_id_limit;
		my ( $from, $to ) = split( '-', $lTmp[1] );
		my $vector_postions = $chr->getVectorByPosition( $from, $to );
		$vector_patient->Intersection( $vector_patient, $vector_postions );
	}
	if ($only_gene) {
		next if ( $chr->id() ne $only_gene->getChromosome->id() );
		my $vector_postions =
		  $chr->getVectorByPosition( ( $only_gene->start - 1000000 ),
			( $only_gene->end + 1000000 ) );
		$vector_patient->Intersection( $vector_patient, $vector_postions );
	}

	$h_chr_vectors->{ $chr->id } = $vector_patient->Clone();
	$h_chr_vectors_counts->{ $chr->id } = $chr->countThisVariants( $h_chr_vectors->{ $chr->id } );
	$j_selected += $h_chr_vectors_counts->{ $chr->id };
	$exists_vector_ratio_10++ if exists $chr->patients_categories->{$patient->name().'_ratio_10'};
	
	if ($h_chr_vectors_counts->{ $chr->id } > 0) {
		$has_regtools_vectors = 1 if exists $chr->global_categories->{'N'};
		$has_regtools_vectors = 1 if exists $chr->global_categories->{'D'};
		$has_regtools_vectors = 1 if exists $chr->global_categories->{'A'};
		$has_regtools_vectors = 1 if exists $chr->global_categories->{'DA'};
		$has_regtools_vectors = 1 if exists $chr->global_categories->{'NDA'};
	}
}
my $cache_id = 'splices_linked_' . $patient->name();
my $no_cache = $patient->get_lmdb_cache("r");
my $h_res    = $no_cache->get_cache($cache_id);
$no_cache->close();
unless ($h_res) {
	$h_var_linked_ids = $h_res;
}

my $cache_html_id = 'splices_HTML_' . $patient->name();
$cache_html_id .= '_maxdv'.$max_dejavu if $max_dejavu;
$cache_html_id .= '_mins'.$min_score if $min_score;
$cache_html_id .= '_only'.$only_gene_name if $only_gene_name;
$cache_html_id .= '_only'.$only_positions if $only_positions;
$cache_html_id .= '_onlydvra10' if $only_dejavu_ratio_10;
$cache_html_id .= '_onlyNDA' if $has_regtools_vectors and $only_junctions_NDA;
$cache_html_id .= '_onlyDA' if $has_regtools_vectors and $only_junctions_DA;
$cache_html_id .= '_onlyA' if $has_regtools_vectors and $only_junctions_A;
$cache_html_id .= '_onlyD' if $has_regtools_vectors and $only_junctions_D;
$cache_html_id .= '_onlyN' if $has_regtools_vectors and $only_junctions_N;
$cache_html_id .= '_'.$j_total;
if (not $only_gene_name and not $only_positions and not $view_polyviewer and not $export_xls) {
	my $no_cache_2 = $patient->get_lmdb_cache("r");
	my $h_html    = $no_cache_2->get_cache($cache_html_id);
	
	#TODO: supress cache html
	#$h_html=undef;
	
	$no_cache_2->close();
	if ($h_html) {
		if ($view_polyviewer) {
			print qq{</div>};
			print $h_html->{html};
			exit(0);
		}
		printJson($h_html);
		exit(0);
	}
}

my ( $is_partial_results, $use_cat, $min_partial_score );
my $no_cache = $patient->get_lmdb_cache("r");
my $cache_vectors_enum_id = $patient->name() . '_' . '_chr_vectors_enum';
my $h_res_v_enum = $no_cache->get_cache($cache_vectors_enum_id);


#warn $j_total;
#die;
#
##$min_partial_score = 0 if ()
#$use_cat = countMinCatToUse($h_res_v_enum);
#print ".$use_cat.";
#foreach my $chr_id ( keys %$h_res_v_enum ) {
#	my $v_filters = $project->getChromosome($chr_id)->getNewVector();
#	$v_filters->from_Enum( $h_res_v_enum->{$chr_id}->{$use_cat} );
#	$h_chr_vectors->{$chr_id}->Intersection( $h_chr_vectors->{$chr_id}, $v_filters );
#}
#$is_partial_results = 1;
#$min_partial_score  = $use_cat;
#$min_partial_score =~ s/min//;
#$min_partial_score -= 4;
$min_partial_score = 0 if $j_total > 20000;

my $percent_dejavu = 0;
if ( defined($use_percent_dejavu) ) {
	$percent_dejavu = $use_percent_dejavu - 90;
}
my $nb_percent_dejavu_value = 90 + $percent_dejavu;
my $nbErrors = 0;


my ($checked_only_dejavu_ratio_10, $checked_regtools_nda, $checked_regtools_d, $checked_regtools_a, $checked_regtools_n, $checked_regtools_da);
$checked_only_dejavu_ratio_10 = qq{checked="checked"} if $only_dejavu_ratio_10;
$checked_regtools_nda = qq{checked="checked"} if $only_junctions_NDA;
$checked_regtools_da = qq{checked="checked"} if $only_junctions_DA;
$checked_regtools_d = qq{checked="checked"} if $only_junctions_D;
$checked_regtools_a = qq{checked="checked"} if $only_junctions_A;
$checked_regtools_n = qq{checked="checked"} if $only_junctions_N;

if (not $view_polyviewer) {
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
	
	my $html_gene_select =
	qq{<input type="text" style="font-size:11px;width:180px;" placeholder="COL4A5, 1:100000-150000" class="form-control" id="input_gene_positions">};
	if ($only_gene_name) {
		$html_gene_select =
	qq{<input type="text" style="font-size:11px;width:180px;" value="$only_gene_name" class="form-control" id="input_gene_positions">};
	}
	elsif ($only_positions) {
		$html_gene_select =
	qq{<input type="text" style="font-size:11px;width:180px;" value="$only_positions" class="form-control" id="input_gene_positions">};
	}
	
	my $html_filters = qq{
		<table style="width:100%;">
			<tr>
	};
	
	if ($has_regtools_vectors) {
		$html_filters .= qq{
			<td style="padding-top:5px;">
				<center>
					<table>
						<tr>
							<td style="padding-left:3px;"> <div class="form-check"><input $checked_regtools_nda class="form-check-input" type="checkbox" value="" id="b_regtools_nda"><label class="form-check-label" for="b_regtools_nda" style="padding-left:10px;font-size:11px;">NDA</label></div> </td>
							<td style="padding-left:3px;"> <div class="form-check"><input $checked_regtools_d class="form-check-input" type="checkbox" value="" id="b_regtools_d"><label class="form-check-label" for="b_regtools_d" style="padding-left:10px;font-size:11px;">D</label></div> </td>
						</tr>
						<tr>
							<td style="padding-left:3px;"> <div class="form-check"><input $checked_regtools_a class="form-check-input" type="checkbox" value="" id="b_regtools_a"><label class="form-check-label" for="b_regtools_a" style="padding-left:10px;font-size:11px;">A</label></div> </td>
							<td style="padding-left:3px;"> <div class="form-check"><input $checked_regtools_n class="form-check-input" type="checkbox" value="" id="b_regtools_n"><label class="form-check-label" for="b_regtools_n" style="padding-left:10px;font-size:11px;">N</label></div> </td>
						</tr>
					</table>
				</center>
			</td>
		};
	}
	
	$html_filters .= qq{
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
	};
	
	if ($has_regtools_vectors) {
		$html_filters .= qq{
			<td style="padding-top:5px;font-size:11px;"><center><b>Filter Category</center></b></td>
		};
	}
	
	$html_filters .= qq{
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
}

my $_tr_lines_by_genes;
my $h_junctions_color;
my $h_dejavu_cnv;
$n = 0;
my ($h_junctions_scores, $h_export_xls);
my $fork     = 5;
my $pm       = new Parallel::ForkManager($fork);
my $nbErrors = 0;
$pm->run_on_finish(
	sub {
		my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $hres ) = @_;
		unless ( defined($hres) or $exit_code > 0 ) {
			$nbErrors++;
			print qq|No message received from child process $exit_code $pid!\n|;
			return;
		}
		
		
		if ($export_xls) {
			my $chr_id = $hres->{chr_id}; 
			delete $hres->{done};
			delete $hres->{chr_id};
			$h_export_xls->{$chr_id} = $hres;
		}
		else {
			foreach my $gene_name ( keys %{ $hres->{genes} } ) {
				foreach my $score ( keys %{$hres->{genes}->{$gene_name}} ) {
					foreach my $html_tr ( @{ $hres->{genes}->{$gene_name}->{$score} } ) {
						push( @{ $_tr_lines_by_genes->{$gene_name}->{$score} }, $html_tr );
					}
				}
			}
	
			foreach my $gene_name ( keys %{ $hres->{score}->{all} } ) {
				foreach my $score ( keys %{ $hres->{score}->{all}->{$gene_name} } )
				{
					$h_junctions_scores->{all}->{$gene_name}->{$score} = undef;
				}
			}
		}
	}
);


print '_save_xls_' if ($export_xls);

foreach my $chr_id ( sort keys %{$h_chr_vectors} ) {
	#next if $chr_id ne '20';
	
	my $chr = $patient->getProject->getChromosome($chr_id);
	next if $chr->countThisVariants( $h_chr_vectors->{$chr_id} ) == 0;

	$patient->getProject->buffer->dbh_deconnect();
	$pm->start and next;
	$patient->getProject->buffer->dbh_reconnect();
	
	if (not $only_gene ) {
		#DV vector
		my $type_vector_dv = 'dejavu';
		my $type_vector_dv_sup = 'dejavu';
		if ($max_dejavu_value >= 90)    {
			$type_vector_dv .= '_90';
		}
		elsif ($max_dejavu_value >= 80) {
			$type_vector_dv .= '_80';
			$type_vector_dv_sup .= '_90';
		}
		elsif ($max_dejavu_value >= 70) {
			$type_vector_dv .= '_70';
			$type_vector_dv_sup .= '_80';
		}
		elsif ($max_dejavu_value >= 60) {
			$type_vector_dv .= '_60';
			$type_vector_dv_sup .= '_70';
		}
		elsif ($max_dejavu_value >= 50) {
			$type_vector_dv .= '_50';
			$type_vector_dv_sup .= '_60';
		}
		elsif ($max_dejavu_value >= 40) {
			$type_vector_dv .= '_40';
			$type_vector_dv_sup .= '_50';
		}
		elsif ($max_dejavu_value >= 30) {
			$type_vector_dv .= '_30';
			$type_vector_dv_sup .= '_40';
		}
		elsif ($max_dejavu_value >= 25) {
			$type_vector_dv .= '_25';
			$type_vector_dv_sup .= '_30';
		}
		elsif ($max_dejavu_value >= 20) {
			$type_vector_dv .= '_20';
			$type_vector_dv_sup .= '_25';
		}
		elsif ($max_dejavu_value >= 15) {
			$type_vector_dv .= '_15';
			$type_vector_dv_sup .= '_20';
		}
		elsif ($max_dejavu_value >= 10) {
			$type_vector_dv .= '_10';
			$type_vector_dv_sup .= '_15';
		}
		else {
			$type_vector_dv .= '_5';
			$type_vector_dv_sup .= '_105';
		}
		$type_vector_dv .= '_r10' if ($only_dejavu_ratio_10);
		$type_vector_dv_sup .= '_r10' if ($only_dejavu_ratio_10);
		if (exists $chr->global_categories->{$type_vector_dv}) {
			$h_chr_vectors->{$chr_id} &= $chr->global_categories->{$type_vector_dv};
#			if ($type_vector_dv_sup and exists $chr->global_categories->{$type_vector_dv_sup}) {
#				my $v_sup = $chr->global_categories->{$type_vector_dv_sup}->Clone();
#				$v_sup -= $chr->global_categories->{$type_vector_dv};
#				$h_chr_vectors->{$chr_id} -= $v_sup;
#			}
		}
	
		# RATIO
		my $type_vector_ratio = $patient->name.'_ratio';
		if ($min_score >= 90)    { $type_vector_ratio .= '_90'; }
		elsif ($min_score >= 80) { $type_vector_ratio .= '_80'; }
		elsif ($min_score >= 70) { $type_vector_ratio .= '_70'; }
		elsif ($min_score >= 60) { $type_vector_ratio .= '_60'; }
		elsif ($min_score >= 50) { $type_vector_ratio .= '_50'; }
		elsif ($min_score >= 40) { $type_vector_ratio .= '_40'; }
		elsif ($min_score >= 30) { $type_vector_ratio .= '_30'; }
		elsif ($min_score >= 20) { $type_vector_ratio .= '_20'; }
		elsif ($min_score >= 10) { $type_vector_ratio .= '_10'; }
		if (exists $chr->patients_categories->{$type_vector_ratio}) {
			$h_chr_vectors->{$chr_id} &= $chr->patients_categories->{$type_vector_ratio};
		}

		if ($has_regtools_vectors) {
			$h_chr_vectors->{$chr_id} -= $chr->global_categories->{'NDA'} if (not $only_junctions_NDA and exists $chr->global_categories->{'NDA'});
			$h_chr_vectors->{$chr_id} -= $chr->global_categories->{'DA'}  if (not $only_junctions_DA and exists $chr->global_categories->{'DA'});
			$h_chr_vectors->{$chr_id} -= $chr->global_categories->{'A'}   if (not $only_junctions_A and exists $chr->global_categories->{'A'});
			$h_chr_vectors->{$chr_id} -= $chr->global_categories->{'D'}   if (not $only_junctions_D and exists $chr->global_categories->{'D'});
			$h_chr_vectors->{$chr_id} -= $chr->global_categories->{'N'}   if (not $only_junctions_N and exists $chr->global_categories->{'N'});
		}
	}
	
	my @lJunctionsChr = @{ $chr->getListVarObjects( $h_chr_vectors->{$chr_id} ) };
	my $size_vector = $h_chr_vectors->{$chr_id}->Size();
	my $hres;
	my $h_td_line;
	
#	my $h_last_j;
	my $h_same_j_description;
	foreach my $junction (@lJunctionsChr) {
		my $pass_same_position;
		my $short_id = $junction->getChromosome->id().'_'.$junction->start().'_'.$junction->end();
		my $min_vid = $junction->vector_id() - 2;
		my $max_vid = $junction->vector_id() + 2;
		foreach my $this_vector_id ($min_vid..$max_vid) {
			next if $this_vector_id < 0;
			next if $this_vector_id >= $size_vector;
			next if $this_vector_id == $junction->vector_id();
			my $this_j = $chr->getVarObject($this_vector_id);
			my $this_short_id = $this_j->getChromosome->id().'_'.$this_j->start().'_'.$this_j->end();
			next if $this_short_id ne $short_id;
			next if $junction->get_nb_new_count($patient) != $this_j->get_nb_new_count($patient);
			if ($junction->get_canonic_count($patient) < $this_j->get_canonic_count($patient)) { $pass_same_position = 1; }
			elsif ($junction->get_canonic_count($patient) == $this_j->get_canonic_count($patient)) {
				if ($this_vector_id > $junction->vector_id()) {
					if ($this_j->isRI($patient)) {
						$h_same_j_description->{$this_vector_id} = $this_j->annex->{$patient->name}->{type_origin_file};
						$h_same_j_description->{$this_vector_id} .= '-'.$this_j->getTypeDescription($patient) if ( $this_j->getTypeDescription($patient) and $this_j->getTypeDescription($patient) ne '---' );
					}
					elsif ($this_j->isSE($patient)) {
						$h_same_j_description->{$this_vector_id} = $this_j->annex->{$patient->name}->{type_origin_file};
						$h_same_j_description->{$this_vector_id} .= '-'.$this_j->getTypeDescription($patient) if ( $this_j->getTypeDescription($patient) and $this_j->getTypeDescription($patient) ne '---' );
					}
					if ($junction->isRI($patient)) {
						$h_same_j_description->{$this_vector_id} .= ' | RI';
						$h_same_j_description->{$this_vector_id} .= '-'.$junction->getTypeDescription($patient) if ( $junction->getTypeDescription($patient) and $junction->getTypeDescription($patient) ne '---' );
					}
					elsif ($junction->isSE($patient)) {
						$h_same_j_description->{$this_vector_id} .= ' | SE';
						$h_same_j_description->{$this_vector_id} .= '-'.$junction->getTypeDescription($patient) if ( $junction->getTypeDescription($patient) and $junction->getTypeDescription($patient) ne '---' );
					}
					$pass_same_position = 1;
				}
			}
		}
		next if $pass_same_position;
		
		$n++;
		print '.' if ( not $only_html_cache and $n % 1000 );
		my $is_junction_linked_filtred;
		next if ( $junction->junction_score_without_dejavu_global($patient) < 0 and not $only_gene);

		next if $junction->start == $junction->end();
		
		my @lGenesNames;
		foreach my $gene ( @{ $junction->getGenes() } ) {
			if ($only_gene) {
				push( @lGenesNames, $gene->id() ) if $only_gene->id() eq $gene->id();
			}
			else  { push( @lGenesNames, $gene->id() ); }
		}
		next unless @lGenesNames;
		
		if ( not $only_gene ) {
			next if ( $junction->get_percent_new_count($patient) < $min_score );

			my $nb_dejavu_pat = 0;
			$use_percent_dejavu = $junction->dejavu_percent_coordinate_similar() unless ($use_percent_dejavu);
		
			if ($only_dejavu_ratio_10) {
				$nb_dejavu_pat =  $junction->dejavu_patients(10,$patient);
			}
			else {
				 $nb_dejavu_pat =  $junction->dejavu_patients("all",$patient);
				#$project->dejavuJunctionsResume->get_nb_junctions($junction->getChromosome->id(), $junction->start(),$junction->end(), $use_percent_dejavu, $patient_name);
			}
			next if ( $nb_dejavu_pat > $max_dejavu_value );
		}

		my $html_validation  = '';
		my $html_to_validate = '';
		my $score = $junction->junction_score($patient, $use_percent_dejavu);

		next if ($min_partial_score and $score < $min_partial_score);
		
		
		if ($export_xls) {
			my @lTypes;
			push( @lTypes, 'RI' ) if $junction->isRI($patient);
			push( @lTypes, 'SE' ) if $junction->isSE($patient);
			my $type_junction_description = join(', ', @lTypes);
			$type_junction_description .= ' - '.$junction->getTypeDescription($patient);
			my $pat_name = $patient->name();
			$hres->{chr_id} = $chr_id;
			$hres->{$junction->id()}->{global}->{'chr'} = $junction->getChromosome->id();
			$hres->{$junction->id()}->{global}->{start} = $junction->start();
			$hres->{$junction->id()}->{global}->{end} = $junction->end();
			$hres->{$junction->id()}->{global}->{type} = $type_junction_description;
			$hres->{$junction->id()}->{global}->{dejavu} = $junction->dejavu_patients('all',$patient);
			$hres->{$junction->id()}->{global}->{dejavu_ratio_10} = $junction->dejavu_patients(10,$patient);
			$hres->{$junction->id()}->{global}->{dejavu_ratio_20} = $junction->dejavu_patients(20,$patient);
			
			$hres->{$junction->id()}->{patients}->{$pat_name}->{name} = $pat_name;
			$hres->{$junction->id()}->{patients}->{$pat_name}->{fam_name} = $patient->getFamily->name();
			$hres->{$junction->id()}->{patients}->{$pat_name}->{sex} = $patient->sex();
			$hres->{$junction->id()}->{patients}->{$pat_name}->{status} = $patient->status();
			$hres->{$junction->id()}->{patients}->{$pat_name}->{nb_new} = $junction->get_nb_new_count($patient);
			$hres->{$junction->id()}->{patients}->{$pat_name}->{nb_canonique} = $junction->get_canonic_count($patient);
			$hres->{$junction->id()}->{patients}->{$pat_name}->{dp} = $junction->get_dp_count($patient);
			$hres->{$junction->id()}->{patients}->{$pat_name}->{ratio} = sprintf( "%.3f", $junction->get_percent_new_count($patient) ) . '%';

			foreach my $gene_name (@lGenesNames) {
				my $g = $project->newGene($gene_name);
				my $this_score = $score;
				my $gscore = int(($g->score/2)+0.5);
				$this_score += $gscore;
				my ($gene_id, $tmp) = split('_', $gene_name);
				$hres->{$junction->id()}->{annotation}->{$gene_name}->{gene_name} = $g->external_name();
				$hres->{$junction->id()}->{annotation}->{$gene_name}->{ensg} = $gene_id;
				$hres->{$junction->id()}->{annotation}->{$gene_name}->{description} = $g->description();
				$hres->{$junction->id()}->{annotation}->{$gene_name}->{phenotypes} = $g->phenotypes();
				$hres->{$junction->id()}->{annotation}->{$gene_name}->{score} = $this_score;
				my $h_exons_introns = $junction->get_hash_exons_introns();
				foreach my $tid ( sort keys %{$h_exons_introns} ) {
					my $t = $patient->getProject->newTranscript($tid);
					$hres->{$junction->id()}->{annotation}->{$gene_name}->{transcripts}->{$tid}->{gene} = $t->gene_external_name;
					$hres->{$junction->id()}->{annotation}->{$gene_name}->{transcripts}->{$tid}->{transcript_name} = $t->external_name;
					$hres->{$junction->id()}->{annotation}->{$gene_name}->{transcripts}->{$tid}->{ccds_name} = $t->ccds_name;
					$hres->{$junction->id()}->{annotation}->{$gene_name}->{transcripts}->{$tid}->{appris_type} = $t->appris_type;
					my @lPos = ( sort keys %{ $h_exons_introns->{$tid}->{by_pos} } );
					my $first_exon_intron = $h_exons_introns->{$tid}->{by_pos}->{ $lPos[0] };
					my $last_exon_intron = $h_exons_introns->{$tid}->{by_pos}->{ $lPos[-1] };
					$hres->{$junction->id()}->{annotation}->{$gene_name}->{transcripts}->{$tid}->{start} = $first_exon_intron;
					$hres->{$junction->id()}->{annotation}->{$gene_name}->{transcripts}->{$tid}->{end} = $last_exon_intron;
				}
			}
			$hres->{done} = 1;
			next;
		}
		
		
		
		my $html_sashimi  = get_sashimi_plot( $junction, $patient );
		my $html_igv      = get_igv_button( $junction, $patient );
		my $html_id       = get_html_id($junction, $h_same_j_description);
		my $html_patients = get_html_patients( $junction, $patient);
		my $html_dv       = get_html_dejavu( $junction, $patient );

		my $jid_tmp = $junction->getChromosome->id().'-'.$junction->start().'-'.$junction->end();
		if ($h_td_line and $h_td_line->{id} eq $jid_tmp) {
			my $max_score = $h_td_line->{max_score};
			foreach my $gene_name (@lGenesNames) {
				my $tmp = pop(@{ $hres->{genes}->{$gene_name}->{$max_score} });
			}
			push (@{$h_td_line->{3}}, "<br>");
			push (@{$h_td_line->{4}}, "<br>");
			push (@{$h_td_line->{7}}, "<br>");
			push (@{$h_td_line->{8}}, "<br>");
			
		}
		else { $h_td_line = undef; }

		$h_td_line->{id} = $jid_tmp;
		
		if (not exists $h_td_line->{1}) { push (@{$h_td_line->{1}}, $html_sashimi); }
		if (not exists $h_td_line->{2}) { push (@{$h_td_line->{2}}, $html_igv); }
#		if (exists $h_td_line->{3}) { push (@{$h_td_line->{3}}, qq{<div>------------</div>}); }
		push (@{$h_td_line->{3}}, $html_id);
#		if (exists $h_td_line->{4}) { push (@{$h_td_line->{4}}, qq{<div>------------</div>}); }
		push (@{$h_td_line->{4}}, $html_patients);
		if (not exists $h_td_line->{5}) { push (@{$h_td_line->{5}}, $html_dv); }

		my $html_push_this_j;
		foreach my $gene_name (@lGenesNames) {
			my $g = $project->newGene($gene_name);
			next unless $g;
			my $this_score = $score;
			
			my $gscore = int(($g->score/2)+0.5);
			$this_score += $gscore;
			
			my $ht = $junction->get_hash_exons_introns();
			
			delete $junction->{get_hash_exons_introns} unless $ht;
			$ht = $junction->get_hash_exons_introns();
			#			next unless $ht;
			#			next if scalar keys %$ht == 0;
			my ( $html_trans, $has_linked_junctions, $tr_found ) = get_html_transcripts( $g, $junction, $patient );
			
			my $score_details_text = get_html_score_details( $junction, $patient, $use_percent_dejavu );
			@{$h_td_line->{6}} = ();
#			@{$h_td_line->{7}} = ();
#			@{$h_td_line->{8}} = ();
			
			push (@{$h_td_line->{6}}, $html_trans);
#			if (not exists $h_td_line->{6}) { push (@{$h_td_line->{6}}, $html_trans);} 
#			if (exists $h_td_line->{7}) { push (@{$h_td_line->{7}}, qq{<div>---</div>}); }
			my $badge_color = '#808080';
			$badge_color = '#92D674' if $this_score >= 0;
			$badge_color = '#FFFF00' if $this_score >= 5;
			$badge_color = '#CC8506' if $this_score >= 8;
			my $badge_color_text = "white";
			$badge_color_text = "black" if $badge_color eq '#FFFF00';
			push (@{$h_td_line->{score}}, $this_score);
			if (not $html_push_this_j) {
				push (@{$h_td_line->{7}}, qq{<span class="badge badge-success badge-xs" style="border-color:$badge_color;background-color:$badge_color;color:$badge_color_text;margin-bottom: 15px;font-size:9px;">$this_score [G: $gscore] </span>});
				push (@{$h_td_line->{8}}, $score_details_text);
			}

			next if not $tr_found;
			
			$hres->{score}->{all}->{$gene_name}->{$this_score} = $junction->id();
			my $html_tr;
			if ($is_junction_linked_filtred) {
				$html_tr .= qq{<tr style="text-align:center;font-size:11px;opacity:0.55;max-width:130px;overflow-y:auto;">};
			}
			else {
				$html_tr .= qq{<tr style="text-align:center;font-size:11px;max-width:130px;overflow-y:auto;">};
			}
			$html_tr .= qq{<td style="width:230px;">@{$h_td_line->{1}}</td>};
			$html_tr .= qq{<td>@{$h_td_line->{2}}</td>};
			$html_tr .= qq{<td>@{$h_td_line->{3}}</td>};
			$html_tr .= qq{<td>@{$h_td_line->{4}}</td>};

			#	$html_tr .= qq{<td></td>};
			$html_tr .= qq{<td>@{$h_td_line->{5}}</td>};

			#	$html_tr .= qq{<td>$html_validation</td>};
			$html_tr .= qq{<td>@{$h_td_line->{6}}</td>};

			#	$html_tr .= qq{<td>$html_to_validate</td>};
			$html_tr .= qq{<td>@{$h_td_line->{7}}</td>};
			$html_tr .= qq{<td>@{$h_td_line->{8}}</td>};
			$html_tr .= qq{</tr>};
			my $max_score = -999;
			foreach my $s (@{$h_td_line->{score}}) { $max_score = $s if $s > $max_score; }
			$h_td_line->{max_score} = $max_score;
			push( @{ $hres->{genes}->{$gene_name}->{$max_score} }, $html_tr );
			
#			$h_last_j = undef;
#			$h_last_j->{$short_id}->{nb_canonique} = $junction->get_canonic_count($patient);
#			$h_last_j->{$short_id}->{nb_new} = $junction->get_nb_new_count($patient);
#			$h_last_j->{$short_id}->{max_score} = $max_score;
#			$h_last_j->{$short_id}->{genes}->{$gene_name} = undef;
			$html_push_this_j++;
		}
	}
	$hres->{done} = 1;
	$pm->finish( 0, $hres );
}
$pm->wait_all_children();

die if $nbErrors > 0;


if ($export_xls) {
	my $session_id = save_export_xls($patient, $h_export_xls);
	my $hash;
	$hash->{session_id} = $session_id;
	printJson($hash);
	exit(0);
}

foreach my $gene_name ( keys %{ $h_junctions_scores->{all} } ) {
	my @lscores = sort { $a <=> $b } keys %{ $h_junctions_scores->{all}->{$gene_name} };
	$h_junctions_scores->{max}->{ $lscores[-1]}->{$gene_name} = undef;
}
my @lGenesNames;
foreach my $score ( sort { $b <=> $a } keys %{ $h_junctions_scores->{max} } ) {
	foreach my $gene_name ( sort keys %{ $h_junctions_scores->{max}->{$score} } ) {
		push( @lGenesNames, $gene_name );
	}
}

if ( scalar(@lGenesNames) == 0 ) {
	$html .= qq{<tr><td><center><b><i>No Result Found...</b></i></center></td></tr>};

	if ($only_gene) {
		my $b_igv = get_igv_button( $only_gene, $patient );
		$html .= qq{<br><tr><td><center><b><i>...but you can view gene <span style="color:orange;">$only_gene_name</span> in IGV :)</b></i></center></td></tr>};
		$html .= qq{<br><tr><td><center><b><i><span class="glyphicon glyphicon-arrow-down"></span></b></i></center></td></tr><br><tr><td><center><b><i>$b_igv</b></i></center></td></tr>};
	}
	else {
		$html .= qq{<br><tr><td><center><b><i><span class="glyphicon glyphicon-exclamation-sign" style="color:red;font-size:18px;"></span> your gene <span style="color:orange;">$only_gene_name</span> is unknown in my database...</b></i></center></td></tr><br><br><br>};
	}

	$html .= qq{</table>};
	my $hash;
	$hash->{html}     = $html;
	$hash->{table_id} = '';
	if ( $release =~ /HG19/ ) {
		$hash->{used_dejavu} = $max_dejavu_value;
		$hash->{used_dejavu} = $max_dejavu if $max_dejavu;
	}
	else {
		$hash->{used_dejavu} = 'not';
	}
	printJson($hash);
	exit(0);
}
elsif ( scalar(@lGenesNames) > 0 and $only_gene_name and not $only_gene ) {
	$html .= qq{<br>} if (not $view_polyviewer);
	$html .= qq{<tr><td><b><i><span class="glyphicon glyphicon-exclamation-sign" style="color:red;font-size:18px;"></span> your gene <span style="color:orange;">$only_gene_name</span> is unknown in my database...</b></i></td></tr><br><br><br>};
}

my @lTablesIds;

my $filters =
qq{data-filter-control="input" data-filter-control-placeholder="Gene Name / Description / Panel"};
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
push( @lTablesIds, 'table_major' );

my ( @l_html, @l_tables_ids );
my $nb_g = 0;
foreach my $gene_name (@lGenesNames) {
	my $class;
	my @lScores = sort { $a <=> $b } keys %{ $h_junctions_scores->{all}->{$gene_name} };
	my $max_score = $lScores[-1];

	my $this_panel_gene_id = 'panel_' . $table_id . '_' . $gene_name;
	my $this_table_id      = $table_id . '_' . $gene_name;
	
	#	push(@lTablesIds, $this_table_id);
	my $g = $project->newGene($gene_name);
	next unless $g;
	my $hgene;
	$hgene->{description}   = $g->description();
	$hgene->{external_name} = $g->external_name();
	$hgene->{max_score}     = $max_score;
	$hgene->{js_id}         = $g->id();
	$hgene->{id}            = $g->id();

	$hgene->{pLI}              = $g->pLI();
	$hgene->{omim_id}          = $g->omim_id();
	$hgene->{omim_inheritance} = $g->omim_inheritance();
	foreach my $panel ( @{ $g->getPanels() } ) {
		$hgene->{panels}->{ $panel->name() } = undef;
	}

	my @lVarHtml;
	foreach my $score (sort {$b <=> $a} keys %{$_tr_lines_by_genes->{$gene_name}}) {
		foreach my $l ( @{ $_tr_lines_by_genes->{$gene_name}->{$score} } ) {
			push(@lVarHtml, $l);
		}
	}
	$hgene->{variants} = \@lVarHtml;

	my $div_id = 'div_' . $this_table_id;
	$hgene->{collapse_with_id} = $div_id;

	my $class_tr_gene->{style} ="background-color:#607D8B;border:1px black solid;padding:0px;white-space: nowrap;height:50px;";

	$html_table_2 .= qq{<tr><td>};
	$html_table_2 .= "<table style='width:100%;background-color:#F3F3F3;'>";
	my $html_gene_panel = update_variant_editor::panel_gene( $hgene, $this_panel_gene_id, $project->name(), $patient );
	$html_gene_panel =~ s/<b>pLI<\/b> Score//;
	$html_gene_panel =~ s/bottom:5px;/bottom:-5px;/;
	$html_table_2 .= $cgi->td( $class_tr_gene, $html_gene_panel );
	$html_table_2 .= "</table>";

	$html_table_2 .= qq{<div class="collapse" id="$div_id">};
	$html_table_2 .=
qq{<table id='$this_table_id' data-sort-name='locus' data-sort-order='desc' data-filter-control='true' data-toggle="table" data-pagination-v-align="both" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="10" data-page-list="[10]" data-resizable='true' class='table table-striped' style='font-size:11px;'>};
	$html_table_2 .= qq{<thead style="text-align:center;">};
	$html_table_2 .=
	  qq{<th data-field="plot"><b><center>Sashimi Plot</center></b></th>};
	$html_table_2 .= qq{<th data-field="igv"><b><center>IGV</center></b></th>};
	$html_table_2 .=
qq{<th data-sortable='true' data-field="var_name"><b><center>Var name</center></b></th>};
	$html_table_2 .=
	  qq{<th data-field="trio"><b><center>Trio</center></b></th>};
	$html_table_2 .=
qq{<th data-field="deja_vu"><b><center>DejaVu (Nb Samples)</center></b></th>};
	$html_table_2 .=
	  qq{<th data-field="transcripts"><b><center>Transcripts</center></b></th>};
	$html_table_2 .=
qq{<th data-sortable='true' data-field="jscore"><b><center>Score</center></b></th>};
	$html_table_2 .=
	  qq{<th data-field="noise"><b><center>Details Score</center></b></th>};
	$html_table_2 .= qq{</thead>};
	$html_table_2 .= qq{<tbody>};


	my $nb_limit = 0;
	foreach my $l (@lVarHtml) {
		next if (not $only_gene_name and $nb_limit == 100);
		$html_table_2 .= $l;
		$nb_limit++;
		if (not $only_gene_name and $nb_limit == 100) {
			$html_table_2 .= "<tr><td colspan='8'><center><span style='color:red'><i> ---------------- Limit junctions view... Search only this gene to see all junctions ---------------- </i></span></td></tr>";
			last;
		}
	}

	$html_table_2 .= qq{</tbody>};
	$html_table_2 .= qq{</table>};
	$html_table_2 .= qq{</div>};
	$html_table_2 .= qq{</td></tr>};
	$nb_g++;
	last if ( $nb_g == 2500 );
}
if ($html_table_2) {
	my $this_html = $html;
	$this_html .= $html_table_1;
	$this_html .= $html_table_2;
	$this_html .= qq{</tbody>};
	$this_html .= qq{</table>};
	push( @l_html, $this_html );
	push( @l_tables_ids, join( ',', @lTablesIds ) );
}

#$no_cnv->close();
$project->dejavuJunctions->close() if ( $release =~ /HG19/ );

my $hash;
$hash->{html}      = $l_html[0];
$hash->{tables_id} = $l_tables_ids[0];
if ( $release =~ /HG19/ ) {
	$hash->{used_dejavu} = $max_dejavu_value;
	$hash->{used_dejavu} = $max_dejavu if $max_dejavu;
}
else {
	$hash->{used_dejavu} = 'not';
}
$hash->{is_partial_results} = $is_partial_results;

if (not $only_gene_name and not $only_positions and not $view_polyviewer) {
	my $no_cache = $patient->get_lmdb_cache("w");
	$no_cache->put_cache_hash( $cache_html_id, $hash );
	$no_cache->close();
}


if ($view_polyviewer) {
	print qq{</div>};
	print $hash->{html};
	exit(0);
}
printJson($hash);
exit(0);

sub transform_results_file {
	my ( $type_analyse_name, $file_name, $path_analysis ) = @_;
	my ( $h_header, $list_res ) = parse_results_file($file_name);

	#my ($list_sort) = sort_results($list_res);
	my ( $table_id, $html ) =
	  add_table_results( $type_analyse_name, $h_header, $list_res,
		$path_analysis );
	return ( $table_id, $html );
}

sub get_sashimi_plot_list_files {
	my ( $path_analysis, $ensg, $patient, $locus, $score ) = @_;
	my $sashimi_plot_file =
	  get_sashimi_plot_file( $path_analysis, $ensg, $patient, $locus, $score );
	my @lFiles;
	if ( -e $sashimi_plot_file ) {
		push( @lFiles, $sashimi_plot_file );
		my $locus_text = $locus;
		$locus_text =~ s/chr//;
		$locus_text =~ s/:/-/;
		my ( $chr_id, $start, $end ) = split( '-', $locus_text );
		my $i = 1;
		while ( $i < 5 ) {
			$start -= ( 1000 * $i );
			$end += ( 1000 * $i );
			my $locus_extended = $chr_id . ':' . $start . '-' . $end;
			my $sashimi_plot_file =
			  get_sashimi_plot_file( $path_analysis, $ensg, $patient,
				$locus_extended, $score );
			push( @lFiles, $sashimi_plot_file );
			$i++;
		}
		return \@lFiles;
	}
	return;
}

sub get_sashimi_plot_file {
	my ( $path_analysis, $ensg, $patient, $locus, $score ) = @_;
	my $locus_text = $locus;
	$locus_text =~ s/chr//;
	$locus_text =~ s/:/-/;
	my $patient_name = $patient->name();
	my $path = $patient->getProject->getProjectPath . '/align/sashimi_plots/';
	my $outfile = $path.'/sashimi_'.$patient_name.'.'.$ensg.'.'.$locus_text.'.pdf';
	return $outfile if ( -e $outfile );
}

sub get_html_dejavu {
	my ( $junction, $patient ) = @_;
	my $color        = 'black';
	my $my_ratio     = $junction->get_percent_new_count($patient);
	my $project_name = $patient->getProject->name();
	my $patient_name = $patient->name();
	my $vector_id = $junction->getChromosome->id() . '-' . $junction->vector_id();
	my $cmd_all = qq{view_deja_vu_rna_junction(\"$project_name\",\"$patient_name\",\"$vector_id\")};
	my $cmd_inthisrun = qq{view_dejavu_nb_int_this_run_patients(\"$project_name\",\"$patient_name\",\"$vector_id\")};
	my $html = $cgi->start_table(
		{
			class => "table table-sm table-striped table-condensed table-bordered table-primary ",
			style => "box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"
		}
	);
	$html .= $cgi->start_Tr();
	$html .= $cgi->th("");
	$html .= $cgi->th("<center><b>Ratio All</b></center>");
	$html .= $cgi->th("<center><b>Ratio >10%</b></center>");
	$html .= $cgi->th("<center><b>Ratio >20%</b></center>");
	$html .= $cgi->end_Tr();
	
	my $dv_other_pat =  $junction->dejavu_patients(undef,$patient);
	my $dv_other_pat_ratio_10 =  $junction->dejavu_patients(10,$patient);
	
	my $dv_other_pat_ratio_20 =  $junction->dejavu_patients(20,$patient);
	$html .= $cgi->start_Tr();
	$html .= $cgi->td("<center><b>DejaVu</b></center>");
	$html .= $cgi->td( obutton( $cmd_all, $dv_other_pat ) );
	$html .= $cgi->td( obutton( $cmd_all, $dv_other_pat_ratio_10));
	$html .= $cgi->td( obutton( $cmd_all, $dv_other_pat_ratio_20 ) );
	$html .= $cgi->end_Tr();
	my $dv_run_other_pat = $junction->in_this_run_patients(0,$patient);
	
	my $dv_run_other_pat_ratio_10 = 0;
	my $dv_run_other_pat_ratio_20 = 0;
	if ( $dv_run_other_pat > 0 ) {
		$dv_run_other_pat_ratio_10 = $junction->in_this_run_patients(10,$patient );
		$dv_run_other_pat_ratio_20 = $junction->in_this_run_patients(20,$patient );
	}
	

	$html .= $cgi->start_Tr();
	$html .= $cgi->td("<center><b>InThisProject</b></center>");
	$html .= $cgi->td( obutton( $cmd_inthisrun, $dv_run_other_pat ));
	$html .= $cgi->td( obutton( $cmd_inthisrun, $dv_run_other_pat_ratio_10 ) );
	$html .= $cgi->td( obutton( $cmd_inthisrun, $dv_run_other_pat_ratio_20 ) );
	$html .= $cgi->end_Tr();
	$html .= $cgi->end_table();
	return $html;
}

sub get_html_transcripts {
	my ( $gene, $junction, $patient ) = @_;
	
	my $has_linked_junctions;
	my $color              = 'black';
	my $project_name       = $patient->getProject->name();
	my $patient_name       = $patient->name();
	my $is_junction_linked = $junction->is_junctions_linked($patient);
	my ( @junctions_ids_linked, $bcolor );
	@junctions_ids_linked =
	  keys %{ $junction->get_hash_junctions_linked_to_me->{ $patient->name() } }
	  if ($is_junction_linked);
	my @l_group_junctions_colors =
	  ( '#5D3EFF', '#FF4571', '#8FFF49', '#FF495F' );
	my $h_exons_introns = $junction->get_hash_exons_introns();
	my $html_tr .= "<div style='max-height:130px;overflow-y:auto;'";
	$html_tr = $cgi->start_table(
		{
			class =>
"table table-sm table-striped table-condensed table-bordered table-primary ",
			style =>
"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px;"
		}
	);
	$html_tr .= "<tr style='background-color:#FFA81E;'>";
	$html_tr .= $cgi->th("<center><b>gene</b></center>");
	$html_tr .= $cgi->th("<center><b>enst</b></center>");
	$html_tr .= $cgi->th("<center><b>nm</b></center>");
	$html_tr .= $cgi->th("<center><b>ccds</b></center>");
	$html_tr .= $cgi->th("<center><b>appris</b></center>");
	$html_tr .= $cgi->th("<center><b>start</b></center>");
	$html_tr .= $cgi->th("<center><b>end</b></center>");
	$html_tr .= $cgi->end_Tr();

	my (@lResMyGene, $hResOthersGenes, $rids, $rids_tr, $tr_found);
	foreach my $tid ( sort keys %{$h_exons_introns} ) {
#		next unless exists $h_g_tr->{$tid};
		my ( $h_junctions_linked, $h_junctions_exons_introns );
		if ( exists $h_junctions_color->{ $junction->id() } ) {
			$bcolor = $h_junctions_color->{ $junction->id() };
		}
		else {
			$bcolor = $l_group_junctions_colors[$nb_group_junctions_colors];
			$nb_group_junctions_colors++;
			$nb_group_junctions_colors = 0
			  if ( $nb_group_junctions_colors + 1 ==
				scalar(@l_group_junctions_colors) );
		}

		my $t = $patient->getProject->newTranscript($tid);
		my $this_html_tr;
		my $rid = "tr_".$gene->external_name.'_'.$t->gene_external_name.'_'.$t->external_name.'_'.time."_".int(rand(50000));
		my $hide = "";
		if ($gene->external_name eq $t->gene_external_name and not $t->isMain) {
			push(@$rids_tr,$rid);
			$hide = "display:none;";
		}
		elsif ($gene->external_name ne $t->gene_external_name) {
			push(@$rids,$rid);
			$hide = "display:none;";
		}
		$this_html_tr .= $cgi->start_Tr({id=>$rid,style=>$hide});
		$this_html_tr .= $cgi->td( "<center>" . $t->gene_external_name . "</center>" );
		$this_html_tr .= $cgi->td("<center>$tid</center>");
		$this_html_tr .= $cgi->td( "<center>" . $t->external_name() . "</center>" );
		$this_html_tr .= $cgi->td( "<center>" . $t->ccds_name() . "</center>" );
		$this_html_tr .= $cgi->td( "<center>" . $t->appris_type() . "</center>" );
		my @lPos = ( sort keys %{ $h_exons_introns->{$tid}->{by_pos} } );
		my $first_exon_intron =
		  $h_exons_introns->{$tid}->{by_pos}->{ $lPos[0] };
		my $last_exon_intron =
		  $h_exons_introns->{$tid}->{by_pos}->{ $lPos[-1] };
		my ( $first_style_color, $last_style_color );

		if ($is_junction_linked) {
			foreach my $other_j_id (@junctions_ids_linked) {
				$first_style_color = "style='background-color:$bcolor;'"
				  if (
					exists $junction->get_hash_junctions_linked_to_me
					->{ $patient->name() }->{$other_j_id}->{$tid}
					->{$first_exon_intron} );
				$last_style_color = "style='background-color:$bcolor;'"
				  if (
					exists $junction->get_hash_junctions_linked_to_me
					->{ $patient->name() }->{$other_j_id}->{$tid}
					->{$last_exon_intron} );
				if ( $first_style_color or $last_style_color ) {

					#NEW
					$h_junctions_linked->{ $junction->id() } =
					  $h_var_linked_ids->{ $junction->id() }->{vector_id};
					foreach my $other_id (
						keys
						%{ $h_var_linked_ids->{ $junction->id() }->{linked_to} }
					  )
					{
						$h_junctions_linked->{$other_j_id} =
						  $h_var_linked_ids->{$other_j_id}->{vector_id};
					}
					$h_junctions_exons_introns->{$tid}->{$first_exon_intron} =
					  undef;
					$h_junctions_exons_introns->{$tid}->{$last_style_color} =
					  undef;

					my $this_chr = $junction->getChromosome();
					my $i        = $junction->vector_id();
					my $min      = $junction->vector_id() - 15;
					$min = 0 if ( $min < 0 );
					while ( $i >= $min ) {
						my $junction2 = $this_chr->getVarObject($i);
						$i--;
						next
						  if (
							not $junction2->get_hash_junctions_linked_to_me() );

						my @lPos = (
							sort keys %{
								$junction2->get_hash_exons_introns->{$tid}
								  ->{by_pos}
							}
						);
						my $first_exon_intron_2 =
						  $junction2->get_hash_exons_introns->{$tid}->{by_pos}
						  ->{ $lPos[0] };
						my $last_exon_intron_2 =
						  $junction2->get_hash_exons_introns->{$tid}->{by_pos}
						  ->{ $lPos[-1] };
						next
						  if (
							not exists $h_junctions_exons_introns->{$tid}
							->{$first_exon_intron_2}
							and not exists $h_junctions_exons_introns->{$tid}
							->{$last_exon_intron_2} );

						foreach my $other_jid ( keys %$h_junctions_linked ) {
							if (
								exists
								$junction2->get_hash_junctions_linked_to_me
								->{ $patient->name() }->{$other_jid}->{$tid} )
							{
								$h_junctions_linked->{ $junction2->id() } =
									$this_chr->id() . '-'
								  . $junction2->vector_id();
								$h_junctions_exons_introns->{$tid}
								  ->{$first_exon_intron_2} = undef;
								$h_junctions_exons_introns->{$tid}
								  ->{$last_exon_intron_2} = undef;
								$min -= 5;
								$min = 0 if ( $min < 0 );
							}
						}
					}
					$i = $junction->vector_id();
					my $max = $junction->vector_id() + 15;
					$max = $this_chr->size_vector()
					  if ( $max >= $this_chr->size_vector() );
					while ( $i < $max ) {
						my $junction2 = $this_chr->getVarObject($i);
						$i++;
						next
						  if (
							not $junction2->get_hash_junctions_linked_to_me() );

						my @lPos = (
							sort keys %{
								$junction2->get_hash_exons_introns->{$tid}
								  ->{by_pos}
							}
						);
						my $first_exon_intron_2 =
						  $junction2->get_hash_exons_introns->{$tid}->{by_pos}
						  ->{ $lPos[0] };
						my $last_exon_intron_2 =
						  $junction2->get_hash_exons_introns->{$tid}->{by_pos}
						  ->{ $lPos[-1] };

						next
						  if (
							not exists $h_junctions_exons_introns->{$tid}
							->{$first_exon_intron_2}
							and not exists $h_junctions_exons_introns->{$tid}
							->{$last_exon_intron_2} );

						foreach my $other_jid ( keys %$h_junctions_linked ) {
							if (
								exists
								$junction2->get_hash_junctions_linked_to_me
								->{ $patient->name() }->{$other_jid}->{$tid} )
							{
								next
								  if (
									scalar keys %{
										$junction2
										  ->get_hash_junctions_linked_to_me
										  ->{ $patient->name() }->{$other_jid}
										  ->{$tid}
									} == 0
								  );
								$h_junctions_linked->{ $junction2->id() } =
									$this_chr->id() . '-'
								  . $junction2->vector_id();
								$h_junctions_exons_introns->{$tid}
								  ->{$first_exon_intron_2} = undef;
								$h_junctions_exons_introns->{$tid}
								  ->{$last_exon_intron_2} = undef;
								$max += 5;
								$max = $this_chr->size_vector()
								  if ( $max >= $this_chr->size_vector() );
							}
						}
					}
				}
			}
		}
		foreach my $jid ( keys %$h_junctions_linked ) {
			$h_junctions_color->{$jid} = $bcolor;
		}
		my $j_linked = join( ',', sort values %$h_junctions_linked );
		$h_junctions_linked = undef;
		my $my_junction_id = $junction->id();
		my $cmd_linked =
qq{view_linked_junctions(\"$patient_name\",\"$tid\",\"$j_linked\",\"$my_junction_id\",\"$min_score\")};
		if ( scalar(@lPos) == 1 ) {
			if ($first_style_color) {
				$this_html_tr .= "<td colspan='2' $first_style_color>"
				  . obutton( $cmd_linked, $first_exon_intron ) . "</td>";
			}
			else {
				$this_html_tr .=
				  "<td colspan='2' $first_style_color>$first_exon_intron</td>";
			}
		}
		else {
			if ($first_style_color) {
				$this_html_tr .= "<td $first_style_color>"
				  . obutton( $cmd_linked, $first_exon_intron ) . "</td>";
			}
			else {
				$this_html_tr .= "<td $first_style_color>$first_exon_intron</td>";
			}
			if ($last_style_color) {
				$this_html_tr .= "<td $last_style_color>"
				  . obutton( $cmd_linked, $last_exon_intron ) . "</td>";
			}
			else { $this_html_tr .= "<td $last_style_color>$last_exon_intron</td>"; }
		}
		$has_linked_junctions = 1 if ( $first_style_color or $last_style_color );
		$this_html_tr .= $cgi->end_Tr();
		
		if ($t->gene_external_name eq $gene->external_name and $t->isMain) {
			$html_tr .= $this_html_tr;
			$tr_found++;
		}
		else {
			$tr_found++ if $t->gene_external_name eq $gene->external_name;
			push(@{$hResOthersGenes->{$t->gene_external_name}}, $this_html_tr);
		} 
	}
	
	if (exists $hResOthersGenes->{$gene->external_name}) {
		my $js = encode_json $rids_tr;
		my $za = "hide_tr_".time."_".int(rand(50000));
		my $nb_skip = scalar(keys @{$hResOthersGenes->{$gene->external_name}});
		$html_tr .=  $cgi->start_Tr({id=>$za});
		$html_tr.= $cgi->td({style=>"box-shadow: 1px 1px 2px #555;background-color:#CECFCE;",colspan=>7,onClick=>qq{showTranscripts($js,"$za");}},qq{<span class="glyphicon glyphicon-plus"></span> }."view $nb_skip ALT transcript(s) from gene ".$gene->external_name);
		$html_tr.= $cgi->end_Tr();
		$html_tr .= join('', @{$hResOthersGenes->{$gene->external_name}});
		delete $hResOthersGenes->{$gene->external_name};
	}
	
	if ($hResOthersGenes and scalar keys %$hResOthersGenes > 0) {
		my $za = "hide_tr_".time."_".int(rand(50000));
		my $nb_skip = scalar(keys %$hResOthersGenes);
		my @l = sort keys %{$hResOthersGenes};
		my $l = join(', ', @l);
		$html_tr .=  $cgi->start_Tr({id=>$za});
		$html_tr.= $cgi->td({style=>"box-shadow: 1px 1px 2px #555;background-color:#CECFCE;",colspan=>7,},"FOUND $nb_skip other(s) gene(s) [$l]");
		$html_tr.= $cgi->end_Tr();
	}
	$html_tr .= qq{</table></div>};
	return ( $html_tr, $has_linked_junctions, $tr_found );
}

sub obutton {
	my ( $onclick, $value ) = @_;
	return
qq{<a class="btn btn-xs btn-primary" onclick=\'$onclick\' target="_blank" style="background-color: #D0D0D0;font-size: 7px;font-family:  Verdana;color:black" role="button">}
	  . $value . "</a>";
}

sub get_html_patients {
	my ( $junction, $patient ) = @_;
	my $h_by_pat;

	my $intspan_junction = Set::IntSpan::Fast::XS->new();
	$intspan_junction->add_range( $junction->start() - 25,
		$junction->start() + 25 );
	$intspan_junction->add_range( $junction->end() - 25,
		$junction->end() + 25 );


	foreach my $pat ( @{ $patient->getProject->getPatients() } ) {
		next if ( not $junction->get_dp_count($pat) );
		next if ( not $junction->get_nb_new_count($pat) and not $junction->isCanonique() );
		my $fam_name = $pat->getFamily->name();
		if ( $pat->isFather() ) {
			if ( $pat->isIll() ) {
				$h_by_pat->{$fam_name}->{ $pat->name() }->{status} =
				  "<center><img src='/icons/Polyicons/male-d.png'></center>";
			}
			else {
				$h_by_pat->{$fam_name}->{ $pat->name() }->{status} =
				  "<center><img src='/icons/Polyicons/male-s.png'></center>";
			}
		}
		if ( $pat->isMother() ) {
			if ( $pat->isIll() ) {
				$h_by_pat->{$fam_name}->{ $pat->name() }->{status} =
				  "<center><img src='/icons/Polyicons/female-d.png'></center>";
			}
			else {
				$h_by_pat->{$fam_name}->{ $pat->name() }->{status} =
				  "<center><img src='/icons/Polyicons/female-s.png'></center>";
			}
		}
		if ( $pat->isChild() ) {
			if ( $pat->sex() eq '1' ) {
				if ( $pat->isIll() ) {
					$h_by_pat->{$fam_name}->{ $pat->name() }->{status} =
"<center><img src='/icons/Polyicons/baby-boy-d.png'></center>";
				}
				else {
					$h_by_pat->{$fam_name}->{ $pat->name() }->{status} =
"<center><img src='/icons/Polyicons/baby-boy-s.png'></center>";
				}
			}
			else {
				if ( $pat->isIll() ) {
					$h_by_pat->{$fam_name}->{ $pat->name() }->{status} =
"<center><img src='/icons/Polyicons/baby-girl-d.png'></center>";
				}
				else {
					$h_by_pat->{$fam_name}->{ $pat->name() }->{status} =
"<center><img src='/icons/Polyicons/baby-girl-s.png'></center>";
				}
			}
		}
		$h_by_pat->{$fam_name}->{ $pat->name() }->{dp} =
		  $junction->get_dp_count($pat);
		$h_by_pat->{$fam_name}->{ $pat->name() }->{nb_new} =
		  $junction->get_nb_new_count($pat);
		$h_by_pat->{$fam_name}->{ $pat->name() }->{nb_normal} =
		  $junction->get_canonic_count($pat);
		$h_by_pat->{$fam_name}->{ $pat->name() }->{percent} =
		  sprintf( "%.3f", $junction->get_percent_new_count($pat) ) . '%';
	}
	my $color         = 'black';
	my $html_patients = $cgi->start_table(
		{
			class =>
"table table-sm table-striped table-condensed table-bordered table-primary ",
			style =>
"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"
		}
	);
	$html_patients .= qq{<thead style="text-align:center;">};

	$html_patients .= qq{<th data-field="famname"><b><center>Fam</center></b></th>};
	$html_patients .=
	  qq{<th data-field="patname"><b><center>Pat</center></b></th>};
	$html_patients .=
	  qq{<th data-field="status"><b><center>Status</center></b></th>};
	$html_patients .=
	  qq{<th data-field="percent"><b><center>Ratio (%)</center></b></th>};
	$html_patients .=
	  qq{<th data-field="nb_new"><b><center>Nb New</center></b></th>};
	$html_patients .=
	  qq{<th data-field="nb_normal"><b><center>Nb Normal</center></b></th>};
	$html_patients .= qq{<th data-field="dp"><b><center>DP</center></b></th>};


	$html_patients .= qq{<th data-field=""><b><center></center></b></th>};
	$html_patients .= qq{</thead>};
	$html_patients .= qq{<tbody>};
	foreach my $fam_name ( sort keys %{$h_by_pat} ) {
		foreach my $pat_name ( sort keys %{ $h_by_pat->{$fam_name} } ) {
			$html_patients .= qq{<tr>};
			$html_patients .= qq{<td>}.$fam_name.qq{</td>};
			$html_patients .= qq{<td>}.$pat_name . qq{</td>};
			$html_patients .= qq{<td>}.$h_by_pat->{$fam_name}->{$pat_name}->{status}.qq{</td>};
			$html_patients .= qq{<td>}.$h_by_pat->{$fam_name}->{$pat_name}->{percent}.qq{</td>};
			$html_patients .= qq{<td>}.$h_by_pat->{$fam_name}->{$pat_name}->{nb_new}.qq{</td>};
			if ($h_by_pat->{$fam_name}->{$pat_name}->{nb_normal} =~ /\./) {
				$html_patients .= qq{<td>}.$h_by_pat->{$fam_name}->{$pat_name}->{nb_normal}.qq{ <br><i>(mean_cov)</i></td>};
			}
			else {
				$html_patients .= qq{<td>}.$h_by_pat->{$fam_name}->{$pat_name}->{nb_normal}.qq{</td>};
			}
			$html_patients .= qq{<td>}.$h_by_pat->{$fam_name}->{$pat_name}->{dp}.qq{</td>};
			$html_patients .= qq{</tr>};
		}
	}
	$html_patients .= qq{</tbody>};
	$html_patients .= qq{</table>};
	return $html_patients;
}

sub get_igv_button {
	my ( $junction, $patient ) = @_;
	my $color    = 'lightgrey';
	my $bam_file = "https://www.polyweb.fr/" . $patient->bamUrl();
	my $list_patients_ctrl =
	  $patient->getPatients_used_control_rna_seq_junctions_analyse();
	if ($list_patients_ctrl) {
		my $nb_control;
		foreach my $other_pat (@$list_patients_ctrl) {
			$bam_file .= ',https://www.polyweb.fr/' . $other_pat->bamUrl();
			$nb_control++;
			last if $nb_control == 3;
		}
	}
	else {
		my $np = 0;
		foreach my $other_pat ( @{ $project->getPatients() } ) {
			next if ( $other_pat->name() eq $patient->name() );
			$bam_file .= ',https://www.polyweb.fr/' . $other_pat->bamUrl();
			$np++;
			last if $np == 3;
		}
	}
	my $gtf = $patient->getProject->get_gtf_genes_annotations_igv();
	$gtf =~ s/\/data-isilon//;
	$gtf = "https://www.polyweb.fr/" . $gtf;
	my $locus =
		$junction->getChromosome->id() . ':'
	  . ( $junction->start() - 100 ) . '-'
	  . ( $junction->end() + 100 );

	my $fasta = "";
	if ( not $patient->getProject->is_human_genome() ) {
		$fasta = $patient->getProject->genomeFasta();
		my $release     = $patient->getProject->getVersion();
		my $fasta_named = $fasta;
		$fasta_named =~ s/all/$release/;
		if ( -e $fasta_named ) { $fasta = $fasta_named; }
		$fasta =~ s/\/data-isilon//;
		$fasta = "https://www.polyweb.fr/" . $fasta;
	}

	my $igv_link =
qq{<button class='igvIcon2' onclick='launch_igv_tool_rna("$fasta", "$bam_file,$gtf","$locus")' style="color:black"></button>};
	return $igv_link;
}

sub get_sashimi_plot {
	my ( $junction, $patient ) = @_;
	my $sashimi_button;
	my $list_sashimi_plot_files = $junction->getListSashimiPlotsPathFiles($patient);
	$sashimi_button .= qq{<center>};
	my @lFiles;
	foreach my $sashimi_plot_file (@$list_sashimi_plot_files) {
		$sashimi_plot_file =~ s/\/\//\//g;
		$sashimi_plot_file =~ s/\/data-isilon\/sequencing\/ngs/\/NGS/;
		push( @lFiles, $sashimi_plot_file );
	}
	my $files = join( ';', @lFiles );
	my $pdf   = $lFiles[0] . '#toolbar=0&embedded=true';
	$sashimi_button .= qq{<button type="button" class="btn btn-default" style="border:2px black double;overflow:hidden;text-align:center;background-color:white;padding-right:20px;padding-left:4px;padding-top:4px;padding-bottom:4px;" onClick="view_pdf_list_files('$files')"><table><td>};
	$sashimi_button .= qq{<image alt="N.A." style="position:relative;width:200px;" loading="lazy" src="$pdf"></image>};
	$sashimi_button .= qq{</td><td style="padding-left:1px;"><span style="writing-mode:vertical-lr !important; font: 12px Verdana, sans-serif;letter-spacing: 1px;">Zoom</span></td></table> </button>};
	$sashimi_button .= qq{</center></};
	return $sashimi_button;
}

sub get_html_id {
	my ($junction, $h_same_j_description)     = @_;
	my $chr_id         = $junction->getChromosome->id();
	my $start          = $junction->start();
	my $end            = $junction->end();
	my $junction_locus = $chr_id . ':' . $start . '-' . $end;
	my $length         = $junction->length();
	my @lTypes;
	push( @lTypes, 'RI' ) if $junction->isRI($patient);
	push( @lTypes, 'SE' ) if $junction->isSE($patient);
	
	
	my ($type_junction, $type_junction_description);
	if (exists $h_same_j_description->{$junction->vector_id()}) {
		$type_junction = $h_same_j_description->{$junction->vector_id()};
		delete $h_same_j_description->{$junction->vector_id()};
	}
	else {
		$type_junction             = join( '+', @lTypes );
		$type_junction_description = $junction->getTypeDescription($patient);
	}
	my $html_id                   = "<center><table>";
	$html_id .= "<tr><td><center><b>$junction_locus</b></center></td></tr>";
	$html_id .= "<tr><td><center>$length nt</td></tr>";
	$html_id .= "<tr><td><center>";
	if ($type_junction) {
		$html_id .= "$type_junction - $type_junction_description" if ( $type_junction_description and $type_junction_description ne '---' );
	}
	else {
		$html_id .= "$type_junction_description" if ( $type_junction_description and $type_junction_description ne '---' );
	}
	$html_id .= "</center></td></tr>";
	$html_id .= "</table></center>";
	return $html_id;
}

sub get_html_score_details {
	my ( $junction, $patient, $use_percent_dejavu ) = @_;
	$junction->dejavu_percent_coordinate_similar($use_percent_dejavu) if $use_percent_dejavu;
	my $noise           = $junction->get_noise_score($patient);
	my $score_pen_ratio = $junction->junction_score_penality_ratio($patient);
	my $score_pen_dp    = $junction->junction_score_penality_dp($patient);
	my $score_pen_new = $junction->junction_score_penality_new_junction($patient);
	my $score_pen_noise = $junction->junction_score_penality_noise($patient);
	my $score_pen_dvrun = $junction->junction_score_penality_dejavu_inthisrun($patient);
	my $score_pen_dv = $junction->junction_score_penality_dejavu($patient);
	my $score_details_text = $cgi->start_table(
		{
			class =>
"table table-sm table-striped table-condensed table-bordered table-primary ",
			style =>
"box-shadow: 1px 1px 6px black;font-size: 7px;font-family:  Verdana;margin-bottom:0px"
		}
	);
 
	if ( $junction->isCanonique() ) {
		$score_details_text .= $cgi->start_Tr();
		$score_details_text .= $cgi->td("<center><b>isCanonique</b></center>");
		$score_details_text .= $cgi->td("<center>- 10</center>");
		$score_details_text .= $cgi->end_Tr();
	}
	if ( $score_pen_ratio > 0 ) {
		$score_details_text .= $cgi->start_Tr();
		$score_details_text .= $cgi->td("<center><b>Ratio</b></center>");
		$score_details_text .= $cgi->td("<center>- $score_pen_ratio</center>");
		$score_details_text .= $cgi->end_Tr();
	}
	if ( $score_pen_dp > 0 ) {
		$score_details_text .= $cgi->start_Tr();
		$score_details_text .= $cgi->td("<center><b>DP</b></center>");
		$score_details_text .= $cgi->td("<center>- $score_pen_dp</center>");
		$score_details_text .= $cgi->end_Tr();
	}
	if ( $score_pen_new > 0 ) {
		$score_details_text .= $cgi->start_Tr();
		$score_details_text .= $cgi->td("<center><b>New Junc</b></center>");
		$score_details_text .= $cgi->td("<center>- $score_pen_new</center>");
		$score_details_text .= $cgi->end_Tr();
	}
	if ( $score_pen_dvrun > 0 ) {
		$score_details_text .= $cgi->start_Tr();
		$score_details_text .= $cgi->td("<center><b>Inthisrun</b></center>");
		$score_details_text .= $cgi->td("<center>- $score_pen_dvrun</center>");
		$score_details_text .= $cgi->end_Tr();
	}
	if ( $score_pen_dv > 0 ) {
		$score_details_text .= $cgi->start_Tr();
		$score_details_text .= $cgi->td("<center><b>DejaVu</b></center>");
		$score_details_text .= $cgi->td("<center>- $score_pen_dv</center>");
		$score_details_text .= $cgi->end_Tr();
	}
	if ( $score_pen_noise > 0 ) {
		$score_details_text .= $cgi->start_Tr();
		$score_details_text .=
		  $cgi->td("<center><b>Noise<br></b>$noise found</center>");
		$score_details_text .= $cgi->td("<center>- $score_pen_noise</center>");
		$score_details_text .= $cgi->end_Tr();
	}
	$score_details_text .= "</table>";
	return $score_details_text;
}


sub countMinCatToUse {
	my ($h_res_v_enum) = @_;
	my $hCount;
	foreach my $chr_id ( keys %$h_res_v_enum ) {
		foreach my $cat ( keys %{ $h_res_v_enum->{$chr_id} } ) {
			my $v_filters = $project->getChromosome($chr_id)->getNewVector();
			$v_filters->from_Enum( $h_res_v_enum->{$chr_id}->{$cat} );
			$hCount->{$cat} += $project->getChromosome($chr_id)->countThisVariants($v_filters);
		}
	}
	my @lCat = sort keys %{$hCount};
	
	foreach my $cat (@lCat) {
		return $cat if $hCount->{$cat} <= 150000;
	}
	return 'min4';
}

sub printJson {
	my ($hashRes) = @_;
	my $json_encode = encode_json $hashRes;
	print ".\",";
	$json_encode =~ s/{//;
	print $json_encode;
	exit(0);
}

sub save_export_xls {
	my ($patient, $hres) = @_;
	$project->buffer->dbh_deconnect();
	$project->buffer->dbh_reconnect();
	my $xls_export = new xls_export();
	$xls_export->title_page('PolySplice_'.$patient->getProject->name().'_'.$patient->name.'.xls');
	$xls_export->{hash_junctions_global} = $hres;
	my $session_id = $xls_export->save();
	return $session_id;
}
