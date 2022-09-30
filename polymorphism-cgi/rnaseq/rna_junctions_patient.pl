#!/usr/bin/perl
$|=1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo";
use lib "$Bin/../../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";

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


$min_score = 20 unless (defined($min_score));


my $buffer = GBuffer->new;
my $project = $buffer->newProjectCache( -name => $project_name );
$project->getChromosomes();
my $patient = $project->getPatient($use_patient);
#$patient->use_not_filtred_junction_files(0) if (not $view_all_junctions);
$patient->use_not_filtred_junction_files(1);

print $cgi->header('text/json-comment-filtered');
print "{\"progress\":\".";

my $patient_name = $patient->name();
my $table_id = 'table_'.$patient_name;
my $html;
#if (not $isChrome) {
#	$html .= qq{<center><span style='color:green;font-size:15px'><u><i><b>TIPS:</b></u> use Chrome web navigator for a best experience :)</i></span></center><br>};
#}
my $release = '<b>Release:</b> '.$project->getVersion();
my $gencode;
$gencode = ' - <b>Gencode version:</b> '.$project->gencode_version() if ($release =~ /HG19/);

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
my @lJunctions;
foreach my $chr (@{$project->getChromosomes()}) {
	my $vector_patient = $patient->getJunctionsVector($chr);
	foreach my $junction (@{$chr->getListVarObjects($vector_patient)}) {
		push(@lJunctions, $junction);
	}
}

my $percent_dejavu = 0;
if (defined($use_percent_dejavu)) {
	$percent_dejavu = $use_percent_dejavu - 90;
}
my $nb_percent_dejavu_value = 90 + $percent_dejavu;

my $html_score = qq{
	<table style="width:100%;">
		<tr style="text-align:right;">
			<td style="padding-left:20px;"><b>Ratio (%):</b></td>
			<td style="padding-left:20px;width:200px;">
				<label for="slider_score" class="form-label" style="font-size:10px;font-weight:300;"><i>Ratio >= <span id="nb_score" style="color:blue;">$min_score%</span></i></label>
				<input type="range" class="form-range" min="0" max="10" value="$score_slider" step="1" id="slider_score" onchange="update_score_span_value()">
			</td>
		</tr>
	</table>};

my $max_dejavu_value = 101;
$max_dejavu_value = $max_dejavu if defined($max_dejavu);
my $html_dejavu = qq{
	<table style="width:100%;">
		<tr style="text-align:right;">
			<td style="padding-left:20px;"><b>DejaVu % Coord Similar:</b></td>
			<td style="padding-left:20px;width:200px;">
				<label for="slider_dejavu_percent" class="form-label" style="font-size:10px;font-weight:300;"><i>if <span id="nb_percent_dejavu" style="color:blue;">$nb_percent_dejavu_value</span> % similar coord</i></label>
				<input type="range" class="form-range" min="0" max="10" value="$percent_dejavu" step="1" id="slider_dejavu_percent" onchange="update_dejavu_percent_span_value()">
			</td>
			<td style="padding-left:20px;"><b>DejaVu Filters:</b></td>
			<td style="padding-left:20px;width:300px;">
				<label for="slider_dejavu" class="form-label" style="font-size:10px;font-weight:300;"><i>in max <span id="nb_max_dejavu_patients" style="color:blue;">$max_dejavu_value</span> Patient(s)</i></label>
				<input type="range" class="form-range" min="0" max="101" value="$max_dejavu_value" step="1" id="slider_dejavu" onchange="update_dejavu_span_value()">
			</td>
		</tr>
	</table>};

my $html_refresh = qq{
	<table style="width:100%;">
		<tr>
			<td style="padding-left:20px;padding-right:20px;width:50px;">
				<center>
				<button type="button" class="btn btn-sm btn-primary" id="b_update_dejavu" onclick="launch_dejavu_span_value()"><span class="glyphicon glyphicon-refresh" aria-hidden="true"></span></button>
				<br><b><span style="color:#363945">REFRESH</span></b>
				</center>
			</td>
		</tr>
	</table>};

#$html .= qq{<center><span style='color:red;font-size:13px'>Patient $patient_name</span>$release$gencode</center><br>};

$html .= qq{<div class="container" style="width:100%;padding-top:10px;padding-bottom:10px;"><div class="row">};
$html .= qq{<div class="col-sm-3" style="border:3px #363945 double;"><center><span style='color:red;font-size:13px'>Patient $patient_name</span><br>$release$gencode</center></div>};
$html .= qq{<div class="col-sm-2">$html_score</div>};
$html .= qq{<div class="col-sm-6" style="text-align:right;padding-right:20px;">$html_dejavu</div>};
$html .= qq{<div class="col-sm-1">$html_refresh</div>};
$html .= qq{</div></div>};

$html .= qq{<table id='$table_id' data-sort-name='locus' data-sort-order='desc' data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="10" data-page-list="[10]" data-resizable='true' class='table table-striped' style='font-size:11px;'>};
$html .= qq{<thead style="text-align:center;">};
$html .= qq{<th data-field="igv"><b><center>IGV</center></b></th>};
#$html .= qq{<th data-field="figs"><b><center>Figure</center></b></th>};
$html .= qq{<th data-field="sashimi"><b><center>Sashimi Plot</center></b></th>};
#$html .= qq{<th data-filter-control='input' data-sortable="true" data-field="ensid"><b><center>ENSID</center></b></th>};
$html .= qq{<th data-filter-control='input' data-sortable="true" data-field="gene"><b><center>Gene</center></b></th>};
$html .= qq{<th data-sortable="true" data-filter-control='input' data-field="locus"><b><center>Locus</center></b></th>};
if ($release =~ /HG19/) {
	if (defined($max_dejavu_value) and $max_dejavu_value == 0) {
		$html .= qq{<th data-sortable="true" data-field="dejavu_junctions"><b><center><div>DV OTHERS Patients<br>-<span style='color:red;'>$nb_percent_dejavu_value% coord</span>-<br>-<span style='color:red;'>only NEW pat.</span>-</div></b></th>};		
	}
	else {
		$html .= qq{<th data-sortable="true" data-field="dejavu_junctions"><b><center><div>DV OTHERS Patients<br>-<span style='color:red;'>$nb_percent_dejavu_value% coord</span>-<br>-<span style='color:red;'>MAX $max_dejavu_value pat.</span>-</div></b></th>};
	}
	$html .= qq{<th data-filter-control='input' $default_filter_dejavu_project data-sortable="true" data-field="dejavu_junctions_project"><b><center>DV InThisRun</b></th>};
}
#$html .= qq{<th data-filter-control='input' data-sortable="true" data-field="chr"><b><center>Chr</center></b></th>};
#$html .= qq{<th data-sortable="true" data-filter-control='input' data-field="start"><b><center>Start</center></b></th>};
#$html .= qq{<th data-sortable="true" data-filter-control='input' data-field="end"><b><center>End</center></b></th>};
$html .= qq{<th data-sortable="true" data-filter-control='input' data-filter-default='>=4' data-field="junc_count"><b><center>Junction Count</center></b></th>};
$html .= qq{<th data-sortable="true" data-filter-control='input' data-field="normal_count"><b><center>Normal Count</center></b></th>};
$html .= qq{<th data-sortable="true" data-filter-control='input' data-filter-default='>=10' data-field="dp_count"><b><center>DP Count</center></b></th>};
$html .= qq{<th data-sortable="true" data-filter-control='input' data-field="ratio"><b><center>Ratio (%)</center></b></th>};
$html .= qq{<th data-sortable="true" data-filter-control='input' data-filter-default='<=10000' data-field="length"><b><center>Length</center></b></th>};
$html .= qq{<th data-filter-control='select' data-sortable="select" data-field="junction_type"><b><center>Junction Type</center></b></th>};
$html .= qq{<th data-filter-control='select' data-sortable="select" data-field="junction_type_description"><b><center>Junction Type Description</center></b></th>};
$html .= qq{</thead>};
$html .= qq{<tbody>};

my $h_dejavu_cnv;
my $n = 0;
foreach my $junction (@lJunctions) {
	$n++;
	print '.' if ($n % 100);
#	next if (not $junction->is_filtred_results($patient));
	
	next if ($junction->isCanonique($patient));
	next if ($junction->get_ratio_new_count($patient) == 1);
	next if ($junction->get_percent_new_count($patient) < $min_score);
	
	$junction->dejavu_percent_coordinate_similar($nb_percent_dejavu_value);
	
#	next if ($junction->length() < 100);
#	next if ($junction->get_score($patient) < 1);
	
	my $nb_dejavu_pat = $junction->dejavu_nb_others_patients();
	next if (not defined($max_dejavu_value));
	next if ($nb_dejavu_pat > $max_dejavu_value);
	
	$junction->getPatients();
	my $ensid = $junction->annex->{$patient->name()}->{ensid};
	my $gene_name = $junction->annex->{$patient->name()}->{gene};
	my $h_exons_introns = $junction->get_hash_exons_introns();
	
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
	my $count_new_junction = $junction->get_nb_new_count($patient);
	my $count_normal_junction = $junction->get_nb_normal_count($patient);
	
	my $dp_count = $junction->get_dp_count($patient);
#	my $dp_count = $junction->raw_coverage($patient);
#	die;
	
	my $score = sprintf("%.3f", $junction->get_percent_new_count($patient));
	my $dv_junctions = '';
	my $jid = $junction->id();
	my $jvectorid = $chr_id.'-'.$junction->vector_id();
	my $dv_text = $nb_dejavu_pat;
	my $color_dejavu = "background-color:#92949A;";
	if ($nb_dejavu_pat <= 5) { $color_dejavu = 'background-color:red;'; } 
	elsif ($nb_dejavu_pat <= 15) { $color_dejavu = 'background-color:orange;'; }
	elsif ($nb_dejavu_pat <= 30) { $color_dejavu = 'background-color:green;'; }
	$dv_junctions = qq{$dv_text<button onclick='view_deja_vu_rna_junction(\"$project_name\",\"$patient_name\",\"$jvectorid\")' class="button button-sm" style="border-color:white;color:white;$color_dejavu;font-size:11px;">$dv_text</span>} ;

	my $dejavu_in_this_run_value_same_ratio = $junction->dejavu_nb_int_this_run_patients($min_score);
	my $dejavu_in_this_run_value = scalar (@{$junction->getPatients()});
	my $dv_run_text = $dejavu_in_this_run_value_same_ratio;
	$dv_run_text = 'F:'.$dejavu_in_this_run_value_same_ratio.' - All:'.$dejavu_in_this_run_value if ($dejavu_in_this_run_value_same_ratio < $dejavu_in_this_run_value);
	my $dejavu_in_this_run = qq{<button onclick='view_dejavu_nb_int_this_run_patients(\"$project_name\",\"$patient_name\",\"$jvectorid\")' class="button button-sm" style="border-color:white;color:grey;font-size:11px;">$dv_run_text</span>} ;
	
	my $dv_junctions_inthisrun = $dejavu_in_this_run_value_same_ratio;
	#$dv_junctions_inthisrun = $dejavu_in_this_run if ($dejavu_in_this_run > $nb_patients);
	$dv_junctions_inthisrun .= $dejavu_in_this_run;
	
	my $color;
	$color = 'lightgrey';
#	$color = 'palegreen' if ($score > 1);
#	$color = 'moccasin' if ($score > 10);
#	$color = 'lightcoral' if ($score > 100);
	
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
	my $gtf = $project->get_gtf_genes_annotations_igv();
	$gtf =~ s/\/data-isilon//;
	$gtf = "https://www.polyweb.fr/".$gtf;
	my $locus = $chr_id.':'.($junction->start()-100).'-'.($junction->end()+100);
	my $igv_link = qq{<button class='igvIcon2' onclick='launch_igv_tool_rna("", "$bam_file,$gtf","$locus")' style="color:black"></button>};
	
#	my $svg_button = 'N.A.';
#	my $svg_patient = $junction->getSvgPlotPath($patient);
#	if (-e $svg_patient) {
#		$svg_patient =~ s/\/\//\//g;
#		$svg_patient =~ s/\/data-isilon\/sequencing\/ngs/\/NGS/;
#		$svg_button = qq{<button type="button" class="btn btn-default" onClick="zoom_file('', '$svg_patient')" style="text-align:center;padding:2px;background-color:$color;max-height:60px;max-width:130px;"><img loading="lazy" src="$svg_patient" alt="Pb" width="100%"></img></button>};
#	}
	
	my $sashimi_button;
	my $list_sashimi_plot_files = $junction->getListSashimiPlotsPathFiles($patient);
	if ($list_sashimi_plot_files and -e $list_sashimi_plot_files->[0]) {
		$sashimi_button .= qq{<center>};
		my @lFiles;
		foreach my $sashimi_plot_file (@$list_sashimi_plot_files) {
			$sashimi_plot_file =~ s/\/\//\//g;
			$sashimi_plot_file =~ s/\/data-isilon\/sequencing\/ngs/\/NGS/;
			push(@lFiles, $sashimi_plot_file);
		}
		my $files = join(';', @lFiles);
		my $pdf = $lFiles[0].'#toolbar=0&embedded=true';
		$sashimi_button .= qq{<button type="button" class="btn btn-default" style="border:2px black double;overflow:hidden;text-align:center;background-color:white;padding-right:20px;padding-left:4px;padding-top:4px;padding-bottom:4px;" onClick="view_pdf_list_files('$files')"><table><td>};
		$sashimi_button .= qq{<image style="position:relative;z-index:2;width:120px;object-position: -50% -50%;transform: scale(1.8) translate(-27px, 10px);" loading="lazy" src="$pdf"></image>};
		$sashimi_button .= qq{</td><td style="padding-left:1px;"><span style="writing-mode:vertical-lr !important; font: 12px Verdana, sans-serif;letter-spacing: 1px;">Zoom</span></td></table> </button>};
		$sashimi_button .= qq{</center></};
	}
	else {
		$sashimi_button .= qq{<center>N.A.</center>};
	}
	
	my @lp;
	foreach my $p (@{$junction->getPatients()}) {
		push(@lp, $p->name);
	}
	
	$html .= qq{<tr style="text-align:center;font-size:11px;">};
	$html .= qq{<td>$igv_link</td>};
#	$html .= qq{<td>$svg_button</td>};
	$html .= qq{<td>$sashimi_button</td>};
	$html .= qq{<td>};
	$html .= qq{<b>[$ensid] $gene_name</b>};
	foreach my $tid (sort keys %{$h_exons_introns}) {
		my @lPos = (sort keys %{$h_exons_introns->{$tid}->{by_pos}});
		if (scalar(@lPos) == 1) { $html .= "<br><i><u>$tid</u></i>:".$h_exons_introns->{$tid}->{by_pos}->{$lPos[0]}; }
		else { $html .= "<br><i><u>$tid</u></i>:".$h_exons_introns->{$tid}->{by_pos}->{$lPos[0]}.'->'.$h_exons_introns->{$tid}->{by_pos}->{$lPos[-1]}; }
	}
	$html .= qq{</td>};
#	$html .= qq{<td>$ensid</td>};
#	$html .= qq{<td>$gene_name</td>};
	$html .= qq{<td>$junction_locus</td>};
	$html .= qq{<td style="max-height:150px;color:white;">$dv_junctions</td>} if ($release =~ /HG19/);
	$html .= qq{<td style="max-height:150px;color:white;">$dv_junctions_inthisrun</td>} if ($release =~ /HG19/);
#	$html .= qq{<td>$chr_id</td>};
#	$html .= qq{<td>$start</td>};
#	$html .= qq{<td>$end</td>};
	$html .= qq{<td>$count_new_junction</td>};
	$html .= qq{<td>$count_normal_junction</td>};
	$html .= qq{<td>$dp_count</td>};
	$html .= qq{<td>$score</td>};
	$html .= qq{<td>$length</td>};
	$html .= qq{<td>$type_junction</td>};
	$html .= qq{<td>$type_junction_description</td>};
	$html .= qq{</tr>};
}
$html .= qq{</tbody>};
$html .= qq{</table>};
#$no_cnv->close();
$project->dejavuJunctions->close() if ($release =~ /HG19/);


my $hash;
$hash->{html} = $html;
$hash->{table_id} = $table_id;
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
	
	if ($cgi->param('update_sashimi')) {
		unless (-d $path) {
			mkdir $path;
			`chmod 755 $path`;
		}
		my $cmd = $patient->getProject->buffer->software('ggsashimi');
		$cmd .= " -b $path_analysis/$ensg/rmdup/$patient_name\_$project_name\_rmdup.bam";
		$cmd .= " -c chr".$locus;
		$cmd .= " -o ".$outfile;
		if ($score and $score >= 100) { $cmd .= " -P /data-isilon/bipd-src/mbras/ggsashimi/ggsashimi-master/colors/red.txt"; }
		elsif ($score and $score >= 10) { $cmd .= " -P /data-isilon/bipd-src/mbras/ggsashimi/ggsashimi-master/colors/orange.txt"; }
		elsif ($score and $score >= 1) { $cmd .= " -P /data-isilon/bipd-src/mbras/ggsashimi/ggsashimi-master/colors/green.txt"; }
		else { $cmd .= " -P /data-isilon/bipd-src/mbras/ggsashimi/ggsashimi-master/colors/black.txt"; }
		$cmd .= " -C 1";
		$cmd .= " --shrink --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18";
		#$cmd .= " -P /data-isilon/bipd-src/mbras/ggsashimi/ggsashimi-master/examples/palette.txt";
		$cmd .= " -g /data-isilon/public-data/repository/HG38/gtf/ensembl/genes.gtf";
		#$cmd .= " -g /data-isilon/public-data/repository/HG38/gtf/ensembl/gencode.v40.annotation.gtf";
		#warn "\n$cmd\n";
		`$cmd`;
		#die;
		return $outfile;
	}
	return;
}


sub printJson {
	my ($hashRes) = @_;
	my $json_encode = encode_json $hashRes;
	print ".\",";
	$json_encode =~ s/{//;
	print $json_encode;
	exit(0);
}
