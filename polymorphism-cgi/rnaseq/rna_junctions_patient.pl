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

my $buffer = GBuffer->new;
my $project = $buffer->newProject( -name => $project_name );
$project->getChromosomes();
my $patient = $project->getPatient($use_patient);
$patient->use_not_filtred_junction_files(0) if (not $view_all_junctions);

print $cgi->header('text/json-comment-filtered');
print "{\"progress\":\".";


my $patient_name = $patient->name();
my $table_id = 'table_'.$patient_name;
my $html;
#if (not $isChrome) {
#	$html .= qq{<center><span style='color:green;font-size:15px'><u><i><b>TIPS:</b></u> use Chrome web navigator for a best experience :)</i></span></center><br>};
#}
my $release = ' - <b>Release:</b> '.$project->getVersion();
my $gencode;
$gencode = ' - <b>Gencode version:</b> '.$project->gencode_version() if ($release =~ /HG19/);

my ($default_filter_dejavu_project, $default_filter_score);
my @lJunctions = @{$patient->getJunctions()};
if (scalar (@lJunctions) > 30) {
	my ($nb_max_patient, $nb_limit);
	my $h_desc = $project->get_hash_patients_description_rna_seq_junction_analyse();
	if ($h_desc) {
		foreach my $pat_name (keys %$h_desc) { $nb_max_patient++ if (exists $h_desc->{$pat_name}->{pat}); }
	}
	else { $nb_max_patient = scalar(@{$project->getPatients()}); }
	if ($nb_max_patient == 1) { $nb_limit = 0; }
	elsif ($nb_max_patient <= 3) { $nb_limit = 1; }
	else { $nb_limit = int($nb_max_patient / 2); }
	$default_filter_dejavu_project = "data-filter-default='<=$nb_limit'" if $nb_limit > 0;
	$default_filter_score = "data-filter-default='>=1'";
}

$html .= qq{<center><span style='color:red;font-size:13px'>Patient $patient_name</span>$release$gencode</center><br>};
$html .= qq{<table id='$table_id' data-sort-name='score' data-sort-order='desc' data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="10" data-page-list="[10]" data-resizable='true' class='table table-striped' style='font-size:11px;'>};
$html .= qq{<thead style="text-align:center;">};
$html .= qq{<th data-field="igv"><b><center>IGV</center></b></th>};
$html .= qq{<th data-field="figs"><b><center>Figure</center></b></th>};
$html .= qq{<th data-field="sashimi"><b><center>Sashimi Plot</center></b></th>};
$html .= qq{<th data-filter-control='select' data-sortable="select" data-field="junction_type"><b><center>Junction Type</center></b></th>};
$html .= qq{<th data-filter-control='select' data-sortable="select" data-field="junction_type_description"><b><center>Junction Type Description</center></b></th>};
$html .= qq{<th data-filter-control='input' data-sortable="true" data-field="ensid"><b><center>ENSID</center></b></th>};
$html .= qq{<th data-filter-control='input' data-sortable="true" data-field="gene"><b><center>Gene</center></b></th>};
$html .= qq{<th data-filter-control='input' data-sortable="true" data-field="chr"><b><center>Chr</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="start"><b><center>Start</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="end"><b><center>End</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="junc_count"><b><center>Junction Count</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="normal_count"><b><center>Normal Count</center></b></th>};
$html .= qq{<th data-filter-control='input' $default_filter_score data-sortable="true" data-field="score"><b><center>Score</center></b></th>};
#$html .= qq{<th data-sortable="true" data-field="dejavu_cnv"><b><center>DejaVu CNVs</th>};
$html .= qq{<th data-filter-control='input' $default_filter_dejavu_project data-sortable="true" data-field="dejavu_junctions_project"><b><center>DejaVu Junctions InThisRun</b><br>-Exact jonction-</th>};
$html .= qq{<th data-sortable="false" data-field="dejavu_junctions"><b><center>DejaVu Junctions</b></th>};
$html .= qq{</thead>};
$html .= qq{<tbody>};

my $h_dejavu_cnv;
my $n = 0;
foreach my $junction (@lJunctions) {
	$n++;
	print '.' if ($n % 100);
	#next if (not $junction->is_filtred_results($patient));
	
	$junction->getPatients();
	my $ensid = $junction->annex->{$patient->name()}->{ensid};
	my $gene_name = $junction->annex->{$patient->name()}->{gene};
	my $chr_id = $junction->getChromosome->id();
	my $start = $junction->start();
	my $end = $junction->end();
	my @lTypes;
	push(@lTypes, 'RI') if $junction->isRI($patient);
	push(@lTypes, 'SE') if $junction->isSE($patient);
	my $type_junction = join('+', @lTypes);
	my $type_junction_description = $junction->getTypeDescription($patient);
	my $count_new_junction = $junction->get_nb_new_count($patient);
	my $count_normal_junction = $junction->get_nb_normal_count($patient);
	my $score = sprintf("%.4f", $junction->get_score($patient));
	my $dv_cnv = '.';
	
	#DEJAVU JUNCTIONS
#	my $dv_junctions_in_this_run = '.';
#	my $junction_dv_id = $type_junction.'_'.$start.'_'.$end;
#	my $nb_pat_run = scalar(@{$junction->getPatients()}) - 1;
#	#$dv_junctions_in_this_run = $nb_pat_run.'/'.($max_patients - 1) if ($nb_pat_run);
#	$dv_junctions_in_this_run = $nb_pat_run;
	
	my $dv_junctions = '';
#	my (@lDvJunctions_resume, $h_tmp_dv);
#	foreach my $h (@{$junction->get_dejavu_list_similar_junctions_resume(98)}) {
#		my $this_id = $h->{id};
#		my $nb_proj = $h->{projects};
#		my $nb_pat = $h->{patients};
#		my $same_as = $h->{same_as};
#		my @lTmp = split('_', $this_id);
#		if ($lTmp[1] == $junction->start() and $lTmp[2] == $junction->end()) {
#			$same_as = '100%';
#		}
#		$same_as =~ s/%//;
##		if ($same_as == 100 and exists $h_junctions_dejavu_run->{$this_id}) {
##			$nb_proj --;
##			$nb_pat -= $h_junctions_dejavu_run->{$this_id};
##		}
#		next if ($nb_proj == 0);
#		my $color = 'green';
#		$color = 'red' if ($same_as eq '100');
#		$color = 'orange' if ($same_as eq '99');
#		my $text = '<b><span style="color:'.$color.';">'.$same_as.'%</span></b>:'.$nb_proj.'/'.$nb_pat;
#		$h_tmp_dv->{$same_as}->{$nb_proj}->{$nb_pat} = $text;
#		if ($same_as == 100) {
#			$dv_junctions_inthisrun = $nb_proj;
#		}
#	}
#	foreach my $same_as (sort {$b <=> $a} keys %$h_tmp_dv) {
#		foreach my $nb_proj (sort {$b <=> $a} keys %{$h_tmp_dv->{$same_as}}) {
#			foreach my $nb_pat (sort {$b <=> $a} keys %{$h_tmp_dv->{$same_as}->{$nb_proj}}) {
#				push(@lDvJunctions_resume, $h_tmp_dv->{$same_as}->{$nb_proj}->{$nb_pat});
#			}
#		}
#	}
#	if (@lDvJunctions_resume) {
		my $jid = $junction->id();
		my $dv_text = $junction->dejavu_nb_projects().'/'.$junction->dejavu_nb_patients();
#		warn $dv_text;
		$dv_junctions = "<button onclick='view_deja_vu_rna_junction(\"$project_name\",\"$patient_name\",\"$jid\")'>$dv_text</button>";
#	}
	my $dv_junctions_inthisrun = $junction->inthisrun_nb_patients();
	
	
	my $color;
	$color = 'lightgrey';
	$color = 'palegreen' if ($score > 1);
	$color = 'moccasin' if ($score > 10);
	$color = 'lightcoral' if ($score > 100);
	
	my $bam_file = "https://www.polyweb.fr/".$patient->bamUrl();
	my $list_patients_ctrl = $patient->getPatients_used_control_rna_seq_junctions_analyse();
	if ($list_patients_ctrl) {
		foreach my $other_pat (@$list_patients_ctrl) {
			$bam_file .= ',https://www.polyweb.fr/'.$other_pat->bamUrl();
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
	
	my $svg_button = 'N.A.';
	my $svg_patient = $junction->getSvgPlotPath($patient);
	if (-e $svg_patient) {
		$svg_patient =~ s/\/\//\//g;
		$svg_patient =~ s/\/data-isilon\/sequencing\/ngs/\/NGS/;
		$svg_button = qq{<button type="button" class="btn btn-default" onClick="zoom_file('', '$svg_patient')" style="text-align:center;padding:2px;background-color:$color;max-height:60px;max-width:130px;"><img loading="lazy" src="$svg_patient" alt="Pb" width="100%"></img></button>};
	}
	
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
		if ($isChrome) {
			my $pdf = $lFiles[0].'#toolbar=0&navpanes=0&scrollbar=0&embedded=true';
			$sashimi_button .= qq{<button type="button" class="btn btn-default" style="text-align:center;padding:2px;background-color:$color;max-height:90px;" onClick="view_pdf_list_files('$files')"><table><td><iframe onClick="view_pdf_list_files('$files')" loading="lazy" loading="lazy" width="200" height="70" src="$pdf" suppressConsole></iframe></td><td><span style="writing-mode:vertical-lr !important; font: 9spx Verdana, sans-serif;padding-left:2px;letter-spacing: 1px;">Zoom</span></td></table> </button>};
		}
		else {
			#Fonctionne mais fait planter la page Web Mozilla si trop de PDF
#			my $pdf = $lFiles[0].'#toolbar=0&embedded=true';
#			$sashimi_button .= qq{<button type="button" class="btn btn-default" style="text-align:center;background-color:$color;" onClick="view_pdf_list_files('$files')"><table><td>};
#			$sashimi_button .= qq{<div class="pdfobject-container">};
#			$sashimi_button .= qq{<embed loading="lazy" class="pdfobject" title="Embedded PDF" src="$pdf" allowfullscreen suppressConsole></embed>};
#			$sashimi_button .= qq{</div>};
#			$sashimi_button .= qq{</td><td><span style="writing-mode:vertical-lr !important; font: 9spx Verdana, sans-serif;padding-right:1px;letter-spacing: 1px;">Zoom</span></td></table> </button>};
			
			$sashimi_button .= qq{<button type="button" class="btn btn-default" onClick="view_pdf_list_files('$files')" style="text-align:center;padding:2px;background-color:$color;"> PDF <span style="top:3px;width:20px;height:20px;" class="glyphicon glyphicon-zoom-in" aria-hidden="true"></span></button>};
		}
		$sashimi_button .= qq{</center></};
	}
	else {
		$sashimi_button .= qq{<center>.</center>};
	}
	
	$html .= qq{<tr style="text-align:center;font-size:11px;">};
	$html .= qq{<td>$igv_link</td>};
	$html .= qq{<td>$svg_button</td>};
	$html .= qq{<td>$sashimi_button</td>};
	$html .= qq{<td>$type_junction</td>};
	$html .= qq{<td>$type_junction_description</td>};
	$html .= qq{<td>$ensid</td>};
	$html .= qq{<td>$gene_name</td>};
	$html .= qq{<td>$chr_id</td>};
	$html .= qq{<td>$start</td>};
	$html .= qq{<td>$end</td>};
	$html .= qq{<td>$count_new_junction</td>};
	$html .= qq{<td>$count_normal_junction</td>};
	$html .= qq{<td>$score</td>};
#	$html .= qq{<td>$dv_cnv</td>};
	$html .= qq{<td style="max-height:150px;">$dv_junctions_inthisrun</td>};
	$html .= qq{<td style="max-height:150px;">$dv_junctions</td>};
	$html .= qq{</tr>};
}
$html .= qq{</tbody>};
$html .= qq{</table>};
#$no_cnv->close();
$project->dejavuJunctions->close();


my $hash;
$hash->{html} = $html;
$hash->{table_id} = $table_id;

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
