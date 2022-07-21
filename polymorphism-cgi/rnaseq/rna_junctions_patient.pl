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
my $isChrome = $cgi->param('is_chrome');

my $buffer = GBuffer->new;
my $project = $buffer->newProject( -name => $project_name );
$project->getChromosomes();
my $patient = $project->getPatient($use_patient);


print $cgi->header('text/json-comment-filtered');
print "{\"progress\":\".";

my $patient_name = $patient->name();
my $table_id = 'table_'.$patient_name;
my $html = qq{<center><span style='color:red;font-size:13px'>Patient $patient_name</span></center>};
$html .= qq{<table id='$table_id' data-sort-name='score' data-sort-order='desc' data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="10" data-page-list="[10, 20, 50, 100, 200, 300]" data-resizable='true' class='table table-striped' style='font-size:11px;'>};
$html .= qq{<thead style="text-align:center;">};
$html .= qq{<th data-field="igv"><b><center>IGV</center></b></th>};
$html .= qq{<th data-field="figs"><b><center>Figure</center></b></th>};
$html .= qq{<th data-field="sashimi"><b><center>Sashimi Plot</center></b></th>};
$html .= qq{<th data-filter-control='select' data-sortable="select" data-field="junction_type"><b><center>Junction Type</center></b></th>};
$html .= qq{<th data-filter-control='input' data-sortable="true" data-field="ensid"><b><center>ENSID</center></b></th>};
$html .= qq{<th data-filter-control='input' data-sortable="true" data-field="gene"><b><center>Gene</center></b></th>};
$html .= qq{<th data-filter-control='input' data-sortable="true" data-field="chr"><b><center>Chr</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="start"><b><center>Start</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="end"><b><center>End</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="junc_count"><b><center>Junction Count</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="normal_count"><b><center>Normal Count</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="score"><b><center>Score</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="dejavu"><b><center>DejaVu</th>};
$html .= qq{<th data-sortable="true" data-field="dejavu_run"><b><center>DejaVu InThisRun</th>};
$html .= qq{</thead>};
$html .= qq{<tbody>};



my $n = 0;
foreach my $junction (@{$patient->getJunctions()}) {
	$n++;
	print '.' if ($n % 100);
	$junction->getPatients();
	
	my $ensid = $junction->annex->{ensid};
	my $gene_name = $junction->annex->{gene};
	my $chr_id = $junction->getChromosome->id();
	my $start = $junction->start();
	my $end = $junction->end();
	my $type_junction;
	$type_junction = 'RI' if $junction->isRI();
	$type_junction = 'SE' if $junction->isSE();
	my $count_new_junction = $junction->nb_new_count();
	my $count_normal_junction = $junction->nb_normal_count();
	my $score = $junction->score();
	
	#my $dv_all = $junction->dejavu_all();
	#my $dv_run = $junction->dejavu_in_this_run();
	my $dv_all = '.';
	my $dv_run = '.';
	
	
	my $color;
	$color = 'lightgrey';
	$color = 'palegreen' if ($score > 1);
	$color = 'moccasin' if ($score > 10);
	$color = 'lightcoral' if ($score > 100);
	
	my $bam_file = "https://www.polyweb.fr/".$patient->bamUrl();
	my $np = 0;
	foreach my $other_pat (@{$project->getPatients()}) {
		next if ($other_pat->name() eq $patient->name());
		$bam_file .= ',https://www.polyweb.fr/'.$other_pat->bamUrl();
		$np++;
		last if $np == 3;
	}
	my $gtf = $project->get_gtf_genes_annotations_igv();
	$gtf =~ s/\/data-isilon//;
	$gtf = "https://www.polyweb.fr/".$gtf;
	my $locus = $chr_id.':'.($junction->start()-100).'-'.($junction->end()+100);
	my $igv_link = qq{<button class='igvIcon2' onclick='launch_igv_tool("", "$bam_file,$gtf","$locus")' style="color:black"></button>};
	
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
			$sashimi_button .= qq{<button type="button" class="btn btn-default" style="text-align:center;padding:2px;background-color:$color;max-height:90px;" onClick="view_pdf_list_files('$files')"><table><td><iframe onClick="view_pdf_list_files('$files')" loading="lazy" loading="lazy" width="200" height="70" src="$pdf"></iframe></td><td><span style="writing-mode:vertical-lr !important; font: 9spx Verdana, sans-serif;padding-left:2px;letter-spacing: 1px;">Zoom</span></td></table> </button>};
		}
		else {
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
	$html .= qq{<td>$ensid</td>};
	$html .= qq{<td>$gene_name</td>};
	$html .= qq{<td>$chr_id</td>};
	$html .= qq{<td>$start</td>};
	$html .= qq{<td>$end</td>};
	$html .= qq{<td>$count_new_junction</td>};
	$html .= qq{<td>$count_normal_junction</td>};
	$html .= qq{<td>$score</td>};
	$html .= qq{<td>$dv_all</td>};
	$html .= qq{<td>$dv_run</td>};
	$html .= qq{</tr>};
}
$html .= qq{</tbody>};
$html .= qq{</table>};


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
