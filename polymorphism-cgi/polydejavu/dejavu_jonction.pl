#!/usr/bin/perl
$|=1;

use strict;
use CGI qw/:standard :html3/;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/";
use lib "$Bin/../../GenBo/lib/obj-nodb/";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
require "$Bin/../../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use Parallel::ForkManager;
use Tabix;
use JSON;
use export_data;
use CGI::Session;
use Storable qw(store retrieve freeze dclone thaw);
use Compress::Snappy;
use Spreadsheet::WriteExcel;
use POSIX qw(strftime);

my $cgi = new CGI;
my $junction_vector_id_search = $cgi->param('vector_id');
my $junction_id_search = $cgi->param('id');
my $project_name = $cgi->param('project');
my $patient_name = $cgi->param('patient');
my $dejavu_percent = $cgi->param('dejavu_percent');
my $is_dejavu_inthis_run = $cgi->param('is_dejavu_inthis_run');
my $dejavu_min_ratio = $cgi->param('min_ratio');


my $buffer = GBuffer->new();
my $project = $buffer->newProjectCache(-name => $project_name);
$project->getChromosomes();
my $patient = $project->getPatient($patient_name);


my $table_id_my_var = "table_dv_".$project_name."_".$patient_name."_".$junction_vector_id_search.'_my_var';
my $html_my_var = qq{<div>};
$html_my_var .= qq{<table id='$table_id_my_var' data-sort-name='identity' data-sort-order='desc' data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-page-size="10" data-resizable='true' class='table table-striped' style='font-size:10px;'>};
$html_my_var .= qq{<thead style="text-align:center;">};
$html_my_var .= qq{<th data-field="pdf"><b><center>Sashimi</center></b></th>};
$html_my_var .= qq{<th data-sortable="true" data-field="identity"><b><center>Identity</center></b></th>};
$html_my_var .= qq{<th data-field="project"><b><center>Project</center></b></th>};
$html_my_var .= qq{<th data-field="patient"><b><center>Patient</center></b></th>};
$html_my_var .= qq{<th data-field="type"><b><center>Type</center></b></th>};
$html_my_var .= qq{<th data-field="start"><b><center>Start</center></b></th>};
$html_my_var .= qq{<th data-field="end"><b><center>End</center></b></th>};
$html_my_var .= qq{<th data-field="count_junctions"><b><center>Count Junctions</center></b></th>};
$html_my_var .= qq{<th data-field="count_normal"><b><center>Count Normal</center></b></th>};
$html_my_var .= qq{<th data-field="ratio"><b><center>Ratio (%)</center></b></th>};
$html_my_var .= qq{</thead>};
$html_my_var .= qq{<tbody>};

my $table_id = "table_dv_".$project_name."_".$patient_name."_".$junction_vector_id_search;
my $html;
$html .= qq{<table id='$table_id' data-sort-name='identity' data-sort-order='desc' data-filter-control='true' data-toggle="table" data-pagination-v-align="bottom" data-toggle="table" data-cache="false" data-pagination-loop="false" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="10" data-resizable='true' class='table table-striped' data-page-list="[10]" style='font-size:10px;'>};
$html .= qq{<thead style="text-align:center;">};
$html .= qq{<th data-field="pdf"><b><center>Sashimi</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="identity"><b><center>Identity</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="project"><b><center>Project</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="patient"><b><center>Patient</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="type"><b><center>Type</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="start"><b><center>Start</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="end"><b><center>End</center></b></th>};
$html .= qq{<th data-sortable="false" data-field="count_junctions"><b><center>Count Junctions</center></b></th>};
$html .= qq{<th data-sortable="false" data-field="count_normal"><b><center>Count Normal</center></b></th>};
$html .= qq{<th data-sortable="false" data-field="ratio"><b><center>Ratio (%)</center></b></th>};
$html .= qq{</thead>};
$html .= qq{<tbody>};

my $sashimi_pdf;
my ($chr_id, $vector_id) = split('-', $junction_vector_id_search);
my $chr = $project->getChromosome($chr_id);
my $vector = $chr->getNewVector();
$vector->Bit_On($vector_id);
my @lJunction = @{$chr->getListVarObjects($vector)};
foreach my $junction (@lJunction) {
	
	my $list_sashimi_plot_files = $junction->getListSashimiPlotsPathFiles($patient);
	if ($list_sashimi_plot_files and -e $list_sashimi_plot_files->[0]) {
		$sashimi_pdf = $list_sashimi_plot_files->[0];
		$sashimi_pdf =~ s/\/\//\//g;
		$sashimi_pdf =~ s/\/data-isilon\/sequencing\/ngs/\/NGS/;
	}
			
	$junction->getPatients();
	
	my $hdv;
	if ($is_dejavu_inthis_run) {
		foreach my $dv_patient (@{$junction->getPatients()}) {
			#next if ($junction->get_percent_new_count($dv_patient) < $dejavu_min_ratio);
			$hdv->{'100%'}->{$dv_patient->getProject->name()}->{$dv_patient->name()}->{start} = $junction->start();
			$hdv->{'100%'}->{$dv_patient->getProject->name()}->{$dv_patient->name()}->{end} = $junction->end();
			$hdv->{'100%'}->{$dv_patient->getProject->name()}->{$dv_patient->name()}->{count_junctions} = $junction->get_nb_new_count($dv_patient);
			$hdv->{'100%'}->{$dv_patient->getProject->name()}->{$dv_patient->name()}->{count_normal} = $junction->get_nb_normal_count($dv_patient);
			$hdv->{'100%'}->{$dv_patient->getProject->name()}->{$dv_patient->name()}->{score} = int($junction->get_percent_new_count($dv_patient));
			$hdv->{'100%'}->{$dv_patient->getProject->name()}->{$dv_patient->name()}->{type} = $junction->getTypeDescription($dv_patient);
		}
	}
	else {
		$dejavu_percent = $junction->dejavu_percent_coordinate_similar() unless ($dejavu_percent);
		$hdv = $junction->get_dejavu_list_similar_junctions($dejavu_percent);
	}
	
	foreach my $identity (reverse sort keys %{$hdv}) {
		foreach my $prn (reverse sort keys %{$hdv->{$identity}}) {
			my $buffer_tmp = GBuffer->new();
			my $project_tmp = $buffer->newProject(-name => $project_name);
			foreach my $ptn (sort keys %{$hdv->{$identity}->{$prn}}) {
				my $type = $hdv->{$identity}->{$prn}->{$ptn}->{type};
				my $count_junctions = $hdv->{$identity}->{$prn}->{$ptn}->{count_junctions};
				my $count_normal = $hdv->{$identity}->{$prn}->{$ptn}->{count_normal};
				my $score = int($hdv->{$identity}->{$prn}->{$ptn}->{score});
				my $start = $hdv->{$identity}->{$prn}->{$ptn}->{start};
				my $end = $hdv->{$identity}->{$prn}->{$ptn}->{end};

				my $color = 'black';
				if ($identity == 100) { $color = 'red' }
				elsif ($identity == 99) { $color = 'orange'; }
				elsif ($identity >= 95) { $color = 'green'; }
				
				my $sashimi_button = qq{<center>N.A.</center>};
				my $locus_text = $chr_id.'-'.($start - 100).'-'.($end + 100);
				my $sashimi_pdf_tmp = $project_tmp->getProjectPath.'/align/sashimi_plots/';
				$sashimi_pdf_tmp =~ s/$project_name/$prn/;
				$sashimi_pdf_tmp .= 'sashimi_'.$ptn.'.'.$junction->annex->{$patient->name()}->{'ensid'}.'.'.$locus_text.'.svg';
				if (-e $sashimi_pdf_tmp) {
					$sashimi_pdf_tmp =~ s/\/\//\//g;
					$sashimi_pdf_tmp =~ s/\/data-isilon\/sequencing\/ngs/\/NGS/;
					$sashimi_pdf_tmp =~ s/\/data-isilon\/sequencing\/ngs/\/NGS/;
					my $files = $sashimi_pdf_tmp;
					my $pdf = $sashimi_pdf_tmp.'#toolbar=0&embedded=true';
					$sashimi_button = qq{<button type="button" class="btn btn-default" style="border:2px black double;overflow:hidden;text-align:center;background-color:white;padding-right:20px;padding-left:4px;padding-top:4px;padding-bottom:4px;" onClick="view_pdf_list_files('$files')"><table><td>};
					$sashimi_button .= qq{<image style="position:relative;z-index:2;width:120px;object-position: -50% -50%;transform: scale(1.8) translate(-27px, 10px);" loading="lazy" src="$pdf"></image>};
					$sashimi_button .= qq{</td><td style="padding-left:1px;"><span style="writing-mode:vertical-lr !important; font: 12px Verdana, sans-serif;letter-spacing: 1px;">Zoom</span></td></table> </button>};
					$sashimi_button .= qq{</center></};
				}
				
				my $this_html;
				if ($dejavu_min_ratio > $score) { $this_html .= qq{<tr style="opacity:0.35;">}; }
				else { $this_html .= qq{<tr>}; }
				$this_html .= qq{<td><center>$sashimi_button</center></td>};
				$this_html .= qq{<td><center><b><span style="color:$color;">$identity%</span></b></center></td>};
				$this_html .= qq{<td><center>$prn</center></td>};
				$this_html .= qq{<td><center>$ptn</center></td>};
				$this_html .= qq{<td><center>$type</center></td>};
				$this_html .= qq{<td><center>$start</center></td>};
				$this_html .= qq{<td><center>$end</center></td>};
				$this_html .= qq{<td><center>$count_junctions</center></td>};
				$this_html .= qq{<td><center>$count_normal</center></td>};
				$this_html .= qq{<td><center>$score</center></td>};
				
				if ($prn eq $project_name and $ptn eq $patient_name) { $html_my_var .= $this_html; }
				else { $html .= $this_html; }
				#$html .= $this_html;
			}
			$project_tmp = undef;
			$buffer_tmp = undef;
		}
	}
}
$html .= qq{</tbody>};
$html .= qq{</table>};

$html_my_var .= qq{</tbody>};
$html_my_var .= qq{</table>};
$html_my_var .= qq{</div>};

my $hashRes;
$hashRes->{html} = $html;
$hashRes->{html_my_variant} = $html_my_var;
$hashRes->{table_id} = $table_id; 
$hashRes->{table_id_my_variant} = $table_id_my_var; 
$hashRes->{pdf_sashimi} = 'none';
#$hashRes->{pdf_sashimi} = $sashimi_pdf if $sashimi_pdf;

my $json_encode = encode_json $hashRes;
print $cgi->header('text/json-comment-filtered');
print $json_encode;
exit(0);

