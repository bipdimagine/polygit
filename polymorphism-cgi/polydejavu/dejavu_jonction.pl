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
my $junction_id_search = $cgi->param('id');
my $project_name = $cgi->param('project');
my $patient_name = $cgi->param('patient');


my $buffer = GBuffer->new();
my $project = $buffer->newProject(-name => $project_name);
$project->getChromosomes();
my $patient = $project->getPatient($patient_name);

my $method_name;
$method_name = 'getFiltredJunctionsRI' if ($junction_id_search =~ /_RI/);
$method_name = 'getFiltredJunctionsSE' if ($junction_id_search =~ /_SE/);


my $table_id_my_var = "table_dv_".$project_name."_".$patient_name."_".$junction_id_search.'_my_var';
my $html_my_var = qq{<div style="border:1px black solid;">};
$html_my_var .= qq{<table id='$table_id_my_var' data-sort-name='score' data-sort-order='desc' data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-page-size="10" data-resizable='true' class='table table-striped' style='font-size:10px;'>};
$html_my_var .= qq{<thead style="text-align:center;">};
$html_my_var .= qq{<th data-field="project"><b><center>Project</center></b></th>};
$html_my_var .= qq{<th data-field="patient"><b><center>Patient</center></b></th>};
$html_my_var .= qq{<th data-field="type"><b><center>Type</center></b></th>};
$html_my_var .= qq{<th data-field="start"><b><center>Start</center></b></th>};
$html_my_var .= qq{<th data-field="end"><b><center>End</center></b></th>};
$html_my_var .= qq{<th data-field="count_junctions"><b><center>Count Junctions</center></b></th>};
$html_my_var .= qq{<th data-field="count_normal"><b><center>Count Normal</center></b></th>};
$html_my_var .= qq{<th data-field="score"><b><center>Score</center></b></th>};
$html_my_var .= qq{</thead>};
$html_my_var .= qq{<tbody>};

my $table_id = "table_dv_".$project_name."_".$patient_name."_".$junction_id_search;
my $html;
$html .= qq{<table id='$table_id' data-sort-name='score' data-sort-order='desc' data-filter-control='true' data-toggle="table" data-cache="false" data-pagination-loop="false" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="10" data-resizable='true' class='table table-striped' style='font-size:10px;'>};
$html .= qq{<thead style="text-align:center;">};
$html .= qq{<th data-sortable="true" data-field="project"><b><center>Project</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="patient"><b><center>Patient</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="type"><b><center>Type</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="start"><b><center>Start</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="end"><b><center>End</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="count_junctions"><b><center>Count Junctions</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="count_normal"><b><center>Count Normal</center></b></th>};
$html .= qq{<th data-sortable="true" data-field="score"><b><center>Score</center></b></th>};
$html .= qq{</thead>};
$html .= qq{<tbody>};

foreach my $junction (@{$patient->$method_name()}) {
	next unless $junction->id() eq $junction_id_search;
	$junction->getPatients();
	my $hdv = $junction->get_dejavu_list_similar_junctions(98);
	foreach my $prn (reverse sort keys %{$hdv}) {
		foreach my $ptn (sort keys %{$hdv->{$prn}}) {
			my $type = $hdv->{$prn}->{$ptn}->{type};
			my $count_junctions = $hdv->{$prn}->{$ptn}->{count_junctions};
			my $count_normal = $hdv->{$prn}->{$ptn}->{count_normal};
			my $score = $hdv->{$prn}->{$ptn}->{score};
			my $identity = $hdv->{$prn}->{$ptn}->{identity};
			my $start = $hdv->{$prn}->{$ptn}->{start};
			my $end = $hdv->{$prn}->{$ptn}->{end};
			
			my $style;
			#$style = "opacity:0.6;" if ($prn eq $project_name);
			$style = "opacity:0.6;color:red;" if ($prn eq $project_name and $ptn eq $patient_name);
			my $this_html;
			$this_html .= qq{<tr style="$style">};
			$this_html .= qq{<td><center>$prn</center></td>};
			$this_html .= qq{<td><center>$ptn</center></td>};
			$this_html .= qq{<td><center>$type</center></td>};
			$this_html .= qq{<td><center>$start</center></td>};
			$this_html .= qq{<td><center>$end</center></td>};
			$this_html .= qq{<td><center>$count_junctions</center></td>};
			$this_html .= qq{<td><center>$count_normal</center></td>};
			$this_html .= qq{<td><center>$score</center></td>};
			
#			if ($prn eq $project_name and $ptn eq $patient_name) { $html_my_var .= $this_html; }
#			else { $html .= $this_html; }
			$html .= $this_html;
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
$hashRes->{table_id_my_variant} = '';
#$hashRes->{table_id_my_variant} = $table_id_my_var; 
my $json_encode = encode_json $hashRes;
print $cgi->header('text/json-comment-filtered');
print $json_encode;
exit(0);

