#!/usr/bin/perl
$|=1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";


use connect;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use JSON;

my $cgi    = new CGI;
my $buffer = GBuffer->new;
my $view_project = $cgi->param('project');

my $html;
my @list_tables_ids;
my @list_projects = @{$buffer->getQuery->getListProjectsControlGIAB()};
foreach my $project_name (reverse sort @list_projects) {
	my $buffer_tmp = GBuffer->new;
	my $project = $buffer_tmp->newProjectCache( -name => $project_name );
	my $json_control_giab = $project->getProjectPath().'/variations/control_giab/'.$project_name.'.giab.json';
	if (-e $json_control_giab) {
		open(FILE, $json_control_giab);
		my $line = <FILE>;
		close(FILE);
		my $hash_giab = decode_json($line);
		$html .= print_resume_project($project,$hash_giab->{stats},$hash_giab->{global});
		$html .= print_table_control_giab($project, $hash_giab);
		$html .= "<br><br>";
		$html .= qq{<center><div style="border:solid 2px grey;width:75%;"></div></center>};
		$html .= "<br><br>";
		push(@list_tables_ids, 'table_'.$project_name);
	}
}
my $hash;
$hash->{html} = $html;
$hash->{tables_ids} = join(';', @list_tables_ids);
my $json_encode = encode_json $hash;
print $cgi->header('text/json-comment-filtered');
print $json_encode;
exit(0);





sub print_resume_project {
	my ($project,$stat,$global) = @_;
	my $project_name = $project->name();
	my $description = $project->description();
	
	my @lCaptures;
	foreach my $c (@{$project->getCaptures()}) { push(@lCaptures, $c->name()); }
	my $captures_names = join(@lCaptures, ', ');
	
	my $html = qq{<table>};
	$html .= qq{<tr><td style="color:blue;"><b><u>$project_name</b></u></td></tr>};
	$html .= qq{<tr><td><i><b>Description:</b> $description</i></td></tr>};
	$html .= qq{<tr><td><i><b>Stats:</b> $stat</i></td></tr>};
	#$html .= qq{<tr><td><i><b>Capture(s):</b> $captures_names</i></td></tr>};
	$html .= qq{</table>};
	$html .= qq{<hr>};
	$html .= qq{<table class="table table-sm table-striped table-condensed table-bordered table-primary" style="width:40%"><tr style="background-color:#363945;color:white">};
	$html .= qq{<th>Type</th>}.qq{<th>False positive</th>}.qq{<th>true negative</th>}.qq{<th>total</th>};
	$html .= qq{</tr><tr>};
	$html.= "<td>Indels</td><td>".$global->{indels}->{fp}."</td><td>".$global->{indels}->{fn}."</td><td>".$global->{indels}->{nb}."</td></tr><tr>";
	$html.= "<td>Substitutions</td><td>".$global->{substitution}->{fp}."</td><td>".$global->{substitution}->{fn}."</td><td>".$global->{substitution}->{nb}."</td></tr></table>";
	$html .= qq{<hr>};
	return $html;
}

sub print_table_control_giab {
	my ($project, $hash_giab) = @_;
	
	my @list;
	my $project_name = $project->name();
	my $patient_name = $hash_giab->{'patient'};
	my $date_analyse = $hash_giab->{'date'};
	foreach my $h (@{$hash_giab->{'errors'}}) {
		my $event = $h->{'event'};
		my $var_id = $h->{'id'};
		my $text_seq = $h->{'text_sequence'};
		my $locus = $h->{'locus'};
		$locus =~ s/chr//;
		my $type = $h->{'vtype'};
		#$type = "-";
		my $control = $h->{'control'};
		my $control_depth = $h->{'control_depth'};
		my $giab = $h->{'giab'};
		my $text_sequence = $h->{'text_sequence'};
		my $gn = $project->getVersion();
		my $v1 = "/";
		my $bam_files = $project->getPatient($patient_name)->bamUrl();
		my $igv_link = qq{<button class='igvIcon2' onclick='launch_web_igv_js("$project_name","$patient_name","$bam_files","$locus","v1","$gn")' style="color:black">IGV</button>};
		my $tr = qq{<tr>};
		$tr .= qq{<td>$event</td>};
		$tr .= qq{<td>$locus</td>};
		$tr .= qq{<td>$type</td>};
		$tr .= qq{<td>$text_seq</td>};
		$tr .= qq{<td>$giab</td>};
		$tr .= qq{<td>$control_depth</td>};
		$tr .= qq{<td>$igv_link</td>};
		$tr .= qq{</tr>};
		push(@list, $tr);
	}
	my $table_id = 'table_'.$project->name();
	my $html = qq{<table data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="10" data-page-list="[10, 25, 50, 100, 200, 300]" data-resizable='true' id='$table_id' class='table table-striped' style='font-size:13px;'>};
	$html .= qq{<thead>};
	$html .= qq{<th data-field="event" data-filter-control="select"><b>Event</b></th>};
	$html .= qq{<th data-field="locus" data-sortable="true"><b>Locus</b></th>};
	$html .= qq{<th data-field="type" data-sortable="true"><b>Type</b></th>};
	$html .= qq{<th data-field="seq" data-sortable="true"><b>Sequence</b></th>};
	$html .= qq{<th data-field="giab" data-sortable="true"><b>GIAB</b></th>};
	$html .= qq{<th data-field="depth" data-sortable="true"><b>Control Depth</b></th>};
	$html .= qq{<th data-field="igv"><b>IGV</b></th>};
	$html .= qq{</thead>};
	$html .= qq{<tbody>};
	foreach my $tr (@list) { $html .= $tr; }
	$html .= qq{</tbody>};
	$html .= qq{</table>};
	return $html;
}