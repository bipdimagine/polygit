#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;

use strict;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../cache_nodb/scripts/";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use GBuffer;
use export_data;
use JSON;
use polyweb_dude;
use VcfMerge;
use GenBoNoSql;
use Set::IntervalTree;
use Spreadsheet::WriteExcel;
use Bit::Vector;
use Bit::Vector::Overload;
use Compress::Snappy;
use Storable qw(store retrieve freeze dclone thaw);
use POSIX qw(strftime);
use List::MoreUtils qw(natatime);
use CGI::Session;
use html; 
use Carp;
use Cache_Commons;
use QueryVectorFilter;
use IO::Handle;
use xls_export;
use session_export;
use MIME::Base64;

require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/html_polygenescout.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";
require  "$Bin/../GenBo/lib/obj-nodb/GenBoDuckDejaVuCnv.pm";
require  "$Bin/../GenBo/lib/obj-nodb/GenBoDuckDejaVuSv.pm";

my $io = IO::Handle->new();
$io->autoflush(1);


my $cgi = new CGI();
my $cnv_id = $cgi->param('id');
my $project_name = $cgi->param('project');


my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => $project_name );



#$type,$chr,$start,$end,$seuil
#warn Dumper $duck->dejavu_details('DUP', 11, 7237519, 7238614, 90);


my @ltmp = split('!', $cnv_id);
my ($type, $chr, $start, $end) = split('_',$ltmp[-1]);

#DEL_5_101278034_101279143
#DUP_11_7237519_7238614

my @lTR;
my ($duck);
if ($type eq 'DUP' or $type eq 'DEL') {
	$duck = GenBoDuckDejaVuCNV->new( project => $project );
	my $h = $duck->dejavu_details($type, $chr, $start, $end, 90);
	
	foreach my $proj_id (sort keys %{$h}) {
		foreach my $pat_id (sort keys %{$h->{$proj_id}}) {
			my $proj_name = $buffer->getProjectNameFromId($proj_id);
			my $b = new GBuffer;
			my $p = $b->newProjectCache( -name => $proj_name );
			foreach my $pat (@{$p->getPatients}) {
				next if $pat->id ne $pat_id;
				my $pat_name = $pat->name();
				my $start = $h->{$proj_id}->{$pat_id}->{start};
				my $end = $h->{$proj_id}->{$pat_id}->{end};
				my $caller_sr = $h->{$proj_id}->{$pat_id}->{caller_sr};
				my $caller_cov = $h->{$proj_id}->{$pat_id}->{caller_coverage};
				my $caller_depth = $h->{$proj_id}->{$pat_id}->{caller_depth};
				my $identity = $h->{$proj_id}->{$pat_id}->{identity}.'%';
				
				my $out = "<td>".$proj_name."</td>";
				$out .= "<td>".$pat_name."</td>";
				$out .= "<td>".$start."</td>";
				$out .= "<td>".$end."</td>";
				$out .= "<td>".$caller_sr."</td>";
				$out .= "<td>".$caller_cov."</td>";
				$out .= "<td>".$caller_depth."</td>";
				$out .= "<td>".$identity."</td>";
				push(@lTR, $out);
			}
		} 
	}
}
else {
	my ($type, $chr1, $start1, $chr2, $start2) = split('_',$ltmp[-1]);
	
	my $duck = GenBoDuckDejaVuSv->new( project => $project );
	my $h = $duck->get_dejavu_details($chr1, $start1, $chr2, $start2, 100);
	
	foreach my $proj_id (sort keys %{$h}) {
		foreach my $pat_id (sort keys %{$h->{$proj_id}}) {
			my $proj_name = $buffer->getProjectNameFromId($proj_id);
			my $b = new GBuffer;
			my $p = $b->newProjectCache( -name => $proj_name );
			foreach my $pat (@{$p->getPatients}) {
				next if $pat->id ne $pat_id;
				my $pat_name = $pat->name();
				
				my $caller_sr = 0;
				my $caller_cov = 0;
				my $caller_depth = 0;
				
				my $start = $h->{$proj_id}->{$pat_id}->{'chr1'}.':'.$h->{$proj_id}->{$pat_id}->{'pos1'};
				my $end = $h->{$proj_id}->{$pat_id}->{'chr2'}.':'.$h->{$proj_id}->{$pat_id}->{'pos2'};
				$caller_sr += $h->{$proj_id}->{$pat_id}->{caller_sr};
				$caller_cov += $h->{$proj_id}->{$pat_id}->{caller_coverage};
				$caller_depth += $h->{$proj_id}->{$pat_id}->{caller_depth};
				my $identity = '-';
				
				my $out = "<td>".$proj_name."</td>";
				$out .= "<td>".$pat_name."</td>";
				$out .= "<td>".$start."</td>";
				$out .= "<td>".$end."</td>";
				$out .= "<td>".$caller_sr."</td>";
				$out .= "<td>".$caller_cov."</td>";
				$out .= "<td>".$caller_depth."</td>";
				$out .= "<td>".$identity."</td>";
				push(@lTR, $out);
			}
		} 
	}
}

my $out2 = $cgi->start_div();
$out2 .= qq{<table data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-v-align="top" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="20" data-page-list="[10, 20]" data-resizable='true' id='table_projects' class='table table-striped' style='font-size:13px;'>};
$out2 .= "<thead>";
$out2 .= $cgi->start_Tr({style=>"background-color:#E9DEFF;font-size:10px"});
$out2 .= qq{<th data-field="name" data-sortable="true" data-filter-control="select" data-filter-control-placeholder="ALL Projects">Name</th>};
$out2 .= qq{<th data-field="samples" data-sortable="true" data-filter-control="input" data-filter-control-placeholder="Patient Name">Sample</th>};
if ($type eq 'DUP' or $type eq 'DEL') {
	$out2 .= qq{<th data-field="start" data-sortable="true">Start</th>};
	$out2 .= qq{<th data-field="end" data-sortable="true">End</th>};
}
else {
	$out2 .= qq{<th data-field="start" data-sortable="true">Coord 1</th>};
	$out2 .= qq{<th data-field="end" data-sortable="true">Coord 2</th>};
}
$out2 .= qq{<th data-field="call_sr" data-sortable="true">Caller SR</th>};
$out2 .= qq{<th data-field="call_cov" data-sortable="true">Caller Coverage</th>};
$out2 .= qq{<th data-field="call_depth" data-sortable="true">Caller Depth</th>};
$out2 .= qq{<th data-field="identity" data-sortable="true">Identity</th>};
$out2 .= $cgi->end_Tr();
$out2 .= "</thead>";
$out2 .= "<tbody>";

foreach my $tr (@lTR) {
	$out2 .= "<tr>".$tr."</tr>";
}

$out2 .= "</tbody>";
$out2 .= "</table>";
$out2 .= "</div>";
$out2 .= "<br><br>";


my $hRes;
$hRes->{html} = $out2;
printJson($hRes);
exit(0);

sub printJson {
	my ($hashRes, $test) = @_;
	my $json_encode = encode_json $hashRes;
	print $cgi->header('text/json-comment-filtered');
	print $json_encode;
	exit(0);
}
