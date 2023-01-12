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
my $junction_vector_ids = $cgi->param('vector_ids');
my $project_name = $cgi->param('project');
my $patient_name = $cgi->param('patient');
my $transcript_name = $cgi->param('transcript');
my $min_ratio = $cgi->param('min_ratio');
my $my_junction_id = $cgi->param('this_junction');



my $buffer = GBuffer->new();
my $project = $buffer->newProjectCache(-name => $project_name);
$project->getChromosomes();
my $patient = $project->getPatient($patient_name);

my ($h_var_linked_ids, $hIds);
my @lIds = split(',', $junction_vector_ids);
foreach my $vid (@lIds) {
	my @lTmp = split('-', $vid);
	$hIds->{$lTmp[0]}->{$lTmp[1]} = undef;
}

my $h_html_pos;

my $chr;
foreach my $chr_id (keys %{$hIds}) {
	$chr = $project->getChromosome($chr_id);
	foreach my $vid (sort {$a <=> $b} keys %{$hIds->{$chr_id}}) {
		my $junction = $chr->getVarObject($vid);
		my $html_j;
		$html_j .= "<div><table style='width:100%;'><td>";
		if ($junction->id() eq $my_junction_id) {
			$html_j .= "<center><span style='color:red;'><b><u>$chr_id:".$junction->start().'-'.$junction->end()."</b></u></span>";
		}
		else {
			$html_j .= "<center><b><u>$chr_id:".$junction->start().'-'.$junction->end()."</b></u></center>";
		}
		$html_j .= "<center>RI Aval</center>" if ($junction->is_ri_aval($patient));
		$html_j .= "<center>RI Amont</center>" if ($junction->is_ri_amont($patient));
		$html_j .= "</td><td><center>".get_html_dejavu($junction, $patient)."</center></td></div></table>";
		$html_j .= "<br>";
		my $h_exons_introns = $junction->get_hash_exons_introns();
		my @lPos = (sort keys %{$h_exons_introns->{$transcript_name}->{by_pos}});
		my $first_exon_intron = $h_exons_introns->{$transcript_name}->{by_pos}->{$lPos[0]};
		$html_j .= get_html_transcripts($junction, $patient, $transcript_name);
		$html_j .= "<br>";
		$html_j .= get_html_patients($junction, $patient);
		$html_j .= "</center>";
		$h_html_pos->{html}->{$first_exon_intron}->{$junction->id()} = $html_j;
		$h_html_pos->{pos}->{$junction->start()}->{$first_exon_intron} = undef;
		$h_html_pos->{score}->{$junction->id()} = $junction->get_percent_new_count($patient);
	}
}

#warn Dumper $h_html_pos; die;


my $html_join = "<td style='padding:5px;margin:5px;'><span class='glyphicon glyphicon-arrow-right' style='font-size:24px;'/><td>";
my ($html, @ltables, @ltablesIndex);

my ($hPos, @lPos);
my $previous_exon_intron = '';
foreach my $start (sort {$a <=> $b} keys %{$h_html_pos->{pos}}) {
	foreach my $exon_intron (keys %{$h_html_pos->{pos}->{$start}}) {
		next if ($previous_exon_intron eq $exon_intron);
		$previous_exon_intron = $exon_intron;
		my $html_pos = "<td style='vertical-align:top;padding:5px;margin:5px;'>";
		$html_pos .= "<table>";
		my $nb_jid;
		foreach my $jid (sort keys %{$h_html_pos->{html}->{$exon_intron}}) {
			$nb_jid++;
			if ($min_ratio and $h_html_pos->{score}->{$jid} < $min_ratio) {
				$html_pos .= "<tr style='border:2px solid black;opacity:0.45;'><td style='padding:15px;'>";
			}
			else {
				$html_pos .= "<tr style='border:2px solid black;'><td style='padding:15px;'>";
			}
			$html_pos .= $h_html_pos->{html}->{$exon_intron}->{$jid};
			$html_pos .= "</td></tr>";
		}
		$html_pos .= "</table>";
		$html_pos .= "</td>";
		push(@ltablesIndex, $exon_intron);
		push(@ltables, $html_pos);
	}
	push(@lPos, $start);
	foreach my $jid (keys %{$h_html_pos->{score}}) {
		my @lTmp = split('_', $jid);	
		my $j_start = $lTmp[1];
		my $j_end = $lTmp[2];
		next if ($j_start ne $start);
		$hPos->{$start}->{start} = $j_start;
		$hPos->{$start}->{end} = $j_end;
	}
}

my ($html_positions, $html_datas);
my @l_datas_cov;
my $i = 0;
foreach my $html_table (@ltables) {
	$html_positions .= "<td><center>(into) <b>".$ltablesIndex[$i]."</b></center><br></td>";
	$i++;
	$html_datas .= $html_table;
	
	next if ($i == scalar(@ltables));
	if ($lPos[0] < $lPos[-1]) {
		$html_datas .= "<td style='vertical-align:top;'><table style='min-width:250px;'><td><div style='width:5px;'/></td><td><span class='glyphicon glyphicon-arrow-right' style='font-size:24px;'></td>";
	}
	else {
		$html_datas .= "<td style='vertical-align:top;'><table style='min-width:250px;'><td><div style='width:5px;'/></td><td><span class='glyphicon glyphicon-arrow-left' style='font-size:24px;'></td>";
	}
	my $start = $hPos->{$lPos[$i-1]}->{end} + 1;
	my $end = $lPos[$i] - 1;
	my ($mean_dp_patient, $length, @l_cov, @l_cov_norm, $z);
	if ($start < $end) {
		$length = $end - $start + 1;
		$mean_dp_patient = sprintf("%.2f",$patient->meanDepth($chr->id(), $start, $end));
		@l_cov = @{$patient->depth($chr->id(), $start, $end)};
		@l_cov_norm = @{$patient->normalize_depth($chr->id(), $start, $end)};
		$z = $start;
		$html_positions .= "<td><center>From $start to $end</center><br></td>";
	}
	else {
		$length = $start - $end + 1;
		$mean_dp_patient = sprintf("%.2f",$patient->meanDepth($chr->id(), $end, $start));
		@l_cov = @{$patient->depth($chr->id(), $end, $start)};
		@l_cov_norm = @{$patient->normalize_depth($chr->id(), $end, $start)};
		$z = $end;
		$html_positions .= "<td><center>From $end to $start</center><br></td>";
	}
	
	my $h_cov_graph;
	my (@X, @Y, @Y2);
	$h_cov_graph->{name} = "From $start to $end - Mean DP: $mean_dp_patient";
	foreach my $cov (@l_cov) {
		$z++;
		push(@X, $z);
		push(@Y, $cov);
	}
	foreach my $cov (@l_cov_norm) {
		push(@Y2, $cov);
	}
	$h_cov_graph->{'patient_name'} = $patient->name();
	$h_cov_graph->{'x_common'} = \@X;
	$h_cov_graph->{'y_common'} = \@Y;
	$h_cov_graph->{'y_norm'} = \@Y2;
	push(@l_datas_cov, $h_cov_graph);
	
	$html_datas .= qq{<td style='vertical-align:top;padding-top:25px;'><center>Length: $length nt<br>Mean DP: $mean_dp_patient<br><button style="margin-top:10px;" id='b_cov_plot$i' onClick="document.getElementById('plot$i').classList.remove('hidden');document.getElementById('b_cov_plot$i').classList.add('hidden')">View Coverage</button><div id='plot$i' class='hidden' style='vertical-align:top !important;width:500px;'></div></center></td>};
	
	if ($lPos[0] < $lPos[-1]) {
		$html_datas .= "<td><span class='glyphicon glyphicon-arrow-right' style='font-size:24px;'></td><td><div style='width:5px;'/></td></table></td>";
	}
	else {
		$html_datas .= "<td><span class='glyphicon glyphicon-arrow-left' style='font-size:24px;'></td><td><div style='width:5px;'/></td></table></td>";
	}
	
}

$html .= "<center><table>";
$html .= "<tr><td><div style='width:25px;'/></td>".$html_positions."<td><div style='width:25px;'/></td></tr>";
$html .= "<tr style='padding-top:15px;'><td><div style='width:25px;'/></td>".$html_datas."<td><div style='width:25px;'/></td></tr>";
$html .= "<tr><td><div style='width:20px;'></td>";
$html .= "<td><div style='width:20px;'></td></tr></table></center><br><br>";


my $hashRes;
$hashRes->{html} = $html;
$hashRes->{list_graphs} = \@l_datas_cov;
my $json_encode = encode_json $hashRes;
print $cgi->header('text/json-comment-filtered');
print $json_encode;
exit(0);



sub get_html_patients {
	my ($junction, $patient) = @_;
	my $h_by_pat;
	foreach my $pat (@{$patient->getFamily->getPatients()}) {
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
		$h_by_pat->{$fam_name}->{$pat->name()}->{nb_new} = $junction->get_nb_new_count($pat);
		$h_by_pat->{$fam_name}->{$pat->name()}->{nb_normal} = $junction->get_nb_normal_count($pat);
		$h_by_pat->{$fam_name}->{$pat->name()}->{percent} = sprintf("%.3f", $junction->get_percent_new_count($pat)).'%';
	}
	my $color = 'black';
	$color = 'red' if ($junction->id() eq $my_junction_id);
	my $html_patients= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
	$html_patients .= qq{<thead style="text-align:center;">};
	#$html_patients .= qq{<th data-field="famname"><b><center>Family</center></b></th>};
	$html_patients .= qq{<th data-field="patname"><b><center>Patient</center></b></th>};
	$html_patients .= qq{<th data-field="status"><b><center>Status</center></b></th>};
	$html_patients .= qq{<th data-field="percent"><b><center>Ratio (%)</center></b></th>};
	$html_patients .= qq{<th data-field="nb_new"><b><center>Nb New</center></b></th>};
	$html_patients .= qq{<th data-field="nb_normal"><b><center>Nb Normal</center></b></th>};
	$html_patients .= qq{<th data-field="dp"><b><center>DP</center></b></th>};
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
			$html_patients .= qq{</tr>};
		}
	}
	$html_patients .= qq{</tbody>};
	$html_patients .= qq{</table>};
	return $html_patients;
}

sub get_html_transcripts {
	my ($junction, $patient) = @_;
	my $color = 'black';
	$color = 'red' if ($junction->id() eq $my_junction_id);
	my $project_name = $patient->getProject->name();
	my $patient_name = $patient->name();
	my $is_junction_linked = $junction->is_junctions_linked($patient);
	my @junctions_ids_linked;
	@junctions_ids_linked = keys %{$junction->get_hash_junctions_linked_to_me->{$patient->name()}} if ($is_junction_linked);
	my $bcolor;
	$bcolor = $h_var_linked_ids->{$junction->id()}->{color} if ($is_junction_linked);
	
	my $h_exons_introns = $junction->get_hash_exons_introns();
	
	my $html_tr = $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
	$html_tr.= "<tr style='background-color:#9A8B4F;'>";
	$html_tr.= $cgi->th("<center><b>enst</b></center>");
	$html_tr.= $cgi->th("<center><b>nm</b></center>");
	$html_tr.= $cgi->th("<center><b>ccds</b></center>");
	$html_tr.= $cgi->th("<center><b>appris</b></center>");
	$html_tr.= $cgi->th("<center><b>start</b></center>");
	$html_tr.= $cgi->th("<center><b>end</b></center>");
	$html_tr.= $cgi->end_Tr();
	
	my $h_junctions_linked;
	foreach my $tid (sort keys %{$h_exons_introns}) {
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
					$h_junctions_linked->{$junction->id()} = $h_var_linked_ids->{$junction->id()}->{vector_id};
					$h_junctions_linked->{$other_j_id} = $h_var_linked_ids->{$other_j_id}->{vector_id};
				};
			}
		}
		my $j_linked = join(',', sort values %$h_junctions_linked);
		my $cmd_linked = qq{view_linked_junctions(\"$patient_name\",\"$j_linked\")};
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
	}
	$html_tr.= $cgi->end_Tr();
	$html_tr.= qq{</table>};
	return $html_tr;
}

sub get_html_dejavu {
	my ($junction, $patient) = @_;
	my $color = 'black';
	$color = 'red' if ($junction->id() eq $my_junction_id);
	my $project_name = $patient->getProject->name();
	my $patient_name = $patient->name();
	my $vector_id = $junction->getChromosome->id().'-'.$junction->vector_id();
	my $cmd_all = qq{view_deja_vu_rna_junction(\"$project_name\",\"$patient_name\",\"$vector_id\")};
	my $cmd_inthisrun = qq{view_dejavu_nb_int_this_run_patients(\"$project_name\",\"$patient_name\",\"$vector_id\")};
	my $html = $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
	$html.= $cgi->start_Tr();
	$html.= $cgi->th("");
	$html.= $cgi->th("<center><b>Samples</b></center>");
	$html.= $cgi->end_Tr();
	$html.= $cgi->start_Tr();
	$html.= $cgi->td("<b>Other<b>");
	$html.= $cgi->td(obutton($cmd_all,$junction->dejavu_nb_others_patients()));
	$html.= $cgi->end_Tr();
	$html.= $cgi->start_Tr();
	$html.= $cgi->td("<b>In this run</b>");
	$html.= $cgi->td(obutton($cmd_inthisrun,scalar(@{$junction->getPatients()})));
	$html.= $cgi->end_Tr();
	$html.=$cgi->end_table();
	return $html;
}

sub obutton {
	my ($onclick,$value) = @_;
	return qq{<a class="btn btn-xs btn-primary" onclick=\'$onclick\' target="_blank" style="background-color: #D0D0D0;font-size: 7px;font-family:  Verdana;color:black" role="button">}.$value."</a>";
}




