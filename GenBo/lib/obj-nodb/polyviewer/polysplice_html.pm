package polysplice_html;
use strict;
use FindBin qw($Bin);
use Moo;
use Data::Dumper;
use JSON::XS;
extends 'polyviewer_html';


my $max_dejavu_value = 51;


has width_default => (is => 'rw');
has view_polyviewer => (is => 'rw');
has value_filter_check_only_dejavu_ratio_10 => (is => 'rw');
has value_filter_percent_similar_dejavu => (is => 'rw');
has value_filter_max_dejavu => (is => 'rw');
has value_filter_min_score => (is => 'rw');
has value_filter_only_gene_name => (is => 'rw');
has value_filter_only_positions => (is => 'rw');
has value_filter_html_gene_positions => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $width_default = $self->width_default();
		return qq{<input type="text" style="font-size:11px;width:$width_default;" placeholder="COL4A5, 1:100000-1500000" class="form-control" id="input_gene_positions">};
	}
);
has background_color_in_polyviewer => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub { return '#FFFFCC'; }
);

sub set_view_polyviewer {
	my ($self, $value) = @_;
	$self->width_default('180px');
	if ($value) {
		$self->width_default('120px');
		return $self->view_polyviewer(1);
	}
}

sub set_filter_check_only_dejavu_ratio_10 {
	my ($self, $value) = @_;
	$self->value_filter_check_only_dejavu_ratio_10(qq{checked="checked"}) if ($value);
}

sub set_filter_only_gene_name {
	my ($self, $value) = @_;
	return unless $value;
	my $width_default = $self->width_default();
	$self->value_filter_html_gene_positions( qq{<input type="text" style="font-size:11px;width:$width_default;" value="$value" class="form-control" id="input_gene_positions">} );
	$self->value_filter_only_gene_name($value);
} 

sub set_filter_only_positions {
	my ($self, $value) = @_;
	return unless $value;
	my $width_default = $self->width_default();
	$self->value_filter_html_gene_positions( qq{<input type="text" style="font-size:11px;width:$width_default;" value="$value" class="form-control" id="input_gene_positions">} );
	$self->value_filter_only_positions($value);
}

sub set_filter_percent_similar_dejavu {
	my ($self, $value) = @_;
	$self->value_filter_percent_similar_dejavu(96);
	$self->value_filter_percent_similar_dejavu($value) if ($value);
}

sub set_filter_max_dejavu {
	my ($self, $value) = @_;
	if (defined($value)) {
		$max_dejavu_value = $value;
		$self->value_filter_max_dejavu($value);
	}
}

sub set_filter_min_score {
	my ($self, $value) = @_;
	$self->value_filter_min_score(20);
	if (defined($value)) {
		$self->value_filter_min_score($value);
	}
}

sub get_html_filters_dejavu {
	my ($self) = @_;
	my $width_default = $self->width_default();
	my $percent_dejavu = $self->value_filter_percent_similar_dejavu();
	my $nb_percent_dejavu_value =  $percent_dejavu - 90;
	my $checked_only_dejavu_ratio_10 = $self->value_filter_check_only_dejavu_ratio_10();
	my $html_dejavu = qq{
		<table style="width:100%;">
			<tr>
				<td style="padding-top:5px;">
					<center>
						<label for="slider_dejavu_percent" class="form-label" style="font-size:10px;font-weight:300;"><nobr><i>if <span id="nb_percent_dejavu" style="color:blue;">$percent_dejavu</span> % similar coord</i></nobr></label>
						<input style="width:$width_default;" type="range" class="form-range" min="0" max="10" value="$nb_percent_dejavu_value" step="1" id="slider_dejavu_percent" onchange="update_dejavu_percent_span_value()">
					</center>
				</td>
				<td style="padding-top:5px;">
					<center>
						<label for="slider_dejavu" class="form-label" style="font-size:10px;font-weight:300;"><i><nobr>in max <span id="nb_max_dejavu_patients" style="color:blue;">$max_dejavu_value</span> Patient(s)</nobr></i></label>
						<input style="width:$width_default;" type="range" class="form-range" min="0" max="51" value="$max_dejavu_value" step="1" id="slider_dejavu" onchange="update_dejavu_span_value()">
					</center>
				</td>
				<td style="padding-top:5px;">
					<center>
						<div class="form-check"><input $checked_only_dejavu_ratio_10 class="form-check-input" type="checkbox" value="" id="b_dejavu_min_ratio_10"><label class="form-check-label" for="b_dejavu_min_ratio_10" style="padding-left:10px;font-size:11px;"><nobr>only ratio >10%</nobr></label></div>
					</center>		
				</td>
						
			</tr>
			<tr>
				<td style="padding-top:5px;font-size:11px;"><center><b>DejaVu Min Similarity</center></b></td>
				<td style="padding-top:5px;font-size:11px;"><center><b>DejaVu Max patients</center></b></td>
				<td style="padding-top:5px;font-size:11px;"><center><b>DejaVu Ratio</center></b></td>
			</tr>
		</table>
	};
	return $html_dejavu;
}

sub get_html_filters_gene_positions_select {
	my ($self) = @_;
	my $width_default = $self->width_default();
	return $self->
	return 
	
}

sub get_html_filters_others {
	my ($self) = @_;
	my $min_score = $self->value_filter_min_score();
	my $score_slider = 0;
	$score_slider = $min_score / 10 if ($min_score and $min_score > 0);
	my $html_gene_select = $self->value_filter_html_gene_positions();
	my $html_filters = qq{
		<table style="width:100%;">
			<tr>
				<td style="padding-top:5px;">
					<center>
						<label for="slider_score" class="form-label" style="font-size:10px;font-weight:300;"><i><nobr>Ratio >= <span id="nb_score" style="color:blue;">$min_score%</span></i></nobr></label>
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
				<td style="padding-top:5px;font-size:11px;"><center><b>Filter Min Ratio</center></b></td>
				<td style="padding-top:5px;font-size:11px;"><center><b>Filter Only Gene / Positions</center></b></td>
			</tr>
		</table>
	};
	return $html_filters;
}

sub get_html_filters_patient_infos {
	my ($self) = @_;
	return if $self->view_polyviewer();
	my $patient_name = $self->patient->name();
	my $release = $self->project->getVersion();
	my $gencode = '-';
	$gencode = $self->project->gencode_version() if ($self->project->getVersion());
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
		</table>
	};
	return $html_patient;
}

sub get_html_filters_refresh {
	my ($self) = shift;
	my $patient_name = $self->patient->name();
	my $html_refresh = qq{
		<table style="width:100%;">
			<tr>
				<td style="padding-top:12px;">
					<center>
						<button style="font-size:25px;" type="button" class="btn btn-sm btn-primary" id="b_update_dejavu" onclick="launch_dejavu_span_value('$patient_name')"><span class="glyphicon glyphicon-refresh" aria-hidden="true"></span></button>
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
		</table>
	};
	return $html_refresh;
}

sub get_html_header_panel_patient_filters {
	my ($self) = @_;
	my $html_dejavu = $self->get_html_filters_dejavu();
	my $html_filters = $self->get_html_filters_others();
	my $html_refresh = $self->get_html_filters_refresh();
	my $html = qq{<div class="container" style="width:100%;height:86px;vertical-align:middle;"><div class="row">};
	if ($self->view_polyviewer()) {
		$html .= qq{<div class="col-sm-7"><div style="height:82px;border:1px solid black;text-align:center;font-size:9px !important;"><center>$html_dejavu</center></div></div>};
		$html .= qq{<div class="col-sm-4"><div style="height:82px;border:1px solid black;text-align:center;font-size:9px !important;"><center>$html_filters</center></div></div>};
		$html .= qq{<div class="col-sm-1"><div style="height:82px;border:1px solid black;text-align:center;font-size:9px !important;"><center>$html_refresh</center></div></div>};
	}
	else {
		my $html_patient = $self->get_html_filters_patient_infos();
		$html .= qq{<div class="col-sm-2"><div style="height:82px;border:3px #363945 double;">$html_patient</div></div>};
		$html .= qq{<div class="col-sm-5"><div style="height:82px;border:1px solid black;text-align:center;"><center>$html_dejavu</center></div></div>};
		$html .= qq{<div class="col-sm-4"><div style="height:82px;border:1px solid black;text-align:center;"><center>$html_filters</center></div></div>};
		$html .= qq{<div class="col-sm-1"><div style="height:82px;border:1px solid black;text-align:center;"><center>$html_refresh</center></div></div>};
	}
	$html .= qq{</div></div>};
	$html .= qq{<br>};
	return $html;
}

sub print_html_header_in_polysplice {
	my ($self) = @_;
	my $html;
	$html .= qq{<thead style="text-align:center;">};
	$html .= qq{<th data-field="plot"><b><center>Sashimi Plot</center></b></th>};
	$html .= qq{<th data-field="igv"><b><center>IGV</center></b></th>};
	$html .= qq{<th data-sortable='true' data-field="var_name"><b><center>Var name</center></b></th>};
	$html .= qq{<th data-field="trio"><b><center>Trio</center></b></th>};
	$html .= qq{<th data-field="deja_vu"><b><center>DejaVu (Nb Samples)</center></b></th>};
	$html .= qq{<th data-field="transcripts"><b><center>Transcripts</center></b></th>};
	$html .= qq{<th data-sortable='true' data-field="jscore"><b><center>Score</center></b></th>};
	$html .= qq{<th data-field="noise"><b><center>Details Score</center></b></th>};
	$html .= qq{</thead>};
	return $html
}

sub print_html_like_header_in_polyviewer {
	my ($self) = @_;
	my $bcolor = $self->background_color_in_polyviewer();
	my $html;
	$html .= qq{<tr>};
	$html .= qq{<td style="background-color:$bcolor;"><b><center>Score</center></b></td>};
	$html .= qq{<td style="background-color:$bcolor;"><b><center>IGV</center></b></td>};
	$html .= qq{<td style="background-color:$bcolor;"><b><center>Details Score</center></b></td>};
	$html .= qq{<td style="background-color:$bcolor;"><b><center>Junction name</center></b></td>};
	$html .= qq{<td style="background-color:$bcolor;"><b><center>Trio</center></b></td>};
	$html .= qq{<td style="background-color:$bcolor;" colspan='3'><b><center>DejaVu (Nb Samples)</center></b></td>};
	$html .= qq{<td style="background-color:$bcolor;"><b><center>Sashimi Plot / Transcripts</center></b></td>};
	$html .= qq{</tr>};
	return $html;
}

sub print_html_line_in_polyviewer_table {
	my ($self, $junction) = @_;
	my $html_sashimi = $self->get_html_sashimi_plot($junction);
	my $html_igv = $self->get_html_igv($junction);
	my $html_id = $self->get_html_id($junction);
	my $html_patients = $self->get_html_patients($junction);
	my $html_dejavu = $self->get_html_dejavu($junction);
	my $html_consequences = $self->get_html_consequences($junction);
	my $html_score_details = $self->get_html_score_details_condensed($junction);
	my $score = $junction->score();
	my $html_score = $self->printBadge($junction->score(),[5,8]);
	$html_score =~ s/#4CAF50/#000000/g if $score < 3;
	my $bcolor = $self->background_color_in_polyviewer();
	
	my $html_tr;
	$html_tr .= qq{<tr style="text-align:center;font-size:11px;">};
	$html_tr .= qq{<td style="border: 2px solid black;border-right: 0px solid black;color:black;font-size:7px;background-color:$bcolor;"><div><b>Score</b><br><br>$html_score</div></td>};
	$html_tr .= qq{<td style="border: 2px solid black;border-right: 0px solid black;border-left: 0px solid black;color:black;font-size:6px;background-color:$bcolor;"><center>$html_score_details</center></td>};
	$html_tr .= qq{<td style="border: 2px solid black;border-right: 0px solid black;border-left: 0px solid black;color:black;background-color:$bcolor;">$html_igv</td>};
	$html_tr .= qq{<td style="border: 2px solid black;border-right: 0px solid black;border-left: 0px solid black;color:black;font-size:8px;background-color:$bcolor;">$html_id</td>};
	$html_tr .= qq{<td style="border: 2px solid black;border-right: 0px solid black;border-left: 0px solid black;color:black;background-color:$bcolor;">$html_patients</td>};
	$html_tr .= qq{<td style="border: 2px solid black;border-right: 0px solid black;border-left: 0px solid black;color:black;background-color:$bcolor;" colspan='3'>$html_dejavu</td>};
	if ($html_sashimi =~ /N.A./) {
		$html_tr .= qq{<td style="border: 2px solid black;border-left: 0px solid black;color:black;background-color:$bcolor;">$html_consequences</td>};
	}
	else {
		my $html = qq{<div class="container" style="width:100%;vertical-align:middle;"><div class="row">};
		$html .= qq{<div class="col-sm-4"><center>$html_sashimi</center></div>};
		$html .= qq{<div class="col-sm-8"><center>$html_consequences</center></div>};
		$html .= qq{</div>};
		$html_tr .= qq{<td style="border: 2px solid black;border-left: 0px solid black;color:black;background-color:$bcolor;">$html</td>};
	}
	$html_tr .= qq{</tr>};
	return $html_tr;
}

sub print_html_line_in_polysplice {
	my ($self, $junction) = @_;
	my $html_sashimi = $self->get_html_sashimi_plot($junction);
	my $html_igv = $self->get_html_igv($junction);
	my $html_id = $self->get_html_id($junction);
	my $html_patients = $self->get_html_patients($junction);
	my $html_dejavu = $self->get_html_dejavu($junction);
	my $html_consequences = $self->get_html_consequences($junction);
	my $html_score = $self->get_html_score($junction);
	my $html_score_details = $self->get_html_score_details($junction);
	
	my $html_tr;
	if ($junction->is_junctions_linked()) { $html_tr .= qq{<tr style="text-align:center;font-size:11px;opacity:0.55;">}; }
	else { $html_tr .= qq{<tr style="text-align:center;font-size:11px;">}; }
	$html_tr .= qq{<td style="width:230px;">$html_sashimi</td>};
	$html_tr .= qq{<td>$html_igv</td>};
	$html_tr .= qq{<td>$html_id</td>};
	$html_tr .= qq{<td>$html_patients</td>};
	$html_tr .= qq{<td>$html_dejavu</td>};
	$html_tr .= qq{<td>$html_consequences</td>};
	$html_tr .= qq{<td>$html_score</td>};
	$html_tr .= qq{<td>$html_score_details</td>};
	$html_tr .= qq{</tr>};
	return $html_tr;
}

sub get_html_score_details {
	my ($self, $junction) = @_;
	my $cgi = $self->cgi;
	my $noise = $junction->nb_noise_junctions();
	my $score_pen_ratio = $junction->score_penality_ratio();
	my $score_pen_dp = $junction->score_penality_dp();
	my $score_pen_new = $junction->score_penality_new_junction();
	my $score_pen_dvrun = $junction->score_penality_dejavu_inthisrun();
	my $score_pen_dv = $junction->score_penality_dejavu();
	my $score_pen_noise = $junction->score_penality_noise();
	my $score_details_text = $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px black;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
	if ($score_pen_ratio > 0) {
		$score_details_text.= $cgi->start_Tr();
		$score_details_text.= $cgi->td("<center><b>Ratio</b></center>");
		$score_details_text.= $cgi->td("<center>- $score_pen_ratio</center>");
		$score_details_text.= $cgi->end_Tr();
	}
	if ($score_pen_dp > 0) {
		$score_details_text.= $cgi->start_Tr();
		$score_details_text.= $cgi->td("<center><b>DP</b></center>");
		$score_details_text.= $cgi->td("<center>- $score_pen_dp</center>");
		$score_details_text.= $cgi->end_Tr();
	}
	if ($score_pen_new > 0) {
		$score_details_text.= $cgi->start_Tr();
		$score_details_text.= $cgi->td("<center><b>New Junc</b></center>");
		$score_details_text.= $cgi->td("<center>- $score_pen_new</center>");
		$score_details_text.= $cgi->end_Tr();
	}
	if ($score_pen_dvrun > 0) {
		$score_details_text.= $cgi->start_Tr();
		$score_details_text.= $cgi->td("<center><b>Inthisrun</b></center>");
		$score_details_text.= $cgi->td("<center>- $score_pen_dvrun</center>");
		$score_details_text.= $cgi->end_Tr();
	}
	if ($score_pen_dv > 0) {
		$score_details_text.= $cgi->start_Tr();
		$score_details_text.= $cgi->td("<center><b>DejaVu</b></center>");
		$score_details_text.= $cgi->td("<center>- $score_pen_dv</center>");
		$score_details_text.= $cgi->end_Tr();
	}
	if ($score_pen_noise > 0) {
		$score_details_text.= $cgi->start_Tr();
		$score_details_text.= $cgi->td("<center><b>Noise<br></b>$noise found</center>");
		$score_details_text.= $cgi->td("<center>- $score_pen_noise</center>");
		$score_details_text.= $cgi->end_Tr();
	}
	$score_details_text .= "</table>";
	return $score_details_text;
}

sub get_html_score_details_condensed {
	my ($self, $junction) = @_;
	my $cgi = $self->cgi;
	my $noise = $junction->nb_noise_junctions();
	my $score_pen_ratio = $junction->score_penality_ratio();
	my $score_pen_dp = $junction->score_penality_dp();
	my $score_pen_new = $junction->score_penality_new_junction();
	my $score_pen_dvrun = $junction->score_penality_dejavu_inthisrun();
	my $score_pen_dv = $junction->score_penality_dejavu();
	my $score_pen_noise = $junction->score_penality_noise();
	my $score_details_text = "<table style='font-size: 7px;font-family: Verdana;margin-bottom:0px'>";
	if ($score_pen_ratio > 0) {
		$score_details_text.= $cgi->start_Tr();
		$score_details_text.= $cgi->td("<center><nobr><b>Ratio</b> -$score_pen_ratio</nobr></center>");
		$score_details_text.= $cgi->end_Tr();
	}
	if ($score_pen_dp > 0) {
		$score_details_text.= $cgi->start_Tr();
		$score_details_text.= $cgi->td("<center><nobr><b>DP</b> -$score_pen_dp</nobr></center>");
		$score_details_text.= $cgi->end_Tr();
	}
	if ($score_pen_new > 0) {
		$score_details_text.= $cgi->start_Tr();
		$score_details_text.= $cgi->td("<center><nobr><b>New Junc</b> -$score_pen_new</nobr></center>");
		$score_details_text.= $cgi->end_Tr();
	}
	if ($score_pen_dvrun > 0) {
		$score_details_text.= $cgi->start_Tr();
		$score_details_text.= $cgi->td("<center><nobr><b>DVrun</b> -$score_pen_dvrun</nobr></center>");
		$score_details_text.= $cgi->end_Tr();
	}
	if ($score_pen_dv > 0) {
		$score_details_text.= $cgi->start_Tr();
		$score_details_text.= $cgi->td("<center><nobr><b>DV</b> -$score_pen_dv</nobr></center>");
		$score_details_text.= $cgi->end_Tr();
	}
	if ($score_pen_noise > 0) {
		$score_details_text.= $cgi->start_Tr();
		$score_details_text.= $cgi->td("<center><nobr><b>Noise</b> ($noise found)</nobr> -$score_pen_noise</center>");
		$score_details_text.= $cgi->end_Tr();
	}
	$score_details_text .= "</table>";
	return $score_details_text;
}

sub get_html_score {
	my ($self, $junction) = @_;
	my $cgi = $self->cgi;
	return "<center>".$junction->score()."</center>";
}

sub get_html_consequences {
	my ($self, $junction) = @_;
	my $color = 'black';
	my $cgi = $self->cgi;
	my $html_tr = $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
	$html_tr.= "<tr style='background-color:#FFA81E;'>";
	$html_tr.= $cgi->th("<center><b>enst</b></center>");
	$html_tr.= $cgi->th("<center><b>nm</b></center>");
	$html_tr.= $cgi->th("<center><b>ccds</b></center>");
	$html_tr.= $cgi->th("<center><b>appris</b></center>");
	$html_tr.= $cgi->th("<center><b>start</b></center>");
	$html_tr.= $cgi->th("<center><b>end</b></center>");
	$html_tr .= $cgi->end_Tr();
	
	
	
	if ($junction->hash_exons_introns()) {
		foreach my $tid (sort keys %{$junction->hash_exons_introns()}) {
			$html_tr .= $cgi->start_Tr();
			$html_tr .= $cgi->td("<center>$tid</center>");
			if ($junction->hash_transcripts_descriptions->{$tid}->{external_name}) {
				$html_tr .= $cgi->td("<center>".$junction->hash_transcripts_descriptions->{$tid}->{external_name}."</center>");
			}
			else { $html_tr .= $cgi->td("<center>-</center>"); }
			if ($junction->hash_transcripts_descriptions->{$tid}->{ccds_name}) {
				$html_tr .= $cgi->td("<center>".$junction->hash_transcripts_descriptions->{$tid}->{ccds_name}."</center>");
			}
			else { $html_tr .= $cgi->td("<center>-</center>"); }
			if ($junction->hash_transcripts_descriptions->{$tid}->{appris_type}) {
				$html_tr .= $cgi->td("<center>".$junction->hash_transcripts_descriptions->{$tid}->{appris_type}."</center>");
			}
			else { $html_tr .= $cgi->td("<center>-</center>"); }
			my @lPos = (sort keys %{$junction->hash_exons_introns->{$tid}->{by_pos}});
			if (@lPos) {
				my $first_exon_intron = $junction->hash_exons_introns->{$tid}->{by_pos}->{$lPos[0]};
				my $last_exon_intron = $junction->hash_exons_introns->{$tid}->{by_pos}->{$lPos[-1]};
				my ($is_first_exon_intron_linked, $is_last_exon_intron_linked);
				foreach my $j_other_id (keys %{$junction->hash_junctions_linked_to_me()}) {
					$is_first_exon_intron_linked = 1 if exists $junction->hash_junctions_linked_to_me->{$j_other_id}->{$tid}->{$first_exon_intron};
					$is_last_exon_intron_linked = 1 if exists $junction->hash_junctions_linked_to_me->{$j_other_id}->{$tid}->{$last_exon_intron};
				}
				
				my $first_style_color = "style='background-color:red;'";
				my $last_style_color = "style='background-color:red;'";
				my $cmd_linked = $junction->{cmd_view_linked_junctions}->{$tid};
				if (scalar(@lPos) == 1) {
					if ($is_first_exon_intron_linked) { $html_tr .= "<td colspan='2' $first_style_color>".$self->obutton($cmd_linked,$first_exon_intron)."</td>"; }
					else { $html_tr .= "<td colspan='2'>$first_exon_intron</td>"; }
				}
				else {
					if ($is_first_exon_intron_linked) { $html_tr .= "<td $first_style_color><center>".$self->obutton($cmd_linked, $first_exon_intron)."</center></td>"; }
					else { $html_tr .= "<td><center>$first_exon_intron</center></td>"; }
					if ($is_last_exon_intron_linked) { $html_tr .= "<td $last_style_color><center>".$self->obutton($cmd_linked, $last_exon_intron)."</center></td>"; }
					else { $html_tr .= "<td><center>$last_exon_intron</center></td>"; }
				}
			}
			else {
				$html_tr .= "<td><center>-</center></td>";
				$html_tr .= "<td><center>-</center></td>";
			}
			$html_tr .= $cgi->start_Tr();
		}
	}
		
		
#		foreach my $tid (sort keys %{$junction->hash_exons_introns()}) {
#			my @lPos = (sort keys %{$junction->hash_exons_introns->{$tid}->{by_pos}});
#			my $first_exon_intron = $junction->hash_exons_introns->{$tid}->{by_pos}->{$lPos[0]};
#			my $last_exon_intron = $junction->hash_exons_introns->{$tid}->{by_pos}->{$lPos[-1]};
#			my $first_style_color = '';
#			my $last_style_color = '';
#			
#			
#			if ($junction->{hash_transcripts_styles_colors}) {
#				warn "\n\n";
#				warn Dumper $junction->{hash_transcripts_styles_colors};
#				warn "\n\n";
#			}
#			
#			$first_style_color = $junction->{hash_transcripts_styles_colors}->{$tid}->{first} if exists $junction->{hash_transcripts_styles_colors}->{$tid}->{first};
#			$last_style_color = $junction->{hash_transcripts_styles_colors}->{$tid}->{last} if exists $junction->{hash_transcripts_styles_colors}->{$tid}->{last};
#			$html_tr .= $cgi->start_Tr();
#			$html_tr .= $cgi->td("<center>$tid</center>");
#			if ($junction->hash_transcripts_descriptions->{$tid}->{external_name}) {
#				$html_tr .= $cgi->td("<center>".$junction->hash_transcripts_descriptions->{$tid}->{external_name}."</center>");
#			}
#			else { $html_tr .= $cgi->td("<center>-</center>"); }
#			if ($junction->hash_transcripts_descriptions->{$tid}->{ccds_name}) {
#				$html_tr .= $cgi->td("<center>".$junction->hash_transcripts_descriptions->{$tid}->{ccds_name}."</center>");
#			}
#			else { $html_tr .= $cgi->td("<center>-</center>"); }
#			if ($junction->hash_transcripts_descriptions->{$tid}->{appris_type}) {
#				$html_tr .= $cgi->td("<center>".$junction->hash_transcripts_descriptions->{$tid}->{appris_type}."</center>");
#			}
#			else { $html_tr .= $cgi->td("<center>-</center>"); }
#			my $cmd_linked = $junction->{cmd_view_linked_junctions}->{$tid};
#			if (scalar(@lPos) == 1) {
#				if ($first_style_color) { $html_tr .= "<td colspan='2' $first_style_color>".$self->obutton($cmd_linked,$first_exon_intron)."</td>"; }
#				else { $html_tr .= "<td colspan='2' $first_style_color>$first_exon_intron</td>"; }
#			}
#			else {
#				if ($first_style_color) { $html_tr .= "<td $first_style_color><center>".$self->obutton($cmd_linked, $first_exon_intron)."</center></td>"; }
#				else { $html_tr .= "<td $first_style_color><center>$first_exon_intron</center></td>"; }
#				if ($last_style_color) { $html_tr .= "<td $last_style_color><center>".$self->obutton($cmd_linked, $last_exon_intron)."</center></td>"; }
#				else { $html_tr .= "<td $last_style_color><center>$last_exon_intron</center></td>"; }
#			}
#		}
#	}
	$html_tr.= $cgi->end_Tr();
	$html_tr.= qq{</table>};
	return $html_tr;
}

sub get_html_dejavu  {
	my ($self, $junction) = @_;
	my $color = 'black';
	my $cgi = $self->cgi;
	my $cmd_all = $junction->cmd_dejavu();
	my $cmd_inthisrun = $junction->cmd_dejavu_inthisrun();
	my $dv_other_pat = $junction->dejavu_nb_other_patients();
	my $dv_other_pat_ratio_10 = $junction->dejavu_nb_other_patients_min_ratio_10();
	my $dv_other_pat_ratio_20 = $junction->dejavu_nb_other_patients_min_ratio_20();
	my $dv_run_other_pat = $junction->dejavu_nb_int_this_run_patients();
	my $dv_run_other_pat_ratio_10 = $junction->dejavu_nb_int_this_run_patients_min_ratio_10();
	my $dv_run_other_pat_ratio_20 = $junction->dejavu_nb_int_this_run_patients_min_ratio_20();
	my $html = $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
	$html .= $cgi->start_Tr();
	$html .= $cgi->th("");
	$html .= $cgi->th("<center><b>Ratio All</b></center>");
	$html .= $cgi->th("<center><b>Ratio >10%</b></center>");
	$html .= $cgi->th("<center><b>Ratio >20%</b></center>");
	$html .= $cgi->end_Tr();
	$html .= $cgi->start_Tr();
	$html .= $cgi->td("<center><b>DejaVu</b></center>");
	$html .= $cgi->td($self->obutton($cmd_all, $dv_other_pat));
	$html .= $cgi->td($self->obutton($cmd_all, $dv_other_pat_ratio_10));
	$html .= $cgi->td($self->obutton($cmd_all, $dv_other_pat_ratio_20));
	$html .= $cgi->end_Tr();
	$html .= $cgi->start_Tr();
	$html .= $cgi->td("<center><b>InThisRun</b></center>");
	$html .= $cgi->td($self->obutton($cmd_inthisrun, $dv_run_other_pat));
	$html .= $cgi->td($self->obutton($cmd_inthisrun, $dv_run_other_pat_ratio_10));
	$html .= $cgi->td($self->obutton($cmd_inthisrun, $dv_run_other_pat_ratio_20));
	$html .= $cgi->end_Tr();
	$html .= $cgi->end_table();
	return $html;
}

sub get_html_patients {
	my ($self, $junction) = @_;
	my $cgi = $self->cgi;
	my $color = 'black';
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
	foreach my $pat_name (sort keys %{$junction->{hash_by_pat_html_status}}) {
		$html_patients .= qq{<tr>};
		#$html_patients .= qq{<td>}.$fam_name.qq{</td>};
		my $ensid = $junction->ensid();
		my $proj_name = $junction->project_name();
		my $cmd_p = qq{view_var_from_proj_gene_pat(\"$proj_name\", \"$ensid\", \"$pat_name\")};
		$html_patients .= "<td>".$self->obutton($cmd_p, $pat_name)."</td>";
		$html_patients .= qq{<td>}.$junction->{hash_by_pat_html_status}->{$pat_name}.qq{</td>};
		$html_patients .= qq{<td>}.$junction->{hash_by_pat_percent}->{$pat_name}.qq{</td>};
		$html_patients .= qq{<td>}.$junction->{hash_by_pat_nb_new}->{$pat_name}.qq{</td>};
		$html_patients .= qq{<td>}.$junction->{hash_by_pat_nb_normal}->{$pat_name}.qq{</td>};
		$html_patients .= qq{<td>}.$junction->{hash_by_pat_dp_count}->{$pat_name}.qq{</td>};
		$html_patients .= qq{</tr>};
	}
	$html_patients .= qq{</tbody>};
	$html_patients .= qq{</table>};
	return $html_patients;
}

sub get_html_id {
	my ($self, $junction) = @_;
	my $chr_id = $junction->chromosome();
	my $start = $junction->start();
	my $end = $junction->end();
	my $junction_locus = $junction->locus();
	my $length = $junction->length();
	my @lTypes;
	push(@lTypes, 'RI') if $junction->isRI();
	push(@lTypes, 'SE') if $junction->isSE();
	my $type_junction = join('+', @lTypes);
	my $type_junction_description = $junction->type_description();
	my $html_id = "<center><table>";
	$html_id .= "<tr><td><center>$length nt</td></tr>";
	$html_id .= "<tr><td><center>";
	$html_id .= "$type_junction";
	$html_id .= " - $type_junction_description" if ($type_junction_description and $type_junction_description ne '---');
	$html_id .= "</center></td></tr>";
	$html_id .= "</table></center>";
	my $region = $junction->chromosome().'-'.$junction->start().'-'.$junction->end();
	return $self->printSimpleBadge(qq{<a href='https://gnomad.broadinstitute.org/region/$region' target = '_blank' style="color:black">$junction_locus</a> }).$html_id;
}


sub get_html_igv {
	my ($self, $junction) = @_;
	my $color = 'lightgrey';
	my $bam_files = $junction->bam_url().','.$junction->bam_controls_urls();
	my $gtf = $junction->gtf_url();
	my $locus = $junction->locus_extended_100nt();
	my $igv_link = qq{<button class='igvIcon2' onclick='launch_igv_tool_rna("", "$bam_files,$gtf","$locus")' style="color:black"></button>};
	return $igv_link;
}

sub get_html_sashimi_plot {
	my ($self, $junction) = @_;
	my $sashimi_button;
	my $list_sashimi_plot_files = $junction->list_sashimi_plots();
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
		$sashimi_button .= qq{<image style="position:relative;width:200px;" loading="lazy" src="$pdf"></image>};
		$sashimi_button .= qq{</td><td style="padding-left:1px;"><span style="writing-mode:vertical-lr !important; font: 12px Verdana, sans-serif;letter-spacing: 1px;">Zoom</span></td></table> </button>};
		$sashimi_button .= qq{</center></};
	}
	else {
		$sashimi_button .= qq{<center>N.A.</center>};
	}
	return $sashimi_button;
}

sub panel_gene {
	my ($self, $hgene,$panel_id,$project_name,$patient) = @_;
	my $cgi  = new CGI();
	my $out;
	my $gene_id = $hgene->{id};
	my $gene;
	$gene = $patient->project->newGene($gene_id) if ($patient);
	my $vval ={};
	my $all_validations;
	($all_validations) = grep {$_ =~ /$gene_id/ } keys %{$patient->validations} if ($patient);
	if ($all_validations) {
		$vval =  $patient->validations->{$all_validations}->[0];
	}
	my $label_id = "label_".$hgene->{uid};
	my $max_score = $hgene->{max_score};
	my $bgcolor2 = "background-color:#607D8B;border-color:#607D8B";
	my $glyph = "";
	my $astyle = "";
	$astyle = "background-color:#E74C3C"  if $vval->{validation} and $vval->{validation}>4;
	my $gene_name = $hgene->{external_name};
	$out .=  $cgi->start_div({class=>" btn-group btn-xs ",style=>'position:relative;float:left;bottom:5px;'});
	my $in;
	if ($gene) { $in = $gene->omim_inheritance; }
	else { $in = $hgene->{omim_inheritance}; }
	$in ="" if $in eq "-";
	$in = "X-linked " if $in =~/X-linked/;
	my $pli = $hgene->{pLI};
	$pli *= 1.0 if $pli ne '-';
	my $bcolor = "grey";
	$bcolor = "green" if $max_score >= 0;
	$bcolor = "yellow" if $max_score >= 5;
	$bcolor = "orange" if $max_score >= 8;
	$bcolor = "coral" if $max_score >= 12;
	$bcolor = "red" if $max_score >= 14;
	my $div_id_gene = "div_splice_".$gene_name;
	$div_id_gene =~ s/\./_/g;
	$div_id_gene =~ s/ //g;
	my $cnv_status = "cnv_none";
	my $this_b_cmd = qq{collapse("$panel_id","$label_id")};
	if (exists $hgene->{specific_cmd}) {
		$this_b_cmd = $hgene->{specific_cmd};
	}
	if (exists $hgene->{collapse_with_id}) {
		my $this_collapse_id = $hgene->{collapse_with_id};
		$out .= qq{<div id="$div_id_gene" class="btn btn-brown btn-xs $bcolor $cnv_status" data-toggle='collapse' data-target="#$this_collapse_id" aria-expanded='false' aria-controls='$this_collapse_id' style="background-color:#4A4F53;border-top: 2px solid #4A4F53;border-bottom: 2px solid #4A4F53;border-right: 4px solid $bcolor;border-left: 4px solid $bcolor;$astyle;font-family: Verdana,Arial,sans-serif; text-shadow:1px 1px 2px black;position:relative;bottom:0px;min-width:150px;" onClick='$this_b_cmd'>  <font style='color:white;'><span id= "$label_id" class="glyphicon glyphicon-triangle-right" style="" aria-hidden="true"  style="float:left;top:4px;"></span> $gene_name<sup>&nbsp;$in</b></font></sup> $glyph}.$cgi->span({class=>"badge1 $bcolor"},$hgene->{max_score}).qq{</div>};
	}
	else {
		$out .= qq{<div id="$div_id_gene" class="btn btn-brown btn-xs $bcolor $cnv_status" style="background-color:#4A4F53;border-top: 2px solid #4A4F53;border-bottom: 2px solid #4A4F53;border-right: 4px solid $bcolor;border-left: 4px solid $bcolor;$astyle;font-family: Verdana,Arial,sans-serif; text-shadow:1px 1px 2px black;position:relative;bottom:0px;min-width:150px;" onClick='$this_b_cmd'>  <font style='color:white;'><span id= "$label_id" class="glyphicon glyphicon-triangle-right" style="" aria-hidden="true"  style="float:left;top:4px;"></span> $gene_name<sup>&nbsp;$in</b></font></sup> $glyph}.$cgi->span({class=>"badge1 $bcolor"},$hgene->{max_score}).qq{</div>};
	}
   	my $nbv = $hgene->{nb};
	my $omim = $hgene->{omim_id};
	if ($omim) { $out .=qq{<a class="btn btn-primary btn-xs" href="http://www.omim.org/entry/$omim" role="button" target="_blank" style="min-width:40px;$bgcolor2;text-shadow:1px 1px 1px black;color:white">Omim</a>}; }
	else { $out .=qq{<a class="btn btn-primary btn-xs"   style="$bgcolor2;min-width:40px"><i class="fa fa-minus"></i></a>}; } 
	#gtex Portal
	my ($gid,$t) = split("_",$hgene->{id});
	$out .=qq{<a class="btn btn-primary btn-xs" href="https://gtexportal.org/home/gene/$gid" role="button" target="_blank" style="min-width:40px;$bgcolor2;text-shadow:1px 1px 1px black;color:white">Gtex</a>};# if $omim ne "";
	my $oid = $hgene->{id};
	my $type ="green";
	#$type = "default" if $pli <= 0.1;
	$type = "orange" if $pli ne '-' and $pli >= 0.75;
	$type = "red" if $pli ne '-' and $pli >= 0.9;
	my $m = $hgene->{max_score};
	#$out .=qq{<a class="btn btn-primary btn-xs" href="https://gnomad.broadinstitute.org/gene/$oid" target="_blank" style="$bgcolor2;min-width:30px;height:22px;padding-top:3px;"><span class="badge" style="color:$type">$pli</span></a>};
 	my ($ensg, $cid) = split('_',$hgene->{id});
	my $b_id_pli = 'b_pli_'.$oid.'_'.$type;
	my $popup_pli = qq{<div data-dojo-type="dijit/Tooltip" data-dojo-props="connectId:'$b_id_pli',position:['above']"><span><b>pLI</b> Score</span></div>};
	$out .=qq{<a class="btn btn-primary btn-xs" href="https://gnomad.broadinstitute.org/gene/$ensg" target="_blank" style="$bgcolor2;min-width:30px"><span id="$b_id_pli" class="badge" style="color:$type">$pli</span>$popup_pli</a>};
 				
	my $panel_name1 = join("-",keys %{$hgene->{panels}});
	my $panel_list;
	foreach my $panel_name (sort keys %{$hgene}) {
 		$panel_list .= '  - '.$panel_name."<br>";
	}
	$panel_name1 = "Panel: ".scalar (keys %{$hgene->{panels}});
	$out .=qq{<a class="btn btn-primary btn-xs" href="#" role="button" style="top:-4px;$bgcolor2" onclick="document.getElementById('span_list_panels').innerHTML='$panel_list';dijit.byId('dialog_list_panels').show();"><p style="font-size:10px;text-shadow:0px 1px 1px #000;position:relative;bottom:-4px">$panel_name1</p></a>} if $panel_name1;
	my ($pheno,$nb_other_terms);
	if ($gene) { ($pheno,$nb_other_terms) = $gene->polyviewer_phentotypes(); }
	elsif (exists $hgene->{phenotypes}->{pheno} and exists $hgene->{phenotypes}->{nb_other_terms}) {
		$pheno = $hgene->{phenotypes}->{pheno};
		$nb_other_terms = $hgene->{phenotypes}->{nb_other_terms};
   	}
	my $color = '';
	$color = qq{ style = "color:#E74C3C"} if $pheno =~/intellectual/ or $pheno =~/mental/ or $pheno =~/retar/;
	my $jtid = 'zz'.time.rand(500);
	my $div_pheno = qq{<a aria-disabled="true" class="btn btn-primary btn-xs" href="#" role="button" style="text-align:left;font-family: proxima-nova, sans-serif;font-style:normal;text-shadow:0px 1px 1px #000000;position:relative;bottom:4px;font-size: 0.9em;color:white;$bgcolor2;">};
	if ($pheno) {
		if (length($pheno) > 70) { $pheno = substr($pheno, 0, 70).'...'; }
		if ($nb_other_terms > 0) { $pheno .= " <span style='color:#5cf0d3'>+ $nb_other_terms terms</span>"; }
		if ($project_name) {
			$div_pheno .= qq{<i class="fa fa-circle fa-xs" $color ></i> <span onclick="update_grid_gene_phenotypes(\'$gene_id\', \'$project_name\')";">$pheno</span>};
		}
		else {
			$div_pheno .= qq{<i class="fa fa-circle fa-xs" $color ></i> <span onclick="update_grid_gene_phenotypes(\'$gene_id\')";">$pheno</span>};
		}
	}
	$div_pheno .= qq{</a>};
	$out .= $div_pheno;
	$out .= $cgi->end_div();
	$out .=  $cgi->start_div({class=>" btn-group btn  ",style=>'position:relative;float:right;bottom:5px;'});
	my $tlabel = "label-grey";
	$tlabel = "label-warning" if exists $hgene->{clinvar_hgmd};
	$tlabel = "label-danger" if exists $hgene->{pathogenic};
	$out .= $cgi->span({class=>"label $tlabel" },qq{<span class="glyphicon glyphicon-alert text-alert " aria-hidden="true" ></span>});
	my $nbv2 = scalar(@{$hgene->{variants}});
	$out .=$cgi->span({class=>"label label-grey"},qq{<span class='badge badge-infos badge-xs ' style="color:#00C851"  >$nbv2 </span> });
	my $style ="";
	$out.= $cgi->end_div(); # end div lavel right 
	return $out;
}


1;