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
my $fork = $cgi->param('fork');
my $fileout = $cgi->param('fileout');
$fork = 1 unless $fork;

my $buffer = GBuffer->new;
my $project = $buffer->newProjectCache( -name => $project_name );


print $cgi->header('text/json-comment-filtered');
print "{\"progress\":\".";


my $h_patient_analisys;
foreach my $patient (@{$project->getPatients}) {
	$h_patient_analisys->{$patient->name()} = undef;
}

my $path_analisys = $project->get_path_rna_seq_root();
opendir my ($dir), $path_analisys;
my @found_files = readdir $dir;
closedir $dir;
foreach my $file (@found_files) {
	next if $file eq '.';
	next if $file eq '..';
	my $tmp_file = $file;
	$tmp_file =~ s/_vs_all//;
	my @lTmp = split('_', $tmp_file);
	if (exists $h_patient_analisys->{$lTmp[-1]}) {
		$h_patient_analisys->{$lTmp[-1]} = $file;
		next;
	}
	$tmp_file =~ s/RNAseqSEA_v1\.0\.8_//;
	my $test_name = join('_', @lTmp);
	if (exists $h_patient_analisys->{$tmp_file}) {
		$h_patient_analisys->{$tmp_file} = $file;
	}
}

my $analyse_name;
if ($use_patient) {
	$analyse_name = $h_patient_analisys->{$use_patient};
}
else {
	foreach my $patient_name (sort keys %$h_patient_analisys) {
		if ($h_patient_analisys->{$patient_name}) {
			$use_patient = $patient_name;
			$analyse_name = $h_patient_analisys->{$patient_name};
			last;
		}
	}
}
#warn $path_analisys;
#warn Dumper $h_patient_analisys;

my $html_table_patients = qq{<table data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="20" data-page-list="[10, 20, 50, 100, 200, 300]" data-resizable='true' id='table_id_patients' class='table table-striped' style='font-size:13px;'>};
$html_table_patients .= qq{<thead>};
$html_table_patients .= qq{<th data-field="name" data-sortable="true"><b>Patient Name</b></th>};
$html_table_patients .= qq{<th data-field="analyse" data-sortable="true"><b>Analyse Name</b></th>};
$html_table_patients .= qq{<th data-field="view"><b>View</b></th>};
$html_table_patients .= qq{<th data-field="bam_igv"><b>Add BAM in IGV</b></th>};
$html_table_patients .= qq{</thead>};
$html_table_patients .= qq{<tbody>};
foreach my $patient_name (sort keys %$h_patient_analisys) {
	my $analyse = $h_patient_analisys->{$patient_name};
	my $tr = qq{<tr>};
	$tr .= qq{<td>$patient_name</td>};
	$tr .= qq{<td>$analyse</td>};
	
	
#	warn "\n\n";
#	warn $path_analisys.'/'.$analyse.'/AllRes/';
#	warn '-> FOUND!' if (-d )
#	warn "\n\n";
	
	if ($patient_name eq $use_patient) { $tr .= qq{<td><button onClick="document.getElementById('a_RI').click();">View</button></td>}; }
	elsif ($analyse and -d $path_analisys.'/'.$analyse.'/AllRes/') { $tr .= qq{<td><button onClick="launch('$patient_name')">View</button></td>}; }
	else { $tr .= qq{<td>N.A.</td>}; }
	my $bam_file = $project->getPatient($patient_name)->bamUrl();
	my $igv_link = qq{<button class='igvIcon2' onclick='add_bam_igv("$bam_file")' style="color:black"></button>};
	$tr .= qq{<td>$igv_link</td>};
	$tr .= qq{</tr>};
	$html_table_patients .= $tr;
}
$html_table_patients .= qq{</tbody>};
$html_table_patients .= qq{</table>};

#die;


my $path_analysis = $project->get_path_rna_seq_analyse($analyse_name);
my $path_analysis_all = $path_analysis.'/AllRes/';

#warn $path_analysis_all;

confess("\n\nERROR: $path_analysis_all not found. Die\n\n") unless (-d $path_analysis_all);
my $file_all_RI = $path_analysis_all.'/AllresRI_f.txt';
confess("\n\nERROR: $file_all_RI not found. Die\n\n") unless (-e $file_all_RI);
my $file_all_SE = $path_analysis_all.'/AllresSE_f.txt';
confess("\n\nERROR: $file_all_SE not found. Die\n\n") unless (-e $file_all_SE);


my $h_header_filters_col;
#$h_header_filters_col->{Junc_RI} = "data-filter-control='select'";
$h_header_filters_col->{ENSID} = "data-filter-control='select'";
$h_header_filters_col->{Chr} = "data-filter-control='select'";
$h_header_filters_col->{Gene} = "data-filter-control='select'";
$h_header_filters_col->{Junc_Normale} = "data-filter-control='select'";
$h_header_filters_col->{Type} = "data-filter-control='select'";
$h_header_filters_col->{Com_ou_Spe} = "data-filter-control='select'";
$h_header_filters_col->{Sample} = "data-filter-control='select'";

my ($table_id_ri, $html_table_ri) = transform_results_file('RI', $file_all_RI, $path_analysis);
my ($table_id_se, $html_table_se) = transform_results_file('SE', $file_all_SE, $path_analysis);

my $html;
$html .= qq{<ul class="nav nav-tabs">};
$html .= qq{<li><a id='a_PATIENTS' data-toggle="tab" href="#PATIENTS"><span style="color:red;">LIST Patients</span></a></li>};
$html .= qq{<li><a id='a_RI' data-toggle="tab" href="#RI">RI Analysis</a></li>};
$html .= qq{<li><a id='a_SE' data-toggle="tab" href="#SE">SE Analysis</a></li>};
$html .= qq{</ul>};

my $onclick_ri = qq{document.getElementById("$table_id_ri").contentWindow.location.reload(true);};
my $onclick_se = qq{document.getElementById("$table_id_se").contentWindow.location.reload(true);};
my $onclick_patients = qq{document.getElementById("table_id_patients").contentWindow.location.reload(true);};

$html .= qq{<div class="tab-content">};
$html .= qq{<div id="PATIENTS" class="tab-pane" onClick=$onclick_patients>};
$html .= qq{$html_table_patients};
$html .= qq{</div>};
$html .= qq{<div id="RI" class="tab-pane" onClick=$onclick_ri>};
$html .= qq{$html_table_ri};
$html .= qq{</div>};
$html .= qq{<div id="SE" class="tab-pane" onClick=$onclick_se>};
$html .= qq{$html_table_se};
$html .= qq{</div>};
$html .= qq{</div>};

my $hash;
$hash->{html} = $html;
$hash->{analyse} = $analyse_name;
$hash->{patient_name} = $use_patient;
$hash->{tables_ids} = $table_id_ri.';'.$table_id_se.';table_id_patients';

#my @list_default_hidden_ri = ("b_ri_ENSID","b_ri_score","b_ri_minMoyNcountParJunc","b_ri_maxMoyNcountParJunc","b_ri_Sample","b_ri_ScoreCorrige","b_ri_pvalFisher","b_ri_Nctrls","b_ri_ratioScores.P_vs_Ctrl.","b_ri_deltaScores.P_vs_Ctrl.","b_ri_pvalCorr","b_ri_MoyScoreCtrlsCorr","b_ri_ratioScores.P_vs_Ctrl.Corr","b_ri_deltaScores.P_vs_Ctrl.Corr");
#$hash->{default_hidden_ri} = join(';', @list_default_hidden_ri);

printJson($hash);
exit(0);


sub transform_results_file {
	my ($type_analyse_name, $file_name, $path_analysis) = @_;
	my ($h_header, $list_res) = parse_results_file($file_name);
	my ($list_sort) = sort_results($list_res);
	my ($table_id, $html) = add_table_results($type_analyse_name, $h_header, $list_sort, $path_analysis);
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

sub add_table_results {
	my ($type_analyse_name, $h_header, $list_res, $path_analysis) = @_;
	my @list;
	if ($cgi->param('update_sashimi')) {
		my $pm = new Parallel::ForkManager($fork);
		foreach my $h_res (@$list_res) {
			my $pid = $pm->start and next;
			print ".";
			my $patient_name = $use_patient;
			my ($start, $end);
			if (exists $h_res->{Junc_RI_Start} and $h_res->{Junc_RI_End}) {
				$start = int($h_res->{Junc_RI_Start}) - 100;
				$end = int($h_res->{Junc_RI_End}) + 100;
			}
			else {
				$start = int($h_res->{Junc_SE_Start}) - 100;
				$end = int($h_res->{Junc_SE_End}) + 100;
			}
			my $locus_extended = $h_res->{Chr}.':'.$start.'-'.$end;
			my $sashimi_plot_file = get_sashimi_plot_file($path_analysis, $h_res->{ENSID}, $project->getPatient($patient_name), $locus_extended, $h_res->{score});
			my $i = 1;
			while ($i < 5) {
				$start -= (1000*$i);
				$end += (1000*$i);
				my $locus_extended = $h_res->{Chr}.':'.$start.'-'.$end;
				my $sashimi_plot_file = get_sashimi_plot_file($path_analysis, $h_res->{ENSID}, $project->getPatient($patient_name), $locus_extended, $h_res->{score});
				$i++;
			}
			$pm->finish();
		}
		$pm->wait_all_children();
	}
	
	foreach my $h_res (@$list_res) {
		print ".";
		my $tr = qq{<tr>};
		my $locus;
		if (exists $h_res->{Junc_RI_Start}) {
			$locus = $h_res->{Chr}.':'.($h_res->{Junc_RI_Start} - 100).'-'.($h_res->{Junc_RI_End} + 100);
		}
		elsif (exists $h_res->{Junc_SE_Start}) {
			$locus = $h_res->{Chr}.':'.($h_res->{Junc_SE_Start} - 100).'-'.($h_res->{Junc_SE_End} + 100);
		}
		my $patient_name = $use_patient;
		my $gn = $project->getVersion();
		my $v1 = "/";
		my $bam_files = $project->getPatient($patient_name)->bamUrl();
		
		my $nb_control = 0;
		foreach my $other_pat (@{$project->getPatients()}) {
			next if $other_pat->name eq $patient_name;
			$bam_files .= ";".$other_pat->bamUrl();
			$nb_control++;
			last if $nb_control == 3;
		}
		
		my $igv_link = qq{<button class='igvIcon2' onclick='launch_web_igv_js("$project_name","$patient_name","$bam_files","$locus","v1","$gn")' style="color:black"></button>};
		$tr .= qq{<td>$igv_link</td>};
		
		my $path_svg = $path_analysis.'/'.$h_res->{ENSID}.'/SpliceRes/Figs/';
		my ($svg_legend, $locus_extended);
		if (exists $h_res->{Junc_RI_Start} and $h_res->{Junc_RI_End}) {
			$locus = $h_res->{Chr}.':'.$h_res->{Junc_RI_Start}.'-'.$h_res->{Junc_RI_End};
			$locus_extended = $h_res->{Chr}.':'.(int($h_res->{Junc_RI_Start})-100).'-'.(int($h_res->{Junc_RI_End})+100);
		}
		else {
			$locus = $h_res->{Chr}.':'.$h_res->{Junc_SE_Start}.'-'.$h_res->{Junc_SE_End};
			$locus_extended = $h_res->{Chr}.':'.(int($h_res->{Junc_SE_Start})-100).'-'.(int($h_res->{Junc_SE_End})+100);
		}
		$svg_legend = 'Patient '.$patient_name.' - '.$h_res->{Gene}.' ('.$h_res->{ENSID}.') - '.$locus;
		my $svg_patient = $path_svg.'/juncPairPos_'.$h_res->{ENSID}.'_'.$h_res->{Chr}.'_'.$patient_name.'_'.$project_name.'.svg';
		my $color;
		$color = 'lightgrey';
		$color = 'palegreen' if ($h_res->{score} > 1);
		$color = 'moccasin' if ($h_res->{score} > 10);
		$color = 'lightcoral' if ($h_res->{score} > 100);
		if (-e $svg_patient) {
			$svg_patient =~ s/\/\//\//g;
			$svg_patient =~ s/\/data-isilon\/sequencing\/ngs/\/NGS/;
			$tr .= qq{<td><button type="button" class="btn btn-default" onClick="zoom_file('$svg_legend', '$svg_patient')" style="text-align:center;padding:2px;background-color:$color;"><img src="$svg_patient" alt="Pb" width="55px"></img><span style="top:3px;width:20px;height:20px;" class="glyphicon glyphicon-zoom-in" aria-hidden="true"></span></button></td>};
		}
		else {
			$tr .= qq{<td>N.A.</td>};
		}
		my $ensg = $h_res->{ENSID};
		
#		$tr .= qq{<td><center><button type="button" class="btn btn-default" onClick="view_sashimi('$analyse_name','$patient_name','$ensg','$locus_extended','$score')" style="text-align:center;padding:2px;background-color:$color;"><span style="top:3px;width:20px;height:20px;" class="glyphicon glyphicon-zoom-in" aria-hidden="true"></span></button></center></td>};
		
		#my $sashimi_plot_file = get_sashimi_plot_file($path_analysis, $h_res->{ENSID}, $project->getPatient($patient_name), $locus_extended, $h_res->{score});
		
		my $list_sashimi_plot_files = get_sashimi_plot_list_files($path_analysis, $h_res->{ENSID}, $project->getPatient($patient_name), $locus_extended, $h_res->{score});
		if ($list_sashimi_plot_files and -e $list_sashimi_plot_files->[0]) {
			$tr .= qq{<td><center>};
			my @lFiles;
			foreach my $sashimi_plot_file (@$list_sashimi_plot_files) {
				$sashimi_plot_file =~ s/\/\//\//g;
				$sashimi_plot_file =~ s/\/data-isilon\/sequencing\/ngs/\/NGS/;
				push(@lFiles, $sashimi_plot_file);
			}
			my $files = join(';', @lFiles);
			$tr .= qq{<button type="button" class="btn btn-default" onClick="view_pdf_list_files('$files')" style="text-align:center;padding:2px;background-color:$color;"><span style="top:3px;width:20px;height:20px;" class="glyphicon glyphicon-zoom-in" aria-hidden="true"></span></button>};
			$tr .= qq{</center></td>};
		}
		else {
			$tr .= qq{<td><center>.</center></td>};
		}
		
		
		foreach my $nb (sort {$a <=> $b} keys %{$h_header}) {
			my $cat = $h_header->{$nb};
			my $value = $h_res->{$cat};
			$tr .= qq{<td>$value</td>};
		}
		
		$tr .= qq{</tr>};
		push(@list, $tr);
	}
	my $table_id = 'table_'.$type_analyse_name;
	
	my $b_filters_id = 'button_filters_'.$type_analyse_name;
	my $html = qq{<h3><div style="text-align: center;">$analyse_name - <span style='color:red'>$type_analyse_name analyse</span></div></h3>};
	$html .= qq{<div style="font-size:10px;">};
	$html .= qq{<button style="background-color:white;" data-toggle="collapse" data-target="#$b_filters_id"><b><i>Show/Hide column(s)</b></button> };
	my @lcolumns;
	my $class_name = "toggle-vis-$type_analyse_name";
	foreach my $nb (sort {$a <=> $b} keys %{$h_header}) {
		my $cat = $h_header->{$nb};
		push(@lcolumns, qq{<a style="" id="b_ri_$class_name" class="$class_name" name="$cat" status="show">$cat</a>});
	}
	$html .= qq{<div id="$b_filters_id" class="collapse">}.join(' | ', @lcolumns).qq{</div>};
	$html .= qq{</i></div>};
	
	$html .= qq{<table data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="20" data-page-list="[10, 20, 50, 100, 200, 300]" data-resizable='true' id='$table_id' class='table table-striped' style='font-size:11px;'>};
	$html .= qq{<thead>};
	$html .= qq{<th data-field="igv"><b>IGV</b></th>};
	$html .= qq{<th data-field="figs"><b>Figure</b></th>};
	$html .= qq{<th data-field="sashimi"><b>Sashimi Plot</b></th>};
	foreach my $nb (sort {$a <=> $b} keys %{$h_header}) {
		my $cat = $h_header->{$nb};
		my $filter;
		$filter = $h_header_filters_col->{$cat} if exists ($h_header_filters_col->{$cat});
		$html .= qq{<th data-field="$cat" $filter data-sortable="true"><b>$cat</b></th>};
	}
	$html .= qq{</thead>};
	$html .= qq{<tbody>};
	foreach my $tr (@list) {
		$html .= $tr;
	}
	$html .= qq{</tbody>};
	$html .= qq{</table>};
	
	return ($table_id, $html);
}

sub sort_results {
	my ($list_res) = @_;
	my ($h, @list_sort);
	foreach my $this_h (@$list_res) {
		my $score = $this_h->{score} * 100;
		my $nb;
		$nb = $this_h->{Junc_RI_Count} if (exists $this_h->{Junc_RI_Count});
		$nb = $this_h->{Junc_SE_Count} if (exists $this_h->{Junc_SE_Count});
		my ($start, $end);
		if (exists $h->{Junc_RI_Start}) {
			$start = $h->{Junc_RI_Start};
			$end = $h->{Junc_RI_End};
		}
		else {
			$start = $h->{Junc_SE_Start};
			$end = $h->{Junc_SE_End};
		}
		my $name = $this_h->{Sample}.'_'.$this_h->{Junc_RI}.'_'.$start.'_'.$end;
		$h->{$score}->{$nb}->{$name} = $this_h;
	}
	foreach my $score (sort {$b <=> $a} keys %$h) {
		foreach my $nb (sort {$b <=> $a} keys %{$h->{$score}}) {
			foreach my $name (sort keys %{$h->{$score}->{$nb}}) {
				push(@list_sort, $h->{$score}->{$nb}->{$name});
			}
		}
	}
	return \@list_sort;
}

sub parse_results_file {
	my ($file_name) = @_;
	open (FILE, $file_name);
	my ($h_header, @l_res);
	my $i = 0;
	while (<FILE>) {
		my $line = $_;
		chomp($_);
		if ($i == 0) { $h_header = parse_header($line); }
		else {
			my $h_res;
			my $nb_col = 0;
			my @l_col = @{parse_line($line)};
			if ($l_col[0] ne 'NA') {
				foreach my $res (@l_col) {
					my $cat = $h_header->{$nb_col};
					confess("error parsing file...") unless ($cat);
					$h_res->{$cat} = $res;
					$nb_col++;
				}
				push(@l_res, $h_res);
			}
		}
		$i++;
	}
	close (FILE);
	return ($h_header, \@l_res);
}

sub parse_line {
	my ($line) = @_;
	chomp($line);
	my @lCol = split(' ', $line);
	return \@lCol;
}

sub parse_header {
	my ($line) = @_;
	my $h_header;
	my $nb_col = 0;
	foreach my $cat (@{parse_line($line)}) {
		$h_header->{$nb_col} = $cat;
		$nb_col++;
	}
	return $h_header;
}

sub printJson {
	my ($hashRes) = @_;
	my $json_encode = encode_json $hashRes;
	print ".\",";
	$json_encode =~ s/{//;
	print $json_encode;
	if ($fileout) {
		open (FILE, ">$fileout");
		print FILE $json_encode;
		close (FILE);
	}
	exit(0);
}
