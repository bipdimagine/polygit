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
use File::Find;
use Getopt::Long;
use Data::Dumper;
use JSON;
use xls_export;
use session_export;

my $cgi = new CGI;
print $cgi->header('text/json-comment-filtered');


my $force = $cgi->param('force');
my $run_name = $cgi->param('run');

my ($h_files, $h_files_date);
my $origin_path = '/data-isilon/sequencing/ngs/demultiplex/';
my ($list_paths_found) = list_html_files_in_dir($origin_path);
my $buffer = new GBuffer;

foreach my $this_path (@$list_paths_found) {
	if ($run_name) {
		next unless ($this_path =~ /$run_name/);
	}
	if (-e $this_path.'/Demultiplex_Stats.csv') {
		mkdir $this_path.'/html/' if (not -d $this_path.'/html/');
		if (-e $this_path.'/Demultiplex_Stats.csv' && -e $this_path.'/Top_Unknown_Barcodes.csv') {
			my @lFiles = ('Demultiplex_Stats.csv', 'Top_Unknown_Barcodes.csv');
			my $json_file = convert_csv_to_json($this_path, \@lFiles);
			my $json_path = $json_file;
			my $use_force;
			$use_force = 1 if (not -e $json_file);
			$json_file =~ s/\/\//\//g;
			$json_file =~ s/.+\/ngs\//https:\/\/www.polyweb.fr\/NGS\//;
			my $name = $json_file;
			if ($name =~ /.+\/NGS\/demultiplex\/(.+)\/json.+/) {
				$name = $1;
			}
			next if ($name eq 'test_xths');
			next if ($name =~ /run1[0-9]+/);
			my $cmd = 'https://www.polyweb.fr/polyweb/demultiplex_view/demultiplex_view_run.html?name='.$name.'&json='.$json_file;
			$cmd .= '&force=1' if ($use_force);
			add_file_json($cmd, $json_path, $name);
		}
#		elsif (-e $this_path.'/Demultiplex_Stats.csv') {
#			convert_csv_to_html($this_path.'/Demultiplex_Stats.csv', $this_path.'/html/laneBarcode.html') if (not -e $this_path.'/html/laneBarcode.html' || $force);
#			add_file_html($this_path.'/html/laneBarcode.html', 'laneBarcode.html');
#		}
#		elsif (-e $this_path.'/Top_Unknown_Barcodes.csv') {
#			convert_csv_to_html($this_path.'/Top_Unknown_Barcodes.csv', $this_path.'/html/top_unknown_barcodes.html') if (not -e $this_path.'/html/top_unknown_barcodes.html' || $force);
#			add_file_html($this_path.'/html/top_unknown_barcodes.html', 'top_unknown_barcodes.html');
#		}
	}
	elsif (-d $this_path.'/html/') {
		opendir my ($dir), $this_path.'/html/';
		my @found_files = readdir $dir;
		closedir $dir;
		my (@lDir, @lFiles);
		foreach my $file (@found_files) {
			next if $file eq '.';
			next if $file eq '..';
			if (-d $this_path.'/html/'.$file) {
				if (-e $this_path.'/html/'.$file.'/all/laneBarcode.html') { add_file_html($this_path.'/html/'.$file.'/all/laneBarcode.html', 'laneBarcode.html'); }
				elsif (-e $this_path.'/html/'.$file.'/all/all/laneBarcode.html') { add_file_html($this_path.'/html/'.$file.'/all/all/laneBarcode.html', 'laneBarcode.html'); }
				elsif (-e $this_path.'/html/'.$file.'/all/all/all/laneBarcode.html') { add_file_html($this_path.'/html/'.$file.'/all/all/all/laneBarcode.html', 'laneBarcode.html'); }
				else { list_html_files_in_dir($this_path.'/html/'); }
			}
		}
	}
	else {
		list_html_files_in_dir($this_path);
	}
}

my $html_table = qq{<table id="table_demultiplex" data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="20" data-page-list="[10, 20, 50, 100, 200, 300]" data-resizable='true' id='table_id_patients' class='table table-striped' style='font-size:13px;'>};
$html_table .= qq{<thead>};
#$html_table .= qq{<th data-field="date" data-sortable="true"><b>Date</b></th>};
$html_table .= qq{<th data-field="path" data-sortable="false" data-filter-control='input'><b>Path</b></th>};
$html_table .= qq{<th data-field="file" data-sortable="false" data-filter-control='input'><b>File</b></th>};
$html_table .= qq{<th data-field="view"><b>View</b></th>};
$html_table .= qq{</thead>};
$html_table .= qq{<tbody>};



foreach my $date (reverse sort keys %$h_files_date) {
	my $file = $h_files_date->{$date};
	my $path_file_origin = $h_files->{$file};
	$path_file_origin =~ s/$origin_path/https:\/\/www.polyweb.fr\/NGS\/demultiplex/;
	my $path_file = $h_files->{$file};
	my ($substring, $color);
	if ($path_file =~ /$origin_path/) {
		$path_file =~ s/$origin_path/DEMULTIPLEX/;
		$path_file =~ s/$file//;
		$path_file =~ s/\/\///g;
		my @l_path = split('/', $path_file);
		$path_file = $l_path[0].'/'.$l_path[1];
		if ($path_file =~ /(run[0-9]+.*$)/) {
			$substring = $1;
			$color = 'red';
		}
		elsif ($path_file =~ /\/(.+)$/) {
			$substring = $1;
			$color = 'green';
		}
		$path_file =~ s/$substring/\/<span style='color:$color;'>$substring<\/span>/;
		$path_file =~ s/\/\//\//;
		$file =~ s/$substring/<span style='color:$color;'>$substring<\/span>/;
		$file =~ s/(.+<\/span>).+-([a-zA-Z]+\.html)$/$1.$2/;
		
		$substring =~ s/_/\./;
		$file =~ s/$substring/<span style='color:$color;'>$substring<\/span>/;
		$file =~ s/(.+<\/span>).+-([a-zA-Z]+\.html)$/$1.$2/;
		
		$file =~ s/\/\//\//;
		$file =~ s/^\.//;
		$file =~ s/^-//;
	}
	else {
		$path_file = qq{DRAGEN/<span style='color:red;'>$file</span>};
		$file = qq{<span style='color:red;'>$file</span>};
	}
	
	
	my $tr = qq{<tr>};
	#$tr .= qq{<td>$date</td>};
	$tr .= qq{<td>$path_file</td>};
	$tr .= qq{<td>$file</td>};
	
	my $name = $file;
	$name =~ s/\.html//;
	my $b_view = qq{ <a title="$name" target="_blank" href="$path_file_origin">VIEW</a> };
	
	#my $b_view = qq{<button onclick='' style="color:black">View</button>};
	$tr .= qq{<td>$b_view</td>};
	$tr .= qq{</tr>};
	$html_table .= $tr;
}
$html_table .= qq{</tbody>};
$html_table .= qq{</table>};

my $hash;
$hash->{html} = $html_table;
$hash->{table_id} = 'table_demultiplex';

my $json_encode = encode_json $hash;
print $json_encode;
exit(0);

sub check_run_id_from_sample_name {
	my $sample_name = shift;
	my $h_runs = $buffer->getQuery->getHashRunIdFromSampleName($sample_name);
	return $h_runs;
}

sub convert_csv_to_json {
	my ($path, $list_files) = @_;
	my $path_out = $path.'/json/';
	unless (-d $path_out) {
		mkdir $path_out;
		`chmod 777 $path_out`;
	}
	my $file_out = "$path_out/demultiplex.json";
	return $file_out if (-e $file_out and not $force);
	return $file_out if (not $run_name);
	my $h;
	
	foreach my $file (sort @$list_files) {
		my $csv_file = $path.'/'.$file;
		my $table_id = 'table_'.$file;
		$table_id =~ s/\.csv//;
		
		my $sorted_col = "data-sort-name='n_reads' data-sort-order='desc'";
		#$sorted_col = "data-sort-name='p_of_unknown_barcodes' data-sort-order='desc'" if (lc($table_id) =~ /top_/);
		
		open (CSV, "$csv_file");
		my $i = 0;
		my $has_RC;
		my $max_nb_header;
		my $nb_samples;
		my $total_reads;
		while (<CSV>) {
			chomp($_);
			my $line = $_;
			my @lCol = split(',', $line);
			if ($i == 0) {
				$h->{$file}->{header} = \@lCol;
				$max_nb_header = scalar(@{$h->{$file}->{header}});
			}
			else {
				my $j = 0;
				my $id = $lCol[0].'_'.$lCol[1];
				if ($id =~ /(.+)_RC/) {
					$has_RC ++ if (exists $h->{$file}->{values}->{$1});
				}
				my $sample_id;
				foreach my $value (@lCol) {
					my $cat = $h->{$file}->{header}[$j];
					$value = $value * 100 if ($cat =~ /\%/);
					$h->{$file}->{values}->{$id}->{$cat} = $value;
					$h->{$file}->{values}->{$id}->{$j} = $value;
					if ($cat eq 'SampleID') {
						$sample_id = $value;
						$nb_samples++ if (not lc($sample_id) eq 'undetermined' and not lc($sample_id) =~ /_rc/);
						my ($run_id, @l_runs_id);
#						if (exists $h->{$file}->{runs_ids} and exists $h->{$file}->{runs_ids}->{$value}) {
#							@l_runs_id = sort keys %{$h->{$file}->{runs_ids}->{$value}};
#						}
#						else {
#							my $h_runs_infos = check_run_id_from_sample_name($value) if (lc($value) ne 'undetermined');
#							@l_runs_id = sort keys %{$h_runs_infos};
#							foreach my $r (@l_runs_id) {
#								$h->{$file}->{runs_ids}->{by_samples}->{$value}->{$r} = undef;
#								$h->{$file}->{runs_ids}->{by_runs}->{$r}->{$value} = undef;
#							}
#						}
						my $runs = join(' ', reverse @l_runs_id);
						$runs = 'no_run_info' unless ($runs);
						
						
						$h->{$file}->{values}->{$id}->{'RunID'} = $runs;
						$h->{$file}->{values}->{$id}->{$max_nb_header} = $runs;
					}
					if ($cat eq '# Reads') {
						$total_reads += int($value)  if (not lc($sample_id) =~ /undetermined/ and not $sample_id =~ /_RC$/);
						$h->{$file}->{values}->{$id}->{'# Reads Norm Only'} = int($value);
						$h->{$file}->{values}->{$id}->{$max_nb_header+1} = int($value);
						if (exists $h->{$file}->{runs_ids}->{by_samples}->{$sample_id}) {
							foreach my $run_id (keys %{$h->{$file}->{runs_ids}->{by_samples}->{$sample_id}}) {
								$h->{$file}->{runs_ids}->{by_samples}->{$sample_id}->{$run_id} = int($value);
							}
						}
					}
					
					$j++;
				}
			}
			$i++;
		}
		close(CSV);
		
		my ($is_demultiplex_file, $is_top_barcodes_file);
		$is_demultiplex_file = 1 if (lc($file) =~ /demultiplex/);
		$is_top_barcodes_file = 1 if (lc($file) =~ /top_/);

		if ($is_demultiplex_file) {
			foreach my $run_id (keys %{$h->{$file}->{runs_ids}->{by_runs}}) {
				my @lSamples = keys %{$h->{$file}->{runs_ids}->{by_runs}->{$run_id}};
				my $nb_pat = 0;
				my $reads_total = 0;
				foreach my $sample (@lSamples) {
					next if (lc($sample) =~ /undetermined/);
					$nb_pat++;
					$reads_total += $h->{$file}->{runs_ids}->{by_samples}->{$sample}->{$run_id};
				}
				$h->{$file}->{runs_ids}->{by_runs}->{$run_id}->{nb_samples} = $nb_pat;
				$h->{$file}->{runs_ids}->{by_runs}->{$run_id}->{mean_reads} = $reads_total / $nb_pat;
			}
			push(@{$h->{$file}->{header}}, 'RunID');
		}
		
		if ($has_RC) {
			push(@{$h->{$file}->{header}}, '# Reads Norm Only');
		}
		if ($is_demultiplex_file) {
			push(@{$h->{$file}->{header}}, 'Seems Ok');
			$sorted_col = "data-sort-name='seems_ok' data-sort-order='desc'";
		}
		
		my $html = qq{<table style="width:100%;" $sorted_col data-filter-control='true'data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="10" data-page-list="[10, 15, 20,30, 50, 100, 200, 300]" data-resizable='true' id='$table_id' class='table table-striped sortable-table' style='font-size:13px;'>};
		$html .= "<thead>\n";
		foreach my $cat (@{$h->{$file}->{header}}) {
			my $cat_name = $cat;
			$cat_name =~ s/# /n_/;
			$cat_name =~ s/% /p_/;
			$cat_name =~ s/ /_/g;
			my $sortable;
			if ($cat =~ /[#%]/) {
				if (lc($table_id) =~ /top_unknown/) {
					if (lc($cat_name) eq 'p_of_unknown_barcodes') {
						$sortable = qq{ data-sortable='true' data-order='desc' class='numeric-sort'};
					}
					else {
						$sortable = qq{ data-sortable='true' class='numeric-sort'};
					}
				}
				if (lc($table_id) =~ /demultiplex/) {
					if ($has_RC) {
						if (lc($cat_name) eq 'n_reads_norm_only') {
							$sortable = qq{ data-sortable='true' data-order='desc' class='numeric-sort'};
						}
					}
					else {
						if (lc($cat_name) eq 'n_reads') {
							$sortable = qq{ data-sortable='true' data-order='desc' class='numeric-sort'};
						}
						else {
							$sortable = qq{ data-sortable='true' class='numeric-sort'};
						}
					}
					
				}
			}
			elsif ($cat eq 'RunID' or $cat eq 'Seems Ok') { $sortable = qq{ data-sortable='true' data-filter-control='select' }; }
			else { $sortable = qq{ data-sortable='true' data-filter-control='input' }; }
			$html .= "<th $sortable data-field='".lc($cat_name)."' ><center><b>$cat</b></center></th>\n";
		}
		
		$html .= "</thead>\n";
		$html .= "<tbody>\n";
		
		my $nb_values = scalar(keys %{$h->{$file}->{values}});
		my ($mean_reads, $limit_errors_reads_min, $limit_errors_reads_max);
		
		foreach my $id (sort keys %{$h->{$file}->{values}}) {
			next if ($id =~ /.+_RC/);
			
			my $probably_run_id;
			my $sample_id = $h->{$file}->{values}->{$id}->{SampleID};
			if (exists $h->{$file}->{runs_ids}->{by_samples}->{$sample_id}) {
				my $h_runs_samples;
				foreach my $run_id (keys %{$h->{$file}->{runs_ids}->{by_samples}->{$sample_id}}) {
					my $this_nb_samples = $h->{$file}->{runs_ids}->{by_runs}->{$run_id}->{nb_samples};
					$h_runs_samples->{$this_nb_samples}->{$run_id} = undef;
				}
				my @lNb_samples = sort {$a <=> $b} keys %{$h_runs_samples};
				my $max_samples = $lNb_samples[-1];
				my @lProb_runs = sort {$a <=> $b} keys %{$h_runs_samples->{$max_samples}};
				$probably_run_id = $lProb_runs[-1];
				$mean_reads = $h->{$file}->{runs_ids}->{by_runs}->{$probably_run_id}->{mean_reads};
			}
			if (exists $h->{$file}->{values}->{$id}->{'RunID'} && $h->{$file}->{values}->{$id}->{'RunID'} eq 'no_run_info') {
				$mean_reads = ($total_reads / $nb_samples);
			}
			$limit_errors_reads_min = $mean_reads * 0.5;
			$limit_errors_reads_max = $mean_reads * 1.5;
			
			my $tr_id = 'tr_'.$file.'_'.$id;
			my $html_tr;
			my $j = 0;
			my $index_barcode_sequence;
			my $with_problem = 0;
			my $with_ambigous = 0;
			while ($j >= 0) {
				my $cat = $h->{$file}->{header}[$j];
				$index_barcode_sequence = $j if (lc($cat) eq 'index');
				my $td_id = 'td_'.$file.'_'.$id.'_'.$j;
				my $value = $h->{$file}->{values}->{$id}->{$j};
				if ($cat ne '# Reads Norm Only' and $cat =~ /[#%]/ && $has_RC) {
					if ($cat eq '# Reads') {
						$h->{$file}->{header}[$j] = '# Reads Only';
					}
					my $value_RC = $h->{$file}->{values}->{$id.'_RC'}->{$j};
					$html_tr .= "<td id=\"$td_id\"><center>";
					#$html_tr .= qq{$total_reads / $nb_samples<br>};
					$html_tr .= "<table><tr><td><b>Norm</b></td><td style='padding-left:5px;'><center>$value</center></td></tr><tr><td><b>RC</b></td><td style='padding-left:5px;'><center>$value_RC</center></td></tr></table>";
					$html_tr .= "</center></td>\n";
					if ($is_demultiplex_file and lc($cat) eq '# reads') {
						$with_problem++ if ($value_RC >= $value);
					}
				}
				else {
					if ($is_demultiplex_file and $has_RC and $cat eq '# Reads Norm Only') {
						if (int($value) >= $limit_errors_reads_max) {
							$with_problem++;
							$h->{$file}->{values}->{$id}->{$j+1} = 'PB +++ reads';
							$h->{$file}->{values}->{$id}->{'Seems Ok'} = 'PB +++ reads';
						}
						elsif (int($value) <= $limit_errors_reads_min) {
							$with_problem++;
							$h->{$file}->{values}->{$id}->{$j+1} = 'PB --- reads';
							$h->{$file}->{values}->{$id}->{'Seems Ok'} = 'PB --- reads';
						}
						else {
							if ($with_problem > 0) {
								$h->{$file}->{values}->{$id}->{$j+1} = 'PB RC reads';
								$h->{$file}->{values}->{$id}->{'Seems Ok'} = 'PB RC reads';
							}
							else {
								$h->{$file}->{values}->{$id}->{$j+1} = 'OK';
								$h->{$file}->{values}->{$id}->{'Seems Ok'} = 'OK';
							}
						}
					}
					if ($is_demultiplex_file and not $has_RC and $cat eq '# Reads') {
						if (int($value) >= $limit_errors_reads_max) {
							$with_problem++;
							$h->{$file}->{values}->{$id}->{'Seems Ok'} = 'PB +++ reads';
						}
						if (int($value) <= $limit_errors_reads_min) {
							$with_problem++;
							$h->{$file}->{values}->{$id}->{'Seems Ok'} = 'PB --- reads';
						}
					}
					
					if ($is_demultiplex_file and $cat eq 'RunID') {
						if ($probably_run_id) {
							$value = $probably_run_id;
						}
						else {
							$value = 'mean_reads: '.int($mean_reads);
						}
						if (not $has_RC) {
							if ($with_problem > 0) {
								$h->{$file}->{values}->{$id}->{$j+1} = $h->{$file}->{values}->{$id}->{'Seems Ok'};
							}
							else {
								$h->{$file}->{values}->{$id}->{$j+1} = 'OK';
								$h->{$file}->{values}->{$id}->{'Seems Ok'} = 'OK';
							}
						}
					}
					
					
					$html_tr .= "<td id=\"$td_id\">".$value."</td>\n";
				}
				if ($is_top_barcodes_file and lc($cat) eq '% of unknown barcodes') {
					if ($value >= 15) {
						if ($h->{$file}->{values}->{$id}->{$index_barcode_sequence} =~ /^A+$/) { $with_ambigous++; }
						elsif ($h->{$file}->{values}->{$id}->{$index_barcode_sequence} =~ /^T+$/) { $with_ambigous++; }
						elsif ($h->{$file}->{values}->{$id}->{$index_barcode_sequence} =~ /^G+$/) { $with_ambigous++; }
						elsif ($h->{$file}->{values}->{$id}->{$index_barcode_sequence} =~ /^C+$/) { $with_ambigous++; }
						else { $with_problem++; }
					}
				}
				
				$j++;
				$j = -1 if (not exists $h->{$file}->{header}[$j] or  not exists $h->{$file}->{values}->{$id}->{$j});
			}
			$html_tr .= "</tr>\n";
			if ($with_problem > 0)     { $html_tr = "<tr style='text-align:center;background-color:red;color:white;' id=\"$tr_id\">\n".$html_tr; }
			elsif ($with_ambigous > 0) { $html_tr = "<tr style='text-align:center;background-color:orange;color:white;' id=\"$tr_id\">\n".$html_tr; }
			elsif (lc($id) =~ /undetermined/) { $html_tr = "<tr style='text-align:center;color:red;' id=\"$tr_id\">\n".$html_tr; }
			else { $html_tr = "<tr style='text-align:center;' id=\"$tr_id\">\n".$html_tr; }
			$html .= $html_tr;
			$with_problem = 0;
			$with_ambigous = 0;
		}
		
		$html .= "</tbody>\n";
		$html .= "</table>\n";
		$h->{$file}->{html} = $html;
		$h->{$file}->{table_id} = $table_id;
	}
	open (JSON, ">$file_out");
	print JSON encode_json $h;
	close (JSON);
	return $file_out;
}

sub convert_csv_to_html {
	my ($csv_file, $html_file) = @_;
	open (HTML, ">$html_file");
	open (CSV, "$csv_file");
	print HTML qq{<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">}."\n";
	print HTML qq{<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">}."\n";
	print HTML qq{<body><table width="100%" class="table table-striped table-hover table-condensed">}."\n";
	my $i = 0;
	while (<CSV>) {
		chomp($_);
		my $line = $_;
		my @lCol = split(',', $line);
		if ($i == 0) {
			print HTML "<thead>\n";
			foreach my $cat (@lCol) {
				print HTML "<th data-field='".lc($cat)."' ><center><b>$cat</b></center></th>\n";
			}
			print HTML "</thead>\n";
			print HTML "<thead>\n";
			print HTML "<tbody>\n";
		}
		else {
			print HTML "<tr>\n";
			foreach my $val (@lCol) {
				print HTML "<td><center>".$val."</center></td>\n";
			}
			print HTML "</tr>\n";
		}
		$i++;
	}
	print HTML "</tbody>\n";
	print HTML "</table>\n";
	close (CSV);
	close(HTML);
}

sub list_html_files_in_dir {
	my ($path) = @_;
	opendir my ($dir), $path;
	my @found_files = readdir $dir;
	closedir $dir;
	my (@lDir, @lFiles);
	foreach my $file (@found_files) {
		next if $file eq '.';
		next if $file eq '..';
		my $path_file = $path.'/'.$file;
		if (-d $path_file) {
			push(@lDir, $path_file);
		}
		elsif (-e $path_file) {
			add_file_html($path_file, $file);
		}
	}
	return (\@lDir);
}

sub add_file_json {
	my ($url, $path, $name) = @_;
	$h_files->{$name} = $url;
	my $date_stat;
	if (-e $path) {
		$date_stat = (stat ($path))[9];
		if (exists $h_files_date->{$date_stat}) {
			my $is_ok = 0;
			while ($is_ok == 0) {
				$date_stat++;
				next if (exists $h_files_date->{$date_stat});
				$is_ok = 1;
			}
		}
	}
	else { $date_stat = '99999999'; }
	$h_files_date->{$date_stat} = $name;
}

sub add_file_html {
	my ($path, $file) = @_;
	next unless ($file =~ /\.html/);
	next if ($file =~ /tree\.html/);
	my $name = $path;
	$name =~ s/$origin_path//;
	$name =~ s/$file//;
	$name =~ s/html//;
	$name =~ s/.all.all.all/.all/;
	$name =~ s/.all/.all/;
	$name =~ s/\/\//\//;
	my @ltmp = split('/', $name);
	$h_files->{join('.', @ltmp).'-'.$file} = $path;
	$h_files_date->{(stat ($path))[9]} = join('.', @ltmp).'-'.$file;
}
