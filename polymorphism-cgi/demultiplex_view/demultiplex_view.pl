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
use Spreadsheet::WriteExcel;
use POSIX;

my $cgi = new CGI;
my $force = $cgi->param('force');
my $run_name = $cgi->param('run');
my $export_xls = $cgi->param('export_xls');

my ($h_files, $h_files_date);
my $origin_path = '/data-isilon/sequencing/ngs/demultiplex/';
my ($list_paths_found) = list_html_files_in_dir($origin_path);
my $buffer = new GBuffer;


if ($run_name and $export_xls) {
	my $this_path = $origin_path.'/'.$run_name;
	if (-e $this_path.'/Demultiplex_Stats.csv' && -e $this_path.'/Top_Unknown_Barcodes.csv') {
		my @lFiles = ('Demultiplex_Stats.csv', 'Top_Unknown_Barcodes.csv');
		convert_csv_to_xls($this_path, \@lFiles);
		exit(0);
	}
}

print $cgi->header('text/json-comment-filtered');

foreach my $this_path (@$list_paths_found) {
	if ($run_name) {
		next if ($this_path ne $run_name);
	}
	if (-e $this_path.'/Demultiplex_Stats.csv') {
		mkdir $this_path.'/html/' if (not -d $this_path.'/html/');
		if (-e $this_path.'/Demultiplex_Stats.csv' && -e $this_path.'/Top_Unknown_Barcodes.csv') {
			my @lFiles = ('Demultiplex_Stats.csv', 'Top_Unknown_Barcodes.csv');
			
#			if ($this_path =~ /NB501645_0664/) {
#				warn "\n\n";
#				warn "$this_path";
#				warn Dumper @lFiles;
#				warn "\n\n";
#				#die;
#			}
			
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
			my $h_run_infos;
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
#die;

my $html_table = qq{<table data-sort-name='date' data-sort-order='desc' id="table_demultiplex" data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="20" data-page-list="[10, 20, 50, 100, 200, 300]" data-resizable='true' id='table_id_patients' class='table table-striped' style='font-size:13px;'>};
$html_table .= qq{<thead>};
$html_table .= qq{<th data-field="date" data-sortable="true" data-filter-control='input'><b>Date</b></th>};
$html_table .= qq{<th data-field="path" data-sortable="true" data-filter-control='input'><b>Path</b></th>};
$html_table .= qq{<th data-field="file" data-sortable="true" data-filter-control='input'><b>Description / File</b></th>};
$html_table .= qq{<th data-field="view"><b>View</b></th>};
$html_table .= qq{</thead>};
$html_table .= qq{<tbody>};



foreach my $date (reverse sort keys %$h_files_date) {
	my $file = $h_files_date->{$date};
	my $path_file_origin = $h_files->{$file};
	
	my $file_to_date = $origin_path.'/'.$file;
	my $found_file;
	if (-d $file_to_date) {
		$file_to_date .= '/Top_Unknown_Barcodes.csv' if -e $file_to_date.'/Top_Unknown_Barcodes.csv';
		$found_file = 1;
	}
	else {
		$file_to_date = $origin_path;
	}
	
	$path_file_origin =~ s/$origin_path/https:\/\/www.polyweb.fr\/NGS\/demultiplex/;
	my $path_file = $h_files->{$file};
	my ($substring, $color, $date_text);
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
		$file =~ s/$substring/<span style='color:black;'>$substring<\/span>/;
		$file =~ s/(.+<\/span>).+-([a-zA-Z]+\.html)$/$1.$2/;
		
		$substring =~ s/_/\./;
		$file =~ s/$substring/<span style='color:black;'>$substring<\/span>/;
		$file =~ s/(.+<\/span>).+-([a-zA-Z]+\.html)$/$1.$2/;
		
		$file =~ s/\/\//\//;
		$file =~ s/^\.//;
		$file =~ s/^-//;
		
		$file = "<i>File: ".$file."</i>";
	}
	else {
		my $description;
		if ($file =~ /^run_[0-9]+$/) {
			my $run_id = $file;
			$run_id =~ s/run_//;
			my $h_run_infos = $buffer->getQuery->getHashRunIdInfos($run_id);
			$description = $h_run_infos->{$run_id}->{description};
			my ($date_tmp, $hour_tmp) = split(' ', $h_run_infos->{$run_id}->{creation_date});
			$date_text = $date_tmp if $date_tmp;
			
		}
		$path_file = qq{DRAGEN/<span style='color:red;'>$file</span>};
		if ($description) { $file = qq{<b>Description:</b> <span style='color:blue;'>$description</span>}; }
		else { $file = qq{<i>File: $file</i>}; }
	}
	
	unless ($found_file) {
		my $file_date = $path_file_origin;
		$file_date =~ s/https:\/\/www.polyweb.fr\/NGS/\/data-isilon\/sequencing\/ngs/;
		$file_to_date = $file_date;
	}
	
	$date_text = POSIX::strftime( "20%y-%m-%d", localtime( ( stat $file_to_date )[9]) ) unless ($date_text);
	
	my $tr = qq{<tr>};
	$tr .= qq{<td>$date_text</td>};
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

sub parse_csv_file {
	my ($file, $h_run_description) = @_;
	my $h;
	$h->{has_RC} = undef;
	
	#warn $file;
	
	open (CSV, "$file");
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
			$h->{header} = \@lCol;
			$max_nb_header = scalar(@{$h->{header}});
		}
		else {
			my $j = 0;
			my $id = $lCol[0].'_'.$lCol[1].'_'.$lCol[2];
			if ($id =~ /(.+)_RC/) {
				$h->{has_RC}++ if (exists $h->{values}->{$1});
			}
			my $sample_id;
			foreach my $value (@lCol) {
				my $cat = $h->{header}[$j];
				$value = $value * 100 if ($cat =~ /\%/);
				$h->{values}->{$id}->{$cat} = $value;
				$h->{values}->{$id}->{$j} = $value;
				if ($cat eq 'SampleID') {
					$sample_id = $value;
					$nb_samples++ if (not lc($sample_id) eq 'undetermined' and not lc($sample_id) =~ /_rc/);
					
					if (not lc($sample_id) =~ /_rc/) {
						my $runs = 'no_run_info';
						if (lc($sample_id) eq 'undetermined') {
							$runs = 'undetermined';
						}
						elsif ($h_run_description and exists $h_run_description->{$sample_id}) {
							$runs = $h_run_description->{$sample_id};
						}
						$h->{runs_ids}->{by_samples}->{$value}->{$runs} = undef;
						$h->{runs_ids}->{by_runs}->{$runs}->{$value} = undef;
						$h->{values}->{$id}->{'RunID'} = $runs;
						$h->{values}->{$id}->{$max_nb_header} = $runs;
					}
				}
				if ($cat eq '# Reads') {
					$total_reads += int($value)  if (not lc($sample_id) =~ /undetermined/ and not $sample_id =~ /_RC$/);
					$h->{values}->{$id}->{'# Reads Norm Only'} = int($value);
					$h->{values}->{$id}->{$max_nb_header+1} = int($value);
					foreach my $run_id (keys %{$h->{runs_ids}->{by_samples}->{$sample_id}}) {
						$h->{runs_ids}->{by_samples}->{$sample_id}->{$run_id} += int($value);
					}
				}
				$j++;
			}
		}
		$i++;
	}
	$h->{total_reads} = $total_reads;
	$h->{nb_samples} = $nb_samples;
	close(CSV);
	return $h;
}

sub convert_csv_to_xls {
	my ($path, $list_files) = @_;
	print "Content-type: application/msexcel\n";
	print "Content-Disposition: attachment;filename=$run_name\_demultiplex.xls\n\n";
	my $workbook = Spreadsheet::WriteExcel->new( \*STDOUT );
	my $formats;
	$formats->{header} = $workbook->add_format( valign => 'vcenter', align => 'center' );
	$formats->{header}->set_color('black');
	$formats->{header}->set_bold();
	$formats->{normal} = $workbook->add_format( valign => 'vcenter', align => 'center' );
	$formats->{normal}->set_color('black');
	foreach my $file (sort @$list_files) {
		my $csv_file = $path.'/'.$file;
		my $table_id = 'table_'.$file;
		$table_id =~ s/\.csv//;
		my $h = parse_csv_file($csv_file);
		my $page_name = $file;
		$page_name =~ s/.csv//;
		my $xls_page = $workbook->add_worksheet($page_name);
		my $i = 0;
		my $j = 0;
		my $max_header = scalar(@{$h->{header}});
		foreach my $cat (@{$h->{header}}) {
			$xls_page->write($i, $j, $cat, $formats->{header});
			$j++;
		}
		$i++;
		$j = 0;
		foreach my $id (sort keys %{$h->{values}}) {
			$j = 0;
			my$continue = 1;
			while ($continue == 1) {
				my $value = $h->{values}->{$id}->{$j};
				$xls_page->write($i, $j, $value, $formats->{normal});
				$j++;
				$continue = undef if (not exists $h->{values}->{$id}->{$j});
				$continue = undef if ($j >= $max_header);
			}
			$i++;
		}
	}
	exit(0);
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
	#return $file_out if (not $run_name);
	my $h;
	
	my $h_run_description;
	my $sample_sheet_file = $path.'/SampleSheet.csv';
	if (-e $sample_sheet_file) {
		open (IN, "$sample_sheet_file");
		my $start_parsing;
		my @l_header;
		while (<IN>) {
			chomp($_);
			my $line = $_;
			if ($line =~ /,Sample_Project,/) {
				@l_header = split (',', $line);
				$start_parsing = 1;
			}
			elsif ($start_parsing == 1) {	
				my $i = 0;
				my @lTmp = split (',', $line);
				my $this_h;
				foreach my $value (@lTmp) {
					my $cat = $l_header[$i];
					$this_h->{$cat} = $value;
					$i++;
				}
				$h_run_description->{$this_h->{Sample_ID}} = $this_h->{Sample_Project};
			}
		}
		close(IN);
	}
	
	foreach my $file (sort @$list_files) {
		
		my $csv_file = $path.'/'.$file;
		my $table_id = 'table_'.$file;
		$table_id =~ s/\.csv//;
		
		
#		if ($csv_file =~ /NB501645_0664/) {
#			warn 'CSV file: '.$csv_file;
#		}
		
		my $sorted_col = "data-sort-name='n_reads' data-sort-order='desc'";
		#$sorted_col = "data-sort-name='p_of_unknown_barcodes' data-sort-order='desc'" if (lc($table_id) =~ /top_/);
		
		$h->{$file} = parse_csv_file($csv_file, $h_run_description);
		
		
#		if ($csv_file =~ /NB501645_0664/ and $csv_file =~ /Top_Unk/) {
#			warn "\n\n";
#			warn Dumper keys %{$h->{$file}->{values}};
#			warn 'CSV file: '.$csv_file;
#			warn "\n\n";
#			die;
#		}
		
		my $has_RC;
		$has_RC = $h->{$file}->{has_RC} if (exists $h->{$file}->{has_RC});
		my $total_reads = $h->{$file}->{total_reads};
		my $nb_samples = $h->{$file}->{nb_samples};
		
		my ($is_demultiplex_file, $is_top_barcodes_file);
		$is_demultiplex_file = 1 if (lc($file) =~ /demultiplex/);
		$is_top_barcodes_file = 1 if (lc($file) =~ /top_/);

		if ($is_demultiplex_file) {
			foreach my $run_id (keys %{$h->{$file}->{runs_ids}->{by_runs}}) {
				my @lSamples = keys %{$h->{$file}->{runs_ids}->{by_runs}->{$run_id}};
				my $nb_pat = 0;
				my $reads_total = 0;
				foreach my $sample (@lSamples) {
					$nb_pat++;
					$reads_total += $h->{$file}->{runs_ids}->{by_samples}->{$sample}->{$run_id};
				}
				$h->{$file}->{runs_ids}->{by_runs}->{$run_id}->{total_reads} = $reads_total;
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
		my $html;
		
		my ($table_id_resume);
		if ($is_demultiplex_file) {
			$html .= qq{<body>};
			
			my ($h_lane_resume, $h_all_patients);
			foreach my $id (sort keys %{$h->{$file}->{values}}) {
				my $lane_id = $h->{$file}->{values}->{$id}->{'Lane'};
				my $sample_id = $h->{$file}->{values}->{$id}->{'SampleID'};
				my $nb_reads = $h->{$file}->{values}->{$id}->{'# Reads'};
				my $nb_perfect_reads = $h->{$file}->{values}->{$id}->{'# Perfect Index Reads'};
				my $nb_one_mismatch_reads = $h->{$file}->{values}->{$id}->{'# One Mismatch Index Reads'};
				my $nb_two_mismatch_reads = $h->{$file}->{values}->{$id}->{'# Two Mismatch Index Reads'};
				
				if (lc($sample_id) eq 'undetermined') { $h_lane_resume->{$lane_id}->{undetermined} += $nb_reads; }
				else { $h_lane_resume->{$lane_id}->{total} += $nb_reads; }
				
				my $run_id = 'ALL_'.$h->{$file}->{values}->{$id}->{'RunID'};
				$h_all_patients->{$sample_id}->{'RunID'} = $run_id;
				$h_all_patients->{$sample_id}->{'Index'} = $h->{$file}->{values}->{$id}->{'Index'};
				$h_all_patients->{$sample_id}->{'lanes'}->{$lane_id} = undef;
				$h_all_patients->{$sample_id}->{'total_reads'} += $nb_reads;
				$h_all_patients->{$sample_id}->{'perfect_reads'} += $nb_perfect_reads;
				$h_all_patients->{$sample_id}->{'one_mismatch_reads'} += $nb_one_mismatch_reads;
				$h_all_patients->{$sample_id}->{'two_mismatch_reads'} += $nb_two_mismatch_reads;
				$h_all_patients->{$sample_id}->{'% Reads'} = '-';
				$h_all_patients->{$sample_id}->{'% Perfect Index Reads'} = '-';
				$h_all_patients->{$sample_id}->{'% One Mismatch Index Reads'} = '-';
				$h_all_patients->{$sample_id}->{'% Two Mismatch Index Reads'} = '-';
				$h_all_patients->{$sample_id}->{'Seems Ok'} = '-';
				
				$h->{$file}->{runs_ids}->{by_runs}->{$run_id}->{total_reads} += $nb_reads;
				$h->{$file}->{runs_ids}->{by_runs}->{$run_id}->{samples}->{$sample_id} = undef
			}
			foreach my $run_id (keys %{$h->{$file}->{runs_ids}->{by_runs}}) {
				next if exists $h->{$file}->{runs_ids}->{by_runs}->{$run_id}->{mean_reads};
				my $reads_total = $h->{$file}->{runs_ids}->{by_runs}->{$run_id}->{total_reads};
				my $nb_pat = scalar keys %{$h->{$file}->{runs_ids}->{by_runs}->{$run_id}->{samples}};
				$h->{$file}->{runs_ids}->{by_runs}->{$run_id}->{nb_samples} = $nb_pat;
				$h->{$file}->{runs_ids}->{by_runs}->{$run_id}->{mean_reads} = $reads_total / $nb_pat;
			}
			
			foreach my $sample_id (sort keys %{$h_all_patients}) {
				my $id = 'all_'.$sample_id;
				$h->{$file}->{values}->{$id}->{'Nb_lanes'} = scalar(keys %{$h_all_patients->{$sample_id}->{'lanes'}});
				$h->{$file}->{values}->{$id}->{'Lane'} = 'ALL';
				$h->{$file}->{values}->{$id}->{'SampleID'} = $sample_id;
				$h->{$file}->{values}->{$id}->{'Index'} = $h_all_patients->{$sample_id}->{'Index'};
				
				my $nb_reads = $h_all_patients->{$sample_id}->{'total_reads'};
				my $nb_perfect_reads = $h_all_patients->{$sample_id}->{'perfect_reads'};
				my $nb_one_mis_reads = $h_all_patients->{$sample_id}->{'one_mismatch_reads'};
				my $nb_two_mis_reads = $h_all_patients->{$sample_id}->{'two_mismatch_reads'};
				$total_reads += 0.0001;
				my $perc_reads = ($nb_reads / $total_reads) * 100;
				my ($perc_perfect_reads, $perc_one_mis_reads, $perc_two_mis_reads) = (0, 0, 0);
				$perc_perfect_reads = ($nb_perfect_reads / $nb_reads) * 100 if ($nb_perfect_reads > 0);
				$perc_one_mis_reads = ($nb_one_mis_reads / $nb_reads) * 100 if ($nb_one_mis_reads > 0);
				$perc_two_mis_reads = ($nb_two_mis_reads / $nb_reads) * 100 if ($nb_two_mis_reads > 0);
				
				$h->{$file}->{values}->{$id}->{'# Reads'} = $h_all_patients->{$sample_id}->{'total_reads'};
				$h->{$file}->{values}->{$id}->{'# Perfect Index Reads'} = $h_all_patients->{$sample_id}->{'perfect_reads'};
				$h->{$file}->{values}->{$id}->{'# One Mismatch Index Reads'} = $h_all_patients->{$sample_id}->{'one_mismatch_reads'};
				$h->{$file}->{values}->{$id}->{'# Two Mismatch Index Reads'} = $h_all_patients->{$sample_id}->{'two_mismatch_reads'};
				$h->{$file}->{values}->{$id}->{'% Reads'} = sprintf("%.2f", $perc_reads);
				$h->{$file}->{values}->{$id}->{'% Perfect Index Reads'} = sprintf("%.2f", $perc_perfect_reads);
				$h->{$file}->{values}->{$id}->{'% One Mismatch Index Reads'} = sprintf("%.2f", $perc_one_mis_reads);
				$h->{$file}->{values}->{$id}->{'% Two Mismatch Index Reads'} = sprintf("%.2f", $perc_two_mis_reads);
				$h->{$file}->{values}->{$id}->{'RunID'} = $h_all_patients->{$sample_id}->{'RunID'};
				my $c = 0;
				foreach my $cat (@{$h->{$file}->{header}}) {
					$h->{$file}->{values}->{$id}->{$c} = $h->{$file}->{values}->{$id}->{$cat};
					$c++;
				}
			}
			
			$html .= qq{<br>};
			$html .= qq{<div class="container" style="width:100%;height:200px;"><div class="row">};
			$html .= qq{<div class="col-sm-12" style="height:200px;">};
			$table_id_resume = $table_id.'_resume';
			$html .= qq{<table id="$table_id_resume" class="table table-striped sortable-table">};
			$html .= "<thead>";
			$html .= "<th data-field='lane_id' ><center><b>Lane ID</b></center></th>";
			$html .= "<th data-field='lane_total' ><center><b>Total Reads</b></center></th>";
			$html .= "<th data-field='lane_perc_total' ><center><b>% Total Reads</b></center></th>";
			$html .= "<th data-field='lane_undetermined' ><center><b>Undetermined</b></center></th>";
			$html .= "<th data-field='lane_perc_undetermined' ><center><b>% Undetermined</b></center></th>";
			$html .= "</thead>";
			$html .= "<tbody>";
			foreach my $lane_id (sort {$a <=> $b} keys %{$h_lane_resume}) {
				my $total_OK = 0;
				$total_OK = $h_lane_resume->{$lane_id}->{total} if exists $h_lane_resume->{$lane_id}->{total};
				my $undetermined = 0;
				$undetermined = $h_lane_resume->{$lane_id}->{undetermined} if exists $h_lane_resume->{$lane_id}->{undetermined};
				my $total = $total_OK + $undetermined+0.0001;
				my $total_perc_OK = sprintf("%.3f", ($total_OK / $total) * 100).'%';
				my $total_perc_undetermined = sprintf("%.3f", ($undetermined / $total) * 100).'%';
				
				$html .= "<tr><td><center>$lane_id</center></td><td><center>$total_OK</center></td><td><center>$total_perc_OK</center></td><td><center>$undetermined</center></td><td><center>$total_perc_undetermined</center></td></tr>";
			
			}
			$html .= "</tbody>";
			$html .= "</table>";
			$html .= qq{</div>};
			$html .= qq{</div>};
			$html .= qq{</div></div>};
		}
		
		
		$html .= qq{<table style="width:100%;" $sorted_col data-filter-control='true'data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="10" data-page-list="[10, 15, 20,30, 50, 100, 200, 300]" data-resizable='true' id='$table_id' class='table table-striped sortable-table' style='font-size:13px;'>};
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
			elsif ($cat eq 'Lane') {
				if ($is_demultiplex_file) { $sortable = qq{ data-sortable='true' data-filter-control='input' data-filter-default="ALL"}; }
				else { $sortable = qq{ data-sortable='true' data-filter-control='input' data-filter-default="1"}; }
			}
			else { $sortable = qq{ data-sortable='true' data-filter-control='input' }; }
			$html .= "<th $sortable data-field='".lc($cat_name)."' ><center><b>$cat</b></center></th>\n";
		}
		
		$html .= "</thead>\n";
		$html .= "<tbody>\n";
		
		my $nb_values = scalar(keys %{$h->{$file}->{values}});
		my ($mean_reads, $limit_errors_reads_min, $limit_errors_reads_max);
		
		foreach my $id (sort keys %{$h->{$file}->{values}}) {
			next if ($id =~ /.+_RC/);
			my $is_all_lanes;
			$is_all_lanes = 1 if ($h->{$file}->{values}->{$id}->{'Lane'} eq 'ALL');
			
			my $probably_run_id;
			my $sample_id = $h->{$file}->{values}->{$id}->{SampleID};
			if ($is_demultiplex_file and $is_all_lanes) {
				$probably_run_id = $h->{$file}->{values}->{$id}->{'RunID'};
				$mean_reads = $h->{$file}->{runs_ids}->{by_runs}->{$probably_run_id}->{mean_reads};
				
			}
			elsif ($is_demultiplex_file and exists $h->{$file}->{runs_ids}->{by_samples}->{$sample_id}) {
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
			if ($is_demultiplex_file and exists $h->{$file}->{values}->{$id}->{'RunID'} && $h->{$file}->{values}->{$id}->{'RunID'} eq 'no_run_info') {
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
						if ($is_all_lanes) {
							$h->{$file}->{values}->{$id}->{$j} = $h->{$file}->{values}->{$id}->{'# Reads'};
							$h->{$file}->{values}->{$id}->{'# Reads Only'} = $h->{$file}->{values}->{$id}->{'# Reads'};
							
						} 
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
						my $use_value = int($value);
						$use_value = $h->{$file}->{values}->{$id}->{'# Reads'} if ($is_all_lanes);
						if ($use_value >= $limit_errors_reads_max) {
							$with_problem++;
							$h->{$file}->{values}->{$id}->{$j+1} = 'PB +++ reads';
							$h->{$file}->{values}->{$id}->{'Seems Ok'} = 'PB +++ reads';
						}
						elsif ($use_value <= $limit_errors_reads_min) {
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
						my $use_value = int($value);
#						if ($is_all_lanes) {
#							$use_value = $use_value / $h->{$file}->{values}->{$id}->{'Nb_lanes'};
#						}
						if ($use_value >= $limit_errors_reads_max) {
							$with_problem++;
							$h->{$file}->{values}->{$id}->{'Seems Ok'} = 'PB +++ reads';
						}
						if ($use_value <= $limit_errors_reads_min) {
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
		$h->{$file}->{table_id_resume} = $table_id_resume;
	}
	open (JSON, ">$file_out");
	print JSON encode_json $h;
	close (JSON);
	return $file_out;
}

sub convert_csv_to_html {
	my ($csv_file, $html_file) = @_;
	my ($html, $html_table_resume, $html_table);
	open (CSV, "$csv_file");
	$html .= qq{<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">};
	$html .= qq{<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">};
	
	my $h_total;
	$html_table .= qq{<body><table width="100%" class="table table-striped table-hover table-condensed">};
	my $i = 0;
	while (<CSV>) {
		chomp($_);
		my $line = $_;
		my @lCol = split(',', $line);
		if ($i == 0) {
			$html_table .= "<thead>";
			foreach my $cat (@lCol) {
				$html_table .= "<th data-field='".lc($cat)."' ><center><b>$cat</b></center></th>";
			}
			$html_table .= "</thead>";
			$html_table .= "<thead>";
			$html_table .= "<tbody>";
		}
		else {
			my $lane_id = $lCol[0];
			my $sample_id = $lCol[1];
			my $nb_reads = $lCol[3];
			if (lc($sample_id) eq 'undetermined') { $h_total->{$lane_id}->{undetermined} += $nb_reads; }
			else { $h_total->{$lane_id}->{total} += $nb_reads; }
			$html_table .= "<tr>";
			foreach my $val (@lCol) {
				$html_table .= "<td><center>".$val."</center></td>";
			}
			$html_table .= "</tr>";
		}
		$i++;
	}
	$html_table .= "</tbody>";
	$html_table .= "</table>";
	close (CSV);
	
	
	$html_table_resume .= qq{<body><table width="100%" class="table table-striped table-hover table-condensed">};
	$html_table_resume .= "<thead>";
	$html_table_resume .= "<th data-field='lane_id' ><center><b>Lane ID</b></center></th>";
	$html_table_resume .= "<th data-field='lane_total' ><center><b>Total Reads</b></center></th>";
	$html_table_resume .= "<th data-field='lane_undetermined' ><center><b>Undetermined</b></center></th>";
	$html_table_resume .= "</thead>";
	$html_table_resume .= "<tbody>";$html_table .= "<tr>";
	
	foreach my $lane_id (sort {$a <=> $b} keys %{$h_total}) {
		$html_table_resume .= "<tr>";
		$html_table_resume .= "<td><center>".$lane_id."</center></td>";
		if (exists $h_total->{$lane_id}->{total}) {
			$html_table_resume .= "<td><center>".$h_total->{$lane_id}->{total}."</center></td>";
		}
		else {
			$html_table_resume .= "<td><center>0</center></td>";
		}
		if (exists $h_total->{$lane_id}->{undetermined}) {
			$html_table_resume .= "<td><center>".$h_total->{$lane_id}->{undetermined}."</center></td>";
		}
		else {
			$html_table_resume .= "<td><center>0</center></td>";
		}
		$html_table_resume .= "</tr>";
	}
	
	$html_table_resume .= "</tbody>";
	$html_table_resume .= "</table>";
	$html_table_resume .= "<br>";
	$html .= $html_table_resume;
	
	$html .= $html_table;
	
	open (HTML, ">$html_file");
	print HTML $html;
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
