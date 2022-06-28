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

my ($h_files, $h_files_date);
my $origin_path = '/data-isilon/sequencing/ngs/demultiplex/';
my ($list_paths_found) = list_html_files_in_dir($origin_path);

foreach my $this_path (@$list_paths_found) {
	if (-e $this_path.'/Demultiplex_Stats.csv') {
		mkdir $this_path.'/html/' if (not -d $this_path.'/html/');
		if (-e $this_path.'/Demultiplex_Stats.csv' && -e $this_path.'/Top_Unknown_Barcodes.csv' ) {
			my @lFiles = ('Demultiplex_Stats.csv', 'Top_Unknown_Barcodes.csv');
			my $json_file = convert_csv_to_json($this_path, \@lFiles);
			my $json_path = $json_file;
			$json_file =~ s/\/\//\//g;
			$json_file =~ s/.+\/ngs\//https:\/\/www.polyweb.fr\/NGS\//;
			my $name = $json_file;
			if ($name =~ /.+\/NGS\/demultiplex\/(.+)\/json.+/) {
				$name = $1;
			}
			next if ($name eq 'test_xths');
			next if ($name =~ /run1[0-9]+/);
			add_file_json('https://www.polyweb.fr/mbras/polyweb/demultiplex_view/demultiplex_view_run.html?name='.$name.'&json='.$json_file, $json_path, $name);
		}
		elsif (-e $this_path.'/Demultiplex_Stats.csv') {
			convert_csv_to_html($this_path.'/Demultiplex_Stats.csv', $this_path.'/html/laneBarcode.html') if (not -e $this_path.'/html/laneBarcode.html');
			add_file_html($this_path.'/html/laneBarcode.html', 'laneBarcode.html');
		}
		elsif (-e $this_path.'/Top_Unknown_Barcodes.csv') {
			convert_csv_to_html($this_path.'/Top_Unknown_Barcodes.csv', $this_path.'/html/top_unknown_barcodes.html') if (not -e $this_path.'/html/top_unknown_barcodes.html');
			add_file_html($this_path.'/html/top_unknown_barcodes.html', 'top_unknown_barcodes.html');
		}
		
		
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



sub convert_csv_to_json {
	my ($path, $list_files) = @_;
	my $path_out = $path.'/json/';
	mkdir $path_out unless (-d $path_out);
	my $file_out = "$path_out/demultiplex.json";
	return $file_out if (-e $file_out);
	my $h;
	foreach my $file (@$list_files) {
		my $csv_file = $path.'/'.$file;
		my $table_id = 'table_'.$file;
		$table_id =~ s/\.csv//;
		my $html = qq{<table style="width:100%;" data-filter-control='true'data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="15" data-page-list="[10, 15, 20,30, 50, 100, 200, 300]" data-resizable='true' id='$table_id' class='table table-striped' style='font-size:13px;'>};
		open (CSV, "$csv_file");
		my $i = 0;
		while (<CSV>) {
			chomp($_);
			my $line = $_;
			my @lCol = split(',', $line);
			if ($i == 0) {
				$h->{$file}->{header} = \@lCol;
				$html .= "<thead>\n";
				foreach my $cat (@lCol) {
					my $cat_name = $cat;
					$cat_name =~ s/# /n_/;
					$cat_name =~ s/% /p_/;
					my $sortable = qq{ data-sortable='true' data-filter-control='input'} if (not $cat =~ /[#%]/);
					$html .= "<th $sortable data-field='".lc($cat_name)."' ><center><b>$cat</b></center></th>\n";
				}
				$html .= "</thead>\n";
				$html .= "<tbody>\n";
			}
			else {
				my $j = 0;
				my $id = $lCol[0].'_'.$lCol[1];
#				my $tr_id = 'tr_'.$file.'_'.$id;
#				$html .= "<tr id=\"$tr_id\">\n";
				foreach my $value (@lCol) {
					my $cat = $h->{$file}->{header}[$j];
					$h->{$file}->{values}->{$id}->{$cat} = $value;
					$h->{$file}->{values}->{$id}->{$j} = $value;
#					my $td_id = 'td_'.$file.'_'.$id.'_'.$j;
#					$html .= "<td id=\"$td_id\"><center>".$value."</center></td>\n";
					$j++;
				}
#				$html .= "</tr>\n";
			}
			$i++;
		}
		close(CSV);
		
		foreach my $id (sort keys %{$h->{$file}->{values}}) {
			next if ($id =~ /.+_RC/);
			my $tr_id = 'tr_'.$file.'_'.$id;
			my $html_tr;
			my $j = 0;
			my $with_problem = 0;
			while ($j >= 0) {
				my $cat = $h->{$file}->{header}[$j];
				my $td_id = 'td_'.$file.'_'.$id.'_'.$j;
				my $value = $h->{$file}->{values}->{$id}->{$j};
				if ($cat =~ /[#%]/ && exists $h->{$file}->{values}->{$id.'_RC'}) {
					my $value_RC = $h->{$file}->{values}->{$id.'_RC'}->{$j};
					$html_tr .= "<td id=\"$td_id\"><center>";
					$html_tr .= "<table><tr><td><b>Norm</b></td><td style='padding-left:5px;'><center>$value</center></td></tr><tr><td><b>RC</b></td><td style='padding-left:5px;'><center>$value_RC</center></td></tr></table>";
					$html_tr .= "</center></td>\n";
					if ($cat eq '# Reads') {
						$with_problem++ if ($value_RC >= $value);
					}
				}
				else {
					$html_tr .= "<td id=\"$td_id\"><center>".$value."</center></td>\n";
				}
				$j++;
				$j = -1 if (not exists $h->{$file}->{values}->{$id}->{$j});
			}
			$html_tr .= "</tr>\n";
			if ($with_problem > 0) { $html_tr = "<tr style='background-color:red;color:white;' id=\"$tr_id\">\n".$html_tr; }
			elsif (lc($id) =~ /undetermined/) { $html_tr = "<tr style='color:red;' id=\"$tr_id\">\n".$html_tr; }
			else { $html_tr = "<tr id=\"$tr_id\">\n".$html_tr; }
			$html .= $html_tr;
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
	my $date_stat = (stat ($path))[9];
	if (exists $h_files_date->{$date_stat}) {
		my $is_ok = 0;
		while ($is_ok == 0) {
			$date_stat++;
			next if (exists $h_files_date->{$date_stat});
			$is_ok = 1;
		}
	}
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