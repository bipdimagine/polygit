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
	if (-d $this_path.'/html/') {
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
#foreach my $file (sort keys %$h_files) {
	my $file = $h_files_date->{$date};
	my $path_file_origin = $h_files->{$file};
	$path_file_origin =~ s/$origin_path/https:\/\/www.polyweb.fr\/NGS\/demultiplex/;
	my $path_file = $h_files->{$file};
	$path_file =~ s/$origin_path/DEMULTIPLEX/;
	$path_file =~ s/$file//;
	$path_file =~ s/\/\///g;
	my @l_path = split('/', $path_file);
	$path_file = $l_path[0].'/'.$l_path[1];
	$path_file =~ s/(run[0-9]+)/<span style='color:red;'>$1<\/span>/;
	
	$file =~ s/(run[0-9]+)/<span style='color:red;'>$1<\/span>/;
	$file =~ s/^\.//;
	
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