#!/usr/bin/perl
$|=1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo";
use lib "$Bin/../../GenBo/lib/obj-nodb";
use Date::Tiny;
use Text::CSV qw( csv );
use connect;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use JSON;
use Carp;
require "$Bin/../../GenBo/lib/obj-nodb/packages/cache/polydiag/utility.pm";

my $cgi    = new CGI;
my $buffer = GBuffer->new;
my $view_project = $cgi->param('project');

my $html;
my @list_tables_ids;
my @list_projects = @{$buffer->getQuery->getListProjectsByType("rhu4")};
my $hash;
foreach my $project_name (reverse sort @list_projects) {
	my $buffer_tmp = GBuffer->new;
	my $project = $buffer_tmp->newProjectCache( -name => $project_name );
	foreach my $patient (@{$project->getPatients}){
		my $name = $patient->name;
		my $id = $name;
		$id =~ s/ATAC-//;
		$id =~ s/ATAC_//;
		my $profile = $patient->getSampleProfile();
		 if ($profile =~/exp/){
		 		$profile = "expression"
		 }
		 elsif ($profile =~/atac/){
		 	$profile = "atac";
		 }
		 else {
		 	next;
		 }
	
		$hash->{$id}->{sort} .= $project->name; 
		$hash->{$id}->{$profile}->{name} = $name;
		$hash->{$id}->{$profile}->{pid} = $patient->id;
		$hash->{$id}->{$profile}->{project} = $project->name;
		$hash->{$id}->{$profile}->{description} = $project->description;
		my $path = $project->getProjectRootPath();
		
		
		if ($profile =~ /atac/){
		my $path1 = $path."/".$patient->name."/outs/cloupe.cloupe";
		my $csv = $path."/".$patient->name."/outs/summary.csv";
		my $aoa = csv (in => "$csv"); 
		my $date_cache =  utility::return_date_from_file($path1);
		my $url = qq{https://www.polyweb.fr/NGS/}.$project->name."/".$patient->name."/outs/web_summary.html";
		$hash->{$id}->{$profile}->{date} = $date_cache;
		#$hash->{$id}->{$profile}->{stats}->[0] = $aoa->[1]->[2];
		$hash->{$id}->{$profile}->{stats}->[1] = $aoa->[1]->[3];
		$hash->{$id}->{$profile}->{stats}->[2] = int($aoa->[1]->[14]);
		$hash->{$id}->{$profile}->{stats}->[3] = ($aoa->[1]->[9]*100);
		$hash->{$id}->{$profile}->{summary_url} = qq{<a href="$url" target="_blank">summary</a>};
		my $url2 = qq{https://www.polyweb.fr/NGS/}.$project->name."/".$patient->name."/outs/cloupe.cloupe";
		$hash->{$id}->{$profile}->{cloupe_url} = qq{<a href="$url2" target="_blank">cloupe</a>};
		}
		elsif  ($profile =~ /expression/){
		my $path1 = $path."/".$patient->name."/outs/cloupe.cloupe";
		my $csv = $path."/".$patient->name."/outs/metrics_summary.csv";
		next unless -e $csv;
		my $aoa = csv (in => "$csv"); 
		my $url = qq{https://www.polyweb.fr/NGS/}.$project->name."/".$patient->name."/outs/web_summary.html";
		my $date_cache =  utility::return_date_from_file($path1);
		$hash->{$id}->{$profile}->{date} = $date_cache;
		$hash->{$id}->{$profile}->{stats}->[0] = $aoa->[1]->[0];
		$hash->{$id}->{$profile}->{stats}->[1] = $aoa->[1]->[1];
		$hash->{$id}->{$profile}->{stats}->[2] = $aoa->[1]->[2];
		
		$hash->{$id}->{$profile}->{date} = $date_cache;
		$hash->{$id}->{$profile}->{summary_url} = qq{<a href="$url" target="_blank">summary</a>};
		my $url2 = qq{https://www.polyweb.fr/NGS/}.$project->name."/".$patient->name."/outs/cloupe.cloupe";
		$hash->{$id}->{$profile}->{cloupe_url} = qq{<a href="$url2" target="_blank">cloupe</a>};
		}
		elsif  ($profile =~ /genome/){
		$hash->{$id}->{$profile}->{date} ="-";
		}
		
		
	}
}
 #print_table($hash);
 #exit(0);
 $html = print_rhu4($hash);
 #print_table($hash);
 #exit(0);
my $hash1;
$hash1->{html} = $html;
$hash1->{tables_ids} = join(';', @list_tables_ids);
my $json_encode = encode_json $hash1;
print $cgi->header('text/json-comment-filtered');
print $json_encode;
exit(0);

sub print_header {
	my $table_id = 'table_'.time;
	push(@list_tables_ids, $table_id);
	my $html = qq{<table data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next"  id='$table_id' class='table table-striped' style='font-size:15px;'>};
	$html .= qq{<thead>};
	$html .="<tr>";
	$html .= qq{<td ></td>};
	$html .= qq{<td style="background-color: #6B5876;font-size:20px;text-align:center;color:white" colspan="6">SingleCell Expression</td>};
	$html .= qq{<td style="background-color: #00758F;font-size:20px;text-align:center;color:white" colspan="6">SingleCell AtacSeq</td>};
	$html .="</tr>";
	
	$html .= qq{<th data-field="event" data-filter-control="search"><b>id</b></th>};
	$html .= qq{<th data-field="Expression"><b>SC Expression Project</b></th>};
	$html .= qq{<th data-field="Expression_desc"><b>Description</b></th>};
	$html .= qq{<th data-field="Expression_date"><b>Analysis date</b></th>};
	$html .= qq{<th data-field="es0" data-sortable="true"><b>Estimated Number of Cells</b></th>};
#	$html .= qq{<th data-field="es2" data-sortable="true"><b>Mean Reads per Cell</b></th>};
#	$html .= qq{<th data-field="es3" data-sortable="true"><b>Median Genes per Cell</b></th>};
	$html .= qq{<th data-field="exp_sum" ><b>10X Summary </b></th>};
	$html .= qq{<th data-field="exp_cloop"><b>cloupe</b></th>};
	$html .= qq{<th data-field="ATAC"><b>SC ATAC-SEQ Project</b></th>};
	$html .= qq{<th data-field="atac_desc"><b>Description</b></th>};
	$html .= qq{<th data-field="Date_atac" data-sortable="true"><b>Date</b></th>};
	$html .= qq{<th data-field="as1" data-sortable="true"><b>Estimated number of cells</b></th>};
#	$html .= qq{<th data-field="as2" data-sortable="true"><b>Median high-quality fragments per cell</b></th>};
#	$html .= qq{<th data-field="as3" data-sortable="true"><b>Fraction of high-quality fragments overlapping peaks</b></th>};
	$html .= qq{<th data-field="atac_sum" "><b>10X Summary </b></th>};
	$html .= qq{<th data-field="atac_cloop"><b>cloupe</b></th>};
	$html .= qq{</thead>};
	$html .= qq{<tbody>};
	return $html;
}
sub print_table {
	my ( $samples) = @_;
	foreach my $id ( sort{$samples->{$a}->{sort} cmp  $samples->{$b}->{sort}} keys %$samples){
			$samples->{$id}->{atac}->{pid} ="-" unless $samples->{$id}->{atac}->{pid};
			$samples->{$id}->{expression}->{pid} ="-" unless $samples->{$id}->{expression}->{pid};
		print $id."\t".$samples->{$id}->{expression}->{pid}."\t".$samples->{$id}->{atac}->{pid}."\n";
	}
}
sub print_rhu4 {
	my ( $samples) = @_;
	
	my $html = print_header();
	
	foreach my $id ( sort{$samples->{$a}->{sort} cmp  $samples->{$b}->{sort}} keys %$samples){
		my $line =  qq{<tr>};
		$line .= qq{<td>}.$id.qq{</td>};
		if (exists $samples->{$id}->{expression}) {
				my $h = $samples->{$id}->{expression};
				$line .= qq{<td>}.$h->{project}.qq{</td>};
				$line .= qq{<td>}.$h->{description}.qq{</td>};
				$line .= qq{<td>}.$h->{date}.qq{</td>};
				$line .= qq{<td>}.$h->{stats}->[0].qq{</td>};
			#	$line .= qq{<td>}.$h->{stats}->[2].qq{</td>};
			#	$line .= qq{<td>}.$h->{stats}->[3].qq{</td>};
				$line .= qq{<td>}.$h->{summary_url}.qq{</td>};
				$line .= qq{<td>}.$h->{cloupe_url}.qq{</td>};
		}
		else {
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
		}
		if (exists $samples->{$id}->{atac}) {
				my $h = $samples->{$id}->{atac};
				$line .= qq{<td>}.$h->{project}.qq{</td>};
				$line .= qq{<td>}.$h->{description}.qq{</td>};
				$line .= qq{<td>}.$h->{date}.qq{</td>};
				$line .= qq{<td>}.$h->{stats}->[2].qq{</td>};
			#	$line .= qq{<td>}.$h->{stats}->[1].qq{</td>};
			#	$line .= qq{<td>}.$h->{stats}->[2].qq{</td>};
				$line .= qq{<td>}.$h->{summary_url}.qq{</td>};
				$line .= qq{<td>}.$h->{cloupe_url}.qq{</td>};
		}
		else {
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
		}
		$html .= $line;
		
	}
	
	
	return $html;
}


