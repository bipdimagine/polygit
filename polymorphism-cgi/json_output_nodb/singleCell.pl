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
my @list_projects = @{$buffer->getQuery->getListProjectsSC()};
my $hash;
foreach my $project_name (reverse sort @list_projects) {
	my $pname = $project_name->{name};
	my $buffer_tmp = GBuffer->new;
	my $project = $buffer_tmp->newProjectCache( -name => $pname );
	foreach my $patient (@{$project->getPatients}){
		my $name = $patient->name;
		my $id = $name;
#		$id =~ s/ATAC-//;
#		$id =~ s/ATAC_//;
		my $profile = $patient->getSampleProfile();
		
		next() if $profile =~/arc/;
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
		$hash->{$id}->{$profile}->{project} = $pname;
		$hash->{$id}->{$profile}->{description} = $project->description;
		my $path = $project->getProjectRootPath();
		
		
		if ($profile =~ /atac/){
		my $path1 = $path."/".$patient->name."/outs/cloupe.cloupe";
		my $csv = $path."/".$patient->name."/outs/summary.csv";
		my $aoa = csv (in => "$csv"); 
		my $date_cache =  utility::return_date_from_file($path1);
		my $url = qq{https://www.polyweb.fr/NGS/}.$pname."/".$patient->name."/outs/web_summary.html";
		$hash->{$id}->{$profile}->{date} = $date_cache;
		#$hash->{$id}->{$profile}->{stats}->[0] = $aoa->[1]->[2];
		$hash->{$id}->{$profile}->{stats}->[1] = $aoa->[1]->[3];
		$hash->{$id}->{$profile}->{stats}->[2] = int($aoa->[1]->[14]);
		$hash->{$id}->{$profile}->{stats}->[3] = ($aoa->[1]->[9]*100);
		$hash->{$id}->{$profile}->{summary_url} = qq{<a href="$url" target="_blank">summary</a>};
		my $url2 = qq{https://www.polyweb.fr/NGS/}.$pname."/".$patient->name."/outs/cloupe.cloupe";
		$hash->{$id}->{$profile}->{cloupe_url} = qq{<a href="$url2" target="_blank">cloupe</a>};
		}
		elsif  ($profile =~ /expression/){
		my $path1 = $path."/".$patient->name."/outs/cloupe.cloupe";
		my $csv = $path."/".$patient->name."/outs/metrics_summary.csv";
		next unless -e $csv;
		my $aoa = csv (in => "$csv"); 
		my $url = qq{https://www.polyweb.fr/NGS/}.$pname."/".$patient->name."/outs/web_summary.html";
		my $date_cache =  utility::return_date_from_file($path1);
		$hash->{$id}->{$profile}->{date} = $date_cache;
		$hash->{$id}->{$profile}->{stats}->[0] = $aoa->[1]->[0];
		$hash->{$id}->{$profile}->{stats}->[1] = $aoa->[1]->[1];
		$hash->{$id}->{$profile}->{stats}->[2] = $aoa->[1]->[2];
		
		$hash->{$id}->{$profile}->{date} = $date_cache;
		$hash->{$id}->{$profile}->{summary_url} = qq{<a href="$url" target="_blank">summary</a>};
		my $url2 = qq{https://www.polyweb.fr/NGS/}.$pname."/".$patient->name."/outs/cloupe.cloupe";
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
 $html .= "</tbody></table>";
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
	
	
	my $html = qq{<table id="table_singlecell" data-toggle="bootstrap-table" data-search="true" data-search-highlight="true" data-filter-control='true' data-toggle="table" data-pagination-v-align="bottom" data-show-extended-pagination="true" data-cache="true" data-pagination-loop="false" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-pagination-successively-size="10" data-page-size="50" data-pagination-parts="['pageInfo', 'pageSize', 'pageList']" class='table table-striped' style='font-size:15px;'>};
					
	
	#my $html = qq{<table data-search="true" data-search-highlight="true" data-filter-control='true' data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next"  id='$table_id' class='table table-striped' style='font-size:15px;'>};
	$html .= qq{<thead>};
	$html .= qq{<tr>};
	$html .= qq{<th colspan="7"><span style="color:#6B5876;"><center><b>SingleCell Expression</b></center></span></th>};
	$html .= qq{<th colspan="7"><span style="color:#00758F;"><center><b>SingleCell AtacSeq</b></center></th>};
	$html .= qq{</tr>};
	$html .= qq{<tr>};
	$html .= qq{<th data-field="event" data-sortable="true"><span style="color:#6B5876;"><center><b>ID</center></span></b></th>};
	$html .= qq{<th data-field="Expression" data-sortable="true"><span style="color:#6B5876;"><center><b>SC Expression Project</center></b></span></th>};
	$html .= qq{<th data-field="Expression_desc" data-sortable="true"><span style="color:#6B5876;"><center><b>Description</center></b></span></th>};
	$html .= qq{<th data-field="Expression_date" data-sortable="true"><span style="color:#6B5876;"><center><b>Analysis date</center></b></span></th>};
	$html .= qq{<th data-field="es0" data-sortable="true"><span style="color:#6B5876;"><center><b>Estimated Number of Cells</center></b></span></th>};
#	$html .= qq{<th data-field="es2" data-sortable="true"><span style="color:#6B5876;"><center><b>Mean Reads per Cell</center></b></span></th>};
#	$html .= qq{<th data-field="es3" data-sortable="true"><span style="color:#6B5876;"><center><b>Median Genes per Cell</center></b></span></th>};
	$html .= qq{<th data-field="exp_sum" data-sortable="true"><span style="color:#6B5876;"><center><b>10X Summary</center></b></span></th>};
	$html .= qq{<th data-field="exp_cloop" data-sortable="true"><span style="color:#6B5876;"><center><b>Cloupe</center></b></span></th>};
	$html .= qq{<th data-field="event_2" data-sortable="true"><span style="color:#00758F;"><center><b>ID</center></b></span></th>};
	$html .= qq{<th data-field="ATAC" data-sortable="true"><span style="color:#00758F;"><center><b>SC ATAC-SEQ Project</center></b></span></th>};
	$html .= qq{<th data-field="atac_desc" data-sortable="true"><span style="color:#00758F;"><center><b>Description</center></b></span></th>};
	$html .= qq{<th data-field="Date_atac" data-sortable="true"><span style="color:#00758F;"><center><b>Date</center></b></span></th>};
	$html .= qq{<th data-field="as1" data-sortable="true"><span style="color:#00758F;"><center><b>Estimated number of cells</center></b></span></th>};
#	$html .= qq{<th data-field="as2" data-sortable="true"><span style="color:#00758F;"><center><b>Median high-quality fragments per cell</center></b></span></th>};
#	$html .= qq{<th data-field="as3" data-sortable="true"><span style="color:#00758F;"><center><b>Fraction of high-quality fragments overlapping peaks</center></b></span></th>};
	$html .= qq{<th data-field="atac_sum"  data-sortable="true"><span style="color:#00758F;"><center><b>10X Summary</center></b></span></th>};
	$html .= qq{<th data-field="atac_cloop" data-sortable="true"><span style="color:#00758F;"><center><b>Cloupe</center></b></span></th>};
	$html .= qq{</tr>};
	$html .= qq{</thead>};
	$html .= qq{<tbody>};
	return $html;
}

sub print_foot_table {
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
	
	
	my $hdone;
	foreach my $id ( sort{$samples->{$a}->{sort} cmp  $samples->{$b}->{sort}} keys %$samples){
		next if exists $hdone->{$id};
		my $line =  qq{<tr>};
		if (exists $samples->{$id}->{expression}) {
			$line .= qq{<td style="background-color:#6b5876;color:white"><b>}.$id.qq{</b></td>};
			my $h = $samples->{$id}->{expression};
			delete $samples->{$id}->{expression};
			
			if (exists $h->{project}) {
				$line .= qq{<td>}.$h->{project}.qq{</td>};
			}
			else {
				$line .= qq{<td>-</td>};
			}
			if (exists $h->{description}) {
				$line .= qq{<td>}.$h->{description}.qq{</td>};
			}
			else {
				$line .= qq{<td>-</td>};
			}
			if (exists $h->{date}) {
				$line .= qq{<td>}.$h->{date}.qq{</td>};
			}
			else {
				$line .= qq{<td>-</td>};
			}
			if (exists $h->{stats}) {
				$line .= qq{<td>}.$h->{stats}->[0].qq{</td>};
			}
			else {
				$line .= qq{<td>-</td>};
			}
			if (exists $h->{summary_url}) {
				$line .= qq{<td>}.$h->{summary_url}.qq{</td>};
			}
			else {
				$line .= qq{<td>-</td>};
			}
			if (exists $h->{cloupe_url}) {
				$line .= qq{<td>}.$h->{cloupe_url}.qq{</td>};
			}
			else {
				$line .= qq{<td>-</td>};
				}
		}
		else {
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
		}
		
		my $combine_datas;
		if (exists $samples->{'ATAC_'.$id}->{atac}) {
			$id = 'ATAC_'.$id;	
			$combine_datas = 1;
		}
		
		if (exists $samples->{$id}->{atac}) {
				$line .= qq{<td style="background-color:#00758F;color:white;"><b>}.$id.qq{</b></td>};
				
				my $h = $samples->{$id}->{atac};
				delete $samples->{$id}->{atac};
				
				if (exists $h->{project}) {
					$line .= qq{<td>}.$h->{project}.qq{</td>};
				}
				else {
					$line .= qq{<td>-</td>};
				}
				if (exists $h->{description}) {
					$line .= qq{<td>}.$h->{description}.qq{</td>};
				}
				else {
					$line .= qq{<td>-</td>};
				}
				if (exists $h->{date}) {
					$line .= qq{<td>}.$h->{date}.qq{</td>};
				}
				else {
					$line .= qq{<td>-</td>};
				}
				if (exists $h->{stats}) {
					$line .= qq{<td>}.$h->{stats}->[2].qq{</td>};
				}
				else {
					$line .= qq{<td>-</td>};
				}
				if (exists $h->{summary_url}) {
					$line .= qq{<td>}.$h->{summary_url}.qq{</td>};
				}
				else {
					$line .= qq{<td>-</td>};
				}
				if (exists $h->{cloupe_url}) {
					$line .= qq{<td>}.$h->{cloupe_url}.qq{</td>};
				}
				else {
					$line .= qq{<td>-</td>};
				}
				$hdone->{$id} = 1 if $combine_datas;
		}
		else {
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
			$line .= qq{<td>-}.qq{</td>};
		}
		$html .= $line;
		
	}
	
	print_foot_table();
	
	return $html;
}


