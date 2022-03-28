#!/usr/bin/perl
$|=1;
use CGI qw/:standard :html3/;
use JSON;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/";
use lib "$Bin/../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../GenBo/lib/obj-nodb/packages/cache/polydiag/";
use Data::Dumper;
use GBuffer;
use QueryVcf;
use Getopt::Long;
use update;
use QueryVectorFilter;

my $cgi = new CGI();
my $user		= $cgi->param('user');
my $pwd			= $cgi->param('pwd');

print_cgi_header_html($cgi);
my $out;
#my $out = qq{<table data-toggle="table" class="display table table-bordered table-hover table-responsive" style="width:100%;height:600px;">};
#$out .= qq{<thead>};
#$out .= qq{<tr>};
#$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;vertical-align : middle;text-align:center;"><b>Project Name</b></th>};
#$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;vertical-align : middle;text-align:center;"><b>Type</b></th>};
#$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;vertical-align : middle;text-align:center;"><b>Description</b></th>};
#$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;vertical-align : middle;text-align:center;"><b>Nb Patients</b></th>};
#$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;vertical-align : middle;text-align:center;"><b>Link</b></th>};
#$out .= qq{</tr>};
#$out .= qq{</thead>};
#$out .= qq{<tbody style="width:100%">};


my $buffer = GBuffer->new();
my $query = $buffer->getQuery();
my $hDone;
foreach my $list_infos (@{$query->listAllProjects($user, $pwd)}) {
	my $project_name = $list_infos->[0];
	my $project_type = $list_infos->[1];
	next unless ($project_name =~ /NGS20/);
	$hDone->{$project_name} = $project_type;
}

my $nb = 0;
my $nb_for_progress = scalar(keys %$hDone) / 100;
print qq{<div hidden>};

foreach my $project_name (reverse sort(keys %{$hDone})) {
	my $project_type = $hDone->{$project_name};
	my $buffer_2 = GBuffer->new();
	my $project = $buffer->newProject( -name => $project_name);
	my $description = $project->description();
	my $nb_patients = scalar @{$query->getPatients($project->id)};
	
	#$out .= print_grid_projects($project_name, $project_type, $description, $nb_patients);
	$out .= print_panels_projects($project_name, $project_type, $description, $nb_patients);
	
	$project = undef;
	$buffer_2 = undef;
	$nb++;
	if ($nb == int($nb_for_progress)) {
		print '@';
		$nb = 0;
	}
}


print qq{</div>};
$out .= '</tbody>';
print $out;



########################



sub print_panels_projects {
	my ($project_name, $project_type, $description, $nb_patients) = @_;
	
	my $collapse_id = 'div_collapse_'.$project_name;
	my $div_container_hgmd = 'div_container_hgmd_'.$project_name;
	my $out = $cgi->start_div( { class => "panel panel-warning", style => "background-color:#ffe7a0;overflow:auto;border: 1px solid black;color-adjust: exact;padding:2px;" } );
	$out .=  $cgi->start_div( { class => "btn-group btn-xs", style => "float:left;min-width:180px;" } );
	$out .=qq{<div class="btn  btn-info btn-xs " style="background-color:#fcb900;position:relative;bottom:1px;min-width:120px;" data-toggle="collapse" data-target="#$collapse_id" onclick="launch('$div_container_hgmd', '$project_name')"> $project_name &nbsp </div>};
	$out .=qq{<div class="btn  btn-info btn-xs " style="background-color:#fcb900;position:relative;bottom:1px;min-width:80px;" data-toggle="collapse" data-target="#$collapse_id" onclick="launch('$div_container_hgmd', '$project_name')"> $project_type &nbsp </div>};
	$out .= $cgi->end_div();
	$out .=  $cgi->start_div( { class => "btn-group btn-xs", style => "float:left;min-width:180px;padding-left:5px;" } );
	$out .=qq{<div class="btn  btn-info btn-xs " style="background-color:#fcb900;position:relative;bottom:1px;min-width:150px;" data-toggle="collapse" data-target="#$collapse_id" onclick="launch('$div_container_hgmd', '$project_name')"> $description &nbsp </div>};
	$out .= $cgi->end_div();
	$out .=  $cgi->start_div( { class => "btn-group btn-xs", style => "float:right;" } );
	$out .=qq{<div class="btn  btn-info btn-xs " style="background-color:#fcb900;position:relative;bottom:1px;min-width:70px;" data-toggle="collapse" data-target="#$collapse_id" onclick="launch('$div_container_hgmd', '$project_name')"> $nb_patients Patient(s)&nbsp </div>};
	$out .= $cgi->end_div();
	$out .= "</div>";
	
	$out .= qq{<div id="$collapse_id" class="collapse" style="width:100%;max-height:640px;">};
	$out .= qq{<div class="container" id="$div_container_hgmd" style="width:100%;max-height:640px;overflow-y:auto;overflow-x:visible;padding-top:5px;"></div>};
	$out .= "<br><br>";
	$out .= "</div>";
	
	
	return $out;
}

sub print_grid_projects {
	my ($project_name, $project_type, $description, $nb_patients) = @_;
	my $url = 'polyhgmd_by_user.html?project='.$project_name;
	my $out = '<tr>';
	$out .= "<td style='vertical-align : middle;text-align:center;'><b>".$project_name.'</b></td>';
	$out .= "<td style='vertical-align : middle;text-align:center;'>".lc($project_type).'</td>';
	$out .= "<td style='vertical-align : middle;text-align:center;'><i>".$description.'</i></td>';
	$out .= "<td style='vertical-align : middle;text-align:center;'>".$nb_patients.'</td>';
	$out .= "<td style='vertical-align : middle;text-align:center;'><a href='$url' target='_blank'>View</a></td>";
	$out .= '</tr>';
	return $out;
}

sub print_cgi_header_html {
	my ($cgi) = @_;
	print $cgi -> header;
    print qq{
    	<script type="text/javascript" src="//maxcdn.bootstrapcdn.com/bootstrap/3.3.2/js/bootstrap.min.js""></script>
    };
	print qq{		
		<style>
			table-mybordered > thead > tr > th, .table-mybordered > thead > tr > th, table-mybordered > tbody > tr > th, .table-mybordered > tbody > tr > th, table-mybordered > tfoot > tr > th, .table-mybordered > tfoot > tr > th, table-mybordered > thead > tr > td, .table-mybordered > thead > tr > td, table-mybordered > tbody > tr > td, .table-mybordered > tbody > tr > td, table-mybordered > tfoot > tr > td, .table-mybordered > tfoot > tr > td {
				border: 1px solid #95A5A6;
		  		vertical-align:middle;
		   		text-align:center; 
		   		padding:2px;
			}
			.bs-callout {
			    padding: 20px;
			    margin: 20px 0;
			    border: 1px solid #eee;
			    border-left-width: 5px;
			    border-radius: 3px;
			}
			.bs-callout h4 {
			    margin-top: 0;
			    margin-bottom: 5px;
			}
			.bs-callout p:last-child {
			    margin-bottom: 0;
			}
			.bs-callout code {
			    border-radius: 3px;
			}
			.bs-callout+.bs-callout {
			    margin-top: -5px;
			}
			.bs-callout-default {
			    border-left-color: #777;
			}
			.badge1 {
			  padding: 1px 9px 2px;
			  font-size: 12.025px;
			  font-weight: bold;
			  white-space: nowrap;
			  color: #ffffff;
			  background-color: #999999;
			  -webkit-border-radius: 9px;
			  -moz-border-radius: 9px;
			  border-radius: 9px;
			}
			
			.hide-bullets {
			    list-style:none;
			    margin-left: -40px;
			    margin-top:20px;
			}
			
			.thumbnail {
			    padding: 0;
			}
			
			.carousel-inner>.item>img, .carousel-inner>.item>a>img {
			    width: 100%;
			}
		</style>	
	}	
}