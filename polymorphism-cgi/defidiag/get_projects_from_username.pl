#!/usr/bin/perl
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";


use connect;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use export_data;
use Carp;
use JSON;

my $cgi    = new CGI;
my $buffer = GBuffer->new;
my $login = $cgi->param('user');
my $pwd   = $cgi->param('pwd');
my $solo_trio = $cgi->param('solo_trio');
getProjectLists($buffer,$login,$pwd);


sub getProjectLists {
	my ( $buffer, $login, $pwd ) = @_;
	#renvoit tous les origin associés au login et au mdp spécifié
	
	my $res = $buffer->getQuery()->getProjectListForUser($login, $pwd );
	my $db_name = $buffer->{config}->{server}->{status};
	
	my $nb_project = scalar(@{$buffer->listProjects()});

	my $h_projects;
	foreach my $project (@$res) {
		next if $project->{validation_db} ne "defidiag";
		next unless ($project->{name} =~ /NGS/);
		warn Dumper $project;
		my $husers = $buffer->getQuery()->getOwnerProject($project->{id});
		my $ph = $buffer->getQuery()->getProjectByName($project->{name});
		$project->{genome_version} = $ph->{version};#join(",",map{$_->{name}}@$release);
		$project->{db_name} = $db_name; 
		if (lc($project->{dbname}) eq "polyexome" || lc($project->{dbname}) eq "polyrock"){
			$project->{dbname} = $project->{dbname}."_".$ph->{version};
		}
		$project->{users} = join("\n",map{$_->{email}} @$husers);		
		my $patients = 	$buffer->getQuery()->getPatients($project->{id});
		
		next if ($solo_trio eq 'solo' and scalar(@$patients) != 1);
		next if ($solo_trio eq 'trio' and scalar(@$patients) != 3);
		
		my %captures_id;
		foreach my $p (@$patients){
			$captures_id{$p->{capture_id}} ++;
		}
		my $find;
		my @c;
		my @t;
		my $nb =0;
		foreach my $cid (keys %captures_id){
			my $capt =  $buffer->getQuery()->getCaptureInfos($cid);
			push(@c,$capt->{name});
			push(@t,$capt->{type});
			$find=1 if $capt->{analyse} ne "exome";
			$nb += split(";",$capt->{transcripts} );
		}
		
		$h_projects->{$project->{name}}->{name}        = $project->{name};
		$h_projects->{$project->{name}}->{description} = $project->{description};
		$h_projects->{$project->{name}}->{users}       = $project->{users};
		$h_projects->{$project->{name}}->{patients}    = scalar(@$patients);
		$h_projects->{$project->{name}}->{captures}    = \@c;
	}
	my @l_proj_names = sort {$b cmp $a} (keys %$h_projects);
	print_cgi_header_html($cgi);
	my $out = qq{<table data-toggle="table" class="display table table-bordered table-hover table-responsive" style="width:100%;background-color:white;">};
	$out .= qq{<thead>};
	$out .= qq{<tr>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;font-size:13px;background-color:white;"><b>Name</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;font-size:13px;background-color:white;"><b>Description</b></th>};
#	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;font-size:13px;background-color:white;"><b>Capture(s)</b></th>};
#	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;font-size:13px;background-color:white;"><b>Pat</b></th>};
#	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;font-size:13px;background-color:white;"><b>User(s)</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;font-size:13px;background-color:white;"><b>PolyQuery<br>PolyDiag</b></th>};
	$out .= qq{</tr>};
	$out .= qq{</thead>};
	$out .= qq{<tbody style="width:100%">};
	foreach my $project_name (@l_proj_names) {
		my $url_polyquery = 'vector/gene.html?project='.$project_name.'&type=ng';
		my $url_polydiag = 'polyviewer.html?project='.$project_name.'&type=ng';
		
		#TODO: TEST ONLY !!!
		if ($project_name eq 'NGS2019_2365' or $h_projects->{$project_name}->{description} =~ /NGS2019_2365/)  {
			$url_polydiag .= '&panel=sysid';
		}
		
		$out .= '<tr>';
		$out .= "<td style='vertical-align : middle;text-align:center;font-size:13px;background-color:white;'>".$h_projects->{$project_name}->{name}.'</td>';
		$out .= "<td style='vertical-align : middle;text-align:center;font-size:13px;background-color:white;'>".$h_projects->{$project_name}->{description}.'</td>';
#		$out .= "<td style='vertical-align : middle;text-align:center;font-size:13px;background-color:white;'>";
#		my $i = 0;
#		foreach my $c (@{$h_projects->{$project_name}->{captures}}) {
#			$out .= "<br>" if ($i > 0);
#			$out .= $c;
#			$i++;
#		}
#		$out .= '</td>';
#		$out .= "<td style='vertical-align : middle;text-align:center;font-size:13px;background-color:white;'>".$h_projects->{$project_name}->{patients}.'</td>';
#		$out .= "<td style='vertical-align : middle;text-align:center;font-size:13px;background-color:white;'>".$h_projects->{$project_name}->{users}.'</td>';
		$out .= "<td style='vertical-align : middle;text-align:center;font-size:13px;background-color:white;'>";
		$out .= "<div class=\"btn-group-vertical\">";
		#$out .= "<div class=\"btn-group btn-group-sm\" role=\"group\">";
		$out .= "<button type=\"button\" class=\"btn btn-success btn-sm\" onclick=\"window.location.href='".$url_polydiag."'\">PolyDiag</button>";
		$out .= "<button type=\"button\" class=\"btn btn-info btn-sm\" onclick=\"window.location.href='".$url_polyquery."'\">PolyQuery</button>";
		$out .= "</div>";
		$out .= "</td>";
		$out .= "</tr>";
	}
	$out .= qq{</table>};
	print $out;
	exit(0);
}

sub print_cgi_header_html {
	my ($cgi) = @_;
	print $cgi -> header;
	print $cgi -> start_html (
						-title => 'PolyWeb - HGMD DB By User',
						-script => {
							-src => "//code.jquery.com/jquery-1.11.2.min.js",
						},
					);
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
