#!/usr/bin/perl

use Carp;
use strict;
use JSON;
use Data::Dumper;
use CGI qw/:standard :html3/;
use Set::IntSpan::Fast::XS;
use Set::IntervalTree;
use List::Util qw[min max];
use FindBin qw($Bin);
use Storable qw(store retrieve freeze);
 
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";

use GBuffer;
use GenBoProject;

my $cgi = new CGI;
my $projects = $cgi->param('projects');


my $html;
my @lProjects = split(',', $projects);

foreach my $projectname (@lProjects) {
	my $buffer = GBuffer->new();	
	my $project = $buffer->newProjectCache( -name => $projectname);
	
	my $hurl;
	foreach my $fam (@{$project->getFamilies()}) {
		my $fam_name = $fam->name();
		foreach my $patient (@{$fam->getPatients()}) {
			my $pat_name = $patient->name();
			if ($patient->sex() == 1) {
				if ($patient->isChild) {
					$hurl->{$fam_name}->{$pat_name}->{sex} = "<img src='/icons/Polyicons/baby-boy.png'>";
				}
				else {
					$hurl->{$fam_name}->{$pat_name}->{sex} = "<img src='/icons/Polyicons/male.png'>";
				}
			}
			elsif ($patient->sex() == 2) {
				if ($patient->isChild) {
					$hurl->{$fam_name}->{$pat_name}->{sex} = "<img src='/icons/Polyicons/baby-girl.png'>";
				}
				else {
					$hurl->{$fam_name}->{$pat_name}->{sex} = "<img src='/icons/Polyicons/female.png'>";
				}
			}
			if ($patient->status() == 1) {
				$hurl->{$fam_name}->{$pat_name}->{status} = "<img src='/icons/Polyicons/bullet_green.png'>";
			}
			elsif ($patient->status() == 2) {
				$hurl->{$fam_name}->{$pat_name}->{status} = "<img src='/icons/Polyicons/pill2.png'>";
			}
			$hurl->{$fam_name}->{$pat_name}->{url_single} = "http://www.polyweb.fr/polyweb/html/manta/Url_allSVEditor_project.html?projectname=$projectname&filename=$pat_name&trio=0";
			if ($patient->isChild()) {
				next unless ($fam->father());
				next unless ($fam->mother());
				$hurl->{$fam_name}->{$pat_name}->{url_trio} = "http://www.polyweb.fr/polyweb/html/manta/Url_allSVEditor_project.html?projectname=$projectname&filename=$pat_name&trio=1";
			}
		}
	}
	
	$html .= "<h3>Project $projectname</h3>";
	$html .= "<br><br>";
	$html .= "<table class='table'>";
	$html .= "<thead>";
	$html .= "<tr>";
	$html .= "<th style='padding-left:10px;' align='left' scope='col'>Family Name</th>";
	$html .= "<th style='padding-left:10px;' align='left' scope='col'>Patient Name</th>";
	$html .= "<th style='padding-left:10px;' align='left' scope='col'>Sex</th>";
	$html .= "<th style='padding-left:10px;' align='left' scope='col'>Status</th>";
	$html .= "<th style='padding-left:10px;' align='left' scope='col'>Link SINGLE</th>";
	$html .= "<th style='padding-left:10px;' align='left' scope='col'>Link TRIO</th>";
	$html .= "</tr>";
	$html .= "</thead>";
	$html .= "<tbody>";
	foreach my $fam_name (sort keys %$hurl) {
		foreach my $pat_name (sort keys %{$hurl->{$fam_name}}) {
			$html .= "<tr>";
			$html .= "<td style='padding-left:10px;' align='left'>".$fam_name."</td>";
			$html .= "<td style='padding-left:10px;' align='left'>".$pat_name."</td>";
			$html .= "<td style='padding-left:10px;' align='left'>".$hurl->{$fam_name}->{$pat_name}->{sex}."</td>";
			$html .= "<td style='padding-left:10px;' align='left'>".$hurl->{$fam_name}->{$pat_name}->{status}."</td>";
			$html .= "<td style='padding-left:10px;' align='left'><a href=\"".$hurl->{$fam_name}->{$pat_name}->{url_single}."\" target=\"_blank\">LINK_SINGLE</a></td>";
			if (exists $hurl->{$fam_name}->{$pat_name}->{url_trio}) {
				$html .= "<td style='padding-left:10px;' align='left'><a href=\"".$hurl->{$fam_name}->{$pat_name}->{url_trio}."\" target=\"_blank\">LINK_TRIO</a></td>";
			}
			else {
				$html .= "<td style='padding-left:10px;' align='left'></td>";
			} 
			$html .= "</tr>";
		}
	}
	$html .= "</table>";
	$html .= "<br><br>";
}
	


print $cgi->header('text/json-comment-filtered');
my $hashRes;
$hashRes->{html} = $html;
my $json_encode = encode_json $hashRes;
print $json_encode;
exit(0);
