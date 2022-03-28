#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
#use lib "/software/polyweb/poly-src/GenBo/lib/obj-nodb/";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag";

use GBuffer;
#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../packages/cache"; 

use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use draw_cnv; 
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/utility.pm";
use html; 
use Carp;
use strict;
use Data::Dumper;
use GenBo;
use Compress::Snappy;
 use File::Temp qw/ tempfile tempdir /;
#require "$Bin/../GenBo/lib/obj-nodb/GBuffer.pm";
use List::MoreUtils qw{ natatime };
use Getopt::Long;
use Carp;
use Set::Intersection;
use Tabix;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Spreadsheet::WriteExcel;
use POSIX;
use validationQuery;
use QueryValidationAcmg;
use Date::Tiny;
use JSON::XS;
use List::MoreUtils qw{part};
#use PDF::API2;
#use PDF::Table;
use constant mm => 25.4 / 72;
use constant in => 1 / 72;
use constant pt => 1;
use Time::HiRes qw ( time alarm sleep );
use Storable qw/thaw freeze/;
 use Digest::MD5 qw(md5 md5_hex md5_base64);
 use File::stat; 
 use LMDB_File qw(:flags :cursor_op :error);
 
$| =1;


my $cgi          = new CGI();
my $values = $cgi->param('values');
my $my_phenotype = $cgi->param('my_phenotype');

my $h;
my @lInfos = split(',', $values);
foreach my $info (@lInfos) {
	my ($project_name, $patient_name, $overlap) = split(':', $info);
	my $buffer = GBuffer->new();
	my $project = $buffer->newProject(-name => $project_name);
	my $patient;
	eval { $patient = $project->getPatient($patient_name); };
	if ($@) {
		$h->{$project->name()}->{$patient_name}->{phenotypes} = '';
		$h->{$project->name()}->{$patient_name}->{family} = '';
		$h->{$project->name()}->{$patient_name}->{sex} = '';
		$h->{$project->name()}->{$patient_name}->{status} = '';
		$h->{$project->name()}->{$patient_name}->{run_id} = '';
		$h->{$project->name()}->{$patient_name}->{machine} = '';
		$h->{$project->name()}->{$patient_name}->{is_child} = '';
		$h->{$project->name()}->{$patient_name}->{capture} = '';
		$h->{$project->name()}->{$patient_name}->{overlap} = '';
		$h->{$project->name()}->{$patient_name}->{contacts} = 'Patient '.$patient_name.' no more in this project. Please update dejavu...';
		$project = undef;
		$buffer = undef;
		next;
	}
	my @lPheno;
	foreach my $pheno_name (@{$project->phenotypes}) {
		if ($my_phenotype eq $pheno_name) { push(@lPheno, "<b><span style='color:red;'>$pheno_name</b></span>"); }
		else { push(@lPheno, $pheno_name); }
	}
	$h->{$project->name()}->{$patient->name()}->{phenotypes} = join(', ', sort @lPheno);
	$h->{$project->name()}->{$patient->name()}->{family} = $patient->getFamily->name();
	$h->{$project->name()}->{$patient->name()}->{sex} = $patient->sex();
	$h->{$project->name()}->{$patient->name()}->{status} = $patient->status();
	$h->{$project->name()}->{$patient->name()}->{run_id} = $patient->getRunId();
	$h->{$project->name()}->{$patient->name()}->{machine} = $patient->getRun->machine();
	$h->{$project->name()}->{$patient->name()}->{is_child} = $patient->isChild();
	my @lCapture;
	foreach my $capture (@{$patient->getCaptures()}) { push(@lCapture, $capture->name()); }
	$h->{$project->name()}->{$patient->name()}->{capture} = join(', ', @lCapture);
	$h->{$project->name()}->{$patient->name()}->{overlap} = $overlap;
	$h->{$project->name()}->{$patient->name()}->{contacts} = join(' - ', @{$project->get_list_emails()});
	
	$project = undef;
	$buffer = undef;
}

my $out;
$out .= qq{<center><div style="width:98%;height:600px;">};
$out .= qq{<table class="table table-bordered">};
$out .= qq{<thead>};
$out .= qq{<tr>};
#$out .= qq{<th scope="col"><center><b>View</b></center></th>};
$out .= qq{<th scope="col"><center><b>Project Name</b></center></th>};
$out .= qq{<th scope="col"><center><b>Phenotype</b></center></th>};
$out .= qq{<th scope="col"><center><b>Family Name</b></center></th>};
$out .= qq{<th scope="col"><center><b>Patient Name</b></center></th>};
$out .= qq{<th scope="col"><center><b>Sex</b></center></th>};
$out .= qq{<th scope="col"><center><b>Status</b></center></th>};
$out .= qq{<th scope="col"><center><b>Run Id</b></center></th>};
$out .= qq{<th scope="col"><center><b>Machine</b></center></th>};
$out .= qq{<th scope="col"><center><b>Capture</b></center></th>};
$out .= qq{<th scope="col"><center><b>Identity</b></center></th>};
$out .= qq{<th scope="col"><center><b>Contacts</b></center></th>};
$out .= qq{</tr>};
$out .= qq{</thead>};
$out .= qq{<tbody>};

my $nb=0;
foreach my $project_name (sort keys %{$h}) {
	foreach my $patient_name (sort keys %{$h->{$project_name}}) {	
		my $phenotypes = $h->{$project_name}->{$patient_name}->{phenotypes};
		my $family_name = $h->{$project_name}->{$patient_name}->{family};
		my $is_child = $h->{$project_name}->{$patient_name}->{is_child};
		my $sex = $h->{$project_name}->{$patient_name}->{sex};
		$nb++;
		if ($is_child eq '1') {
			if ($sex eq '1') { 
				$sex = "<img src='/icons/Polyicons/baby-boy.png'>"; 
			}
			elsif ($sex eq '2') { $sex = "<img src='/icons/Polyicons/baby-girl.png'>"; 
				
			}
		}
		else {
			if ($sex eq '1') { $sex = "<img src='/icons/Polyicons/male.png'>"; }
			elsif ($sex eq '2') { $sex = "<img src='/icons/Polyicons/female.png'>"; }
		}
		my $status = $h->{$project_name}->{$patient_name}->{status};
		if ($status eq '1') { 
			$status = "<img src='/icons/Polyicons/bullet_green.png'>"; 
		}
		elsif ($status eq '2') { 
			$status = "<img src='/icons/Polyicons/pill2.png'>"; 
		}
		my $run_id = $h->{$project_name}->{$patient_name}->{run_id};
		my $machine = $h->{$project_name}->{$patient_name}->{machine};
		my $capture = $h->{$project_name}->{$patient_name}->{capture};
		my $overlap = $h->{$project_name}->{$patient_name}->{overlap};
		my $contacts = $h->{$project_name}->{$patient_name}->{contacts};
		
		$out .= qq{<tr>};
		#$out .= qq{<td><center><div><button  id="btn.$nb" onclick="LaunchPolycyto('$project_name','$patient_name');"> <span style='font-size:15px;'>&#128270;</span></button></div></center></td>};
		$out .= qq{<td><center>$project_name</center></td>};
		$out .= qq{<td><center>$phenotypes</center></td>};
		$out .= qq{<td><center>$family_name</center></td>};
		$out .= qq{<td><center>$patient_name</td>};
		$out .= qq{<td><center>$sex</center></td>};
		$out .= qq{<td><center>$status</center></td>};
		$out .= qq{<td><center>$run_id</center></td>};
		$out .= qq{<td><center>$machine</center></td>};
		$out .= qq{<td><center>$capture</center></td>};
		$out .= qq{<td style="margin:0px;width:5%;">};
		my @tabres = split(/ /,$overlap);
		foreach my $res (@tabres)
		{
			my ($p,$c,$d,$f) = split(/_/,$res);   # pourcentage caller debut fin
			my $color="purple";
			if ($c=~ m/w/){$color="green";}
			elsif($c=~ m/c/){$color="orange";}
			
			$out .= qq{<div style="margin:0px;"><center><font style="font-size: 11px;color:$color">$c : $p</font></div>};
		}
		$out .= qq{</td>};		
		#$out .= qq{<td><center><font style="font-size: 10px;">$overlap</font></center></td>};
		$out .= qq{<td><center>$contacts</center></td>};
		$out .= qq{</tr>};
	}
}
$out .= qq{</tbody>};
$out .= qq{</table>};
$out .= qq{</div></center>};


my $hashRes;
$hashRes->{html_table} = $out;
print $cgi->header('text/json-comment-filtered');
my $json_encode = encode_json $hashRes;
print $json_encode;
exit(0);

