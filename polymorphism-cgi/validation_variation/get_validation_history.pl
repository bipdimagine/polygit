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
my $project_name = $cgi->param('project');
my $id = $cgi->param('vid');

my $my_phenotype = $cgi->param('my_phenotype');
my $no_css = $cgi->param('no_css');
my $buffer = GBuffer->new();
my $project = $buffer->newProject(-name => $project_name);
my $vquery = $project->validations_query(1);

#my $vid = $vquery->getIdFromVariation($id);

my $h1 = $vquery->getAllValidationsForVariation($id);

#die();
#my $values;
#my $h;
#my @lInfos = split(',', $values);
my $h;
my ($tab) = values %$h1;
my $out;
$out .= qq{<center><div style="width:98%;height:600px;">};
$out .= qq{<table class="table table-bordered" style="font-size:1.5em;">};
$out .= qq{<thead>};
$out .= qq{<tr>};
$out .= qq{<th scope="col"><center><b>Project Name</b></center></th>};
$out .= qq{<th scope="col"><center><b>Sample</b></center></th>};
$out .= qq{<th scope="col"><center><b>status</b></center></th>};
$out .= qq{<th scope="col"><center><b>Validation</b></center></th>};
$out .= qq{<th scope="col"><center><b>User</b></center></th>};
$out .= qq{<th scope="col"><center><b>Date</b></center></th>};
$out .= qq{</tr>};
$out .= qq{</thead>};
$out .= qq{<tbody>};
foreach my $info (@$tab) {
	
	my $project_name = $info->{project_name};
	my $patient_name  =$info->{sample_name};
	
	#my ($project_name, $patient_name) = split(':', $info);
	my $buffer = GBuffer->new();
	my $project = $buffer->newProject(-name => $project_name);
	my $patient;
	eval { $patient = $project->getPatient($patient_name); };
	if ($@) {
		
		next;
	}
	
	$out.= $cgi->start_Tr();
	$out.= $cgi->td($info->{project_name});
	$out.= $cgi->td($info->{sample_name});
	$out.= $cgi->td($patient->return_icon());
	$out.= $cgi->td($info->{term});
	$out.= $cgi->td($info->{user_name});
	$out.= $cgi->td($info->{modification_date});
	
	$out.= $cgi->td($info->{phenotypes_name});
	$out.= $cgi->end_Tr();
	
	$project = undef;
	$buffer = undef;
}


$out .= qq{</tbody>};
$out .= qq{</table>};
$out .= qq{</div></center>};

if ($no_css) { print $cgi->header(); }
else { html::print_cgi_header($cgi); }
print $out;
exit(0);
my $hashRes;
$hashRes->{html_table} = $out;
print $cgi->header('text/json-comment-filtered');
my $json_encode = encode_json $hashRes;
print $json_encode;
exit(0);

