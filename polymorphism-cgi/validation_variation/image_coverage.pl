#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../GenBo/lib/kyoto";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";
use Carp;
use export_data;
use strict;
use Set::IntSpan::Fast::XS;
use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Carp;
use Set::Intersection;
use Tabix;
use Storable qw/store thaw retrieve/;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use coverage;
use Spreadsheet::WriteExcel;
use POSIX;
use validationQuery;
use Date::Tiny;
use lib "$Bin/../packages/validation_variation"; 
use draw_cnv; 
use infos_coverage_exons;
use GD;
use image_coverage;

my $buffer = GBuffer->new();
my $cgi          = new CGI();
print $cgi->header("image/png");
binmode STDOUT;

my $project_name = $cgi->param('project');
my $transcript_name = $cgi->param('transcript');

my $utr =$cgi->param('utr')+0;
my $intronic = $cgi->param('intronic')+0;
my $limit = $cgi->param('limit');

my $padding = $cgi->param('span');
my $patient_names =  $cgi->param('patients');

my $project = $buffer->newProject(-name=>$project_name);
$limit= 15 if $project->validation_db() eq 'glucogen' ;
#my $f = $project->coverageImagesCoverageFile;
#my $kyoto_id = join("_",("all",$transcript_name,$utr,$intronic,$limit,$padding));
my $tr1 = $project->newTranscript($transcript_name);
##my $image;
#eval {
#	if (-e $f){
#	my $db = $buffer->open_kyoto_db($f,"r");
#	$image =  $db->get($kyoto_id);
#	
#	}
#};
#
#if ($image && $patient_names eq 'all'){
#	print $image;
#	exit(0);
#}
my $patients =  $project->get_list_patients($patient_names,",");
my $ret;
if ($patients->[0]->isNoSqlDepth) {
	$ret  = image_coverage::image_depth_lmdb ($patients, $tr1,$intronic,$utr, $padding, $limit );	
}
else {
	$ret  = image_coverage::image ($patients, $tr1,$intronic,$utr, $padding, $limit );
}
print $ret->{image}->png;


exit(0);
