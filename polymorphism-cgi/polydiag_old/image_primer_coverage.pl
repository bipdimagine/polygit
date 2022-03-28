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
use GenBoNoSql;

use image_primer_coverage;

my $buffer = GBuffer->new();
my $cgi          = new CGI();
print $cgi->header("image/png");
binmode STDOUT;

my $project_name = $cgi->param('project');
my $capture_name = $cgi->param('capture');
my $multiplex_name = $cgi->param('multiplex');
my $utr =$cgi->param('utr')+0;
my $intronic = $cgi->param('intronic')+0;
my $limit = $cgi->param('limit');
my $padding = $cgi->param('span');

my $project = $buffer->newProject(-name=>$project_name);
my $dir = $project->getCacheDir();
warn $dir;
my $no =  GenBoNoSql->new(dir=>$dir."/images",mode=>"w");

#my $f = $project->coverageImagesCoverageFile;
#my $kyoto_id = join("_",("all",$transcript_name,$utr,$intronic,$limit,$padding));
my $multi;
my $primers;
my $multi = $cgi->param("multiplex");
#$no->put("multipex","test","arg");
#my $ret = $no->get("multiplex",$multi);
#unless ($ret){
my $capture = $project->getCapture($capture_name);
 my $ret  = image_primer_coverage::image($project->getPatients, $capture,$multiplex_name, $limit );
#$no->put("multiplex",$multi,$ret);
#}

print $ret->{image}->png;
exit(0);
#my $tr1 = $project->newTranscript($transcript_name);
#my $image;
#eval {
#	if (-e $f){
#	my $db = $buffer->open_kyoto_db($f,"r");
#	$image =  $db->get($kyoto_id);
#	
#	}
#};
#
#my $patient_name = "all";
#my $patients = $project->getPatients();
#my $ret  = image_coverage::image ($patients, $tr1,$intronic,$utr, $padding, $limit );
#
#print $ret->{image}->png;
#
#
#exit(0);
