#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/kyoto";

#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
#use Set::;
use Storable qw/store thaw retrieve/;
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
use Storable qw/thaw/;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use coverage;
use Spreadsheet::WriteExcel;
use POSIX;
use validationQuery;
use Date::Tiny;
use lib "$Bin/../packages/validation_variation"; 
use draw_cnv; 
use infos_coverage_exons;
use polyweb_coverage;
use html;


my $buffer = GBuffer->new();

my $cgi          = new CGI();
#print $cgi -> header;
#print $CSS_TABLE."\n";
	html::print_header_polydiag($cgi);
my $project_name = $cgi->param('project');

my $project = $buffer->newProject(-name=>$project_name);

my $patient_name = $cgi->param('patients');
$patient_name ="all" unless $patient_name;
my $patients = $project->get_list_patients($patient_name,",");


my $cgi_transcript =  $cgi->param('transcripts');
my $transcript  = $project->newTranscript($cgi_transcript);

warn $transcript;
warn $transcript->id;

my $padding = $cgi->param('span');
my $limit = $cgi->param('limit');

my $utr  = $cgi->param('utr');


                                                                                     #(Any :$patients!,        Any :$transcript!,             Int :$utr!,  Int :$limit!, Int :$padding) {
	my $coverage = polyweb_coverage->new(patients=>$patients,transcript=>$transcript,utr=>$utr,padding=>$padding,limit=>$limit,intronic=>0);

print $coverage->html_table();
#print image_coverage::table_depth_lmdb_from_matrix( $transcript,$patients,$matrix,$limit);



