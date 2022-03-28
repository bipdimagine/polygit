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
use lib "$Bin/../packages/validation_variation"; 
use html; 

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
use JSON;
use validationQuery;
use Date::Tiny;
use lib "$Bin/../packages/validation_variation"; 
use draw_cnv; 
use infos_coverage_exons;
use preload_coverage;
$|=1;
                
my $cgi          = new CGI();
my $project_name = $cgi->param('project');


my $buffer = GBuffer->new();
my $project = $buffer->newProject(-name=>$project_name);

my $hRes;
$hRes->{isDude} = 'yes';

my $cache_dir = $project->transcriptsDudeDir();
foreach my $patient (@{$project->getPatients()}) {
	my $file = $cache_dir.'/'.$patient->name().'.dude.transcripts';
	unless (-e $file) {
		$hRes->{isDude} = 'no';
		last;
	}
}

print $cgi->header('text/json-comment-filtered');
print encode_json $hRes;
exit(0);