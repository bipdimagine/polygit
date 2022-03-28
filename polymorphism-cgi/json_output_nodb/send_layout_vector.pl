#!/usr/bin/perl
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/layout";
use connect;
use GBuffer;
use Data::Dumper;
use JSON::XS;
use layout_vector;


my %types = ( variations => "variations" );

my $cgi    = new CGI();
my $buffer = new GBuffer;

my $project_name = $cgi->param('project');
my $type_layout = $cgi->param('type');



my $project = $buffer->newProject( -name => $project_name );
die( "unknown project" . $project_name ) unless $project;

my $layoutVariation = layout_vector::returnLayout_nodb($type_layout,$project);

print $cgi->header('text/json-comment-filtered');
my %all;
$all{identifier} = "name";
$all{label}      = "name";
$all{items}      = $layoutVariation;
	
print encode_json \%all;
exit(0);
