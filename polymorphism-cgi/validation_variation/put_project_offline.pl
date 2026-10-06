#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);

use lib "$Bin/../../GenBo";
use lib "$Bin/../../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use GBuffer;
use GenBoProject;
use export_data;


my $cgi    = new CGI;
my $project_name = $cgi->param('project');
my $buffer  = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );

$buffer->getQuery->setProjectOffline($project->id());
my $h;
$h->{done} = 1;
export_data::print_simpleJson($cgi,[$h]);
exit(0);
