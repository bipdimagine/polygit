#!/usr/bin/perl
$|=1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo";
use lib "$Bin/../../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";


use connect;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use JSON;
use xls_export;
use session_export;
use Parallel::ForkManager;

my $cgi    = new CGI;
my $project_name = $cgi->param('project');
my $fork = $cgi->param('fork');
$fork = 1 unless ($fork);

my $buffer = GBuffer->new;
my $project = $buffer->newProject( -name => $project_name );
my @l_cmds;
foreach my $patient (@{$project->getPatients}) {
	my $patient_name = $patient->name();
	print "\n\n### PATIENT $patient_name ###\n\n";
	my $cmd = "$Bin/rna_view_analyse.pl project=$project_name patient=$patient_name update_sashimi=1 fork=$fork"; 
	`$cmd`;
}





