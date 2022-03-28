#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
#use lib "/software/polyweb/poly-src/GenBo/lib/obj-nodb/";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use GBuffer;
 use File::stat; 
 use LMDB_File qw(:flags :cursor_op :error);
 
 my $cgi          = new CGI();
 my $project_name = $cgi->param('project');
 my $patient_name =  $cgi->param('patient');
 my $buffer = GBuffer->new();
my $project = $buffer->newProjectCache( -name 			=> $project_name);
my $patient = $project->getPatient($patient_name);

my $fork = 3;
my $pm   = new Parallel::ForkManager($fork);
$| =1;
foreach my $chr ( @{ $project->getChromosomes } ) {
		my $pid = $pm->start and next;
		#print "!";
		my $no       = $project->getChromosome($chr->name)->lmdb_polyviewer_variants( $patient, "r" );
		$no->test(1);
		$no->lmdb($no->name);
		$pm->finish();
		#$no->lmdb($chr->name);
	}
	$pm->wait_all_children();
exit(0);