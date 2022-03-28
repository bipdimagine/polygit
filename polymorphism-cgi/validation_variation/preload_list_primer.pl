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
use validationQuery;
use Date::Tiny;
use lib "$Bin/../packages/validation_variation"; 
use draw_cnv; 
use infos_coverage_exons;
use preload_coverage;

use JSON::XS;
$|=1;
my $buffer = GBuffer->new();




my $cgi          = new CGI();
my $project_name = $cgi->param('project');
my $project = $buffer->newProject(-name=>$project_name);
$project->print_waiting(1);
my $no = $project->noSqlCoverage();
 print $cgi->header('text/json-comment-filtered');
  print "{\"bundle\":\"";
foreach my $chr (@{$project->getChromosomes}){
	if ($no->exists("primers",$chr->name)){
		print"@ @";
	}
	
	my $primers = $chr->getPrimers();
	#warn scalar(@$primers);
	print "@";
	map {delete $_->{project}; delete $_->{buffer}} @$primers;
	eval{
	$no->put("primers",$chr->name,$primers);
	}
	
}
foreach my $capture (@{$project->getCaptures}){
	unless ($no->exists($project->name(),$capture->id)){
		my $list = $capture->getListPrimers();
		$no->put($project->name(),$capture->id,$list);
	}
	#
}
print"\"}";
#response({status=>"ok"});

exit(0);

sub response {
	my ($rep) = @_;
	print $cgi->header('text/json-comment-filtered');
	print encode_json $rep;
	exit(0);
}