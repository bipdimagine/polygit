#!/usr/bin/perl
use CGI;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-lite";
use lib "$Bin/../GenBo/lib/GenBoDB/writeDB";
use lib "$Bin/../packages/export";
use GBuffer;
use Data::Dumper;
use GenBoFilter;
use export_data;
my $cgi = new CGI;

my $buffer = GBuffer->new();
my $project_name = $cgi->param('project');
my $project = $buffer->newProject( -name => $project_name );
die( "unknown project" . $project_name ) unless $project;
my $filter_name= $cgi->param('filter_name');
my $user_name= $cgi->param('user');
die();
if ($filter_name) {
	my $user_id = GenBoFilter::get_user_id($buffer->dbh,$user_name);	
	die("no user") unless $user_id;
	my $filter_id =GenBoFilter::getid($buffer->dbh,$filter_name,$project->id,$user_id);	
	my $params = GenBoFilter::getParam($buffer->dbh,$filter_id);
	
	my $regions_text  =  $params->{filter_region}->{value}; 
	warn $regions_text;
	my @data;
	my $zz = 1;
	foreach my $r (split(" ",$regions_text)){
		my @rt = split(":",$r);
		my %region;
		$region{id} = $region{name} = $zz++;
		$region{chromosome} = $rt[0];
		$region{start} = $rt[1];
		$region{end} = $rt[2];
		$region{include} = 0;
		push(@data,\%region);
	}
	export_data::print($project,$cgi,\@data);
}
exit(0);