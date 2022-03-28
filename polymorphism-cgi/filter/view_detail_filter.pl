#!/usr/bin/perl
use CGI;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/GenBoDB/writeDB";
use lib "$Bin/../packages/export";
use GBuffer;
use Data::Dumper;
use GenBoFilter;
use export_data;

my $cgi = new CGI;
 

my @params = $cgi->param();
my $buffer = GBuffer->new();

my $project_name = $cgi->param("project");
my $user_name = $cgi->param("user");
my $table = $cgi->param("table");
my $query = $buffer->queryFilter();
my $user_id =$query->get_user_id(user_name=>$user_name);

my $project = $buffer->newProject( -name => $project_name );
die( "unknown project" . $project_name ) unless $project;

my $filter_name  = $cgi->param("filter_name");
unless ($user_id) {
	my $out;
	my %item;
	$item{name} = "problem no user";
	$item{id} = "1";
	$item{pp} = "zz";
	push(@$out,\%item);
	export_data::print($project,$cgi,$out);
}




if ($filter_name) {

die("no user") unless $user_id;
my $filter_id =$query->get_filter_id(filter_name=>$filter_name,project_id=>$project->id, user_id=> $user_id);	
my $params = $query->getParam(filter_id=>$filter_id);
my $out;
my $nb = 1;
foreach my $name (keys %$params) {
	#next if ($name =~ /region/);
	$params->{$name}->{value} =~ s/\+/ /g;
	my %d;
	$d{name} = $name;
	$d{value} = $params->{$name}->{value};
	$d{id} = $nb++;
	push(@$out,\%d);
	
}
export_data::print($project,$cgi,$out);
}






exit(0);

