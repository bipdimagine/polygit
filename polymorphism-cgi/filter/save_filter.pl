#!/usr/bin/perl
use CGI;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/GenBoDB/writeDB";
use GBuffer;
use Data::Dumper;

my $cgi = new CGI;
print $cgi->header();

my @params = $cgi->param();
my $buffer = GBuffer->new();


my $filter_name = $cgi->param("name");
my $project_name = $cgi->param("project");
my $user_name = $cgi->param("user");
my $query = $buffer->queryFilter();
my $user_id =$query->get_user_id(user_name=>$user_name);
die("no user") unless $user_id;

my $project = $buffer->newProject( -name => $project_name );
die( "unknown project" . $project_name ) unless $project;

my $buffer2 = GBuffer->new();
my $filter_id = $query->newFilter(filter_name=>$filter_name,project_id=>$project->id,user_id=>$user_id);
die("no filter") unless $filter_id;

foreach my $param (@params){
	next if $param !~ /filter/;
	my $value = $cgi->param($param);
	next if $value eq "";
	next if $value eq "all";
	warn $param." ".$cgi->param($param);
	$query->addParam(filter_id=>$filter_id,param_name=>$param,param_value=>$value);
}
sendOK("coucou");
exit(0);

sub sendOK {
	my ($text) = @_;
	my $resp;
	$resp->{status}  = "OK";
	$resp->{message} = $text;
	response($resp);
}

sub sendError {
	my ($text) = @_;

	my $resp;
	$resp->{status}  = "error";
	$resp->{message} = $text;
	response($resp);
}

sub response {
	my ($rep) = @_;
	print qq{<textarea>};
	print to_json($rep);
	print qq{</textarea>};
	exit(0);
}
