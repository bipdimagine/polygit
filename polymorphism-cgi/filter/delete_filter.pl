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
use GenBoFilter;
  use JSON::XS;
my $cgi = new CGI;
print $cgi->header();

my @params = $cgi->param();
my $buffer = GBuffer->new();
my $query = $buffer->queryFilter();

my $filter_ids = $cgi->param("filter_id");
my $project_name = $cgi->param("project");
my $user_name = $cgi->param("user");

my $user_id = $query->get_user_id(user_name=>$user_name);
die("no user") unless $user_id;

my $project = $buffer->newProject( -name => $project_name );
die( "unknown project" . $project_name ) unless $project;

foreach my $filter_id (split(";",$filter_ids)){
	#my $filter_id = $query->newFilter(filter_name=>$filter_name,project_id=>$project->id,user_id=>$user_id);
	#next unless $filter_id;
	$query->delete_filter(filter_id=>$filter_id,user_id=>$user_id);
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
	print encode_json($rep);
	print qq{</textarea>};
	exit(0);
}
