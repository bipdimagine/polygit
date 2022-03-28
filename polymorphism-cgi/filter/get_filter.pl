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
use lib "$Bin/../packages/layout";
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
#warn $user_id;

my $project = $buffer->newProject( -name => $project_name );
die( "unknown project" . $project_name ) unless $project;


unless ($user_id) {
	my $out;
	my %item;
	$item{name} = "problem no user";
	$item{id} = "1";
	$item{pp} = "zz";
	push(@$out,\%item);
	export_data::print($project,$cgi,$out);
}



my $all_filter = $query->getAllFilterName(project_id=>$project->id,user_id=>$user_id); 

my $out;
unless ($table) {

foreach my $v (sort {$a->{id} <=> $b->{id} } values %$all_filter){
	my ($date,$h) = split(" ",$v->{date});
	$v->{date} = $date;
	push(@$out,$v);
}

unless (keys %$all_filter) {
	my %item;
	$item{name} = "no filter";
	$item{id} = "zz";
	$item{pp} = "zz";
	$item{date} = "...";
	push(@$out,\%item);
}
}
export_data::print($project,$cgi,$out);

exit(0);

