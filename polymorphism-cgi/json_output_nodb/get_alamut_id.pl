#!/usr/bin/perl
$|=1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";


use connect;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use JSON;

my $cgi    = new CGI;
my $project_name = $cgi->param('project');
my $var_id = $cgi->param('var');

#$var_id = '14_58427532_AAAAT_A';

my $buffer = GBuffer->new;
my $project = $buffer->newProject( -name => "$project_name" );

my $hash;
$hash->{init} = $var_id;
my @ltmp = split('_', $var_id);
if ($ltmp[2] eq '-') {
	my $chr = $project->getChromosome($ltmp[0]);
	my $start = $ltmp[1] - 1;
	my $ref_all = $chr->sequence($start, $start);
	$var_id = $ltmp[0].'_'.$start.'_'.$ref_all.'_'.$ref_all.$ltmp[3];
}
elsif ($ltmp[2] eq $ltmp[3]) {
	my $chr = $project->getChromosome($ltmp[0]);
	my $start = $ltmp[1] - 1;
	my $ref_all = $chr->sequence($start, $start);
	$var_id = $ltmp[0].'_'.$start.'_'.$ref_all.$ltmp[2].'_'.$ref_all;
}

my $v = $project->_newVariant($var_id);
my $alamut_id = $v->alamut_id();


$hash->{alamut} = $alamut_id;
$hash->{id} = $var_id;
my $json_encode = encode_json $hash;
print $cgi->header('text/json-comment-filtered');
print $json_encode;
exit(0);

