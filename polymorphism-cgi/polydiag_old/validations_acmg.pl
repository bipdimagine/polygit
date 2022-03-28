#!/usr/bin/perl
use CGI;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo";
use lib "$Bin/../../GenBo/lib/GenBoDB";
use lib "$Bin/../../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/validations";
#use lib "$Bin/../GenBo/lib/obj";

use Getopt::Long;
use Data::Dumper;
use Carp;
use JSON::XS;
use QueryValidationAcmg;
use Try::Tiny;
use GBuffer;

my $cgi    = new CGI;

#chargement du buffer 



eval {
my $buffer = GBuffer->new;
$buffer->dbh->{'AutoCommit'} = 0;
my $project_name = $cgi->param('project');

my $sample = $cgi->param('sample');
my $value = $cgi->param('value');
my $vid = $cgi->param('vid');
my $user = $cgi->param('user');
my $gene_id = $cgi->param('gene');
my $project = $buffer->newProjectCache(-name=>$project_name);
die() unless $project;
my $vquery = QueryValidationAcmg->new(dbh=>$buffer->dbh,database=>$project->validation_db());

my $patient = $project->getPatient($sample);

my $variation = $project->_newVariant($vid);

my $bam_lines2 ="";
my $gene = $project->newGene($gene_id);
my $id = $vquery->createVariation($variation);

die() unless $id;
$vquery->createValidation($id,$variation,$patient,$gene,$user,$value);
$buffer->dbh->commit() if $id;

sendOK("ok");
1;
}
or do {
   sendError("Insertion Problem");
   
};


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
	print $cgi->header('text/json-comment-filtered');
	print encode_json $rep;
	exit(0);
}
