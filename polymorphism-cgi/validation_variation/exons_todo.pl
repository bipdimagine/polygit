#!/usr/bin/perl
use CGI;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo";
use lib "$Bin/../../GenBo/lib/GenBoDB";

use lib "$Bin/../../GenBo/lib/obj-nodb";
#use lib "$Bin/../GenBo/lib/obj";
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use Carp;
use JSON::XS;
use GenBoFilter;
use validationQuery;

my $cgi    = new CGI;

#chargement du buffer 
#!/usr/bin/perl




my $buffer = GBuffer->new;
my $project_name = $cgi->param('project');
my $project = $buffer->newProject(-name=>$project_name);
die() unless $project;
my $capture = $project->getCapture();
my $vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>$capture->validation_db());

my $user_name = $cgi->param("user_name");

my @ids_exons = split("!",$cgi->param('exons'));



foreach my $id_exon (@ids_exons){
	warn  $id_exon;
	my $hres;
	$hres->{id} = $id_exon;
	$hres->{project_name} = $project_name;
	$hres->{user_name} = $user_name;
	$hres->{todo} = 1;
	 ($hres->{sample_name},$hres->{gene},$hres->{transcript},$hres->{name},$hres->{chromosome},$hres->{start},$hres->{end}) = split(":",$id_exon);
	 $vquery->exon_todo(%$hres);
} 

sendOK("ok");
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
