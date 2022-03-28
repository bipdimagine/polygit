#!/usr/bin/perl
use CGI;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/kyoto";
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
my $types = {
							he=>2,
							ho=>3,
							td=>-3,
							nc=>-5
							
						};
					
my $validations_variations = $cgi->param("variations");
my $validations_exons = $cgi->param("exons");


foreach my $validation (split(",",$validations_variations)){
	my ($vid,$st_heho) = split(":",$validation);
	die() unless exists $types->{$st_heho};

	$vquery->update_variation_validation(validation_id=>$vid,heho=>$types->{$st_heho},method=>"sanger");
	
}
foreach my $str  (split(",",$validations_exons)){
	my ($exon_id,$value) = split(/@/,$str);
	$vquery->update_exon_validation(exon_id=>$exon_id,value=>$value);
	
}


sendOK("ok");
exit(0);	


sub response {
	my ($rep) = @_;
	print $cgi->header('text/json-comment-filtered');
	print encode_json $rep;
	exit(0);
}

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
