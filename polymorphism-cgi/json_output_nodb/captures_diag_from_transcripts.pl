#!/usr/bin/perl
use CGI qw/:standard :html3/;
use JSON;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo/lib/obj-nodb/";
use GBuffer;
use Data::Dumper;


my $cgi    = new CGI();
my $buffer = new GBuffer;
my $tr_ids = $cgi->param('tr_ids');

my $hCaptures;
foreach my $tr_id (split(',', $tr_ids)) {
	foreach my $capture_name (@{$buffer->getListCaptureDiagFromTransId($tr_id)}) {
		$hCaptures->{$capture_name} = undef;
	}
}

my $hashRes;
$hashRes->{'label'} = 'name';
# print only genes found in annotation
my @lCaptOk = keys %{$hCaptures};
$hashRes->{'items'} = \@lCaptOk;
my $json_encode = encode_json $hashRes;

my $cgi = new CGI();
print $cgi->header('text/json-comment-filtered');
print $json_encode;
exit(0);
