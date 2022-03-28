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
my $buffer = GBuffer->new;
my $user_name = $cgi->param('user_name');

my $hash = $buffer->get_alamut_api_key_from_user_name($user_name);
my $json_encode = encode_json $hash;
print $cgi->header('text/json-comment-filtered');
print $json_encode;
exit(0);

