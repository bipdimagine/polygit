#!/usr/bin/perl
$|=1;
use CGI qw/:standard :html3/;
use JSON;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/";
use lib "$Bin/../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../GenBo/lib/obj-nodb/packages/cache/polydiag/";
use Data::Dumper;
use GBuffer;




my $cgi = new CGI();
my $buffer = GBuffer->new();
my $h;
$h->{'hgmd_last_version'} = $buffer->queryHgmd->database();

print $cgi->header('text/json-comment-filtered');
print encode_json $h;
exit(0);