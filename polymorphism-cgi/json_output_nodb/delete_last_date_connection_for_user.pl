#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../cache_nodb/scripts/";

require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";

use connect;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use export_data;
use Carp;
use JSON;

my $cgi    = new CGI;
my $buffer = GBuffer->new;
my $user_name = $cgi->param('user_name');

$buffer->getQuery->delteLastConnectionUser($user_name);
my $hashRes;
$hashRes->{ok} = 'ok';
my $json_encode = encode_json $hashRes;
print $cgi->header('text/json-comment-filtered');
print $json_encode;
exit(0);

