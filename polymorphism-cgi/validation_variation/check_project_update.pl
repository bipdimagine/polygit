#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/kyoto";
use lib "$Bin/../packages/export";
use Carp;
use export_data;
use strict;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use JSON::XS;
use validationQuery;


my $buffer = GBuffer->new();
my $cgi    = new CGI();
my $project_name = $cgi->param('project');
my $project = $buffer->newProjectCache(-name=>$project_name);


my @out;
my $item;
if ($project->isUpdate()) {
	$item->{is_updated} = 'yes';
}
else {
	$item->{is_updated} = 'no';
}
push(@out,$item);
export_data::print($project,$cgi,\@out);
exit(0);