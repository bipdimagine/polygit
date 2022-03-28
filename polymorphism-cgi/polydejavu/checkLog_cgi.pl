#!/usr/bin/perl

use strict;
use CGI qw/:standard :html3/;use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/";
use lib "$Bin/../../GenBo/lib/obj-nodb/";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use Parallel::ForkManager;
use Tabix;
use JSON;
use export_data;
use QueryMooseNoDb;

my $cgi = new CGI;
my $user = $cgi->param('login');
my $pass = $cgi->param('pwd');

my @lHashRes;
printJson(@lHashRes) unless ($user and $pass);

my $buffer = GBuffer->new();
my $query = $buffer->getQuery();
my @list = @{$query->getProjectListForUser($user, $pass)};
unless (scalar(@list) > 0) { printJson(\@lHashRes); }
foreach my $hash (@{$query->getProjectsList($user, $pass)}) {
	my $hashRes->{'projName'} = $hash->{'name'};
	$hashRes->{'pat'} = $hash->{'name'};
	unless ($hash->{'description'}) { $hashRes->{'description'} = ' '; }
	else { $hashRes->{'description'} = $hash->{'description'}; }
	my @lMails;
	foreach my $hashOwners (@{$query->getOwnerProject_byName($hash->{'name'})}) { push(@lMails, $hashOwners->{'email'});  }
	$hashRes->{'owner'} = join('; ', @lMails); 
	push(@lHashRes, $hashRes);
}
printJson(\@lHashRes);


sub printJson {
	my ($listHash) = @_;
	my $hash;
	$hash->{'identifier'} = 'projName';
	$hash->{'label'} = 'projName';
	$hash->{'items'} = $listHash;
	print $cgi->header('text/json-comment-filtered');
	print encode_json $hash;
	exit(0);
}