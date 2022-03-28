#!/usr/bin/perl
use CGI qw/:standard :html3/;
use Data::Dumper;

use strict;
use CGI::Session;
use JSON;

my $cgi = new CGI();

my $h;
my $h1;
my $name = `hostname`;
chomp($name);
my $res2 = `printenv`;
chomp($res2);
my $host;
my @lEnv = split ("\n", $res2);
foreach my $line (@lEnv) {
	my @lTmp = split ('=', $line);
	if ($lTmp[0] eq 'SERVER_NAME'){
		my $ip_complete = $lTmp[-1];
		my @lIpInt = split ('\.', $ip_complete);
		$host = 'SERVER '.$lIpInt[-1];
	}
}

if ($name =~ /\d+/){
	$h1->{'name'} = $host;
}
else {
	$h1->{'name'} = $name.' SERVER';
}

my @lHash = ($h1);

$h->{'label'} = 'name';
$h->{'items'} = \@lHash;


print $cgi->header('text/json-comment-filtered');
print encode_json $h;