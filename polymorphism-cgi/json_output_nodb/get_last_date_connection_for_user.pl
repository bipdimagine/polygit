#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../cache_nodb/scripts/";


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
my $is_new_polybtf = $cgi->param('is_new_polybtf');


my $hashRes;
my $date = $buffer->getQuery->getLastConnectionUser($user_name);

if ($is_new_polybtf) {
	my $buffer = new GBuffer;
	my $h_polybtf_infos = $buffer->polybtf_infos();
	if ($h_polybtf_infos->{date_last_days} <= 30) { $hashRes->{$user_name} = 1; }
	elsif ($date) {
		my @lTmp = split(' ', $date);
		my ($day_user, $month_user, $year_user) = reverse split('-', $lTmp[0]);
		my ($day_polybtf, $month_polybtf, $year_polybtf) = split('/', $buffer->polybtf_infos->{date_release});
		$year_user =~ s/20//;
		$hashRes->{$user_name} = 1;
		if ($year_polybtf < $year_user) { $hashRes->{$user_name} = 0; }
		elsif ($month_polybtf < $month_user) { $hashRes->{$user_name} = 0; }
		elsif ($day_polybtf < $day_user) { $hashRes->{$user_name} = 0; }
	}
	else { $hashRes->{not} = undef; }
}

else {
	if ($date) { $hashRes->{$user_name} = $date; }
	else { $hashRes->{not} = undef; }
}

my $json_encode = encode_json $hashRes;
print $cgi->header('text/json-comment-filtered');
print $json_encode;
exit(0);

