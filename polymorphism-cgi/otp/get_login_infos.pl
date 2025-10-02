#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";

use GBuffer;
use Getopt::Long;
use Data::Dumper;
use Carp;
use JSON;

my $cgi    = new CGI;
my $login = $cgi->param('login');
my $show_global_infos = $cgi->param('show_global_infos');
my $generate_otp_key = $cgi->param('generate_otp_key');


my $buffer = GBuffer->new;
my $dbh = $buffer->getQuery();

$login = lc($login);

if ($show_global_infos) {
	my $h_user_infos = $buffer->getQuery->getUserInfos($login);
	if (defined($h_user_infos->{$login}->{'uKey'}) and $h_user_infos->{$login}->{'uKey'} == 1) {
		$h_user_infos->{$login}->{'otp'} = "OTP Auth: YES";
		$h_user_infos->{$login}->{'otp'} .= "<br>key:".$h_user_infos->{$login}->{'Key'};
	}
	else {
		my $cmd = qq{generate_otp_key('$login')};
		$h_user_infos->{$login}->{'otp'} = "<center>";
		$h_user_infos->{$login}->{'otp'} .= "OTP Auth: NO";
		if (defined($h_user_infos->{$login}->{'Key'}) and $h_user_infos->{$login}->{'Key'} =~ /[a-z0-9]+/) {
			$h_user_infos->{$login}->{'otp'} .= "<br>key:".$h_user_infos->{$login}->{'Key'};
		}
		else {
			$h_user_infos->{$login}->{'otp'} .= "<br><button id='b_generate_key_otp' onClick=\"$cmd\">Generate KEY</button>";
		}
		$h_user_infos->{$login}->{'otp'} .= "</center>";
	}
	
	my $h_alamut = $buffer->get_alamut_api_key_from_user_name($login);
	if ($h_alamut) {
		$h_user_infos->{$login}->{'alamut'} = "Licence: ".$h_alamut->{licence_alamut};
		$h_user_infos->{$login}->{'alamut'} .= "<br>Institution Key: ".$h_alamut->{institution_key};
		$h_user_infos->{$login}->{'alamut'} .= "<br>API Key: ".$h_alamut->{api_key};
	}
	else {
		$h_user_infos->{$login}->{'alamut'} = 'No Licence';
	}
	
	my $json_encode = encode_json $h_user_infos->{$login};
	print $cgi->header('text/json-comment-filtered');
	print $json_encode;
	exit(0);
}

if ($generate_otp_key) {
	my $hres;
	my $new_key = $buffer->google_auth->generate_secret32();
	$hres->{'otp'} = "<center>";
	$hres->{'otp'} .= "OTP Auth: NO";
	$hres->{'otp'} .= "<br>";
	$hres->{'otp'} .= "<br>KEY: ".$new_key.'<br><b><i>now in DB!</b></i>';
	my $dbh = $buffer->getQuery->getDbh();
	my $sql = qq{
		UPDATE `bipd_users`.`USER` SET `Key`='$new_key' WHERE `LOGIN`=?;
	};
	my $sth = $dbh->prepare($sql);
	$sth->execute($login);
	$hres->{'otp'} .= "</center>";
	my $json_encode = encode_json $hres;
	print $cgi->header('text/json-comment-filtered');
	print $json_encode;
	exit(0);
}


exit(0);
 