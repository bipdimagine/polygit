#!/usr/bin/perl
$|=1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo";
use lib "$Bin/../../GenBo/lib/obj-nodb";

use Data::Dumper;
use JSON;
use GBuffer;

my $cgi = new CGI;
print $cgi->header('text/json-comment-filtered');
print "{\"progress\":\".";

my $user = $cgi->param('user');
unless ($user) {
	return_json();
	exit(0);
}

my $project_name = $cgi->param('project');
unless ($project_name) {
	return_json();
	exit(0);
}

my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name );
my $specie = $project->getVersion();
unless ($specie) {
	return_json();
	exit(0);
}

my $pathTest = $project->get_path_rna_seq_polyrna_root();
if (not $pathTest =~ /\/$specie\//) {
	my @lPotentialRelease = ('HG19', 'HG19_CNG', 'HG19_MT', 'HG38', 'HG38-ERCC', 'MM38', 'MM39');
	foreach my $rel2 (@lPotentialRelease) {
		if ($pathTest =~ /\/$rel2\//) {
			$specie = $rel2;
			last;
		}
	}
}

my $url_link;
my $found_port = read_log_file();
if ($found_port) {
	$url_link = $buffer->get_base_url_polyrna().':'.$found_port;
	return_json($url_link);
	exit(0);
}

my $file_docker_to_server = $buffer->get_polyrna_file_docker_to_server();
open (FIFO, ">>$file_docker_to_server");
print FIFO $user.','.$project_name.','.$specie."\n";
close (FIFO);

my $i = 0;
my $please_wait = 1;
while ($please_wait == 1) {
	print '.';
	$i++;
	sleep 5;
	$found_port = read_log_file();
	if ($found_port) {
		$please_wait = undef;
		$url_link = $buffer->get_base_url_polyrna().':'.$found_port;
		return_json($url_link);
		exit(0);
	}
	if ($i == 15) {
		$please_wait = undef;
		return_json();
		exit(0);
	}
}
return_json();
exit(0);




sub return_json {
	my $url = shift;
	my $h;
	if ($url) {
		$h->{url_docker} = $url;
		$h->{error} = '';
	}
	else {
		$h->{url_docker} = '';
		$h->{error} = 'error';
	}
	my $json_encode = encode_json $h;
	print ".\",";
	$json_encode =~ s/{//;
	print $json_encode;
	exit(0);
}

sub read_log_file {
	my $found_port;
	my $file_server_to_docker = $buffer->get_polyrna_file_server_to_docker();
	open (LOG, $file_server_to_docker);
	while (<LOG>) {
		chomp($_);
		my $line = $_;
		my ($log_user_name, $log_project, $log_port) = split(',', $line);
		if ($log_user_name eq $user and $log_project eq $project_name) {
			$found_port = $log_port;
		}
	}
	close(LOG);
	return $found_port;
}
