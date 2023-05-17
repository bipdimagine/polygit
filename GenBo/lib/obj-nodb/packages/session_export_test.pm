package session_export_test;

use strict;
use Moo;

use Data::Dumper;
use Config::Std;
use POSIX qw(strftime);
extends "session_export";




has tmp_dir => (
	is		=> 'rw',
	lazy 	=> 1,
	default => sub {
		my $self = shift;
		my $buffer = new GBuffer;
		my $tmp_dir = $buffer->config->{project_data}->{archives}.'/session_test/';
		unless ( -d $tmp_dir ) {
			`mkdir $tmp_dir`;
			`chmod 777 $tmp_dir`;
		}
		return $tmp_dir;
	}
);

has session_expiration => (
	is		=> 'rw',
	lazy 	=> 1,
	default => "1y",
);

sub get_config_file_name {
	my ($self, $name) = @_;
	my $file_cfg = $self->tmp_dir.'/'.$name.'.cfg';
	return $file_cfg;
}

sub get_session_id_from_config {
	my ($self, $name) = @_;
	my $file_cfg = $self->get_config_file_name($name);
	confess("\n\nERROR: no session PROD for $name. Die.\n\n") unless (-e $file_cfg);
	read_config $file_cfg => my %cfg;
	my $cgi_old_session_id = $cfg{$name}->{session_id};
	return $cgi_old_session_id;
}

sub update_config_session {
	my ($self, $name) = @_;
	my $file_cfg = $self->get_config_file_name($name);
	if (-e $file_cfg) {
		my $cgi_old_file = $self->tmp_dir.'/cgisess_'.$self->get_session_id_from_config($name);
		`rm $cgi_old_file` if (-e $cgi_old_file);
		`rm $file_cfg`;
	}
	$self->write_config_session($name);
}

sub write_config_session {
	my ($self, $name) = @_;
	my $file_cfg = $self->get_config_file_name($name);
	confess("\n\nERROR: $file_cfg already exists. Die.\n\n") if (-e $file_cfg);
	my $date = strftime "%d/%m/%Y", localtime;
	open(FILE, ">$file_cfg");
	print FILE "[$name]\n";
	print FILE "session_id:".$self->{session_id}."\n";
	print FILE "date:".$date."\n";
	close(FILE);
}

1;