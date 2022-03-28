package session_export;
use strict;
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../";
use lib "$Bin/../../";
use List::Util qw( max min sum);
use Carp qw(confess croak);
use Data::Dumper;
use Moose;
use MooseX::Method::Signatures;
use Spreadsheet::WriteExcel;
use CGI::Session;
use Compress::Snappy;
use Storable qw(store retrieve freeze dclone thaw);

has tmp_dir => (
	is		=> 'rw',
	lazy 	=> 1,
	default => sub {
		my $self = shift;
		my $buffer = new GBuffer;
		my $tmp_dir = $buffer->config->{project_data}->{global_search};
		#my $tmp_dir = '/data-isilon/bipd-src/mbras/tmp_session/';
		unless ( -d $tmp_dir ) {
			`mkdir $tmp_dir`;
			`chmod 777 $tmp_dir`;
		}
		return $tmp_dir;
	}
);

has session_id => (
	is		=> 'rw',
	lazy 	=> 1,
	default => undef,
);

has session => (
	is		=> 'rw',
	lazy 	=> 1,
	default => sub {
		my $self = shift;
		my $cgi = new CGI();
		my $session = new CGI::Session( undef, $cgi, { Directory => $self->tmp_dir() } );
		$session->expire("1d"); 
		$self->{session_id} = $session->id();
		return $session;
	}
);

sub id() {
	my ($self) = @_;
	return $self->session_id();
}

sub load_session {
	my ($self, $my_sid) = @_;
	$self->{session_id} = $my_sid;
	$self->{session} = new CGI::Session( undef, $my_sid, { Directory => $self->tmp_dir() } );
}

sub save {
	my ($self, $key, $value) = @_;
	$self->session->param( $key, $value );
}

sub save_compress {
	my ($self, $key, $value) = @_;
	my $value2 = compress( freeze $value );
	$self->save($key, $value2);
}

sub load {
	my ($self, $key) = @_;
	confess ("\n\nNo session found... die.\n\n") unless ($self->session());
	my $value;
	eval{ $value = $self->session->param($key); };
	return $value;
}

sub load_compress {
	my ($self, $key) = @_;
	return thaw( decompress $self->load($key) );
}

sub check_if_expired {
	my $self = shift;
	if ( $self->session->is_expired ) {
    	$self->session->delete();
    	return 1;
	}
	if ( $self->session->is_empty ) {
	    $self->session->delete();
	    return 1;
	}
	return;
}

1;