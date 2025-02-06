package GenBoNoSqlLmdbSereal;
use strict;
use warnings;
use DBD::SQLite;
use Moo;

use Data::Dumper;
 use POSIX;
 use Sereal  qw(sereal_encode_with_object sereal_decode_with_object );
 
extends "GenBoNoSqlLmdb";



has sereal_encoder => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		#return Sereal::Encoder->new();
		return Sereal::Encoder->new({compress=>Sereal::SRL_ZLIB});
		return 0;
	},
);

has sereal_decoder => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return Sereal::Decoder->new({compress=>Sereal::SRL_ZLIB});
		return 0;
	},
);


sub decode {
	my ($self,$code) = @_;
	return sereal_encode_with_object( $self->sereal_encoder, $value );
	
}


sub encode {
		my ($self,$code) = @_;
		return sereal_decode_with_object( $self->sereal_decoder, $value );
}