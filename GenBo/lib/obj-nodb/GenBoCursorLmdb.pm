package GenBoCursorLmdb;

use Moo;
use strict;
use warnings;
use Data::Dumper;


use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use LMDB_File qw(:flags :cursor_op :error);
 use Digest::MD5 qw(md5 md5_hex md5_base64);
 use Compress::Snappy;
use Storable qw/thaw freeze/;


has lmdb_index => (
	is		=> 'ro',
	required=> 1,
);


has start => (
	is		=> 'rw',
	required=> 1,
);


has end => (
	is		=> 'ro',
	required=> 1,
);

has current => (
	is		=> 'rw',
);
has current_data => (
	is		=> 'rw',
);
has last_id => (
	is		=> 'rw',
);
has last_data => (
	is		=> 'rw',
);


has cursor_index =>(
is		=> 'ro',
	lazy =>1,
default => sub {
		my $self = shift;
		my $cursor =  $self->lmdb_index->Cursor();
		my $data;
		my $z;
		$cursor->get($z, $data, MDB_LAST);
		$self->last_data($data);
		$self->last_id($z);

		if ($self->start == 0){
				my $t ;
				$cursor->get($t, $data, MDB_FIRST);
				$self->current($t);
		
		}
		else {
		
		$self->current($self->start);
		$cursor->get($self->start, $data, MDB_SET);
		}
		$self->current_data($data);
		return $cursor;
}
);

sub  current_index {
	my $self = shift;
	confess() unless $self->current();
	return $self->current() -1;
}

sub  next_key {
	my $self = shift;
	my $cursor = $self->cursor_index();
	confess() unless defined $self->current();
	return (undef) if $self->current() > $self->end();
	
	my $id1 = $self->current_data();
	my $key;
	my $value;
	if ($self->current ne $self->last_id){
		$cursor->get($key,$value,MDB_NEXT);
		$self->current($key);
		$self->current_data($value);
	}
	else {
		my $z = $self->end+1;
		$self->current($z);
		$self->current_data(undef);
	}
	
	
	 return ($id1);
	
}

1;