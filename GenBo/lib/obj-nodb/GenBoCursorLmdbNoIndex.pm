package GenBoCursorLmdbNoIndex;

use Moo;
use strict;
use warnings;
use Data::Dumper;
use LMDB_File qw(:flags :cursor_op :error);

 use Compress::Snappy;
use Storable qw/thaw freeze/;


has lmdb => (
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
has first_id => (
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
		my $cursor =  $self->lmdb->Cursor();
		if  ($self->lmdb->stat->{entries} == 0){
			$self->current(-99);
			return $cursor;
		}
		
		my $z;
		my $data;		
		$cursor->get($z, $data, MDB_LAST);
		$self->last_id($z);
		$self->current(-1);
		return $cursor;
}
);

sub  current_index {
	my ($self) =@_;
	confess() unless $self->current();
	return $self->current() -1;
}

sub  next_key {
	my ($self) =@_;
	my $cursor = $self->cursor_index();
	return (undef) if $self->current() eq -99;
		my $z;
		my $data;		
	if ($self->current() eq -1){
			$cursor->get($z, $data, MDB_FIRST);	
		
	}
	elsif  ($self->current() eq 1){
			$cursor->get($z, $data, MDB_NEXT);	
			$self->current(1) ;
		
	}
	
		if ($z ne $self->last_id){
			$self->current(1);
		}
		else {
				$self->current(-99);
		}
	
	return ($z,$data);
	return (undef) if $self->current() eq $self->end();
	
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