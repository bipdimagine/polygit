package GenBoNoSqlLmdbPipeline;
use strict;
use warnings;
use DBD::SQLite;
use Moose;
use MooseX::Method::Signatures;
use Data::Dumper;
 use JSON::XS;
 use POSIX;
 use LMDB_File qw(:flags :cursor_op :error);
extends "GenBoNoSqlLmdb";

has dir_prod => (
	is		=> 'ro',
	required=> 1,
);

has filename_prod =>(
	is		=> 'rw',
	lazy=>1,
default => sub {
		my $self =shift;
		my $filename = $self->dir_prod."/".$self->name;
		return $filename;
		
}
);
has filename_prod_index =>(
	is		=> 'rw',
	lazy=>1,
default => sub {
		my $self =shift;
		my $filename = $self->dir_prod."/".$self->name."_index";
		return $filename;
		
}
);

sub clean_files {
	my ($self) = @_;
	die() unless $self->mode eq "c";
	unlink $self->filename  if -e $self->filename;
	unlink $self->filename_index  if -e $self->filename_index;
	unlink $self->filename_prod  if -e $self->filename_prod;
	unlink $self->filename_prod_index  if -e $self->filename_prod_index;
}
sub close {
	my ($self) = @_;
	 $self->SUPER::close();
	 $self->_put_in_prod();
}
sub _put_in_prod {
	my ($self) = @_;
	warn $self->dir_prod.'-'.$self->filename;
	system('mkdir -p '.$self->dir_prod.'&& chmod a+rwx '.self->dir_prod ) unless -e $self->dir_prod;
	system("rsync -av ".$self->filename.' '.$self->dir_prod.'/' )  if -e $self->filename;
	warn "rsync -av ".$self->filename.' '.$self->dir_prod.'/';
	system("rsync -av ".$self->filename_index.' '.$self->dir_prod.'/' )  if -e $self->filename_index;

	unlink $self->filename;
	unlink $self->filename_index  if -e $self->filename_index;
	
}
1;