package GenBoCoverageLmdb;
use strict;
use Moose;
extends "GenBoCoverage";
use Data::Dumper;
use Config::Std;
use Storable qw(store retrieve freeze thaw fd_retrieve dclone);
use Clone 'clone';
use List::Util qw(first max maxstr min minstr reduce shuffle sum0);

has chromosome => (
	is		=> 'rw',
	required =>1,
);

has patient => (
	is		=> 'rw',
	required =>1,
);

 
has array => (
	is		=> 'ro',
		lazy=>1,
	default => sub {
	my $self = shift;
		return $self->patient->depth($self->chromosome->name,$self->start,$self->end);
	#	my $gc  = GenBoCoverage->new(start=>$start,end=>$end,array=>\@v);
	},

);

sub add_value {
	my ($self,$pos_array,$pos_genomic) = @_;
	confess();
}
1;