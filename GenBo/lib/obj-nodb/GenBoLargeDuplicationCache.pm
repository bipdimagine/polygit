package GenBoLargeDuplicationCache;
use strict;
use Moo;

use Data::Dumper;
extends 'GenBoVariantCache','GenBoLargeDuplication';




has allele_length => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return abs($self->start - $self->end)+1;
	},
);

1;