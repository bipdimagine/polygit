package GenBoSVDup;

use strict;
use Moo;
use Data::Dumper;
use Config::Std;
extends "GenBoStructuralVariant";

has isSVDuplication => (
	is		=> 'ro',
	default	=> 1
);


has type_object => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "svduplications_object";
	},
);

1;