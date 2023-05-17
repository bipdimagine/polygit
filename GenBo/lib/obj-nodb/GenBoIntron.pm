package GenBoIntron;

use strict;
use Moo;

extends "GenBoExon";

has isIntron => (
	is		=> 'ro',
	default	=> 1,
);
has isExon => (
	is		=> 'ro',
	default	=> 0,
);
1;