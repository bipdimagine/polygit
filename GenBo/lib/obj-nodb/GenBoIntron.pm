package GenBoIntron;

use strict;
use Moose;

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