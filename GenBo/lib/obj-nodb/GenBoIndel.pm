package GenBoIndel;
use strict;
use Moo;
extends "GenBoVariation";

has type_object => (
	is		=> 'ro',
	default	=> "indels_object",
);

1;