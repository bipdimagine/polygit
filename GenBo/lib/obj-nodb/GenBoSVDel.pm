package GenBoSVDel;

use strict;
use Vcf;
use Moose;
use Data::Dumper;
use Config::Std;
use Position;
extends "GenBoStructuralVariant";

has isSVDeletion => (
	is		=> 'ro',
	default	=> 1
);

has type_object => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
	
		return "svdeletions_object";
	},
);
1;