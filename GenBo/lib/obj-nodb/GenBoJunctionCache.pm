package GenBoJunctionCache;
use strict;
use Moo;
use Carp;
use Data::Dumper;
extends  'GenBoJunction','GenBoVariantCache';



sub setPatients {
	my $self = shift;
	my $h;
	foreach my $patient (@{$self->getProject->getPatients()}) {
		$h->{$patient->id} = undef if exists $self->annex->{$patient->name()};
	}
	return $h;
}


1;