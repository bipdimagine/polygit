package GenBoJunctionCache;
use strict;
use Moose;
use MooseX::Method::Signatures;
use Data::Dumper;
extends  'GenBoVariantCache','GenBoJunction';



sub setPatients {
	my $self = shift;
	my $h;
	foreach my $patient (@{$self->getProject->getPatients()}) {
		$h->{$patient->id} = undef if exists $self->annex->{$patient->name()};
	}
	return $h;
}

1;