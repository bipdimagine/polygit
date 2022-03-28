package GenBoSomaticGroup;

use strict;
use Moose;
use Data::Dumper;
extends "GenBoSampleGroup";



has somatic_details => (
	is		 => 'rw',
	required => 1,
);



##### METHODS ######



sub getPatientsSomatic {
	my $self = shift;
	my @lPatients;
	foreach my $patient (@{$self->getPatients()}) {
		push(@lPatients, $patient) if ($patient->isSomatic());
	}
	return \@lPatients;
}

sub getPatientsGerminal {
	my $self = shift;
	my @lPatients;
	foreach my $patient (@{$self->getPatients()}) {
		push(@lPatients, $patient) if ($patient->isGerminal());
	}
	return \@lPatients;
}



1;