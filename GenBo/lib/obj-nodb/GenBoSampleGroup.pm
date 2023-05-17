package GenBoSampleGroup;

use strict;
use Moo;
use Data::Dumper;
extends "GenBo";




has hash_samples => (
	is		 => 'rw',
	required => 1,
);

has hash_fictif_samples => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my %h;
		foreach my $sample_name (keys %{$self->hash_samples()}) { $h{$sample_name} = undef; }
		foreach my $sample (@{$self->project->getPatients()})   { delete $h{$sample->name()}; }
		return \%h;
	},
);


##### METHODS #####



sub getPatients {
	my $self = shift;
	my @lPat;
	foreach my $pat_name (keys %{$self->hash_samples()}) {
		unless (exists $self->hash_fictif_samples->{$pat_name}) {
			push(@lPat, $self->project->getPatient($pat_name));
		}
	}
	return \@lPat;
}

sub getPatientsHealthy() {
	my $self = shift;
	my @lPat;
	foreach my $patient (@{$self->getPatients()}) {
		push(@lPat, $patient) if ($patient->isHealthy());
	}
	return \@lPat;
}

sub getPatientsIll() {
	my $self = shift;
	my @lPat;
	foreach my $patient (@{$self->getPatients()}) {
		push(@lPat, $patient) if ($patient->isIll());
	}
	return \@lPat;
}


1;