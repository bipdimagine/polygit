package GenBoPrimerCache;
use strict;
use Moose;
use Data::Dumper;
use base 'Exporter';
extends 'GenBoPrimer', 'GenBoCache';

#sub sd {
#	my ( $self , $patient) =  @_;
#	return $self->{sd}->{$patient->id} if exists $self->{sd}->{$patient->id};
#	return -1;
#	confess() if exists $self->runs_object->{$patient->getRun->id()};
#	
#		return -1;
#	
#	
#}

sub level {
	my ( $self , $patient,$compute) =  @_;
		my $run_id = $patient->getRun->id();
		return -1 unless exists $self->{runs_object}->{$run_id};
		return $self->{level}->{$patient->id} if exists $self->{level}->{$patient->id};
		return -1;
		warn $patient->getRun->id();
		warn $self->id;
		confess();
}

sub cnv_score {
	my ($self,$patient,$compute) = @_;
	my $run_id = $patient->getRun->id();
	return -1 unless exists $self->{runs_object}->{$run_id};
	return $self->cnv->{$patient->id} if exists $self->cnv->{$patient->id};
	return -1;
		foreach my $z  (@{$self->getRuns}){
			warn $z->id;
		}
		warn $patient->getRun->id();
		warn $self->id;
		warn $patient->id."  ======== ";
	confess();
}
sub zscore {
	my ($self,$patient,$compute) = @_;
	
	confess() unless $patient;
	my $run_id = $patient->getRun->id();
	return 0 unless exists $self->{runs_object}->{$run_id};
	return $self->{zscore}->{$patient->id} if exists $self->{zscore}->{$patient->id};
	return 0;
	confess();
}

1;

