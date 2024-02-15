package GenBoReference;

use strict;
use Vcf;
use Moo;
use Carp;
use Data::Dumper;
use Config::Std;
extends "GenBoGenomic";


has deja_set => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	default	=> sub {
		return;
	}
);

has parse_hgmd => (
	is	 => 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		#return 1 if ($self->getProject->is_hgmd_project());
		return;
	},
);

my $set_variants;

sub setStructuralVariants {
	my ($self, $typeVar) = @_;
	my @objs;
	foreach my $p (@{$self->getProject()->getPatients}){
		foreach my $o (@{$p->setStructuralVariantsForReference($self)}){
			$o->{references_object}->{$self->id} = undef;
			$self->{$o->type_object}->{$o->id} = undef;
		}
	}
}

sub setVariants {
	my ($self, $typeVar) = @_;
	return if $self->deja_set;
	$self->deja_set (1);
	#$set_variants = 1;
	my @objs;
	foreach my $p (@{$self->getProject()->getPatients}){
		foreach my $o (@{$p->setVariantsForReference($self, $typeVar)}){
			$o->{references_object}->{$self->id} = undef;
			$self->{$o->type_object}->{$o->id} = undef;
		}
	}
	
}
sub setVariations {
	my $self = shift;
	$self->setVariants('variations');
	return $self->{variations_object} ;
}

sub setDeletions {
	my $self = shift;
	$self->setVariants('deletions');
	return $self->{deletions_object} ;
}

sub setInsertions {
	my $self = shift;
	$self->setVariants('insertions');
	return $self->{insertions_object} ;
}
sub setMnps {
	my $self = shift;
	$self->setVariants('mnps');
	return $self->{mnps_object} ;
}
sub setLargeDeletions {
	my $self = shift;
	$self->setVariants('large_deletions');
	return $self->{large_deletions_object} ;
}
sub setInversions {
	my $self = shift;
	$self->setVariants('inversions');
	return $self->{inversions_object} ;
}
sub setBoundaries {
	my $self = shift;
	$self->setVariants('boundaries');
	return $self->{boundaries_object} ;
}
sub setLargeInsertions {
	my $self = shift;
	$self->setVariants('large_insertions');
	return $self->{large_insertions_object} ;
}
sub setLargeDuplications {
	my $self = shift;
	$self->setVariants('large_duplications');
	return $self->{large_duplications_object} ;
}
sub refreshObjects{
	my $self = shift;
	$self->setVariations();
	$self->setDeletions();
	$self->setInsertions();
	$self->setLargeDuplications();
	$self->setLargeDeletions();
	$self->setMnps();
}

1;