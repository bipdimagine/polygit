package GenBoInversion;

use strict;
use Moose;
use MooseX::Method::Signatures;
extends 'GenBoCnv';


has isInversion => (
	is		=> 'ro',
	default	=> 1,
);


has isCnv => (
	is		=> 'ro',
	default	=> undef,
);

has type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		return "inversion";
	},
);

has sv_type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		return "INV";
	},
);

sub limit_cnv_value {
	my ($self,$value) = @_;
	
	return undef;
}

sub limit_cnv_value_ho {
	my ($self,$value) = @_;
	return undef;
}

has string_dejavu => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		#my $h = $self->project->get_deja_vu_from_overlap($self->getChromosome->name,$self->start,$self->start+$self->length,"DL");
		#return $h->[1] if $h;
		return "";
	},
);


has type_public_db => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "invertions";
	},
);

has type_object => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "inversions_object";
	},
);


has alleles => (
	is		=> 'ro',
	lazy	=> 1,
	default=> sub {
		my $self = shift;
		return 'del '.$self->start().' to '.$self->end().' (length:'.$self->length().')';
	},
);

sub constructNomenclature {
	my ( $self, $transcript, $debug ) = @_;
	confess("none transcript for nomenclature") unless $transcript;
	 return $self->name;
	
	#coding 
}


sub delete_sequence {
	my ($self,$obj) = @_;
	confess() if $obj;
	my $chr = $self->getChromosome();
	return $chr->sequence($self->position($chr)->start,$self->position($chr)->end);
}






1;


1;