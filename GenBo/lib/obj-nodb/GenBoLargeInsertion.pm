package GenBoLargeInsertion;
use strict;
use Moose;
use MooseX::Method::Signatures;
use Data::Dumper;
use Config::Std;
use GenBoCapture;
use Position;
extends "GenBoCnv";


has imprecise => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		return undef;
	},
);

has isLargeInsertion => (
	is		=> 'rw',
	default	=> 1,
);
has isCnv => (
	is		=> 'ro',
	default	=> undef,
);

has length => (
	is		=> 'ro',
	#isa		=> 'Int',
	default	=> sub {
		my $self = shift;
		return 1;
	},
);


has string_dejavu => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		my $h = $self->project->get_deja_vu_from_overlap($self->getChromosome->name,$self->start,$self->start+$self->length,"IL");
		return $h->[1] if $h;
		return "";
	},
);
has sv_type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "INS";
	},
);

sub limit_cnv_value {
	my ($self,$value) = @_;
	return 1 if ($value >= 1.4);
	return;
}
sub limit_cnv_value_ho {
	my ($self,$value) = @_;
	return 1 if ($value >= 1.8);
	return;
}






has type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "large_insertion";
	},
);


has structural_type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "l_ins";
	},
);

has type_public_db => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "insertions";
	},
);

has type_object => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "large_insertions_object";
	},
);


has alleles => (
	is		=> 'ro',
	lazy	=> 1,
	default=> sub {
		my $self = shift;
		return self->getChromosome->id().'-'.$self->start().'-ins-?' if $self->imprecise;
		return self->getChromosome->id().'-'.$self->start().'-ins-'.$self->sequence();
	},
);





sub constructNomenclature {
	my ( $self, $transcript, $debug ) = @_;
	return $self->name();
	return;
	#coding 
}



sub return_interval_tree {
	my $self = shift;
	confess();
	return $self->getChromosome->large_duplication_interval_tree();
}
has allele_length => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return -1;
	},
);


1;