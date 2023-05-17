package GenBoLargeDuplication;

use strict;
use Moo;
use Data::Dumper;
use Config::Std;
use GenBoCapture;
use Position;
extends "GenBoCnv";




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
		return "DUP";
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






has length => (
is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->allele_length;
	},
);


has allele_length => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return abs($self->start - $self->end)+1;
	},
);


has isLargeDuplication => (
	is		=> 'rw',
	default	=> 1,
);

has type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "large_duplication";
	},
);

has isCnv => (
	is		=> 'rw',
	default	=> 1,
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
		return "large_duplications_object";
	},
);

has spliceAI => (
	is   => 'rw',
	lazy =>1,
	default => sub {  return "-" ;}
);

has alleles => (
	is		=> 'ro',
	lazy	=> 1,
	default=> sub {
		my $self = shift;
		return 'dup-'.$self->start().' to '.$self->end().' (length:'.$self->length().')';
	},
);




sub constructNomenclature {
	my ( $self, $transcript, $debug ) = @_;
	return $self->name();
	#coding 
}


sub return_interval_tree {
	my $self = shift;
	return $self->getChromosome->large_duplication_interval_tree();
}



1;