package GenBoLargeDeletion;

use strict;
use Moo;
use Data::Dumper;
use Config::Std;
use GenBoCapture;
use Position;
use Carp;
extends 'GenBoCnv';


has rocksdb_id => (
	is		=> 'ro',
	lazy=>1,
	default => sub {
		my $self = shift;
	my ($chr,$pos,$ref,$alt) = split("-",$self->gnomad_id);
	 $pos  = sprintf("%010d", $self->start());
		return  ($pos."!".$self->length);
	},
	
);

has name => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->getChromosome->id().'-'.$self->start().'-'.$self->nomenclatureType.'-'.$self->length;
	},
);

has gnomad_id => (
	is		=> 'ro',
	lazy=> 1,
	default=> sub {
	my $self = shift;
	return $self->name;
	}
	);
	
has isLargeDeletion => (
	is		=> 'ro',
	default	=> 1,
);



has type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		return "large_deletion";
	},
);

has sv_type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		return "DEL";
	},
);

sub limit_cnv_value {
	my ($self,$value) = @_;
	return 1 if  $value <= 0.6;
	return undef;
}

sub limit_cnv_value_ho {
	my ($self,$value) = @_;
	return 1 if  $value <= 0.15;
	return undef;
}



has structural_type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "l_del";
	},
);

has string_dejavu => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		my $h = $self->project->get_deja_vu_from_overlap($self->getChromosome->name,$self->start,$self->start+$self->length,"DL");
		return $h->[1] if $h;
		return "";
	},
);


has type_public_db => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "deletions";
	},
);

has type_object => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "large_deletions_object";
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


sub setLargeDeletions {
	my $self = shift;
	confess();

}

sub return_interval_tree {
	my $self = shift;
	return $self->getChromosome->large_deletion_interval_tree();
}


1;