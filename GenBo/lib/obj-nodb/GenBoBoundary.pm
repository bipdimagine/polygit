package GenBoBoundary;

use strict;
use Moose;
use MooseX::Method::Signatures;
extends 'GenBoCnv';


has isBoundary => (
	is		=> 'ro',
	default	=> 1,
);

has name => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $al = $self->allele_length;
		
		#$al = "*" if $al == -1;
		my $id = $self->getChromosome->name().'-'.$self->start().'-trans-'. $self->mate_chr.":".$self->mate_pos;
		return $id;
	},
);

has isCnv => (
	is		=> 'ro',
	default	=> undef,
);

has mate_chr =>(
is		=> 'ro',
	default	=> "XX",
);
has mate_pos =>(
is		=> 'ro',
	default	=> "XX",
);
has mate_id => (
	is		=> 'ro',
		default	=> sub {
		my $self = shift;
		return $self->mate_chr."_".$self->mate_pos."_".$self->ref_allele."_bnd-".$self->getChromosome->name.":".$$self->position;
	},
		
);

has event_id => (
	is		=> 'ro',
);

has type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		return "boundary";
	},
);

has sv_type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		return "BND";
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
		return "boundaries";
	},
);

has type_object => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "boundaries_object";
	},
);


has alleles => (
	is		=> 'ro',
	lazy	=> 1,
	default=> sub {
		my $self = shift;
		return "*";
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