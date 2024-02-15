package GenBoLargeInsertion;
use strict;
use Moo;
use Data::Dumper;
use Config::Std;
use GenBoCapture;
use Carp;
use Position;
extends "GenBoInsertion";

has imprecise => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		return undef;
	},
);

has rocksdb_id => (
	is		=> 'ro',
	lazy=> 1,
	default => sub {
	my ($self) = @_;
	#my ($chr,$pos,$ref,$alt) = split("-",$self->gnomad_id);
	my $pos  = sprintf("%010d", $self->start());
	return  ($pos."!+".$self->sequence);
	},
);

has name => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->getChromosome->id().'-'.$self->start().'-'.$self->nomenclatureType.'-'.$self->sequence;
	},
);

#has gnomad_id => (
#	is		=> 'ro',
#	lazy=> 1,
#	default=> sub {
#	my $self = shift;
#	warn $self->vcf_id;
#	die();
#	}
#	);
	
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
		return $self->sequence();
	},
);






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
		return "?" if $self->imprecise;
		return length($self->sequence);
	},
);

sub annotation_coding {
	my ( $self, $tr, $annot,$span ) = @_;
	my $protid  = $tr->protein();
	my $gid   = $tr->getGene->id();
	my $trid = $tr->id;
	my $project = $self->getProject();
	my $l = length($self->sequence()) % 3;
	my $start = $self->position($self->getChromosome())->start;
	my $consequence =  $tr->codonsConsequenceForLargeInsertion($self);
	$annot->{ $tr->id }->{coding}->{sequences} = $consequence;
	$annot->{ $protid }->{coding}->{sequences} = $consequence;
	$annot->{ $tr->id }->{coding}->{frameshift}++;
	$annot->{ $protid }->{coding}->{frameshift}++;
	$annot->{$gid}->{coding}->{frameshift}++;
	$annot->{all}->{coding}->{frameshift}++;
	$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("frameshift");
}


sub constructNomenclature {
	my ( $self, $transcript, $debug ) = @_;
	confess("none transcript for nomenclature") unless $transcript;
	#return;
	#$self->annotation() unless exists $self->{annotation};
	my $id = $transcript->id();
	my $text;
	if ($self->isCoding($transcript)) {
		my $cons = $self->{annotation}->{ $id }->{coding}->{sequences};
		return "c.".($cons->{orf_position}-1)."_ins";
		
	}
	elsif ($self->isUtr($transcript)){
		#return;
		#warn $transcript->name();
		return "UTR-ins-".$self->length;
	}
	else {
		return "ins-".$self->length;
	}	
	#coding 
}


1;