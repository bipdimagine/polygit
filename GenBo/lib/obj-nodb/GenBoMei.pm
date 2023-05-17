package GenBoMei;

use strict;
use Moo;

use Data::Dumper;
use Config::Std;
use GenBoCapture;
use Position;
extends "GenBoInsertion";


has name => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->getChromosome->name().'-'.$self->start().'ins-'.$self->mei_type;
	},
);

has isMei => (
	is		=> 'ro',
	default	=> sub {
		my $self = shift;
		return 1;
	},
);
has mei_type => (
	is		=> 'ro',

	default	=> sub {
		my $self = shift;
		return "ALU";
	},
);
has gnomad_id => (
	is		=> 'ro',
	lazy=> 1,
	default=> sub {
	my $self = shift;
	return $self->getChromosome->id().'-'.$self->start().'ins-'.$self->mei_type;
	}
	);

has alleles => (
	is		=> 'ro',
	lazy	=> 1,
	default=> sub {
		my $self = shift;
		return $self->getChromosome->id().'-'.$self->start().'-ins-'.$self->mei_type;
	},
);

has alamut_id => (
	is		=> 'ro',
	lazy	=> 1,
	default=> sub {
		my $self = shift;
		my $seq = $self->sequence();
		my $len = length($seq) - 1;
		my $alamut_id = 'chr'.$self->getChromosome->id().':'.$self->start().'_'.($self->start() + $len).'ins'.$seq;
		return $alamut_id;
	},
);




has length => (
	is		=> 'ro',
	#isa		=> 'Int',
	default	=> 1
);




sub protein_nomenclature {
	my ( $self, $prot ) = @_;
		confess() unless $prot->isProtein();
		my $pos = $self->getProteinPosition($prot);
		return "p.".$self->changeAA($prot).$pos."fs".$self->mei_type;
}

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

sub getSequence{
	my $self = shift;
	return $self->mei_type;
}

sub nomenclatureType {
	my ($self) = @_;
	return "insertion";
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
		return "c.".($cons->{orf_position}-1)."_MEI";
		
	}
	elsif ($self->isUtr($transcript)){
		#return;
		#warn $transcript->name();
		return "UTR-ins-".$self->mei_type;
	}
	else {
		return "ins-".$self->mei_type;
	}	
	#coding 
}

##### METHODS #####



1;