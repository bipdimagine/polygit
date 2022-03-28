package GenBoLargeDuplication;

use strict;
use Moose;
use MooseX::Method::Signatures;
use Data::Dumper;
use Config::Std;
use GenBoCapture;
use Position;
extends "GenBoCnv";




has name => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $id = $self->getChromosome->id().'-'.$self->start().'-cnv-dup-'.abs($self->allele_length()-1);
		return $id;
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
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->allele_length() if $self->allele_length();
		return abs($self->start - $self->end);
	},
);

has allele_length => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return length($self->var_allele())-1;
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

has structural_type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "l_dup";
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
		return 'dup '.$self->start().' to '.$self->end().' (length:'.$self->length().')';
	},
);





sub annotation_coding {
	my ( $self, $tr, $annot,$span ) = @_;
	my $protid  = $tr->protein();
	my $gid   = $tr->getGene->id();
	my $trid = $tr->id;
	my $project = $self->getProject();
	warn $self->start." ".$self->end unless $self->sequence();
	confess(Dumper $self->annex) unless $self->sequence();
	my $l = length($self->sequence()) % 3;
	my $start = $self->position($self->getChromosome())->start;
	my $consequence =  $tr->codonsConsequenceForDuplication($self);
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
		#warn $self->{annotation}->{ $id }->{orf_position};
		#return "c.".$cons->{orf_position}.$cons->{seq_orf}.">".$cons->{seq_mut};	
		
		
		return "c.".($cons->{orf_position}-1)."_".($cons->{orf_end}-1)."ins".uc($cons->{seq_orf})	
		
	}
	elsif ($self->isUtr($transcript)){
		#return;
		#warn $transcript->name();
		return $self->getNomenclatureForUtr($transcript);
	}
		else {
			my $r = $transcript->find_exon_intron($self->start,$self->end);
			return "?" unless $r;
		if ($r->{type} eq "exon"){
		my $pos_transcript = $transcript->translate_position($self->start);
		 my $seqv = $self->sequence();
		 my $seqr = $self->getChromosome()->sequence($self->start,$self->start);
		 die() unless  $seqr;
		if ($transcript->strand() == -1 ){
				 $seqv = BioTools::complement_sequence($seqv);
				 $seqr = BioTools::complement_sequence($seqr);
			}
			return "n.".$pos_transcript."ins".$seqv;
		#return n
	}	
		my ($dist,$pos) = $self->getPositionForNomenclatureIntronic($transcript);
		 my $st_dist = $dist;
		confess($dist) unless defined $dist;			 
		 $st_dist="+".$dist if $dist >0;
		 #return;
		  my $seqv = $self->sequence();
		  #warn $self->getSequence()." ".$seqv;
		 # die();
		 	if ($transcript->strand() == -1 ){
				 $seqv = BioTools::complement_sequence($seqv);
			}
		 return"c.$pos$st_dist"."dup".length($seqv);
	}
	return;
	#coding 
}

has intspan => (
	is		=> 'rw',
	reader	=> 'getGenomicSpan',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $intSpan = Set::IntSpan::Fast::XS->new();
    	$intSpan->add_range($self->start(), ($self->start() + $self->allele_length()));
    	return $intSpan;
	},
);

sub getReference {
	my $self = shift;
	return $self->getChromosome->getReferences($self->start()-1, ($self->start() + $self->allele_length() +1))->[0];
}

sub setLargeDuplications {
	my $self = shift;
	foreach my $o (@{$self->getChromosome->getLargeDuplications()}){
		next unless ($self->getGenomicSpan()->intersection($o->getGenomicSpan()));
#		next if ($self->id() eq $o->id());
		next if ($self->name() eq $o->name());
		my $areSameCNV = $self->getChromosome->areSameCNV($self->start(), ($self->start() + $self->allele_length()), $o->start(), ($o->start() + $o->allele_length()));
		next unless ($areSameCNV);
		$o->{references_object}->{$self->id} = undef;
		$self->{$o->type_object}->{$o->id} = undef;
	}
	return $self->{large_duplications_object} ;
}

sub return_interval_tree {
	my $self = shift;
	return $self->getChromosome->large_duplication_interval_tree();
}



1;