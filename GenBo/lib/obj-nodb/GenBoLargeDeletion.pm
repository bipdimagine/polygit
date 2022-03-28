package GenBoLargeDeletion;

use strict;
use Moose;
use MooseX::Method::Signatures;
use Data::Dumper;
use Config::Std;
use GenBoCapture;
use Position;
extends 'GenBoCnv';



has name => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->getChromosome->id().'-'.$self->start().'-cnv-del-'.abs($self->length()-1);
	},
);

has isLargeDeletion => (
	is		=> 'ro',
	default	=> 1,
);
has isCnv => (
	is		=> 'rw',
	default	=> 1,
);
has type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "large_deletion";
	},
);

has sv_type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
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

has spliceAI => (
	is   => 'rw',
	lazy =>1,
	default => sub {   return "-" ; }
);

has alleles => (
	is		=> 'ro',
	lazy	=> 1,
	default=> sub {
		my $self = shift;
		return 'del '.$self->start().' to '.$self->end().' (length:'.$self->length().')';
	},
);




sub annotation_coding {
	my ( $self, $tr, $annot ) = @_;
	my $span =  $self->getGenomicSpan()->intersection( $tr->getSpanCoding );
	my $prot  = $tr->getProtein();
	my $gid   = $tr->getGene->id();
	my $trid = $tr->id;
	my $namep = $prot->name();
	my $pos   = $self->start();
	my $seq   = $self->sequence();
	my @array = $span->as_array();
	my $project = $self->getProject();
	my $start_tr = $array[0];
	my $end_tr = $array[-1];
	my $tres = scalar(@array) % 3;
	my $pos_transcript = $tr->translate_position($pos);
	$end_tr = $start_tr unless $end_tr; 
	my $consequence =  $tr->codonsConsequenceForVariations($self,$start_tr,$end_tr);
	$annot->{ $tr->id }->{coding}->{sequences} = $consequence;
	$annot->{ $prot->id }->{coding}->{sequences} = $consequence;
	$annot->{ $tr->id }->{coding}->{frameshift}++;
	$annot->{ $prot->id }->{coding}->{frameshift}++;
	$annot->{$gid}->{coding}->{frameshift}++;
	$annot->{all}->{coding}->{frameshift}++;
	$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("frameshift");
}

sub constructNomenclature {
	my ( $self, $transcript, $debug ) = @_;
	confess("none transcript for nomenclature") unless $transcript;
	my $id = $transcript->id();
	my $text;

		if ($self->isCoding($transcript) && $transcript->getProtein) {
		my $cons = $self->annotation()->{ $id }->{coding}->{sequences};
		unless ($cons){
		return "c.largedeletion".$self->length;
		}
		my $seq =$cons->{codon} ;
	
		#warn $self->{annotation}->{ $id }->{orf_position};
		#return "c.".$cons->{orf_position}.$cons->{seq_orf}.">".$cons->{seq_mut};	
		return "c.".$cons->{orf_position}."del".uc($cons->{seq_orf}) if $self->length == 1;	
		
		return "c.".$cons->{orf_position}."_".$cons->{orf_end}."del".uc($cons->{seq_orf})	;
		
	}
	elsif ($self->isUtr($transcript)){
		return $self->getNomenclatureForUtr($transcript);
	}
		else {
				my $r = $transcript->find_exon_intron($self->start,$self->end);
		if ($r and $r->{type} eq "exon"){
		my $pos_transcript = $transcript->translate_position($self->start);
		 my $seqv = $self->sequence();
		 my $seqr = $self->getChromosome()->sequence($self->start,$self->start);
		 die() unless  $seqr;
		if ($transcript->strand() == -1 ){
				 $seqv = BioTools::complement_sequence($seqv);
				 $seqr = BioTools::complement_sequence($seqr);
			}
			return "n.".$pos_transcript."+".length($seqr)."del";
	}
		my ($dist,$pos) = $self->getPositionForNomenclatureIntronic($transcript);
		 my $st_dist = $dist;
		 

		 my $seqv = $self->delete_sequence();
		 	if ($transcript->strand() == -1 ){
				 $seqv = BioTools::complement_sequence($seqv);
			}
		 return"c.$pos$st_dist"."del".length($seqv);
	}
	return;
	#coding 
}

sub delete_sequence {
	my ($self,$obj) = @_;
	if ($obj){
		my $id = $obj->id;
		if (exists $self->{annotation}->{ $id }->{coding}->{sequences}){
			my $cons = $self->{annotation}->{ $id }->{coding}->{sequences};
			return uc($cons->{seq_orf});
		}
	}
	my $chr = $self->getChromosome();
	return $chr->sequence($self->position($chr)->start,$self->position($chr)->end);
}

sub setLargeDeletions {
	my $self = shift;
	foreach my $o (@{$self->getChromosome->getLargeDeletions()}){
		next unless ($self->getGenomicSpan()->intersection($o->getGenomicSpan()));
#		next if ($self->id() eq $o->id());
		next if ($self->name() eq $o->name());
		my $areSameCNV = $self->getChromosome->areSameCNV($self->start(), $self->end(), $o->start(), $o->end());
		next unless ($areSameCNV);
		$o->{references_object}->{$self->id} = undef;
		$self->{$o->type_object}->{$o->id} = undef;
	}
	return $self->{large_deletions_object} ;
}

sub return_interval_tree {
	my $self = shift;
	return $self->getChromosome->large_deletion_interval_tree();
}



1;