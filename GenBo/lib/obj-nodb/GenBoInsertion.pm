package GenBoInsertion;

use strict;
use Moo;
use Data::Dumper;
use Config::Std;
use GenBoCapture;
use Position;
extends "GenBoVariant";


	
has name => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		#if ($self->rs_name() and $self->rs_name() ne '.') { return $self->rs_name(); }
		my $len = length($self->var_allele())-1;
		if ($len > 15) {
			my $id = $self->getChromosome->id().'-'.$self->start().'-'.$self->nomenclatureType.'-'.$len;
			return $id;
		}
		else {
			return $self->gnomad_id;
		}
		if ($self->getChromosome->name eq "MT"){
			
		}
		return $self->id();
	},
);
has alleles => (
	is		=> 'ro',
	lazy=> 1,
	default=> sub {
		my $self = shift;
		my $previous_nt = $self->getChromosome->sequence($self->start()-1, $self->start()-1);
		return $previous_nt."/".$previous_nt.$self->sequence();
	},
);

has rocksdb_id => (
	is		=> 'ro',
	lazy=> 1,
	default => sub {
	my ($self) = @_;
	my ($chr,$pos,$ref,$alt) = split("-",$self->gnomad_id);
	$pos =~ s/ins//;
	$pos  = sprintf("%010d", ($pos));
	$alt = "X".$ref unless $alt; 
	my $seqid = $alt;
	$seqid = "+".substr($alt, 1);
	return  (sprintf("%010d", ($self->start-1))."!".$seqid);
	},
	
);


has kyoto_id => (
	is		=> 'rw',
	#isa		=> 'Str',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $chr = $self->getChromosome();
		my $type;
		 $type = 'Ins'; 
		return join("_", ($chr->name, $self->start(), $type));
		
	},
);
has type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "insertion";
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


has new_name => (
	is		=> 'rw',
	
	default	=> "0",
		lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hash = $self->getChromosome()->get_lmdb_public_db("insertions")->get($self->start);
		return $hash->{$self->sequence()}->{name} if $hash && exists $hash->{$self->sequence()}->{name};
		return 0;
		
	},
);

has type_object => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "insertions_object";
	},
);

 has public_db  =>(
	is		=> 'ro',
	lazy=> 1,
	default=> sub {
		my $self = shift;
		return  $self->project->lite_public_insertions();
	},

);
has kyoto_id_2 => (
	is		=> 'rw',
	#isa		=> 'Str',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $chr = $self->getChromosome();
		my $type;
		 $type = 'Ins'; 
		return join("_", ($chr->name, ($self->start()-1), $type));
		
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

has dbsnp => (
	is		=> 'ro',
	lazy=> 1,
	default=> sub {
		my $self = shift;
		my $res =  $self->getProject->public_data->dbsnp->get_indel(chr=>$self->getChromosome()->name ,id=>$self->kyoto_id,sequence=>$self->getSequence);
		unless ($res->{dbname}){
			 $res =  $self->getProject->public_data->dbsnp->get_indel(chr=>$self->getChromosome()->name ,id=>$self->kyoto_id_2,sequence=>$self->getSequence);
			 
		}
		return unless $res->{dbname};
		$res->{real_freq} = $res->{freq} ; 
		if ($res->{dbname} eq "pheno_snp") {$self->{clinical} = 1; }
		if (($res->{freq} == 0) && ($res->{dbname} eq "dbsnp")) {$res->{freq} = -1; };
		return  $res;
	}
);


has isInsertion => (
	is		=> 'ro',
	lazy =>1,
	default	=> sub {
		return 1;
	}
);

has length => (
	is		=> 'ro',
	#isa		=> 'Int',
	default	=> 1
);

sub polyphenStatus {
	my ( $self, $obj ) = @_;
	if (not $obj) {
		foreach my $t (@{$self->getTranscripts()}) {
			return 3 if ($self->isHighImpact($t));
		}
	}
	elsif ($obj->isGene()) {
		foreach my $t (@{$obj->getTranscripts()}) {
			return 3 if ($self->isHighImpact($t));
		}
	}
	elsif ($obj->isTranscript()) {
		return 3 if ($self->isHighImpact($obj));
	}
	return 0;
}

sub polyphenStatusText {
	my ( $self, $obj ) = @_;
	return;
}

sub siftStatus {
	my ( $self, $obj ) = @_;
	if (not $obj) {
		foreach my $t (@{$self->getTranscripts()}) {
			return 2 if ($self->isHighImpact($t));
		}
	}
	elsif ($obj->isGene()) {
		foreach my $t (@{$obj->getTranscripts()}) {
			return 2 if ($self->isHighImpact($t));
		}
	}
	elsif ($obj->isTranscript()) {
		return 2 if ($self->isHighImpact($obj));
	}
	return 0;
}

sub siftStatusText {
	my ( $self, $obj ) = @_;
	return;
}

sub polyphenScore {
	return '-';
}

sub siftScore {
	return '-';
}



sub protein_nomenclature {
	my ( $self, $prot ) = @_;
		confess() unless $prot->isProtein();
		my $pos = $self->getProteinPosition($prot);
		
			if ($self->isFrameshift($prot->getTranscripts->[0])){
			return "p.".$self->changeAA($prot).$pos."fs";
		}
		
		my $aa1 = $prot->getSequence($pos);
		my $aa2 = $prot->getSequence($pos+1);
		$aa2 ="" unless $aa2; 
		return "p.$aa1".$pos."_".$aa2.($pos+1)."ins".$self->changeAA($prot);
}

sub annotation_coding {
	my ( $self, $tr, $annot,$span ) = @_;
	my $protid  = $tr->getProtein->id;
	my $gid   = $tr->getGene->id();
	my $trid = $tr->id;
	my $project = $self->getProject();
	my $l = length($self->sequence()) % 3;
	my $start = $self->position($self->getChromosome())->start;
	my $consequence =  $tr->codonsConsequenceForVariations($self);
	$annot->{ $tr->id }->{coding}->{sequences} = $consequence;
	$annot->{ $protid }->{coding}->{sequences} = $consequence;
	
	if ( $tr->getSpanCodonStartEnd()->contains($start)){
		$annot->{$trid}->{mask} = $annot->{$trid}->{mask} |$project->getMaskCoding("phase");
		}
		
	elsif ($l > 0){
	
		
	#	$annot->{ $tr->id }->{coding}->{frameshift}++;
	#	$annot->{ $protid }->{coding}->{frameshift}++;
	#	$annot->{$gid}->{coding}->{frameshift}++;
	#	$annot->{all}->{coding}->{frameshift}++;
		$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("frameshift");
	}
	else {
		#$annot->{ $tr->id }->{coding}->{nonframeshift}++;
		#$annot->{ $protid }->{coding}->{nonframeshift}++;
		#$annot->{$gid}->{coding}->{nonframeshift}++;
		#$annot->{all}->{coding}->{nonframeshift}++;
		$annot->{$trid}->{mask} = $annot->{$trid}->{mask} |$project->getMaskCoding("nonframeshift");
	}
	
	
		
}

sub getSequence{
	my $self = shift;
	return $self->{var_allele};
}
sub nomenclatureType {
	my ($self) = @_;
	my $seq1 = $self->sequence();
	my $seq2 = $self->getChromosome()->getSequence($self->start-2,$self->start+$self->length+5);
	 #warn "dup".$self->id." ".$seq1." ".$self->length.$seq2  if $seq2 =~ /$seq1/;
	return "dup" if $seq2 =~ /$seq1/;
	return "ins";
}
sub constructNomenclature {
	my ( $self, $transcript, $debug ) = @_;
	confess("none transcript for nomenclature") unless $transcript;
	#return;
	#$self->annotation() unless exists $self->{annotation};
	my $id = $transcript->id();
	my $text;
	my $type = $self->nomenclatureType();
	if ($self->isCoding($transcript)) {
		my $cons = $self->{annotation}->{ $id }->{coding}->{sequences};
		#warn $self->{annotation}->{ $id }->{orf_position};
		#return "c.".$cons->{orf_position}.$cons->{seq_orf}.">".$cons->{seq_mut};	
		return $self->string_nomenclature("c.".($cons->{orf_position}-1)."_".($cons->{orf_end}-1)."$type",uc($cons->{seq_orf}));	
		
	}
	elsif ($self->isUtr($transcript)){
		#return;
		#warn $transcript->name();
		return $self->getNomenclatureForUtr($transcript);
	}
		else {
			my $r = $transcript->find_exon_intron($self->start,$self->end);
		if ($r && $r->{type} eq "exon"){
		my $pos_transcript = $transcript->translate_position($self->start);
		 my $seqv = $self->sequence();
		 my $seqr = $self->getChromosome()->sequence($self->start,$self->start);
		 die() unless  $seqr;
		if ($transcript->strand() == -1 ){
				 $seqv = BioTools::complement_sequence($seqv);
				 $seqr = BioTools::complement_sequence($seqr);
			}
			return $self->string_nomenclature("n.".$pos_transcript."$type",$seqv);
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
			
		 return $self->string_nomenclature("c.$pos$st_dist"."$type",$seqv);
	}
	return;
	#coding 
	
	
	
}



##### METHODS #####



1;