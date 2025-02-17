package GenBoDeletion;

use strict;
use Vcf;
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
		if ($self->length() > 20) {
			return $self->getChromosome->id().'-'.$self->start().'-del-'.$self->length();
		}
		else {
			return $self->gnomad_id;
		}
		#if ($self->getChromosome->name eq "MT"){
		#	return $self->start.$self->sequence();
		#}
		#return $self->id();
	},
);

has gnomad_id => (
	is		=> 'ro',
	lazy=> 1,
	default=> sub {
	my $self = shift;
	my $b =  $self->getChromosome->sequence($self->start-1,$self->start-1);
	my $name = $self->getChromosome->name."-".($self->start()-1)."-". $b.$self->getChromosome->sequence($self->start,$self->end)."-".$b;
	return $name;
	my $vn=$self->vcf_id;
	#confess unless $vn;
	$vn =~ s/_/-/g;
	$vn=~ s/chr//;
	return $vn;
	}
	);


has alleles => (
	is		=> 'ro',
	lazy	=> 1,
	default=> sub {
		my $self = shift;
		return $self->getChromosome->id().'-'.$self->start().'-del-'.$self->length() if ($self->length() > 20);
		return $self->getChromosome->sequence($self->start(), $self->end())."/".$self->sequence();
	},
);
sub nomenclatureType {
	my ($self) = @_;
	return "del";
}

has rocksdb_id => (
	is		=> 'ro',
	lazy=>1,
	default => sub {
		my $self = shift;
		my ($chr,$pos,$ref,$alt) = split("-",$self->gnomad_id);
		return sprintf("%010d", ($pos))."!".$self->length();
	},
	
);


has kyoto_id => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $chr = $self->getChromosome();
		my $type;
	 	$type = 'Del'; 
		return  join("_", ($chr->name, $self->start(), $type));
		
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

sub alternate_allele {
	my ($self) = @_;
	return $self->delete_sequence();
} 

has alamut_id => (
	is		=> 'ro',
	lazy	=> 1,
	default=> sub {
		my $self = shift;
		my $len = length($self->ref_allele()) - 1;
		my $alamut_id = 'chr'.$self->getChromosome->id().':'.$self->start().'_'.($self->start() + $len).'del';
		return $alamut_id;
	},
);




 has public_db  =>(
	is		=> 'ro',
	lazy=> 1,
	default=> sub {
		my $self = shift;
		return  $self->project->lite_public_deletions();
	},

);

has type_object => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "deletions_object";
	},
);

has type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "deletion";
	},
);

has isDeletion => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> 1
);



sub theoric_vcf_id {
	my ($self) = @_;
	return $self->{theoric_vcf} if exists $self->{theoric_vcf} ;
	my $start = $self->start;
	my $pos = $self->start -1;
	my $ref = $self->getChromosome()->sequence($start-1,$start-1);
	my $alt = $ref.$self->alternate_allele;
	
	$self->{theoric_vcf}->{chromsome} = $self->getChromosome()->fasta_name;
	$self->{theoric_vcf}->{ref} = $ref;
	$self->{theoric_vcf}->{alt} = $alt;
	$self->{theoric_vcf}->{position} = $pos;
	
	
	
	return $self->{theoric_vcf};
}


sub checkLargeDeletion_newVariant {
	my $self = shift;
	return 1 if (length($self->ref_allele()) >= 30);
	return;
}

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
	unless ( $self->getGenomicSpan()->intersection($tr->getSpanCodonStartEnd())->is_empty){
		$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("phase");
	}
	if ($tres > 0 ){
		#$annot->{ $tr->id }->{coding}->{frameshift}++;
		#$annot->{ $prot->id }->{coding}->{frameshift}++;
	#	$annot->{$gid}->{coding}->{frameshift}++;
	#	$annot->{all}->{coding}->{frameshift}++;
		$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("frameshift");
	}
	else {
	#	$annot->{ $tr->id }->{coding}->{nonframeshift}++;
	#	$annot->{ $prot->id }->{coding}->{nonframeshift}++;
	#	$annot->{$gid}->{coding}->{nonframeshift}++;
	#	$annot->{all}->{coding}->{nonframeshift}++;
		$annot->{$trid}->{mask} = $annot->{$trid}->{mask} |$project->getMaskCoding("nonframeshift");
	}
}


##### METHODS #####
sub constructNomenclature {
	my ( $self, $transcript, $debug ) = @_;
	confess("none transcript for nomenclature") unless $transcript;
	my $id = $transcript->id();
	my $text;

		if ($self->isCoding($transcript)) {
		my $cons = $self->annotation()->{ $id }->{coding}->{sequences};
		my $seq =$cons->{codon} ;
	
		#warn $self->{annotation}->{ $id }->{orf_position};
		#return "c.".$cons->{orf_position}.$cons->{seq_orf}.">".$cons->{seq_mut};	
		return $self->string_nomenclature("c.".$cons->{orf_position}."del",uc($cons->{seq_orf})) if $self->length == 1;	
		
		return $self->string_nomenclature("c.".$cons->{orf_position}."_".$cons->{orf_end}."del",uc($cons->{seq_orf}));	
		
	}
	elsif ($self->isUtr($transcript)){
		return $self->getNomenclatureForUtr($transcript);
	}
		else {
				my $r = $transcript->find_exon_intron($self->start,$self->end);

		if (exists $r->{type} && $r->{type} eq "exon"){
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
		$st_dist="+".$dist if $dist >0;		 


		 my $seqv = $self->delete_sequence();
		 	if ($transcript->strand() == -1 ){
				 $seqv = BioTools::complement_sequence($seqv);
			}
		 return $self->string_nomenclature("c.$pos$st_dist"."del",$seqv);
	}
	return;
	#coding 
}


sub protein_nomenclature {
	my ( $self, $prot ) = @_;
		confess() unless $prot->isProtein();
		my $pos = $self->getProteinPosition($prot);
		my $seq1 = $self->annotation->{ $prot->id }->{coding}->{sequences}->{aa};
		my @seqa = split("",$seq1);

			
			
			if ($self->isFrameshift($prot->getTranscripts->[0])){
				
			return "p.". $seqa[0].$pos."fs";
		}
		
		my $aa1 = $seqa[0];
		my $aa2 = $seqa[-1];
		my $p1  = $pos+length($seq1);
		return "p.$aa1".$pos."_".$aa2.$p1."del";
		
		 
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


1;