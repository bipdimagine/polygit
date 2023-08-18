package GenBoVariation;
use strict;
use Moo;

use Data::Dumper;
use Config::Std;
use GenBoCapture;
use Position;
extends "GenBoVariant";


has isVariation => (
	is		=> 'ro',
	default	=> 1,
);


#has spliceAI => (
#is              => 'rw',
#lazy =>1,
#default => sub {
#        my $self = shift;
#        my $v =   $self->getChromosome->score_spliceAI($self->start,$self->alternate_allele);
#        return "-" unless defined $v;
#        return $v;
#
#
#}
#);

has type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "substitution";
	}
	);


has type_public_db => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "snps";
	},
);

has type_object => (
	is		=> 'ro',
	default	=> "variations_object",
);

has kyoto_id => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $chr = $self->getChromosome();
		my $type;
		 $type = $self->var_allele(); 
		return join("_", ($chr->name, $self->start(), $type));
		
	},
);

has rocksdb_id => (
	is		=> 'ro',
	lazy=>1,
	default => sub {
	my ($self) = @_;
	my ($chr,$pos,$ref,$alt) = split("-",$self->gnomad_id);
	 $pos  = sprintf("%010d", $pos);
	return  ($pos."!".$alt);
	},
	
);

#has cadd_score => (
#	is		=> 'rw',
#	lazy	=> 1,
#	default => sub {
#		my $self = shift;
#		my $sc =  $self->getChromosome()->get_lmdb_score("cadd",$self);
#		return "-" unless $sc;
#		return $sc;
#		
#	},
#);

has ncboost_score => (
	is              => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return if ($self->getChromosome->id() eq 'MT');
		my $h =  $self->getChromosome()->get_lmdb_score("ncboost",$self);
		
		return unless ($h);
		return $h->{'b'};
	}    
);

has ncboost_category => (
	is              => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return unless ($self->ncboost_score());
		return if ($self->ncboost_score() eq '-');
		
		my $find_cat;
		my $ncboost = $self->ncboost_score() * 100;
		foreach my $sort_cat (sort keys %{$self->project->buffer->config->{scaled_priority_ncboost}}) {
			my $cat = $self->project->buffer->config->{scaled_priority_ncboost}->{$sort_cat};
			my $ncboost_score = $self->project->buffer->config->{scaled_score_ncboost}->{$cat} * 100;
			if ($ncboost >= $ncboost_score) {
				$find_cat = $cat;
			}
		}
		return $find_cat;
	},
);

has revel_score => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $sc =  $self->getChromosome->get_lmdb_score("revel",$self);
		return "-" unless $sc;
		return $sc;
		
	},
);
has dbscsnv_ada => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $chr = $self->getChromosome();
		my $sc =  $chr->get_lmdb_score("dbscsnv_ada",$self);
		return "-" unless $sc;
		return $sc;
		
	},
);
has dbscsnv_rf => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $chr = $self->getChromosome();
		my $sc =  $chr->get_lmdb_score("dbscsnv_rf",$self);
		return "-" unless $sc;
		return $sc;
		
	},
);



sub getPredictions {
	my ($self) = @_;
	$self->{predictions} = undef;
	$self->{predictions}->{all}->{mask} =0;
	my $transcripts = $self->getTranscripts();
	my $max_score;
	my $max_pred;
	my %methods;
	foreach my $tr (@$transcripts){
		next unless $self->isCoding($tr);
		warn $tr->id unless $tr->getProtein();
		#my $protid = $tr->getProtein()->id;
		my ($protid, $chr_id) = split('_', $tr->getProtein()->id());
		my $gene_id = $tr->getGene()->id;
		$self->{predictions}->{$protid}->{mask} = 0;
		if ($self->isSilent($tr) ) {
			$self->{predictions}->{$protid}->{mask} = $self->getProject->getMaskPrediction("completely_benign");	 	
		}
		elsif ($self->isStop($tr) || $self->isPhase($tr) || $self->isEssentialSplicing($tr) || $self->isMatureMiRNA($tr) || $self->isFrameshift($tr) ) {
			$self->{predictions}->{$protid}->{mask} = $self->getProject()->getMaskPrediction("completely_damaging");
			$self->{predictions}->{$gene_id}->{mask} = $self->getProject()->getMaskPrediction("completely_damaging");
			$self->{predictions}->{all}->{mask} = $self->getProject()->getMaskPrediction("completely_damaging");
		}
		else {
			my $pos = $self->position($tr->getProtein)->start();
			my $aa  = $self->changeAA($tr->getProtein);
			#warn "\t".$protid." $pos ".$aa."\n";
			$self->{predictions}->{$tr->getProtein()->id()} = $self->getProject()->getPredictions($self->getChromosome(),$protid,$pos,$aa);
			$self->{predictions}->{$protid} = $self->getProject()->getPredictions($self->getChromosome(),$protid,$pos,$aa);
			$self->{predictions}->{$tr->id}=$self->{predictions}->{$protid};
			#warn Dumper $self->{predictions}->{$protid};
			 foreach my $method (keys %{$self->{predictions}->{$protid}}){
				next if $method eq "mask";
			#	if ($self->{predictions}->{$protid}->{$method}->{pred} > $self->{predictions}->{all}->{$method}->{pred}){
					$self->{predictions}->{$tr->getGene->id}->{$method}->{score} = $self->{predictions}->{$protid}->{$method}->{score} ;
					$self->{predictions}->{$tr->getGene->id}->{$method}->{pred} = $self->{predictions}->{$protid}->{$method}->{pred};
			#	}
			 }
		}
		$self->{predictions}->{$tr->id}->{mask} = $self->{predictions}->{$protid}->{mask};
		$self->{predictions}->{$gene_id}->{mask} =  0 unless exists $self->{predictions}->{$gene_id}->{mask};
		$self->{predictions}->{$gene_id}->{mask} = $self->{predictions}->{$gene_id}->{mask}  | $self->{predictions}->{$protid}->{mask};
		$self->{predictions}->{all}->{mask} =0 unless exists $self->{predictions}->{all}->{mask};
		$self->{predictions}->{all}->{mask} = $self->{predictions}->{all}->{mask} | $self->{predictions}->{$protid}->{mask};
	}
}
sub getPredictionMask {
	my ($self,$obj) =  @_;
	$self->getPredictions() unless exists $self->{predictions};
	
	my $id = "all";
	if ($obj){
		die($obj->name ." is not a protein or a gene") unless $obj->isProtein() || $obj->isGene() || $obj->isTranscript;
		$id = $obj->id();
	}
	return $self->{predictions}->{$id}->{mask};
}


sub polyphenScore {
	my ( $self, $obj ) = @_;
	$self->getPredictions() unless exists $self->{predictions};
	die($obj->name) unless $obj->isProtein()||  $obj->isTranscript();
	$self->{predictions}->{$obj->id}->{polyphen_humvar}->{score} = "-" unless $self->{predictions}->{$obj->id}->{polyphen_humvar}->{score};
	return $self->{predictions}->{$obj->id}->{polyphen_humvar}->{score};
}
sub polyphenScore2 {
	my ( $self, $obj ) = @_;
	$self->getPredictions() unless exists $self->{predictions};
	
	my $chr_name = $self->getChromosome()->name();
	my $h = $self->project->liteScore->get($chr_name, $self->start);

	return unless $h;
	my $sc = $h->{$self->sequence}->{$obj->name}->{Polyphen2_HDIV_score};
	warn Dumper $h if $sc; 
	return $sc;
	#warn Dumper $h;
	#die($obj->name) unless $obj->isProtein()||  $obj->isTranscript();
	#$self->{predictions}->{$obj->id}->{polyphen_humvar}->{score} = "-" unless $self->{predictions}->{$obj->id}->{polyphen_humvar}->{score};
	#return $self->{predictions}->{$obj->id}->{polyphen_humvar}->{score};
}

sub siftScore {
	my ( $self, $obj ) = @_;
	$self->getPredictions() unless exists $self->{predictions};
	die() unless $obj->isProtein() ||  $obj->isTranscript();
	return '-' unless ($self->{predictions}->{$obj->id}->{sift});
	return '-' unless (exists $self->{predictions}->{$obj->id}->{sift}->{score});
 	return 0 if (exists $self->{predictions}->{$obj->id}->{sift}->{score} and $self->{predictions}->{$obj->id}->{sift}->{score} == 0);
	return $self->{predictions}->{$obj->id}->{sift}->{score}+0;
}

sub polyphenStatusText {
	my ( $self, $obj ) = @_;
	if ($self->polyphenStatus($obj) == 3) { return "Probably damaging"; }
	if ($self->polyphenStatus($obj) == 2) { return "Possibly damaging"; }
	if ($self->polyphenStatus($obj) == 1) { return "Benign damaging"; }
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
	my $mask = $self->getPredictionMask($obj);
	$mask = 0 unless $mask;
	return 3 if $mask & $self->getProject()->getMaskPrediction('polyphen_probably damaging');
	return 2 if $mask & $self->getProject()->getMaskPrediction('polyphen_possibly damaging');
	return 1 if $mask & $self->getProject()->getMaskPrediction( 'polyphen_benign');
	return 0 ;
}

sub siftStatusText {
	my ( $self, $obj ) = @_;
	if ($self->siftStatus($obj) == 2) { return "Deleterious"; }
	if ($self->siftStatus($obj) == 1) { return "Tolerated"; }
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
	my $mask = $self->getPredictionMask($obj);
	return 0 unless $mask;
	#'unknown' => 0,
    #    'tolerated'     => 1,
    #   'deleterious'   => 2,
	#return 0 if $mask & $self->buffer()->getMaskPrediction('sift_unknown');
	return 2 if $mask & $self->getProject->getMaskPrediction( 'sift_deleterious');
	return 1 if $mask & $self->getProject()->getMaskPrediction('sift_tolerated');
	
	return 0;
	
}




###
# divers cosequences method
#####
sub consequence {
	my ( $self, $prot ) = @_;
	confess() unless $prot;
	
	return "" unless  exists $self->annotation->{$prot->id};
	my $hseq = $self-> annotation->{ $prot->id }->{coding}->{sequences};
	return (  $hseq->{aa}. '/' . $hseq->{aa_mut});
}



##### METHODS #####
=head2 getproteinPosition
	Title   : getproteinAA
 	
=cut


sub annotation_coding {
	my ( $self, $tr, $annot ) = @_;
	my $project = $self->getProject();
	 if ($tr->{protein} =~ /ENST/){
	 	warn $self->id;
		warn $tr->name();
	 }
	return if $tr->{protein} =~ /ENST/;
	
	
	#return unless $tr->getProtein;
	my $prot  = $tr->getProtein->id;
	my $gid   = $tr->getGene->id();
	my $trid = $tr->id();

	my $consequence =  $tr->codonsConsequenceForVariations($self);
	$annot->{ $tr->id }->{coding}->{sequences} = $consequence;
	$annot->{ $prot}->{coding}->{sequences} = $consequence;
	if ( $consequence->{aa} eq $consequence->{aa_mut} ) {
		$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("silent");
	}
	elsif ( $consequence->{aa_mut}  eq "*"  ) {
		$annot->{$trid}->{mask} = $annot->{$trid}->{mask} |$project->getMaskCoding("stop");
		

	}
	elsif ( $consequence->{aa} eq "*" || $consequence->{orf_position} <= 3 ) {
		
		$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("phase");
	}

	else {
		$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("nonsynonymous");
	}
}



sub constructNomenclature {
	my ( $self, $transcript, $debug ) = @_;
	confess("none transcript for nomenclature") unless $transcript;
	my $id = $transcript->id();
	my $text;
		#my $cons = $self->annotation()->{ $id }->{coding}->{sequences};
	if ($self->isCoding($transcript)) {
		#return;
		my $cons = $self->annotation()->{ $id }->{coding}->{sequences};
		warn $self->return_mask_testing($transcript,"coding")  unless $cons;
		return "c.".$cons->{orf_position}.$cons->{seq_orf}.">".$cons->{seq_mut} ;
	}
	elsif ($self->isUtr($transcript)){
		#return;
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
			return "n.".$pos_transcript."".$seqr.">".$seqv ;#
		#return n
	}
			#return $self->getNomenclatureForUtr($transcript);
		my ($dist,$pos) = $self->getPositionForNomenclatureIntronic($transcript);
		
	#	warn "cuicui";
		 my $st_dist = $dist;
		 $st_dist="+".$dist if $dist >0;
		# warn $self->{$self->{var_allele};};
		 my $seqv = $self->{var_allele};#"";#$self->{sequence};
		#warn $self->{ref_allele};
		 my $seqr = $self->{ref_allele};;
		  #return;
			if ($transcript->strand() == -1 ){
				 $seqv = BioTools::complement($seqv);
				  $seqr = BioTools::complement($seqr);
			}
		 return"c.$pos$st_dist".$seqr.">".$seqv;
	}
	return;
	#coding 
}




1;