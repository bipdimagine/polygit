package GenBoTranscriptCache;
use strict;
use Storable qw(retrieve);
use Moo;
use Bit::Vector;
use Bit::Vector::Overload;
use Data::Dumper;
extends 'GenBoTranscript', 'GenBoCache';



has gene_object => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> undef,
);



##### METHODS #####

sub getPrimers {
	my ($self) = @_;
	return $self->project->getPrimersByObjects($self);
}

sub get_score_impact {
	my ($self, $score) = @_;
	my $v = $self->getChromosome->getVectorScore($self->id."_".$score);
	warn "---->".$self->id." $score" unless defined $v;
	warn $self->getChromosome->getVectorScore($self->id) unless defined $v;
	confess() unless defined $v;
	return $v;
}

sub get_consequences {
 	my ($self,$consequences) = @_;
 	confess() unless $consequences;
 	my $s;
 	foreach my $c (@$consequences){
 		$self->project->getMaskCoding($c);
 		my $v = $self ->getChromosome()->getVectorScore($self->id."_".$c);
 		unless ($s) {
 			$s = $v;
 		}
 		else {
 			$s |= $v;
 		}
 	}
	return $s;
}

has vector => (
	is		=> 'ro',
	#reader	=> 'getUcscName',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		if ($self->getChromosome->isVectorScore){
			return $self->get_score_impact("all");
		}
		else {
		my $gene = $self->getGene;
		my $ids = $gene->variants_tree()->fetch($self->start,$self->end);
		my $vector = $self->getChromosome->getNewVector();
			foreach my $id (@$ids){
			$vector->Bit_On($id); 
			}
		return $vector;
		}
	},
);





sub get_score_ratio {
	my ($self, $patient,$score) = @_;
	return $self->getChromosome->getVectorScore($patient->id."_ratio_".$score);
}

has vector_dm => (
	is		=> 'ro',
	#reader	=> 'getUcscName',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->vector & $self->getChromosome->getVectorScore("dm") ;
	},
);


sub getVectorDM{
	my ($self,$patient) = @_;
	return $self->vector_dm unless $patient;
	my $v = $self->vector_dm;
	$v &= $patient->getVectorOrigin($self->getChromosome);
	return $v;
	return ($self->vector_dm &&  $patient->getVectorOrigin($self->getChromosome));
}

has vector_clinvar => (
	is		=> 'ro',
	#reader	=> 'getUcscName',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->vector & $self->getChromosome->getVectorScore("clinvar_pathogenic") ;
	},
);



sub getVectorClinvar{
	my ($self,$patient) = @_;
	return $self->vector_clinvar unless $patient;
	return ($self->vector_clinvar &  $patient->getVectorOrigin($self->getChromosome));
}

sub getVectorPatients {
	my ($self,$patient) = @_;
	return $self->{hvector}->{$patient->id} if exists  $self->{hvector}->{$patient->id};
	#warn $self->getChromosome()->name." ".$self->get_score_impact("all")->Size." patient size :".$patient->getVectorOrigin($self->getChromosome)->Size." ".$self->getChromosome->getVectorScore("clinvar_pathogenic")->Size()." ==> ".$self->getChromosome->getNewVector()->Size;
	$self->{hvector}->{$patient->id} =  $self->vector && $patient->getVectorOrigin($self->getChromosome);
	#$self->{vector}->{$patient->id} &= $patient->getVectorOrigin($self->getChromosome);
	return  $self->{hvector}->{$patient->id};
}

#TODO: strict-denovo
#TODO: recesiif + compound
sub getVariants {
	my ($self,$filters) = @_;
	my $chr = $self->getChromosome();

	my @lVar;
	my $gene = $self->getGene();
	my $v =  $self->vector();
	return [] if $v->is_empty;
	return $self->getChromosome->getListVarObjects($v) unless $filters;
	my $patient = $filters->{patient};
	die() unless $patient;
	$filters->{score_impact} = "all" unless $filters->{score_impact} ;
	return [] if $v->is_empty;
	$v &= $self->getVectorPatients($patient);
	
	if ($filters->{patient}) {
			if (exists $filters->{score_ratio} && $filters->{score_ratio} >0){
				$v &= $chr->getVectorScore($patient->id."_ratio_".$filters->{score_ratio});
			}
			else {
				$v &= $patient->getVectorOrigin($self->getChromosome);
			}
		
	}


	#warn $v;
	
	
	
	return [] if $v->is_empty;
	if (exists $filters->{score_impact} && $filters->{score_impact} ne "all"){
	 	$v &= $self->get_score_impact($filters->{score_impact});
	}
	return [] if $v->is_empty;
	
	if (exists $filters->{consequences}){
	 	$v &= $self->get_consequences($filters->{consequences});
	}
	return [] if $v->is_empty;
	if (exists $filters->{score_frequency} && $filters->{score_frequency} ne "all"){
	 $v &= $chr->getVectorScore($filters->{score_frequency});
	}
	if (exists $filters->{score_gnomad_ac} && $filters->{score_gnomad_ac} ne "all"){
	 $v &= $chr->getVectorScore($filters->{score_gnomad_ac});
	}
	### GENETICS FAMILIAL MODELS
	if ($filters->{strict_denovo}) {
		$v->Intersection($v, $patient->getFamily->getModelVector_fam_strict_denovo($chr));
	}
	if ($filters->{denovo}) {
		$v->Intersection($v, $patient->getFamily->getModelVector_fam_denovo($chr));
	}
	if ($filters->{dominant}) {
		$v->Intersection($v, $patient->getFamily->getModelVector_fam_dominant($chr));
	}
	if ($filters->{recessif}) {
		$v->Intersection($v, $patient->getFamily->getModelVector_fam_recessif($chr));
	}
	if ($filters->{mosaique}) {
		$v->Intersection($v, $patient->getFamily->getModelVector_fam_mosaique($chr));
	}
	if ($filters->{mosaique}) {
		$v->Intersection($v, $patient->getFamily->getModelVector_fam_mosaique($chr));
	}
	if ($filters->{uniparental_disomy}) {
		$v->Intersection($v, $patient->getFamily->getModelVector_fam_uniparental_disomy($chr));
	}
	
	
	if ($filters->{score_frequency} ){
		my $c = $filters->{score_frequency};
		$v &=  $self->getChromosome->lmdb_score_impact->get($c) if  $self->getChromosome->lmdb_score_impact &&  $self->getChromosome->lmdb_score_impact->get($c) ;
	}
	$self->getVectorDM()->to_Enum;
	$v |= $self->getVectorDM($patient);
	$v |= $self->getVectorClinvar($patient);
	
	
	if ($filters->{vector}) {
		$v->Intersection($v, $filters->{vector});
	}
	
	return $self->getChromosome->getListVarObjects($v);
	

}
 
sub getVariations {
	my $self = shift;
	return $self->getVariants('getVariations');
}

sub getInsertions {
	my $self = shift;
	return $self->getVariants('getInsertions');
}

sub getDeletions {
	my $self = shift;
	return $self->getVariants('getDeletions');
}

sub getLargeDeletions {
	my $self = shift;
	return $self->getVariants('getLargeDeletions');
}

sub getIndels{
	my $self = shift;
	my @lRes;
	push(@lRes,@{$self->getInsertions()});
	push(@lRes,@{$self->getDeletions()});
    return \@lRes;
}

sub getStructuralVariations {
	my $self = shift;
	my @lRes;
	push(@lRes,@{$self->getVariations()});
	push(@lRes,@{$self->getIndels()});
    return \@lRes;
}

1;