package GenBoVariantCache;
use strict;
use Moo;
use Carp;
use Data::Dumper;
extends 'GenBoCache';


# objet GenBoChromosomeCache
has chromosome => (
	is  => 'ro',
);

# id in vector
has vector_id => (
	is	=> 'ro',
);



sub dejaVuInfosForDiag2 {
	my ($self,$key) = @_;
	if (exists $self->{dejaVuInfosForDiag2}){
		my $hash = delete $self->{dejaVuInfosForDiag2};
		$self->{array_dejavu} = $self->buffer->hash_to_array_dejavu($hash);
	}
	if (exists $self->{cad}){
		my $a = delete $self->{cad}; #aka compress dejavu for diag
		my @aa = unpack("w*",$a);
		$self->{array_dejavu} = \@aa;
	}
	return $self->{array_dejavu} unless $key;
	my $index = $self->buffer->index_dejavu($key);
	return  $self->{array_dejavu}->[$index];
}

has global_vector_id => (
	is	=> 'ro',
	lazy => 1,
	default => sub {
		my $self = shift;
		return $self->getChromosome->name()."!".$self->vector_id;
	}
);


sub isHomozygote {
	my ($self, $patObj) = @_;
	return $patObj->getVectorOriginHo($self->getChromosome)->contains($self->vector_id);
}
sub isHeterozygote {
	my ($self, $patObj) = @_;
	return $patObj->getVectorOriginHe($self->getChromosome)->contains($self->vector_id);
}

# quick hash annotation for global annotation (from vectors cache)
has vector_infos => (
	is => 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		my $hAnnot;
		my $project = $self->project();
		foreach my $gene (@{$self->getGenes()}) {
			foreach my $annot (keys %{$project->ensembl_annotations()}) {
				next unless (exists $gene->categories->{$annot});
				if ($gene->categories->{$annot}->contains($self->vector_id())) {
					my $text_annot = $annot;
					$text_annot = 'splice_site' if ($annot eq 'splicing');
					$text_annot = 'nonsynonymous' if ($annot eq 'coding');
					$text_annot = 'nonframeshift' if ($annot eq 'non-frameshift');
					$hAnnot->{annot}->{$text_annot} = undef;
				}
			}
		}
		foreach my $annot (keys %{$hAnnot->{annot}}) {
			$hAnnot->{annot}->{all}->{mask} += $self->project->getMaskCoding($annot);
		}
		confess("\n\nERROR: impossible to annot ".$self->id()." (vector id: ".$self->vector_id()."). Die.\n\n") unless ($hAnnot);
		return $hAnnot;
	},	
);

# hash annotation for all genes / transcripts / proteins of this variant



##### METHODS #####



#sub return_mask_testing {
#	my ($self, $obj, $type) = @_;
#	my $id = "all";
#	return unless ($obj);
#	unless ($obj) {
#		return $self->vector_infos->{annot}->{all}->{mask} & $self->project()->getMaskCoding($type);
#	}
#	confess("\n\nERROR: call only on transcript. Die.\n\n") unless $obj->isTranscript();
#	$id = $obj->id ;
#	unless (exists $self->annotation()->{$id}->{mask}) {
#		return;
##		confess();
##		warn "\n";
##		warn 'Self: '.ref($self);
##		warn 'ID: '.$id;
##		warn 'Type: '.$type;
#	}
#	return $self->annotation()->{$id}->{mask} & $self->getProject()->getMaskCoding($type);
#}



sub setPatients {
	my $self = shift;
	confess();
	my $h;
	my $chr = $self->project->getChromosome($self->chromosome_name());
	my $v_id = $self->vector_id();
	my $type = 'substitution';
	$type = 'insertion' if ($self->isInsertion());
	$type = 'deletion' if ($self->isDeletion());
	foreach my $pat (@{$chr->getPatients()}) {
		$pat->getVariantsVector($chr);
		$h->{$pat->id()} = undef if ($pat->global_categories->{$chr->id()}->{$type}->contains($v_id));
	}
	return $h;
}

 


sub isCompoundTransmission {
	my ($self,$fam,$child,$gene,$father_mother) = @_;
	confess("\n\ndon't forget the familly man\n\n") unless $fam;
	confess("\n\ndon't forget the child man\n\n") unless $child;
	confess("\n\ndon't forget the gene man\n\n") unless $gene;
	confess("\n\ndon't forget the father_mother dude\n\n") unless $father_mother;
	my $proj_name = $fam->getProject->name();
	my $v_parent = $gene->getVariantsVector()->Clone();
	$v_parent->Intersection($v_parent, $child->getVariantsVector($gene->getChromosome()));
	if ($father_mother eq 'father' and $fam->getFather()) {
		$v_parent->Intersection($v_parent, $fam->getFather->getVariantsVector($gene->getChromosome()));
		$v_parent->AndNot($v_parent, $fam->getMother->getVariantsVector($gene->getChromosome()));
	}
	if ($father_mother eq 'mother' and $fam->getMother()) {
		$v_parent->Intersection($v_parent, $fam->getMother->getVariantsVector($gene->getChromosome()));
		$v_parent->AndNot($v_parent, $fam->getFather->getVariantsVector($gene->getChromosome()));
	}
	return 0 if ($v_parent->is_empty());
	return 1;
}

sub isUniparentalDisomyTransmission {
	my ($self,$fam,$child) = @_;
	my $vector = $fam->getVector_individual_uniparental_disomy($self->getChromosome,$child);
	return $vector->contains($self->vector_id);
#	#warn '-';
#	return 0  unless 	$self->getPourcentAllele($child) ;
#	#warn '-';
#	return 0  if 	$self->getPourcentAllele($child) eq "-" ;
#	#warn '-';
#	#return 0 if 	$self->getPourcentAllele($child) > 75;
#	#warn '-';
#	#die();
#	die("don't forget the familly dude") unless $fam;
#	return 0 if ($child->getVectorHe($self->getChromosome())->contains($self->vector_id));
#	#die();
#	if ($self->isFatherTransmission($fam, $child)) {
#		return 1 if ($fam->getFather->getVectorHe($self->getChromosome())->contains($self->vector_id));
#	}
#	elsif ($self->isMotherTransmission($fam, $child)) {
#		return 1 if ($fam->getMother->getVectorHe($self->getChromosome())->contains($self->vector_id));
#	} 
#	return 0;
}




sub isBothTransmission{
		my ($self,$fam,$child,$debug) = @_;
	my $vector = $fam->getVectorBothTransmission($self->getChromosome,$child);
	#my $vector2 = $fam->getVectorRecessiveTransmission($self->getChromosome,$child);
	#$vector -= $vector2;
	#warn "\t\t\t vector : ".$self->vector_id." ".$vector->contains($self->vector_id) if $debug;
	#warn "\t\t\t vector2 : ".$self->vector_id." ".$vector2->contains($self->vector_id) if $debug;
	#die() if $debug;
	return $vector->contains($self->vector_id);
	
}

sub isDominantTransmission {
	my ($self,$fam,$child) = @_;

	return undef unless $fam->isDominant();
	die("don't forget the familly dude".' '.$fam) unless $fam;
	my $vector;
	if($child) {
		$vector = $fam->getVector_individual_dominant($self->getChromosome,$child);
	}
	else {
		$vector = $fam->getVector_family_dominant($self->getChromosome);
	}
	eval { return $vector->contains($self->vector_id); };
	
	eval { return $vector->contains($self->vector_id); };
	 if ($@){
	 	 warn $@;
	 	 warn ref $vector;
	 	 warn $vector;
	 	confess();
	 }
	return $vector->contains($self->vector_id);
}

sub isRecessiveTransmission {
	my ($self,$fam,$child) = @_;
	die("don't forget the familly dude".' '.$fam) unless $fam;
	my $vector;
	if($child) {
		$vector = $fam->getVector_individual_recessive($self->getChromosome,$child);
	}
	else {
		$vector = $fam->getVector_family_recessive($self->getChromosome);
	}
	eval { return $vector->contains($self->vector_id); };
	
	eval { return $vector->contains($self->vector_id); };
	 if ($@){
	 	 warn $@;
	 	 warn ref $vector;
	 	 warn $vector->Size;
	 	confess();
	 }
	return $vector->contains($self->vector_id);
}

sub isFatherTransmission {
	my ($self,$fam,$child) = @_;
	my $vector;
	if($child) {
		$vector = $fam->getVector_individual_father($self->getChromosome,$child);
	}
	else {
		$vector = $fam->getVector_family_father($self->getChromosome);
	}
	eval { return $vector->contains($self->vector_id); };
	 if ($@){
	 	 warn $@;
	 	 warn ref $vector;
	 	 warn $vector;
	 	confess();
	 }
	return $vector->contains($self->vector_id);
		

}

sub isMotherTransmission {
	my ($self,$fam,$child,$debug) = @_;
	my $vector;
	if($child) {
		warn "coucou   ::: ".$self->vector_id if $debug;
		
		$vector = $fam->getVector_individual_mother($self->getChromosome,$child,$debug);
	}
	else {
		$vector = $fam->getVector_family_mother($self->getChromosome);
	}
	eval { return $vector->contains($self->vector_id); };
	 if ($@){
	 	 warn $@;
	 	 warn ref $vector;
	 	 warn $vector;
	 	confess();
	 }
	return $vector->contains($self->vector_id);
	

}


sub isDenovoTransmission {
	my ($self,$fam,$child) = @_;
	#my $vector = $fam->getModelVector_fam_denovo($self->getChromosome);
	my $vector;
	if($child) {
		$vector = $fam->getVector_individual_denovo($self->getChromosome,$child);
	}
	else {
		$vector = $fam->getVector_family_denovo($self->getChromosome);
	}
	eval { return $vector->contains($self->vector_id); };
	 if ($@){
	 	 warn $@;
	 	 warn ref $vector;
	 	 warn $vector;
	 	confess();
	 }
	return $vector->contains($self->vector_id);
}


sub isStrictDenovoTransmission{
	my ($self,$fam,$child) = @_;
	my $vector;
	if($child) {
		$vector = $fam->getVector_individual_strict_denovo($self->getChromosome,$child);
	}
	else {
		$vector = $fam->getVector_family_strict_denovo($self->getChromosome);
	}
	eval { return $vector->contains($self->vector_id); };
	 if ($@){
	 	 warn $@;
	 	 warn ref $vector;
	 	 warn $vector;
	 	confess();
	 }
	return $vector->contains($self->vector_id);
}


sub isMosaicTransmission {
	my ($self,$fam,$child) = @_;
	warn $self->id;
#	my $vector = $fam->getModelVector_fam_mosaique($self->getChromosome);
#	return $vector->contains($self->vector_id);
	my $var_id = $self->id();
	my $chr = $self->getChromosome();
	return unless ($child->isIll());
	return unless ($child->getHe($chr)->contains($self->vector_id()));
	warn $child->getHe($chr)->contains($self->vector_id()) ;
	warn $child->getHo($chr)->contains($self->vector_id()) ;
	warn $self->getSequencingGenotype($child);
	return if $self->isBothTransmission($fam,$child); 
	my @lOk;
	my $pc = $self->getPourcentAllele($child);
	return if $pc eq "-";
	my $sample_var_1 = $self->text_heho($child);
	warn $sample_var_1;
	unless ($sample_var_1) {
			warn "\n\nERROR: pb with $var_id and patient ".$child->name()."... Die\n\n";
			confess ();
		}
	foreach my $parent (@{$fam->getParents()}) {
		#if ($chr->ploidy()
		next unless ($parent->getHe($chr)->contains($self->vector_id()));
		my $pcp = $self->getPourcentAllele($parent);
		warn $pcp;
		next if $pcp > 20;
		return if $pcp eq  "-";
		next if ($pc - $pcp ) < 20;
		#->{he_ho_details};
		
		my $sample_var_2 = $self->text_heho($parent);#->{he_ho_details};
		unless ($sample_var_2) {
			warn "\n\nERROR: pb with $var_id (parent: ".$parent->name().")... Die\n\n".$self->vector_id()."\n";
			confess ($self->id);
		}
		
#		next if ($pc - $pcp ) < 20;
		my ($p_fisher, $res_fisher) = $chr->project->buffer->test_fisher($sample_var_2, $sample_var_1, $self->project->mosaique_min_p());
		if ($res_fisher == 1) {
			push (@lOk, 'father') if ($parent->sex() == 1);
			push (@lOk, 'mother') if ($parent->sex() == 2);
		}
	}
	return join(' and ', @lOk) if (@lOk);
	return;
}
sub sequencing_details {
	my ($self,$patient) = @_;
	
	confess();
}
#
sub annex {
	confess();
}
#sub getSequencingInfos {
#        my ($self,$patient) = @_;
#        my $pid = $patient->id;
#        return $self->{seq}->{$pid}->{text} if exists  $self->{seq}->{$pid}->{text};
#        warn $pid;
#        warn Dumper($self->{seq}->{$pid});
#       die();
#}
#
#

#sub getRatio {
#        my ($self,$patient,$method) = @_;
#        confess() if $method;
#        my $pid = $patient->id;
#       
#        return $self->{seq}->{$pid}->{ratio} if exists  $self->{seq}->{$pid}->{ratio};
#         warn Dumper $self->{seq};
#        confess($patient->name." ".$patient->id);
#
#        
#}

sub getTransmissionModelType {
		my ($self,$fam,$child,$gene) = @_;
		return "-" unless (exists $self->{patients_object}->{$child->id});
		#my $is_mosaic = $self->isMosaicTransmission($fam,$child);
		if($child->isMother ){
			return "-m" unless exists $self->{patients_object}->{$child->id};
			return "+m"; 
		}
		if ($child->isFather ){
			return "-f" unless exists $self->{patients_object}->{$child->id};
			return "+f"; 
		}
		if ($self->isRecessiveTransmission($fam,$child) == 1 ){
			return "recessive";
		}
		elsif ($self->isUniparentalDisomyTransmission($fam,$child) == 1  && $self->getPourcentAllele($child) > 75 ){
			
			return "Uniparental Disomy";	
		}
		elsif  ($self->isMosaicTransmission($fam,$child)  ){
			return "Mosaic";
		}
		elsif ($self->isBothTransmission($fam,$child) == 1  ){
			return "Both"
		}
		elsif ($self->isFatherTransmission($fam,$child) ==1){
			if ($gene and $self->isCompoundTransmission($fam,$child,$gene,'mother') == 1){
				return "Father_c";
			}
			return "Father";
		}
		elsif ($self->isMotherTransmission($fam,$child) == 1){
			if ($gene and $self->isCompoundTransmission($fam,$child,$gene,'father') == 1){
				return "Mother_c";
			}
			return "Mother";
		}
		elsif  ($self->isStrictDenovoTransmission($fam,$child) == 1  ){
			return "strict_denovo";
		}
		elsif  ($self->isDenovoTransmission($fam,$child) == 1  ){
			
			if ($fam->father && $fam->mother){
				return "denovo";
			}
			else {
				return "denovo/?";
			}
			
		}
		
		else {
			return "?";
		}
		confess("code it by your self")
}


sub getTransmissionModel {
		my ($self,$fam,$child,$gene,$debug) = @_;
		#my $is_mosaic = $self->isMosaicTransmission($fam,$child);
		if ($self->isRecessiveTransmission($fam,$child) == 1 ){
			return "recessive";
		}
		elsif ($self->isUniparentalDisomyTransmission($fam,$child) == 1  && $self->getPourcentAllele($child) > 75 ){
			
			return "Uniparental Disomy";	
		}
		elsif  ($self->isMosaicTransmission($fam,$child)  ){
			return "Mosaic";
		}
		elsif ($self->isBothTransmission($fam,$child,$debug) == 1  ){
			return "Both"
		}
		elsif ($self->isFatherTransmission($fam,$child) ==1){
			if ($gene and $self->isCompoundTransmission($fam,$child,$gene,'mother') == 1){
				return "Father_c";
			}
			return "Father";
		}
		elsif ($self->isMotherTransmission($fam,$child) == 1){
			if ($gene and $self->isCompoundTransmission($fam,$child,$gene,'father') == 1){
				return "Mother_c";
			}
			return "Mother";
		}
		elsif  ($self->isStrictDenovoTransmission($fam,$child) == 1  ){
			return "strict_denovo";
		}
		elsif  ($self->isDenovoTransmission($fam,$child) == 1  ){
			
			if ($fam->father && $fam->mother){
				return "denovo";
			}
			else {
				return "denovo/?";
			}
			
		}
		else {
			return "?";
		}
		confess("code it by your self")
}

# return list genes not filtred by interface filters
sub setFiltredGenes {
	my $self = shift;
	my ($hGenesIds, @lGenes);
	foreach my $id (@{$self->getGenesId()}) {
		$hGenesIds->{$id} = undef;
	}
	my $hGenesId;
	my $chr = $self->project->getChromosome($self->chromosome_name());
	foreach my $gene (@{$chr->getGenes()}) {
		if (exists $hGenesIds->{$gene->ensg().'_'.$chr->id()}) {
			push(@lGenes, $gene);
			$hGenesId->{$gene->id()} = undef;
		}
	}
	return $hGenesId;
}
sub is_model {
	my ($self, $ind_fam_som, $model, $patient) = @_;
	my $chr = $self->project->getChromosome($self->chromosome_name());
	if ($patient and ref($patient) eq 'GenBoPatientCache') {
		return 1 if ($chr->models_categories->{$ind_fam_som}->{'saved_model'}->{$model}->{$patient->name()}->contains($self->vector_id()));
	}
	else {
		return 1 if ($chr->models_categories->{$ind_fam_som}->{$model}->contains($self->vector_id()));
	}
	return;
}

sub is_individual_recessif {
	my ($self, $patient) = @_;
	return $self->is_model('individual', 'recessif', $patient);
}

sub is_individual_compound {
	my ($self, $patient) = @_;
	return $self->is_model('individual', 'compound', $patient);
}

sub is_familial_denovo {
	my ($self, $patient) = @_;
	return $self->is_model('familial', 'denovo', $patient);
}

sub is_familial_strict_denovo {
	my ($self, $patient) = @_;
	return $self->is_model('familial', 'strict_denovo', $patient);
}

sub is_familial_dominant {
	my ($self, $patient) = @_;
	return $self->is_model('familial', 'dominant', $patient);
}

sub is_familial_mosaic {
	my ($self, $patient) = @_;
	return $self->is_model('familial', 'mosaic', $patient);
}

sub is_familial_recessif {
	my ($self, $patient) = @_;
	return $self->is_model('familial', 'recessif', $patient);
}

sub is_familial_compound {
	my ($self, $patient) = @_;
	return $self->is_model('familial', 'compound', $patient);
}

sub is_familial_recessif_compound {
	my ($self, $patient) = @_;
	return $self->is_model('familial', 'recessif_compound', $patient);
}

sub is_somatic_loh {
	my ($self, $patient) = @_;
	return $self->is_model('somatic', 'loh', $patient);
}

sub is_somatic_dbl_evt {
	my ($self, $patient) = @_;
	return $self->is_model('somatic', 'dbl_evt', $patient);
}

1;