package GenBo;

use strict;
use Moose;
use Data::Dumper;
use Config::Std;
use Array::IntSpan;
use Storable qw(store retrieve freeze thaw fd_retrieve);
use Digest::MD5 qw(md5 md5_hex md5_base64);




my $hashObjectType = {
	'GenBoChromosome'	=> 'chromosomes',
 	'GenBoCapture'		=> 'captures',
 	'GenBoVariation'	=> 'variations',
 	'GenBoDeletion'		=> 'deletions',
 	'GenBoInsertion'	=> 'insertions',
 	'GenBoLargeDeletion'	=> 'large_deletions',
 	'GenBoInversion'	=> 'inversions',
 	'GenBoBoundary'	=> 'boudaries',
 	'GenBoLargeDuplication'	=> 'large_duplications',
 	'GenBoGene'			=> 'genes',
 	'GenBoReference'	=> 'references',
 	'Position'			=> 'positions',
};

has name => (
	is		=> 'ro',
	required=> 1,
);


has id => (
	is		=> 'ro',
	#isa		=> 'Str',
	required=> 1,
	reader	=> 'id',
);

has md5_id => (
	is		=> 'ro',
	#isa		=> 'Str',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return md5_hex($self->id);;
	},
);

has project => (
	is		=> 'rw',
	reader	=> 'getProject',
);

has type 	=> (
	is		=> 'ro',
	lazy	=> 1,
	reader	=> 'getType',
	default => sub {
		my $self = shift;
		my $ref = ref $self;
		return $hashObjectType->{$ref};
	},
);

has maxProc => (
	is		=> 'ro',
	#isa		=> 'Int',
	reader	=> 'getMaxProc',
	writer	=> 'setMaxProc',
	lazy	=> 1,
	default => 1,
);

has isPatient => (
	is		=> 'ro',
	#isa		=> 'Bool',
	reader	=> 'isPatient',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $type = ref $self;
		if ($type eq 'GenBoPatient') { return 1; }
		return 0;
	},
);
has isCnv => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> undef,
);
has isChromosome => (
	is		=> 'ro',
	#isa		=> 'Bool',
	reader	=> 'isChromosome',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $type = ref $self;
		if ($type eq 'GenBoChromosome') { return 1; }
		return 0;
	},
);
has isCache => (
	is		=> 'ro',
	#isa		=> 'Bool',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $type = ref $self;
		if ($type  =~ /Cache/) { return 1; }
		return 0;
	},
);
has isCapture => (
	is		=> 'ro',
	#isa		=> 'Bool',
	reader	=> 'isCapture',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $type = ref $self;
		if ($type eq 'GenBoCapture') { return 1; }
		return 0;
	},
);

has isGene => (
	is		=> 'ro',
	default	=> 0,
);

has isTranscript => (
	is		=> 'ro',
	default	=> 0,
);

has isProtein => (
	is		=> 'ro',
	default	=> 0,
);

has isVariation => (
	is		=> 'ro',
	default	=> 0,
);

has isMnp => (
	is		=> 'ro',
	default	=> 0
);

has isDeletion => (
	is		=> 'ro',
	default	=> 0
);

has isComplex => (
	is		=> 'ro',
	default	=> 0
);

has isInsertion => (
	is		=> 'ro',
	default	=> 0
);
has isInvertion => (
	is		=> 'ro',
	default	=> 0
);

has isVariant => (
	is		=> 'ro',
	default	=> 0
);

has isLargeDeletion => (
	is		=> 'ro',
	default	=> 0
);

has isLargeDuplication => (
	is		=> 'ro',
	default	=> 0
);

has isPanel => (
	is		=> 'ro',
	default	=> 0
);


###
##### SET/MOOSE OBJECTS #####
###

### GenBoPhenotype
#
has phenotypes_object => (
	is		=> 'ro',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->setPhenotypes();	
	}
);

sub setPhenotypes {
	my $self = shift;
	$self->_errorMethod('setPhenotypes');
}

sub getPhenotypes {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->phenotypes_object(), "phenotypes");
}



### GenBoPanel
#
has panels_object => (
	is		=> 'ro',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->setPanels();	
	}
);

sub setPanels {
	my $self = shift;
	$self->_errorMethod('setPanels');
}

sub getPanels {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->panels_object(), "panels");
}

### GenBoBundle
#
has bundles_object => (
	is		=> 'ro',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->setBundles();	
	}
);

sub setBundles {
	my $self = shift;
	$self->_errorMethod('setBundles');
}

sub getBundles {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->bundles_object(), "bundles");
}

### GenBoRun
#
has runs_object => (
	is		=> 'ro',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->setRuns();	
	}
);

sub setRuns {
	my $self = shift;
	$self->_errorMethod('setRun');
}

sub getRuns {
	my $self = shift;
#	warn $self->id unless $self->getProject();
	return $self->getProject()->myflushobjects($self->runs_object(), "runs");
}



sub getRun {
	my ($self,$name) = @_;
	my $toto = $self->getRuns();
	return $self->_getOneObjectByName($toto,$name);
	
}

### GenBoSomaticGroup
#
has somatic_groups_object => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->setSomaticGroups();	
	}
);

sub setSomaticGroups {
	my $self = shift;
	$self->_errorMethod('setSomaticGroups');
}

sub getSomaticGroups {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->somatic_groups_object(), "somatic_groups");
}

### GenBoFamily
#
has families_object => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->setFamilies();	
	}
);

sub setFamilies {
	my $self = shift;
	$self->_errorMethod('setFamilies');
}

sub getFamilies {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->families_object(), "families");
}

### GenBoPatient
#
has patients_object => (
	is		=> 'ro',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $toto = $self->setPatients();	
		return $toto;	
	}
);

sub setPatients {
	my $self = shift;
	$self->_errorMethod('setPatients');
}

sub getPatients {
	my $self = shift;
	my @lPat;
	foreach my $patient (@{$self->getProject()->myflushobjects($self->patients_object(), "patients")}) {
		next if $patient->is_control();
		push (@lPat, $patient);
	}
	return \@lPat;
}

### GenBoChromosome
#
has chromosomes_object => (
	is		=> 'ro',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->setChromosomes();	
	}
);

sub setChromosomes {
	my $self = shift;
	$self->_errorMethod('setChromosomes');
}

sub getChromosomes {
	my $self = shift;
	confess($self->name) unless $self->getProject();
	#warn $self->project->name;
	#my @toto = grep {$_->karyotypeId =~ /KI/}   @{$self->getProject()->myflushobjects($self->chromosomes_object(), "chromosomes")};
	#confess($self->name) if scalar(@toto) ;
	my @chrs = sort {$a->karyotypeId <=> $b->karyotypeId} @{$self->getProject()->myflushobjects($self->chromosomes_object(), "chromosomes")};

	#my @chrs = sort @{$self->getProject()->myflushobjects($self->chromosomes_object(), "chromosomes")};
	return \@chrs;
}

### GenBoVariation
#

has variations_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hRes = $self->setVariations();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;
	}
);



sub setVariations {
	my $self = shift;
	$self->setVariants('variations');
	return $self->{variations_object} ;
}


sub getVariations {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->variations_object(), "variations");
}


########################
### StructuralVariants
#######################

#SVDELetions 

has isSVDeletion => (
	is		=> 'ro',
	default	=> 0
);


has svdeletions_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		confess();
		my $hRes = $self->setSVDeletions();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;
	}
);

sub getSVDeletions {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->svdeletions_object(), "svdeletions");
}

#SV Duplications
has isSVDuplication => (
	is		=> 'ro',
	default	=> 0
);
sub setSVDuplications {
	my $self = shift;
	confess();
}

has svduplications_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		confess();
		my $hRes = $self->setSVDuplications();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;
	}
);

sub getSVDuplications {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->svduplications_object(), "svduplications");
}


########################
### Inversions
#######################

sub setInversions {
	my $self = shift;
	$self->setVariants('inversions');
	return $self->{inversions_object} ;
}

has inversions_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hRes = $self->setInversions();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;
	}
);

sub getInversions {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->inversions_object(), "inversions");
}

########################
### Boudaries
#######################

sub setBoundaries {
	my $self = shift;
	$self->setVariants('boundaries');
	return $self->{inversions_object} ;
}

has boundaries_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hRes = $self->setBoundaries();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;
	}
);

sub getBoundaries {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->boundaries_object(), "boundaries");
}



########################
### LargeDeletions
#######################

sub setLargeDeletions {
	my $self = shift;
	$self->setVariants('large_deletions');
	return $self->{large_deletions_object} ;
}

has large_deletions_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hRes = $self->setLargeDeletions();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;
	}
);

sub getLargeDeletions {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->large_deletions_object(), "large_deletions");
}

########################
### LargeDuplications
#######################

sub setLargeDuplications {
	my $self = shift;
	$self->setVariants('large_duplications');
	return $self->{large_duplications_object} ;
}

has large_duplications_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hRes = $self->setLargeDuplications();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;
	}
);

sub getLargeDuplications {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->large_duplications_object(), "large_duplications");
}

########################
### LargeInsertions
#######################

sub setLargeInsertions {
	my $self = shift;
	$self->setVariants('large_insertions');
	return $self->{large_insertions_object} ;
}

has large_insertions_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hRes = $self->setLargeInsertions();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;
	}
);

sub getLargeInsertions {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->large_insertions_object(), "large_insertions");
}

########################
### Mnps
#######################


sub setMnps {
	my $self = shift;
	$self->setVariants('mnps');
	return $self->{mnps_object} ;
}


sub getMnps {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->mnps_object(), "mnps");
}


has mnps_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hRes = $self->setMnps();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;
	}
);
### GenBoComplex
#
has complex_object => (
	is		=> 'ro',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hRes = $self->setComplex();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;
	}
);



sub setComplex {
	my $self = shift;
	$self->setVariants('complex');
	return $self->{complex_object} ;
}

sub getComplex {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->complex_object(), "complex");
}

### GenBoInsertions
#
has insertions_object => (
	is		=> 'ro',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hRes = $self->setInsertions();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;
	}
);



sub setInsertions {
	my $self = shift;
	$self->setVariants('insertions');
	return $self->{insertions_object} ;
}


sub getInsertions {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->insertions_object(), "insertions");
}



### GenBoDeletions
#
has deletions_object => (
	is		=> 'ro',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hRes = $self->setDeletions();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;
	}
);
sub setDeletions {
	my $self = shift;
	$self->setVariants('deletions');
	return $self->{deletions_object} ;
}
sub getDeletions {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->deletions_object(), "deletions");
}

## GenBoDeletions
#

sub setVariants {
	my $self = shift;
	$self->_errorMethod('setVariants');
}


sub getIndels{
	my $self  =shift;
	my @lRes;
		push(@lRes,@{$self->getInsertions()});
		push(@lRes,@{$self->getDeletions()});
	    return \@lRes;
}

sub getCnvs{
	my $self = shift;
	my @lRes;
	push(@lRes,@{$self->getLargeDuplications()});
	push(@lRes,@{$self->getLargeInsertions()});
	push(@lRes,@{$self->getLargeDeletions()});
    return \@lRes;
}
sub getSV{
	my $self = shift;
	my @lRes;
	push(@lRes,@{$self->getLargeDuplications()});
	push(@lRes,@{$self->getLargeInsertions()});
	push(@lRes,@{$self->getLargeDeletions()});
	push(@lRes,@{$self->getInversions()});
	push(@lRes,@{$self->getMeis()});
    return \@lRes;
}
#has indels => (
#	is		=> 'ro',
#	#isa		=> 'HashRef',
#	reader	=> 'getIndels',
#	lazy	=> 1,
#	default	=> sub {
#		my $self = shift;
#		my @lRes;
#		push(@lRes,@{$self->getInsertions()});
#		push(@lRes,@{$self->getDeletions()});
#	    return \@lRes;
#	}
#);
sub getStructuralVariations {
		my $self  =shift;
		my @lRes;
		push(@lRes,@{$self->getVariations()});
		push(@lRes,@{$self->getIndels()});
		push(@lRes,@{$self->getCnvs()});	
		#push(@lRes,@{$self->getMnps()});	
	#	die() if @{$self->getMnps()};
	    return \@lRes;
}

#has structural_variations => (
#	is		=> 'ro',
#	#isa		=> 'HashRef',
#	reader	=> 'getStructuralVariations',
#	lazy	=> 1,
#	default	=> sub {
#		my $self = shift;
#		my @lRes;
#		push(@lRes,@{$self->getVariations()});
#		push(@lRes,@{$self->getIndels()});
#	#	push(@lRes,@{$self->getMnps()});	
#	#	die() if @{$self->getMnps()};
#	    return \@lRes;
#	}
#);



### GenBoCapture
#
has captures_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->setCaptures();	
	}
);

sub setCaptures {
	my $self = shift;
	$self->_errorMethod('setCaptures');
}

sub getCaptures {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->captures_object(), "captures");
}

sub getCapture {
	my ($self, $name) = @_;
	return $self->_getOneObjectByName($self->getCaptures(),$name);
}
### Position
#
has positions_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->setPositions();	
	}
);

sub setPositions {
	my ($self, $objRef) = @_;
	$self->_errorMethod('setPositions');
}

sub getPosition {
	my ($self, $objRef) = @_;
	return $self->getProject()->myflushobjects($self->positions_object($objRef), "positions");
}

### GenBoReference
#
has references_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->setReferences();	
	}
);

sub setReferences {
	my $self = shift;
	my $hRefIds;
	foreach my $chr (@{$self->getProject()->getChromosomes}) {
		foreach my $ref (@{$chr->getReferences()}) {
			$hRefIds->{$ref->id()} = undef;
		}
	}
	return $hRefIds;
}

sub getReferences {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->references_object(), "references");
}
### GenBoRegulatoryRegion
#
has isRegulatoryRegion => (
	is		=> 'ro',
	default	=> 0,
);

has regulatory_region_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		
		my $hRes = $self->setRegulatoryRegions();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;	
	}
);

sub setRegulatoryRegions {
	my $self = shift;
	$self->_errorMethod('setGenes');
}

sub getRegulatoryRegions {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->regulatory_region_object(), "regulatory_regions");
}

sub getRegulatoryRegion {
	my ($self, $name) = @_;
	return $self->_getOneObjectByName($self->getRegulatoryRegions(),$name);
}

### GenBoGene
#
has genes_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hRes = $self->setGenes();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;	
	}
);

sub setGenes {
	my $self = shift;
	$self->_errorMethod('setGenes');
}

sub getGenes {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->genes_object(), "genes");
}

### GenBoTranscript
#
has transcripts_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hRes = $self->setTranscripts();
		
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;
	}
);

sub setTranscripts {
	my $self = shift;
	$self->_errorMethod('setTranscripts');
}

sub getTranscripts {
	my $self = shift;
	my %hash;
	#foreach my $id (keys %{$self->transcripts_object()}) {
	#	my $count = () = $id =~ /\Q_/g;
	#	delete $self->transcripts_object()->{$id} if $count >1;
	#}
	
	return $self->getProject()->myflushobjects($self->transcripts_object(), "transcripts");
}



### GenBoExons
#
has exons_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hRes = $self->setExons();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;
	}
);

sub setExons {
	my $self = shift;
	$self->_errorMethod('setExons');
}

sub getExons {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->exons_object(), "exons");
}

### GenBoIntrons
#
has introns_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hRes = $self->setIntrons();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;
	}
);

sub setIntrons {
	my $self = shift;
	$self->_errorMethod('setIntrons');
}

sub getIntrons {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->introns_object(), "introns");
}

### GenBoPrimers
#
has primers_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hRes = $self->setPrimers();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;
	}
);

sub setPrimers {
	my $self = shift;
	die();
	$self->_errorMethod('setPrimers');
}

sub getPrimers {
	my $self = shift;
	return $self->getProject()->myflushobjects($self->primers_object(), "primers");
}

sub getPrimer {
	my ($self, $primer_name,$notest) = @_;


	my $p = $self->_getOneObjectByName($self->getPrimers(), $primer_name,1);
	
	return $p if $p;
	
	confess("can't find $primer_name") unless $notest;
	#return $self->_getOneObjectByName($self->getPatients(), $patientIdName,$notest);
}

sub getOrderedPrimers {
		my $self = shift;
		return $self->{ordered_primers} if exists $self->{ordered_primers};
		
		my $foo = [];
		my @objs =  sort {$a->start <=>$b->start or $a->end <=>$b->end} @{$self->getProject()->myflushobjects($self->primers_object(), "primers")};
		foreach my $o (@objs) {
			push(@$foo,{start=>$o->start,end=>$o->end,id=>$o->id});
		#	$foo->set_range($o->start, $o->end, $o->id);
		}
		$self->{array_ordered} = $foo;

		 $self->{ordered_primers} = \@objs;
		 return   $self->{ordered_primers} ;
		
}


sub find_primers {
	my ($self,$start,$end) = @_;
#	confess();
	return [] unless @{$self->{array_ordered}};
	my $array = $self->{array_ordered};
	my $id= $self->_find_pos($start);
	my @primers_id;
	$id -- if $id >0;
	for (my $i=$id;$i<@$array;$i++){
		my $p = $array->[$i];
		next if $p->{end} <= $start;
		last if ($p->{start} > $end);
		push(@primers_id,$p->{id});
	}
	return \@primers_id;
}
sub _find_pos {
  my $self = shift;
  my $val  = shift;
  my $low  = shift || 0;
my $array = $self->{array_ordered};
  my $high = scalar( @$array ) -1;

  while ( $low < $high ) {
  	
    my $mid = int( ( $low + $high ) / 2 );
   # warn $val." :: $mid :: $high :: ".$array->[$mid]->{start};
    if ( $val < $array->[$mid]->{end}) {
      $high = $mid;
    }
    elsif ( $val >  $array->[$mid]->{end}) {
      $low = $mid + 1;
    }
    else {
      return $mid;
    }
  }

  return $low;
}

### GenBoProtein
#
has proteins_object => (
	is		=> 'rw',
	#isa		=> 'HashRef',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hRes = $self->setProteins();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;
	}
);

sub setProteins {
	my $self = shift;
	$self->_errorMethod('setTranscripts');
}

sub getProteins {
	my $self = shift;

	my %hash;
	foreach my $id (keys %{$self->proteins_object()}) {
		my $count = () = $id =~ /\Q_/g;
		delete $self->proteins_object()->{$id} if $count >1;
	}
		my $toto = $self->proteins_object();

	return $self->getProject()->myflushobjects($self->proteins_object(), "proteins");
}

sub getProtein {
	my ($self, $id) = @_;
	my $prots = $self->getProteins();
	if (scalar @$prots > 1){
		confess ("more than one gene") unless $id;
		return $self->getObject('proteins', $id);
	}
	
	return $prots->[0];
}

###
##### METHODS #####
###



sub getObjects {
	my ($self, $type) = @_;
	my $refObject = ref $self;
	my @lObjects;
	my $project = $self->getProject();
	if (exists ($project->{objects}->{$type})) {
		 @lObjects = values(%{$project->{objects}->{$type}}); 
	}
	return \@lObjects
}

sub getObject {
	my ($self, $type, $id) = @_;
	my $project = $self->getProject();
	if (exists $project->{objects}->{$type}->{$id}) { return $project->{objects}->{$type}->{$id}; }
	confess("\n\nERROR: No $type found with id $id !!\n");
}

sub _errorMethod {
	my ($self, $methodName) = @_;
	my $refObject = ref $self;
	confess("\n\nERROR: No reason to call '$methodName' with '$refObject' object. Exit.\n\n");
}

sub _getOneObjectByName {
	my ($self, $objs, $name,$notest) = @_;
	if (scalar @$objs == 1 && !(defined $name)) { return $objs->[0]; }
	confess($self->name) unless $name;
	my (@find) = grep {$_->name() eq $name} @$objs;
	return $find[0] if scalar(@find) == 1 ;
	return if scalar(@find) == 0 && $notest ;
	confess("no object name -$name- found in ".join(";",map{$_->name} @$objs)." ==> ".$self->id);
}

sub getPanel {
	my ($self, $name) = @_;
	return $self->_getOneObjectByName($self->getPanels(), $name);
}

sub getBundle {
	my ($self, $name) = @_;
	return $self->_getOneObjectByName($self->getBundles(), $name);
}

sub getFamily {
	my ($self, $name) = @_;
	return $self->_getOneObjectByName($self->getFamilies(), $name);
}

sub getSomaticGroup {
	my ($self, $name) = @_;
	return $self->_getOneObjectByName($self->getSomaticGroups(), $name);
}

sub getPatient {
	my ($self, $patientIdName,$notest) = @_;
	$self->getPatients();
	$self->patients_object;
	unless ($patientIdName) { $self->_errorMethod('getPatient'); }

	my $p = $self->_getOneObjectByName($self->getPatients(), $patientIdName,1);
	
	return $p if $p;
	if (exists($self->getProject()->{objects}->{patients}->{$patientIdName})) { 
		confess();
	#	confess("problem") if  $self->getProject()->{objects}->{patients}->{$patientIdName}->name ne $patientIdName ; 
		return $self->getProject()->{objects}->{patients}->{$patientIdName};
		 }
	confess("can't find $patientIdName") unless $notest;
	return;
	#return $self->_getOneObjectByName($self->getPatients(), $patientIdName,$notest);
}
sub getExon {
	my ($self, $name) = @_;
	return $self->_getOneObjectByName($self->getExons(),$name);
}
sub getChromosome {
	my ($self, $name,$nodie) = @_;
	if ($name){
	$name =~ s/chr//;
	$name = "MT" if $name eq "M";
	}
	return $self->_getOneObjectByName($self->getChromosomes(),$name,$nodie);
}

sub getReference {
	my ($self, $name) = @_;
	my $toto = $self->getReferences();

	return $self->_getOneObjectByName($self->getReferences(),$name);
}

sub getVariation {
	my ($self, $name) = @_;
	return $self->_getOneObjectByName($self->getVariations(),$name);
}

sub getInsertion {
	my ($self, $name) = @_;
	return $self->_getOneObjectByName($self->getInsertions(),$name);
}

sub getDeletion {
	my ($self, $name) = @_;
	return $self->_getOneObjectByName($self->getDeletions(),$name);
}

sub getIndel {
	my ($self, $id) = @_;
	confess();
	if (exists $self->{objects}->{'insertions'}->{$id}) { return $self->{objects}->{'insertions'}->{$id}; }
	elsif (exists $self->{objects}->{'deletions'}->{$id}) { return $self->{objects}->{'deletions'}->{$id}; }
	confess("\n\nERROR: No indel found in deletions and insertions objects !!\n")
}

sub getIntSpan {
	my $self = shift;
	my $refObject = ref $self;
	confess("\n\nERROR: no IntSpan for " . $refObject . " object !!\n\n");
}

sub getGene {
	my ($self, $name) = @_;
	return $self->_getOneObjectByName($self->getGenes(),$name);
}

sub getTranscript {
	my ($self, $name) = @_;
	return $self->_getOneObjectByName($self->getTranscripts(),$name);
}

sub disconnect {
		my $self = shift;
		$self->buffer->disconnect();
}

sub buffer {
	my $self = shift;
	confess() unless $self->getProject();
	return $self->getProject()->buffer();
}
sub DESTROY {
	my $self = shift;
}

1;