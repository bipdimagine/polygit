package GenBoGenomic;

use strict;
use Vcf;
use Moo;
use Data::Dumper;
 use GenBoCoverageSamtools;
 use GenBoCoverageTabix;
 use GenBoCoverageLmdb;
 use Carp;
 #use GenBoCoverageLmdb;
 use Storable qw(store retrieve freeze thaw);
  use Compress::Snappy;

extends "GenBo";

has kyotoId => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { my $self =shift; $self->{genbo_id} },
#	required=> 1,
);
	

has chromosomes_object => (
	is		=> 'ro',
	#required=> 1,
);

has chromosome_object => (
	is		=> 'rw',
	lazy    => 1,
	default => undef,
);

has chromosome_name =>(
	is		=> 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $chr_name;
		foreach my $chr (keys %{$self->chromosomes_object()}) {
			$chr_name = $chr;
		}
		return $chr_name;
	}
	#required=> 1,
	
);
has start => (
	is		=> 'rw',
	#required=> 1,
);

has end => (
	is		=> 'rw',
	#required=> 1,
);

has strand => (
	is		=> 'rw',
	default => 1,
);

has length => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $length = $self->end() - $self->start() + 1;
		return $length;
	},
);


has intspan => (
	is		=> 'rw',
	reader	=> 'getGenomicSpan',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $intSpan = Set::IntSpan::Fast::XS->new();
    	$intSpan->add_range($self->start(), $self->end());
    	return $intSpan;
	},
);

has span_position_genomic => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $intSpan = Set::IntSpan::Fast::XS->new();
    	$intSpan->add_range($self->start(), $self->end());
    	return $intSpan;
	},
);

has genesId => (
	is		=> 'rw',
	lazy	=> 1,
	reader	=> 'getGenesId',
	default => sub {
		my $self = shift;
		confess();
		my $padding = 200;
		$padding = 0 if $self->getChromosome->name eq "MT";  
		my $t = $self->getChromosome()->genesIntervalTree->fetch($self->start-$padding,$self->end+$padding+1);
		return $t;
	},
);

has isMT => (
	is		=> 'ro',
	lazy		=> 1,
	default => sub {
		my $self = shift;
		my $mt;
		 $mt=1 if $self->getChromosome->name eq 'MT'; 
		return $mt;
	}
);

has _transcriptsId => (
	is		=> 'rw',
	lazy	=> 1,
	#reader	=> 'getTranscriptsId',
	default => sub {
		my $self = shift;
	
		return $self->getChromosome()->getFastTranscriptsIds($self);
	},
);




has annot_only_transcritpts=> (
	is		=> 'rw',
	lazy    => 1,
	default => undef,
);


###### SET OBJECTS #####


sub test {
	my ($self,$tt)=@_;
	confess($tt);
}

sub inside {
	my ($self,$obj) = @_;
	my $set = $self->span_position_genomic->intersection($obj->span_position_genomic);
	return scalar($set->as_array);
	# return !($set->is_empty());
}

sub setReferences {
	my ($self) = @_;
	my $chr = $self->getChromosome();
	my $nb_ref;
	my $return_hash;
	my @refs;
	my $max_coverage = -1;
	my $aspan = $chr->findReferences($self);
	confess() unless scalar(@$aspan);
	foreach my $span (@$aspan){
		$return_hash->{$span->[2]} = undef;
	}
	return $return_hash;
}

sub objectsInside{
	my ($self,$objs,$hash) = @_;
	my $h;
	foreach my $o (@$objs){
		$h->{$o->id} = undef if $self->inside($o);
	}
	return $h;
}

sub setVariations {
	my $self = shift;
	my $hash;
	foreach my $ref (@{$self->getReferences()}){
		$hash =  $self->objectsInside($ref->getVariations());
	}
	return $hash;
}
sub setMnps {
	my $self = shift;
	my $hash;
	foreach my $ref (@{$self->getReferences()}){
		$hash =  $self->objectsInside($ref->getMnps());
	}
	return $hash;
}

sub setLargeDeletions {
	my $self = shift;
	my $hash;
	foreach my $ref (@{$self->getReferences()}){
		$hash =  $self->objectsInside($ref->getLargeDeletions());
	}
	return $hash;
}
sub setInversions {
	my $self = shift;
	my $hash;
	foreach my $ref (@{$self->getReferences()}){
		$hash =  $self->objectsInside($ref->getInversions());
	}
	return $hash;
}
sub setBoundaries {
	my $self = shift;
	my $hash;
	foreach my $ref (@{$self->getReferences()}){
		$hash =  $self->objectsInside($ref->getBoundaries());
	}
	return $hash;
}
sub setLargeDuplications {
	my $self = shift;
	my $hash;
	foreach my $ref (@{$self->getReferences()}){
		$hash =  $self->objectsInside($ref->getLargeDuplications());
	}
	return $hash;
}
sub setLargeInsertions {
	my $self = shift;
	my $hash;
	
	foreach my $ref (@{$self->getReferences()}){
		
		$hash =  $self->objectsInside($ref->getLargeInsertions());
	}
	return $hash;
}
sub setDeletions {
	my $self = shift;
	my $hash;
	foreach my $ref (@{$self->getReferences()}){
		 $hash = $self->objectsInside($ref->getDeletions());
	}
	return $hash;
}

sub setInsertions {
	my $self = shift;
	my $hash;
	foreach my $ref (@{$self->getReferences()}){
		$hash = $self->objectsInside($ref->getInsertions());
	}
	return $hash;
}

sub setGenes {
	my $self = shift;
#	return [keys %{$self->vannot_chr->{genes}}];

	my $hGenesid;
	my $upstream = $self->buffer->config->{definitions}->{upstream};
	my $downstream = $self->buffer->config->{definitions}->{downstream};
	my $t = $self->getChromosome()->genesIntervalTree->fetch($self->start,$self->end+1);
	$self->vannot({}) unless $t;
	return $t;
}



has vannot => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $h ={};
		foreach my $g (@{$self->getGenes}){
			my $hash = $g->getTanscriptsAnnotations($self->start,$self->end);
			$h = {%$h,%$hash};
		}
		
		return $h;
	},
);

sub setTranscripts {
	my $self = shift;
	my $hreturnTranscriptsId ={};#= $self->_transcriptsId();
	return $self->vannot;
	#return [keys %{$self->vannot}];
	
}





sub setProteins {
	my $self = shift;
	die("yeah it's a problem !!! you have to check that man");
	my $t = $self->getChromosome->getFastProteinsIds($self);
	
	my $hreturns = {};#= $self->_transcriptsId();
	
	foreach my $tr (@{$self->getTranscripts()}){
		my $prot = $tr->getProtein();
		next unless $prot;
		$hreturns->{$prot->id} = undef if $self->inside($prot);
	}
	return $hreturns;
}


sub isInside {
	my ($self,$obj) = @_;
	return undef if $self->getChromosome()->name ne $obj->getChromosome()->name ;
	my $span = $self->getGenomicSpan()->intersection($obj->getGenomicSpan());
	return !($span->isEmpty);
}

sub position {
	my ($self,$obj) = @_;
	confess() unless $obj;
	return $self->{position}->{$obj->id} if exists $self->{position}->{$obj->id};
	if ($obj->isChromosome){
		my $start = $self->start();
		my $end = $self->end();
		my $strand = $self->strand();
		$self->{position}->{$obj->id} = new Position({start=>$start,end=>$end,strand=>1});
		return $self->{position}->{$obj->id};
	}
	else {
		confess();
	}
}


sub setPrimers {
	my ($self) = @_;
	my $chr = $self->getChromosome();
	#my $primers = $chr->getPrimers();
	my $tree = $chr->tree_primers();
	my $results = $tree->fetch($self->start,$self->end);
	my $hreturns;
	foreach my $r (@$results){
		$hreturns->{$r} =undef
	}
	return $hreturns;
	

	}


### #################
# methods for coverage
####################

has coverage_start => (
is		=> 'rw',
	lazy => 1,
	default => sub {
		my $self =shift;
		my $x = $self->start - 100;
		return 0 if $x <0;
		return $x;
		
	},
);

has coverage_end => (
is		=> 'rw',
	lazy => 1,
	default => sub {
		my $self =shift;
		return ($self->end + 100);
		
	},
);



sub get_coverage {
	my ($self,$p) = @_;
	return  $self->{coverage_obj}->{$p->id} if exists $self->{coverage_obj}->{$p->id};
	$self->_set_coverage($p) unless  exists $self->{coverage_obj}->{$p->id};
	
	return $self->{coverage_obj}->{$p->id};
}


sub set_coverage {
	my ($self,$p,$v) = @_;
	$self->{coverage_obj}->{$p} = $v;
		return;
}

sub _set_coverage {
	my ($self,$p) = @_;
	# $self->set_coverage($p->id, $self->return_raw_coverage_obj($p));
	 #return;
	  if ($p->isNoSqlDepth){
	 	 my $gc =  GenBoCoverageLmdb->new(chromosome=>$self->getChromosome, patient=>$p, start=>$self->coverage_start, end=>$self->coverage_end);
	 	 $self->set_coverage($p->id, $gc);
	 }
	 my $no =  $self->project->noSqlCoverage();
	 my $gc1 = $no->get($p->name,$self->id);
	 
	if ($gc1){
		$gc1->{patient} = $p;
		$gc1->chromosome($self->getChromosome);
		 $self->set_coverage($p->id, $gc1);
		
		 return;
	}
	else {
		    my $gc = $self->return_raw_coverage_obj($p);
			#die($self->name);
		    #$gc->coverage($self->start,$self->end);
	 		$self->set_coverage($p->id, $gc);
	}
}

sub return_raw_coverage_obj{
	my ($self,$p) = @_;
	my $gc;
	 if ($p->tabix_coverage){
		 $gc =  GenBoCoverageTabix->new(chromosome=>$self->getChromosome, patient=>$p, start=>$self->coverage_start, end=>$self->coverage_end);
	 }
	 else {
	 	$gc = GenBoCoverageSamtools->new(chromosome=>$self->getChromosome, patient=>$p, start=>$self->coverage_start, end=>$self->coverage_end);
	 }
	return $gc;
	
}

sub getChromosome {
	my ($self, $name) = @_;
	return $self->chromosome_object() if ($self->chromosome_object());
	if ($name){
		$name =~ s/chr//;
		$name = "MT" if $name eq "M";
	}
	return $self->_getOneObjectByName($self->getChromosomes(),$name);
}

1;