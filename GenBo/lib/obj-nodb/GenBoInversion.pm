package GenBoInversion;

use strict;
use Moo;
use Carp;
use Data::Dumper;

extends 'GenBoVariant';


has isInversion => (
	is		=> 'ro',
	default	=> 1,
);

has isSrPr => (
	is		=> 'rw',
	default	=> undef,
);

has name => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->getChromosome->id().'-'.$self->start().'-'.$self->nomenclatureType.'-'.$self->length;
	},
);

has gnomad_id => (
	is		=> 'ro',
	lazy=> 1,
	default=> sub {
	my $self = shift;
	return $self->name;
	}
	);


has alamut_id => (
	is		=> 'ro',
	lazy	=> 1,
	default=> sub {
		my $self = shift;
		my $seq = $self->sequence();
		my $len = length($seq) - 1;
		my $build = 'GRCh37';
		$build = 'GRCh38' if ($self->project->getVersion() =~ /HG38/);
		my $alamut_id = 'chr'.$self->getChromosome->id().':'.$build.':'.$self->start().'_'.($self->start() + $len).'inv'.$seq;
		return $alamut_id;
	},
);

has rocksdb_id => (
	is		=> 'ro',
	lazy=>1,
	default => sub {
		my $self = shift;
		my $pos  = sprintf("%010d", $self->start());
		return $pos."!"."@".$self->length;
	},
	
);

has dejavu => (
	is		=> 'ro',
	lazy => 1,
	default => sub {
		my $self = shift;
		my $dvlite = $self->getProject->dejavuDuckSV->get_dejavu($self->getChromosome->id(), $self->start(), $self->getChromosome->id(), $self->end(), 90);
		return $dvlite;
	}
);

sub other_projects  {
	my ($self) = @_;
	return $self->dejavu->{"nb_projects"} if $self->dejavu->{"nb_projects"};
	return 0;
}

sub other_projects_ho  {
	my ($self) = @_;
	return 0;
}

sub other_patients{
	my $self = shift;
	return $self->dejavu->{"nb_patients"} if $self->dejavu->{"nb_patients"};
	return 0;
}

sub other_patients_ho {
	my $self = shift;
	return 0;
}
	
	
sub nomenclatureType {
	my ($self) = @_;
	return "inv";
}

has isLargeDeletion => (
	is		=> 'ro',
	default	=> 1,
);


has type_object => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "inversions_object";
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
has type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		return "inversion";
	},
);

has sv_type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		return "INV";
	},
);



has structural_type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "inv";
	},
);


has alleles => (
	is		=> 'ro',
	lazy	=> 1,
	default=> sub {
		my $self = shift;
		return 'inv '.$self->length();
	},
);




sub protein_nomenclature {
	my ( $self, $prot ) = @_;
		confess() unless $prot->isProtein();
		my $pos = $self->getProteinPosition($prot);
		return "p.inv".$self->lengthe;
}

sub annotation_coding {
	my ( $self, $tr, $annot,$span ) = @_;
	my $protid  = $tr->protein();
	my $gid   = $tr->getGene->id();
	my $trid = $tr->id;
	my $project = $self->getProject();
	my $l = length($self->sequence()) % 3;
	my $start = $self->position($self->getChromosome())->start;
	my $consequence =  $tr->codonsConsequenceForInversion($self);
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



sub constructNomenclature {
	my ( $self, $transcript, $debug ) = @_;
	confess("none transcript for nomenclature") unless $transcript;
	#return;
	#$self->annotation() unless exists $self->{annotation};
#	my $id = $transcript->id();
#	my $text;
#	if ($self->isCoding($transcript)) {
#		my $cons = $self->{annotation}->{ $id }->{coding}->{sequences};
#		return "c.".($cons->{orf_position}-1)."_INV";
#		
#	}
#	elsif ($self->isUtr($transcript)){
#		#return;
#		#warn $transcript->name();
#		return "UTR-inv-".$self->mei_type;
#	}
#	else {
#		return "ins-".$self->mei_type;
#	}	
	return "INVERSION";
	#coding 
}

##### METHODS #####

sub polyphenStatus {
	my ( $self, $obj ) = @_;
	return 0;
}

sub polyphenStatusText {
	my ( $self, $obj ) = @_;
	return;
}

sub siftStatus {
	my ( $self, $obj ) = @_;
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


1;



