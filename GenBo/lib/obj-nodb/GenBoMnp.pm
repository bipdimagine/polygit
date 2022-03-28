package GenBoMnp;
use strict;
use Moose;
use MooseX::Method::Signatures;
use Data::Dumper;
use Config::Std;
use GenBoCapture;
use Position;
extends "GenBoVariant";

has isMnp => (
	is		=> 'ro',
	default	=> 1
);
has isVariation => (
	is		=> 'ro',
	default	=> 1,
);

has kyoto_id => (
	is		=> 'rw',
	#isa		=> 'Str',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $chr = $self->getChromosome();
		my $type;
		 $type = 'mnp'; 
		return join("_", ($chr->name, $self->start(), $type));
		
	},
);

has type_object => (
	is		=> 'ro',
	default	=> "mnps_object",
);


has kyoto_id_2 => (
	is		=> 'rw',
	#isa		=> 'Str',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $chr = $self->getChromosome();
		my $type;
		 $type = 'mnp'; 
		return join("_", ($chr->name, ($self->start()-1), $type));
		
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



sub polyphenScore {
	return 0;
}
sub polyphenStatus {
	return 0;
}
sub siftStatus {
	return 0;
}
sub siftScore {
	return 0;
}


sub annotation_coding {
	my ( $self, $tr, $annot ) = @_;
	my $project = $self->getProject();
	my $prot  = $tr->protein;
	my $gid   = $tr->getGene->id();
	my $trid = $tr->id();

	my $consequence =  $tr->codonsConsequenceForVariations($self);
	$annot->{ $tr->id }->{coding}->{sequences} = $consequence;
	$annot->{ $prot}->{coding}->{sequences} = $consequence;
	if ( $consequence->{aa} eq $consequence->{aa_mut} ) {
		$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("silent");
	}
	elsif ( $consequence->{aa_mut}  =~ /\*/  ) {
		$annot->{$trid}->{mask} = $annot->{$trid}->{mask} |$project->getMaskCoding("stop");
		

	}
	elsif ( $consequence->{aa} =~ /\*/  || $consequence->{orf_position} <= 3 ) {
		
		$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("phase");
	}

	else {
		$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("nonsynonymous");
	}
	
}

sub getNomenclature {
	my ( $self, $transcript, $debug ) = @_;
	confess("none transcript for nomenclature") unless $transcript;
	
	my $id = $transcript->id();
	my $text;
	
	if ($self->isCoding($transcript)) {
		my $cons = $self->annotation()->{ $id }->{coding}->{sequences};
		#warn Dumper $self->annotation()->{$id}->{coding}  unless $cons;
		warn $self->return_mask_testing($transcript,"coding")  unless $cons;
		return "c.".$cons->{orf_position}.$cons->{seq_orf}.">".$cons->{seq_mut} ;
	}
	return;
	#coding 
}


1;