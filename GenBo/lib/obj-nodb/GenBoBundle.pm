
package GenBoBundle;

use Moo;
use Data::Dumper;
use Config::Std;
use FindBin qw($Bin);
use lib "$Bin/";
use lib "$Bin/../GenBoDB";
use Set::IntSpan::Fast::XS;
use Storable qw(store retrieve freeze thaw);
use List::Util qw( shuffle sum max min);
extends "GenBo";



has isBundle => (
	is		=> 'ro',
	default	=> 1
);

has name => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->infos->{name};
	},
);

has description => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->infos->{description};
	},
);

has version => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->infos->{version};
	},
);

has infos => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $query = $self->getProject->{buffer}->queryPanel();
		my $hash = $query->getBundleInfos($self->id);
		return $hash->{$self->id()};
	}
);

has genes_name => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $query = $self->getProject->buffer->queryPanel();
		my %hash;
		map { $hash{$_}++ } @{$query->getBundleGenesByEnsg($self->id())};
	#	warn Dumper %hash;
		#warn Dumper %hash;
		
		return \%hash;
	},
);

has genes_id => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $hash ={};
		foreach my $gene_name (keys %{$self->genes_name()}) {

			next unless $gene_name;
			my $gene = $self->getProject->newGene($gene_name);
#			warn $gene_name unless $gene;
			next unless $gene;
			$hash->{$gene->id()} = $gene->getChromosome->id();

		}
		return $hash;
	},
);
sub setGenes {
	my $self = shift;
	return $self->genes_id();
}


sub setTranscripts {
	my $self = shift;
	my $hRes;
	foreach my $gene (@{$self->getGenes()}) {
		my $hTranscriptsId = $gene->setTranscripts();
		foreach my $id (keys(%$hTranscriptsId)) { $hRes->{$id} = undef; }
	}
	return $hRes;
}



1;
