package GenBoPanel;

use Moose;
use Data::Dumper;
use Config::Std;
use FindBin qw($Bin);
use lib "$Bin/";
use lib "$Bin/../GenBoDB";
use Set::IntSpan::Fast::XS;
use Storable qw(store retrieve freeze thaw);
use List::Util qw( shuffle sum max min);
extends "GenBo";



has isPanel => (
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

has infos => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $query = $self->getProject->{buffer}->queryPanel();
		my $hash = $query->getPanelInfos($self->id);
		return $hash->{$self->id()};
	}
);

has genes_id => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my %h;
		foreach my $bundle (@{$self->getBundles}) {
			foreach my $g ( keys %{$bundle->genes_id()}){
				$h{$g}++;
				
			}
			
		}
		return \%h; 
	},
);

has genes_name => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my %h;
		foreach my $bundle (@{$self->getBundles}) {
			map { $h{$_}++ } keys %{$bundle->genes_name()};
		}
		return \%h; 
	},
);

sub setBundles {
	my $self = shift;
	my $query = $self->getProject->{buffer}->queryPanel();
	my %hIds;
	map { $hIds{$_}++ } @{$query->getBundlesIdsFromPanelId($self->id())};
	return \%hIds;
}

sub setPhenotypes {
	my $self = shift;
	my %hids;
	my $query = $self->buffer->queryPhenotype();
	
	return $query->getPhenotypesIdFromPanelId($self->id);

}


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

has intspan => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $hint;
		foreach my $chr (@{$self->project->getChromosomes}){
			$hint->{$chr->name} = Set::IntSpan::Fast::XS->new();# unless exists ;
		}
		foreach my $g (@{$self->getGenes}){
			my $chr= $g->getChromosome();
			$hint->{$chr->name}->add_range($g->start,$g->end);
		}
		return $hint;
	},
);

sub getIntspan {
	my ($self,$chr) = @_;
	return $self->intspan()->{$chr->name};
}

has main_transcripts => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my @trs;
		my $t = time; 
		foreach my $gene (@{$self->getGenes}){
			push(@trs, map{$_->id} @{$gene->getMainTranscripts})
		}
		return \@trs; 
	},
);
sub getMainTranscripts {
	my ($self) = @_;
	return $self->project->myflushobjects($self->main_transcripts,"transcripts");
}

1;
