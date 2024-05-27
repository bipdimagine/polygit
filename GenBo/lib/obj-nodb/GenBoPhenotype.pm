package GenBoPhenotype;

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



has isPanel => (
	is		=> 'ro',
	default	=> 1
);

has name => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		confess() unless $self->infos;
		return $self->infos->{name};
	}
);

has short_name => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->infos->{short_name};
	}
);

has concept => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->infos->{concept};
	}
);


has infos => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $query = $self->buffer->queryPhenotype();
		my $hash = $query->getPhenotypeInfos($self->id);
		return $hash->{$self->id()};
	}
);

has keywords => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $query = $self->buffer->queryPhenotype();
		return $query->getKeywords($self->id);
	}
);
has projects_id => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $query = $self->buffer->queryPhenotype();
		return $query->getProjectsId($self->id);
	}
);
has projects_name => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $ids = $self->projects_id;
		my $query = $self->buffer->getQuery();
		my $names ;
		foreach my $id (@$ids){
			push(@$names,$query->getProjectNameFromId($id));
		}
		
		return $names;
	}
);
sub setPanels {
	my $self = shift;
	my %hids;
	my $query = $self->buffer->queryPhenotype();
	return $query->getPanelsId($self->id);

}
has nb_panels =>(
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return scalar(@{$self->getPanels});
	}
);


has statistic_genes =>(
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
			my $hids = {};
			#
		foreach my $panel (@{$self->getPanels}){
		my $genes =  $panel->genes_name();
		next unless $genes;
			foreach my $g (keys %{$genes}){
				$hids->{$g} ++;
				}
		}
	return $hids;
	}
);

sub setGenes {
	my $self = shift;
	my $hids;
	foreach my $panel (@{$self->getPanels}){
		my $genes =  $panel->getGenes();
		foreach my $g (@{$genes}){
			$self->{statistics_genes}->{$g->id} ++;
			#$hids->{$g->id} ++ ;
		}
	}
	return $hids;
}


1;