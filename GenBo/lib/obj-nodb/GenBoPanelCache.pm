package GenBoPanelCache;

use Moose;
use Data::Dumper;
use Config::Std;
use FindBin qw($Bin);
use lib "$Bin/";
use lib "$Bin/../GenBoDB";
use Set::IntSpan::Fast::XS;
use Storable qw(store retrieve freeze thaw);
use List::Util qw( shuffle sum max min);
extends "GenBoPanel";




has hash_variants_vector => (
	is      => 'rw',
	lazy    => 1,
	default => sub { {} },
);

sub getVariantsVector {
	my ($self, $chr) = @_;
	confess("\n\nERROR: need GenBoChromomsomeCache argument in GenBoBundleCache::getVariants. Die.\n\n") unless($chr);
	unless (exists $self->hash_variants_vector->{$chr->id()}) { 
		my $vector = $chr->getNewVector();
		foreach my $gene_name (keys %{$self->genes_name()}) {
			my $gene = $self->getProject->newGene($gene_name);
			next if ($gene->getChromosome->id() ne $chr->id());
			my $v_region = $chr->getVariantsVector_from_coord($gene->start(), $gene->end());
			if ($chr->countThisVariants($v_region) > 1) {
				my @lIds = @{$chr->getListVarVectorIds($v_region)};
			}
			$vector += $v_region;
		}
		$self->hash_variants_vector->{$chr->id()} = $vector;
	}
	return $self->hash_variants_vector->{$chr->id()};
}

has vector => (
	is      => 'rw',
	lazy    => 1,
	default => sub { 
		my $self = shift;
		my $hv;
		foreach my $chr (@{$self->project->getChromosomes}){
			$hv->{$chr->name} = $chr->getNewVector();
		}
		foreach my $tr (@{$self->getMainTranscripts}){
			my $chr = $tr->getChromosome();
			$hv->{$chr->name} |= $tr->vector;
		}
		return $hv;
		 },
);
sub getVectorDM {
	my ($self, $chr) = @_;
	confess("\n\nERROR: need GenBoChromomsomeCache argument in GenBoBundleCache::getVariants. Die.\n\n") unless($chr);
	my $v = $self->getVector($chr)->Clone();
	$v &= $chr->vectorDM();
	return $v;
}
sub getVectorClinvarPathogenic {
	my ($self, $chr) = @_;
	confess("\n\nERROR: need GenBoChromomsomeCache argument in GenBoBundleCache::getVariants. Die.\n\n") unless($chr);
	return $self->getVector($chr) & $chr->vectorClinvarPathogenic();
}
sub getVector {
	my ($self, $chr) = @_;
	confess("\n\nERROR: need GenBoChromomsomeCache argument in GenBoBundleCache::getVariants. Die.\n\n") unless($chr);
	
	return $self->vector->{$chr->name};	
}


1;