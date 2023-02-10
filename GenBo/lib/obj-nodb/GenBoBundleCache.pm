package GenBoBundleCache;

use Moose;
use Data::Dumper;
use Config::Std;
use FindBin qw($Bin);
use lib "$Bin/";
use lib "$Bin/../GenBoDB";
use Set::IntSpan::Fast::XS;
use Storable qw(store retrieve freeze thaw);
use List::Util qw( shuffle sum max min);
extends "GenBoBundle";




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
		foreach my $gene (@{$self->getGenes($chr)}) {
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



1;