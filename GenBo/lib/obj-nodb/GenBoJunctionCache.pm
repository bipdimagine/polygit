package GenBoJunctionCache;
use strict;
use Moo;
use Data::Dumper;
extends  'GenBoJunction','GenBoVariantCache';



sub setPatients {
	my $self = shift;
	my $h;
	foreach my $patient (@{$self->getProject->getPatients()}) {
		$h->{$patient->id} = undef if exists $self->annex->{$patient->name()};
	}
	return $h;
}


sub get_hash_noise {
	my ($self, $patient) = @_;
	return $self->{h_noise}->{$patient->name()} if (exists $self->{h_noise}->{$patient->name()});
	my ($h_noise, @lJunctions, $my_id_in_list);
	my $h_exons_introns = $self->get_hash_exons_introns();
	my $intspan_self = $self->getStrictGenomicSpan();
	my @lVectorIds = sort {$a <=> $b} @{$self->getChromosome->getListVarVectorIds($patient->getVectorOrigin($self->getChromosome()))};
	my $j = 0;
	foreach my $v_id (@lVectorIds) {
		if ($self->vector_id() eq $v_id) {
			$my_id_in_list = $j;
			last;
		}
		$j++;
	}
	unless (defined($my_id_in_list)) {
		warn join(', ', @lVectorIds);
		warn $my_id_in_list;
		confess();
	}
	my $min = $my_id_in_list - 15;
	$min = 0 if ($min < 0);
	my $max = $my_id_in_list + 15;
	$max = scalar(@lVectorIds) - 1 if ($max > scalar(@lVectorIds) - 1);
	my $i = $min;
	
	while ($i <= $max) {
		my $vid = $lVectorIds[$i];
		if ($i == $my_id_in_list) {
			$i++;
			next;
		}
		my $junction2 = $self->getChromosome->getVarObject($vid);
		my $h_exons_introns_2 = $junction2->get_hash_exons_introns;
		if ($h_exons_introns_2) {
			foreach my $tid (keys %{$h_exons_introns_2}) {
				next if $tid eq 'by_pos';
				next unless (exists $h_exons_introns->{$tid});
				my @lPos = (sort keys %{$junction2->get_hash_exons_introns->{$tid}->{by_pos}});
				next if (@lPos == 0);
				foreach my $pos (keys %{$junction2->get_hash_exons_introns->{$tid}->{by_pos}}) {
					if (exists $h_exons_introns->{$tid}->{by_pos}->{$pos}) {
						my $jid2 = $junction2->id();
						$jid2 =~ s/_RI//;
						$jid2 =~ s/_SE//;
						$h_noise->{$tid}->{$pos}->{$jid2} = undef;
						$h_noise->{all}->{$jid2} = undef;
					}
				}
			}
		}
		my $intspan_this = $junction2->getStrictGenomicSpan();
		my $inter1 = $intspan_self->intersection( $intspan_this );
		if (not $inter1->is_empty()) {
			my $jid2 = $junction2->id();
			$jid2 =~ s/_RI//;
			$jid2 =~ s/_SE//;
			$h_noise->{all}->{$jid2} = undef;
		}
		$i++;
	}
	$self->{h_noise}->{$patient->name()} = $h_noise;
	return $h_noise;
}


1;