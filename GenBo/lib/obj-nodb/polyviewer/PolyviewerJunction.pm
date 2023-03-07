package PolyviewerJunction;
use strict;
use FindBin qw($Bin);
use Moose;
use Data::Dumper;
use CGI qw/:standard :html3/;
extends 'PolyviewerVariant';

my $cgi = new CGI();



has project_name => (is => 'rw');

has min_score_used => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub { return 20; }
);

has vector_id => (is => 'rw');

has ensid => (is => 'rw');

has isRI => (is	=> 'rw');

has isSE => (is => 'rw');

has locus_extended_100nt => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my ($self) = @_;
		return $self->chromosome.":".($self->start()-100)."-".($self->end()+100);
	}
);

has hash_by_pat_dp_count => (is => 'rw');

has hash_by_pat_nb_new => (is => 'rw');

has hash_by_pat_nb_normal => (is => 'rw');

has hash_by_pat_percent => (is => 'rw');

has hash_by_pat_family_name => (is => 'rw');

has hash_by_pat_status => (is => 'rw');

has hash_by_pat_html_status => (is => 'rw');

has type_description => (is => 'rw');

has list_sashimi_plots => (is => 'rw');

has bam_url => (is => 'rw');

has bam_controls_urls => (is => 'rw');

has gtf_url => (is => 'rw');

has cmd_dejavu => (is => 'rw');

has cmd_dejavu_inthisrun => (is => 'rw');

has dejavu_nb_other_patients => (is => 'rw');

has dejavu_nb_other_patients_min_ratio_10 => (is => 'rw');

has dejavu_nb_other_patients_min_ratio_20 => (is => 'rw');

has dejavu_nb_int_this_run_patients => (is => 'rw');

has dejavu_nb_int_this_run_patients_min_ratio_10 => (is => 'rw');

has dejavu_nb_int_this_run_patients_min_ratio_20 => (is => 'rw');

has is_junctions_linked => (is => 'rw');

has hash_exons_introns => (is => 'rw');

has hash_junctions_linked_to_me => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub { return {}; }
);

has list_colors_linked_junctions => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my ($self) = @_;
		my @l_group_junctions_colors = ('#5D3EFF', '#FF4571', '#8FFF49', '#FF495F');
		return \@l_group_junctions_colors;
	}
);

has index_group_junction_colors => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub { return 0; }
);

has hash_transcripts_descriptions => (is => 'rw');

has hash_transcripts_styles_colors => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub { return {}; }
);

has hash_junctions_linked_colors => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub { return {}; }
);

has cmd_view_linked_junctions => (is => 'rw');

has score => (is => 'rw');

has score_penality_ratio => (is => 'rw');

has score_penality_dp => (is => 'rw');

has score_penality_new_junction => (is => 'rw');

has score_penality_noise => (is => 'rw');

has nb_noise_junctions => (is => 'rw');

has score_penality_dejavu => (is => 'rw');

has score_penality_dejavu_inthisrun => (is => 'rw');



sub setJunction {
	my ($self, $junction, $patient, $min_score_used) = @_;
	$self->min_score_used($min_score_used) if (defined($min_score_used));
	my $project_name = $patient->getProject->name();
	my $patient_name = $patient->name();
	$self->project_name($project_name);
	$self->ensid($junction->annex->{$patient->name()}->{ensid});
	$self->id($junction->id());
	$self->vector_id($junction->vector_id());
	$self->start($junction->start());
	$self->end($junction->end());
	$self->chromosome($junction->getChromosome->id());
	$self->isCnv(0);
	$self->isJunction(1);
	$self->isRI(1) if $junction->isRI($patient);
	$self->isSE(1) if $junction->isSE($patient);
	$self->type_description($junction->getTypeDescription($patient));
	$self->type($junction->type());
	$self->name($junction->id());
	$self->list_sashimi_plots($junction->getListSashimiPlotsPathFiles($patient));
	$self->bam_url("https://www.polyweb.fr/".$patient->bamUrl());
	$self->bam_controls_urls( join(',', @{$self->get_list_bam_controls_urls($junction, $patient)}) );
	$self->gtf_url($self->get_gtf_url($junction->getProject()));
	
	foreach my $patient_family (@{$patient->getFamily->getPatients()}) {
		$self->{hash_by_pat_family_name}->{$patient_family->name()} = $patient_family->getFamily->name();
		$self->{hash_by_pat_status}->{$patient_family->name()} = $patient_family->status();
		if ($patient_family->isFather()) {
			if ($patient_family->isIll()) { $self->{hash_by_pat_html_status}->{$patient_family->name()} = "<center><img src='/icons/Polyicons/male-d.png'></center>"; }
			else { $self->{hash_by_pat_html_status}->{$patient_family->name()} = "<center><img src='/icons/Polyicons/male-s.png'></center>"; }
		}
		if ($patient_family->isMother()) {
			if ($patient_family->isIll()) { $self->{hash_by_pat_html_status}->{$patient_family->name()} = "<center><img src='/icons/Polyicons/female-d.png'></center>"; }
			else { $self->{hash_by_pat_html_status}->{$patient_family->name()} = "<center><img src='/icons/Polyicons/female-s.png'></center>"; }
		}
		if ($patient_family->isChild()) {
			if ($patient_family->sex() eq '1') { 
				if ($patient_family->isIll()) { $self->{hash_by_pat_html_status}->{$patient_family->name()} = "<center><img src='/icons/Polyicons/baby-boy-d.png'></center>"; }
				else { $self->{hash_by_pat_html_status}->{$patient_family->name()} = "<center><img src='/icons/Polyicons/baby-boy-s.png'></center>"; }
			}
			else {
				if ($patient_family->isIll()) { $self->{hash_by_pat_html_status}->{$patient_family->name()} = "<center><img src='/icons/Polyicons/baby-girl-d.png'></center>"; }
				else { $self->{hash_by_pat_html_status}->{$patient_family->name()} = "<center><img src='/icons/Polyicons/baby-girl-s.png'></center>"; }
			}
		}
		$self->{hash_by_pat_dp_count}->{$patient_family->name()} = $junction->get_dp_count($patient_family);
		$self->{hash_by_pat_nb_new}->{$patient_family->name()} = $junction->get_nb_new_count($patient_family);
		$self->{hash_by_pat_nb_normal}->{$patient_family->name()} = $junction->get_nb_normal_count($patient_family);
		$self->{hash_by_pat_percent}->{$patient_family->name()} = sprintf("%.3f", $junction->get_percent_new_count($patient_family)).'%';
	}
	
	my $vector_dv_id = $junction->getChromosome->id().'-'.$junction->vector_id();
	$self->cmd_dejavu( qq{view_deja_vu_rna_junction(\"$project_name\",\"$patient_name\",\"$vector_dv_id\")} );
	$self->cmd_dejavu_inthisrun( qq{view_dejavu_nb_int_this_run_patients(\"$project_name\",\"$patient_name\",\"$vector_dv_id\")} );
	$self->update_dejavu_values($junction, $patient);
	
	$self->is_junctions_linked($junction->is_junctions_linked($patient));
	$self->hash_exons_introns($junction->get_hash_exons_introns());
	my $h_junctions_linked_to_me = $junction->get_hash_junctions_linked_to_me->{$patient->name()} if ($junction->get_hash_junctions_linked_to_me() and exists $junction->get_hash_junctions_linked_to_me->{$patient->name()});
	$self->hash_junctions_linked_to_me($h_junctions_linked_to_me) if $h_junctions_linked_to_me;
	
	if ($self->hash_exons_introns()) {
		foreach my $tid (keys %{$self->hash_exons_introns()}) {
			my $t;
			eval { $t = $patient->getProject->newTranscript($tid); };
			if ($@) {
				$self->{hash_transcripts_descriptions}->{$tid}->{external_name} = 'pb...';
				$self->{hash_transcripts_descriptions}->{$tid}->{ccds_name} = 'pb...';
				$self->{hash_transcripts_descriptions}->{$tid}->{appris_type} = 'pb...';
			}
			else {
				$self->{hash_transcripts_descriptions}->{$tid}->{external_name} = $t->external_name();
				$self->{hash_transcripts_descriptions}->{$tid}->{ccds_name} = $t->ccds_name();
				$self->{hash_transcripts_descriptions}->{$tid}->{appris_type} = $t->appris_type();
			}
		}
		$self->prepare_transcripts_consequence($patient);
	}
	
	$self->nb_noise_junctions( $junction->get_noise_score($patient) );
	$self->score_penality_noise( $junction->junction_score_penality_noise($patient) );
	$self->score_penality_ratio( $junction->junction_score_penality_ratio($patient) );
	$self->score_penality_dp( $junction->junction_score_penality_dp($patient) );
	$self->score_penality_new_junction( $junction->junction_score_penality_new_junction($patient) );
	$self->score_penality_dejavu_inthisrun( $junction->junction_score_penality_dejavu_inthisrun($patient) );
	$self->update_junction_score_value($junction, $patient);
	return $self;
}

sub prepare_transcripts_consequence {
	my ($self, $patient) = @_;
	my $patient_name = $patient->name();
	my $project_name = $self->project_name();
	my (@junctions_ids_linked, $bcolor);
	@junctions_ids_linked = keys %{$self->hash_junctions_linked_to_me()} if $self->is_junctions_linked();
	my @l_group_junctions_colors = @{$self->list_colors_linked_junctions()};
	
	if (@junctions_ids_linked and $self->is_junctions_linked()) {
		foreach my $tid (sort keys %{$self->hash_exons_introns()}) {
			my ($h_junctions_linked, $h_junctions_exons_introns);
			
			my @lPos = (sort keys %{$self->hash_exons_introns->{$tid}->{by_pos}});
			next unless @lPos;
			my $first_exon_intron = $self->hash_exons_introns->{$tid}->{by_pos}->{$lPos[0]};
			my $last_exon_intron = $self->hash_exons_introns->{$tid}->{by_pos}->{$lPos[-1]};
			
			my $h_j_by_tr;
			my ($first_found, $last_found);
			foreach my $other_j_id (@junctions_ids_linked) {
				$first_found = 1 if (exists $self->hash_junctions_linked_to_me->{$other_j_id}->{$tid}->{$first_exon_intron});
				$last_found = 1 if (exists $self->hash_junctions_linked_to_me->{$other_j_id}->{$tid}->{$last_exon_intron});
				
				if ($first_found or $last_found) {
					$h_junctions_exons_introns->{$tid}->{$first_exon_intron} = undef;
					my $this_chr = $patient->getProject->getChromosome($self->chromosome());
					my $id_in_list;
					my $v_j_p = $patient->getJunctionsVector($this_chr);
					my @lVectorIds = @{$this_chr->getListVarVectorIds($v_j_p)};
					my $i = 0;
					foreach my $v_id (@lVectorIds) {
						if ($v_id == $self->vector_id()) {
							$id_in_list = $i;
							last;
						}
						$i++;
					}
					my $min = $id_in_list - 15;
					$min = 0 if $min < 0;
					my $max = $id_in_list + 15;
					$max = scalar(@lVectorIds) - 1 if $max >= scalar(@lVectorIds);
					
					$i = $id_in_list - 1;
					while ($i >= $min) {
						my $junction2 = $this_chr->getVarObject($lVectorIds[$i]);
						$i--;
						next if (not $junction2->get_hash_junctions_linked_to_me());
						my @lPos = (sort keys %{$junction2->get_hash_exons_introns->{$tid}->{by_pos}});
						next unless (@lPos);
						my $first_exon_intron_2 = $junction2->get_hash_exons_introns->{$tid}->{by_pos}->{$lPos[0]};
						my $last_exon_intron_2 = $junction2->get_hash_exons_introns->{$tid}->{by_pos}->{$lPos[-1]};
						
						$h_j_by_tr->{$first_exon_intron_2}->{$this_chr->id().'-'.$junction2->vector_id()} = undef;
						$h_j_by_tr->{$last_exon_intron_2}->{$this_chr->id().'-'.$junction2->vector_id()} = undef;
					}
					$i = $id_in_list + 1;
					while ($i <= $max) {
						my $junction2 = $this_chr->getVarObject($lVectorIds[$i]);
						$i++;
						next if (not $junction2->get_hash_junctions_linked_to_me());
						my @lPos = (sort keys %{$junction2->get_hash_exons_introns->{$tid}->{by_pos}});
						next unless (@lPos);
						my $first_exon_intron_2 = $junction2->get_hash_exons_introns->{$tid}->{by_pos}->{$lPos[0]};
						my $last_exon_intron_2 = $junction2->get_hash_exons_introns->{$tid}->{by_pos}->{$lPos[-1]};
						
						$h_j_by_tr->{$first_exon_intron_2}->{$this_chr->id().'-'.$junction2->vector_id()} = undef;
						$h_j_by_tr->{$last_exon_intron_2}->{$this_chr->id().'-'.$junction2->vector_id()} = undef;
					}
				}
			}
			
			my $h_all;
			if (exists $h_j_by_tr->{$first_exon_intron}) {
				$h_all->{$self->chromosome().'-'.$self->vector_id()} = undef;
				foreach my $vid (keys %{$h_j_by_tr->{$first_exon_intron}}) { $h_all->{$vid} = undef; }
			}
			if (exists $h_j_by_tr->{$last_exon_intron}) {
				$h_all->{$self->chromosome().'-'.$self->vector_id()} = undef;
				foreach my $vid (keys %{$h_j_by_tr->{$last_exon_intron}}) { $h_all->{$vid} = undef; }
			}
			my $j_linked = join(',', sort keys %$h_all);
			
			my $min_score_used = $self->min_score_used();
			my $my_junction_id = $self->id();
			my $cmd_linked = qq{view_linked_junctions(\"$patient_name\",\"$tid\",\"$j_linked\",\"$my_junction_id\",\"$min_score_used\")};
			$self->{cmd_view_linked_junctions}->{$tid} = $cmd_linked;
		}
	}
	return;
}

sub update_dejavu_values {
	my ($self, $junction, $patient) = @_; 
	$self->dejavu_nb_other_patients($junction->dejavu_nb_other_patients($patient));
	if ($self->dejavu_nb_other_patients() > 0) {
		$self->dejavu_nb_other_patients_min_ratio_10($junction->dejavu_nb_other_patients_min_ratio_10($patient));
		$self->dejavu_nb_other_patients_min_ratio_20($junction->dejavu_nb_other_patients_min_ratio_20($patient));
	}
	else {
		$self->dejavu_nb_other_patients_min_ratio_10(0);
		$self->dejavu_nb_other_patients_min_ratio_20(0);
	}
	
	$self->dejavu_nb_int_this_run_patients($junction->dejavu_nb_int_this_run_patients($patient));
	if ($self->dejavu_nb_int_this_run_patients() > 0) {
		$self->dejavu_nb_int_this_run_patients_min_ratio_10($junction->dejavu_nb_int_this_run_patients($patient, 10));
		$self->dejavu_nb_int_this_run_patients_min_ratio_20($junction->dejavu_nb_int_this_run_patients($patient, 20));
	}
	else {
		$self->dejavu_nb_int_this_run_patients_min_ratio_10(0);
		$self->dejavu_nb_int_this_run_patients_min_ratio_20(0);
	}
}

sub update_junction_score_value {
	my ($self, $junction, $patient) = @_; 
	$self->score_penality_dejavu( $junction->junction_score_penality_dejavu($patient) );
	$self->score($junction->junction_score($patient));
}

sub get_gtf_url {
	my ($self, $project) = @_; 
	my $gtf = $project->get_gtf_genes_annotations_igv();
	$gtf =~ s/\/data-isilon//;
	$gtf = "https://www.polyweb.fr/".$gtf;
	return $gtf;
}

sub get_list_bam_controls_urls {
	my ($self, $junction, $patient) = @_; 
	my @lBams;
	my $list_patients_ctrl = $patient->getPatients_used_control_rna_seq_junctions_analyse();
	if ($list_patients_ctrl) {
		my $nb_control;
		foreach my $other_pat (@$list_patients_ctrl) {
			push (@lBams, 'https://www.polyweb.fr/'.$other_pat->bamUrl());
			$nb_control++;
			last if $nb_control == 3;
		}
	}
	else {
		my $np = 0;
		foreach my $other_pat (@{$patient->getProject->getPatients()}) {
			next if ($other_pat->name() eq $patient->name());
			push (@lBams, 'https://www.polyweb.fr/'.$other_pat->bamUrl());
			$np++;
			last if $np == 3;
		}
	}
	return \@lBams;
}

1;
