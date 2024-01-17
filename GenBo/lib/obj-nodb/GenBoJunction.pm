package GenBoJunction;

use strict;
use Moo;
use Data::Dumper;
use GenBoCapture;
use Carp;
use List::Util qw[min max];
extends "GenBoGenomic";


has type => (
	is		=> 'ro',
	default	=> "junction",
);

has type_public_db => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return "junctions";
	},
);

has type_object => (
	is		=> 'ro',
	default	=> "junctions_object",
);

has name => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->id();
	},
);

has isJunction => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> 1,
);

#has coverage_start => (
#	is		=> 'ro',
#	lazy 	=> 1,
#	default	=> sub {
#		my $self = shift;
#		return $self->start();
#	},
#);
#
#has coverage_end => (
#	is		=> 'ro',
#	lazy 	=> 1,
#	default	=> sub {
#		my $self = shift;
#		return $self->end();
#	},
#);
#
#sub raw_coverage {
#	my ($self, $patient) = @_;
#	
#	warn "\n\n";
#	warn $self->id;
#	
#	my $gc_start = GenBoCoverageSamtools->new(chromosome=>$self->getChromosome, patient=>$patient, start=>($self->start()-5), end=>$self->start());
#	
#	warn "\n";
#	warn ($self->start()-5).' - '.$self->start();
#	warn Dumper $gc_start->{raw};
#	
#	my $gc_end = GenBoCoverageSamtools->new(chromosome=>$self->getChromosome, patient=>$patient, start=>$self->end(), end=>($self->end()+5));
#	warn "\n";
#	warn $self->end().' - '.($self->end()+5);
#	warn Dumper $gc_end->{raw};
#	
#	die;
#}

has get_hash_exons_introns => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		my @lPat = @{$self->getPatients()};
		my $ensid;
		$ensid = $self->annex->{$lPat[0]->name()}->{ensid} if (exists $self->annex->{$lPat[0]->name()}->{ensid});
		my $h_exons_introns;
		my (@lExons, @lIntrons);
		my $intspan_junction = $self->getStrictGenomicSpan();
		my $i = 0;
		foreach my $g (@{$self->getGenes()}) {
			my @ltmp = split('_', $g->id());
			next if ($ensid and $ltmp[0] ne $ensid);
			foreach my $t (@{$g->getTranscripts()}) {
				next unless $t->isMain();
				my $intspan_t = $t->getGenomicSpan();
				my $inter1 = $intspan_junction->intersection( $intspan_t );
				next if $inter1->is_empty();
				my @ltmp2 = split('_', $t->id());
				my $t_id_all = $t->id();
				my $t_id = $ltmp2[0];
				foreach my $e (@{$t->getExons()}) {
					my $intspan_e = $e->getGenomicSpan();
					my $inter2 = $intspan_junction->intersection( $intspan_e );
					next if ($inter2->is_empty());
					my $e_id = $e->id();
					$e_id =~ s/$t_id_all//;
					push(@lExons, $e_id);
					$h_exons_introns->{$t_id}->{exons}->{$e_id} = undef;
					$h_exons_introns->{$t_id}->{by_pos}->{$e->start} = $e_id;
					$i++;
				}
				foreach my $i (@{$t->getIntrons()}) {
					my $intspan_i = $i->getGenomicSpan();
					my $inter2 = $intspan_junction->intersection( $intspan_i );
					next if ($inter2->is_empty());
					my $i_id = $i->id();
					$i_id =~ s/$t_id_all//;
					push(@lIntrons, $i_id);
					$h_exons_introns->{$t_id}->{introns}->{$i_id} = undef;
					$h_exons_introns->{$t_id}->{by_pos}->{$i->start} = $i_id;
					$i++;
				}
			}
		}
		return $h_exons_introns;
	},
);

has get_hash_junctions_linked_to_me => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> undef,
);

sub is_junctions_linked {
	my ($self, $patient) = @_;
	confess() unless ($patient);
	return unless ($self->{get_hash_junctions_linked_to_me});
	return unless (exists $self->get_hash_junctions_linked_to_me->{$patient->name()});
	return 1; 
}

sub is_ri_aval {
	my ($self, $patient) = @_;
	confess() unless ($patient);
	my $type = lc($self->getTypeDescription($patient)); 
	return 1 if ($type =~ /ri_aval/);
	return;
}

sub is_ri_amont {
	my ($self, $patient) = @_;
	confess() unless ($patient);
	my $type = lc($self->getTypeDescription($patient)); 
	return 1 if ($type =~ /ri_amont/);
	return;
}

has is_dragen => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		foreach my $patient (@{$self->getPatients}) {
			return 1 if lc($self->getTypeDescription($patient)) eq 'dragen';
		}
		return;
	},
);

has is_star => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		foreach my $patient (@{$self->getPatients}) {
			return 1 if lc($self->getTypeDescription($patient)) eq 'star';
		}
		return;
	},
);

has isCanonique => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return if not -d $self->getChromosome()->get_lmdb_junctions_canoniques('r')->{dir}.'/'.$self->getChromosome->id();
		if ($self->is_dragen() or $self->is_star()) {
			my $id = $self->getChromosome->id().'_'.$self->start().'_'.$self->end();
			return 1 if $self->getChromosome()->get_lmdb_junctions_canoniques('r')->get($id);
		}
		else {
			my $id2 = $self->getChromosome->id().'_'.($self->start() + 1).'_'.($self->end() - 1);
			return 1 if $self->getChromosome()->get_lmdb_junctions_canoniques('r')->get($id2);
		}
		return;
	},
);

sub getTypeDescription {
	my ($self, $patient) = @_;
	confess() unless ($patient);
	return $self->annex->{$patient->name()}->{type}; 
}

sub isRI {
	my ($self, $patient) = @_;
	return $self->{isRi}->{$patient->name} if (exists $self->{isRi}->{$patient->name});
	my $isRi;
	$isRi = 1 if exists $self->annex->{$patient->name}->{type_origin_file} and $self->annex->{$patient->name}->{type_origin_file} eq 'RI';
	$isRi = 1 if $self->getTypeDescription($patient) and $self->getTypeDescription($patient) =~ /RI/; 
	$self->{isRi}->{$patient->name} = $isRi;
	return $isRi;
}

sub isSE {
	my ($self, $patient) = @_;
	return $self->{isSe}->{$patient->name} if (exists $self->{isSe}->{$patient->name});
	my $isSe;
	$isSe = 1 if exists $self->annex->{$patient->name}->{type_origin_file} and $self->annex->{$patient->name}->{type_origin_file} eq 'SE';
	$isSe = 1 if $self->getTypeDescription($patient) and $self->getTypeDescription($patient) =~ /SE/; 
	$self->{isSe}->{$patient->name} = $isSe;
	return $isSe;
}

has isCnv => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> 0,
);

has isLarge => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> 0,
);

has intspan => (
	is		=> 'ro',
	reader	=> 'getGenomicSpan',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $intSpan = Set::IntSpan::Fast::XS->new();
		$intSpan->add_range($self->start()-100, $self->start() + $self->length() +100);
    	return $intSpan;
	},
);

has strict_intspan => (
	is		=> 'ro',
	reader	=> 'getStrictGenomicSpan',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $intSpan = Set::IntSpan::Fast::XS->new();
		$intSpan->add_range($self->start(), $self->start() + $self->length());
    	return $intSpan;
	},
);

has annex => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> 0,
);

sub is_filtred_results {
	my ($self, $patient) = @_;
	confess() unless $patient;
	return 1 if (exists $self->annex->{$patient->name()}->{filtred_result} and $self->annex->{$patient->name()}->{filtred_result} == 1);
	return;
}

sub get_score {
	my ($self, $patient) = @_;
	confess() unless $patient;
	return $self->annex->{$patient->name()}->{score} if (exists $self->annex->{$patient->name()}->{score});
	return;
}

sub get_nb_new_count {
	my ($self, $patient) = @_;
	confess() unless $patient;
	return $self->{nb_new_count}->{$patient->name()} if (exists $self->{nb_new_count}->{$patient->name()});
	my $count = 0;
	if (exists $self->annex->{$patient->name()}->{junc_new_count}) {
		$count  = $self->annex->{$patient->name()}->{junc_new_count};
	}
	else {
		$count += $self->annex->{$patient->name()}->{junc_ri_count} if ($self->isRI($patient) and exists $self->annex->{$patient->name()}->{junc_ri_count});
		$count += $self->annex->{$patient->name()}->{junc_se_count} if ($self->isSE($patient) and exists $self->annex->{$patient->name()}->{junc_se_count});
	}
	$self->{nb_new_count}->{$patient->name()} = $count;
	return $count;
}

sub get_nb_normal_count {
	my ($self, $patient) = @_;
	confess() unless $patient;
	return $self->annex->{$patient->name()}->{junc_normale_count} if (exists $self->annex->{$patient->name()}->{junc_normale_count});
	return;
}

sub get_dp_count {
	my ($self, $patient) = @_;
	confess() unless $patient;
	return $self->{dp_count}->{$patient->name()} if (exists $self->{dp_count}->{$patient->name()});
	my $new_count = $self->get_nb_new_count($patient);
	$new_count = 0 if $self->get_nb_new_count($patient) eq '---';
	my $normal_count = $self->get_nb_normal_count($patient);
	$normal_count = 0 if $self->get_nb_normal_count($patient) eq '---';
	$self->{dp_count}->{$patient->name()} = ($new_count + $normal_count);
	return $self->{dp_count}->{$patient->name()};
}

sub get_ratio_new_count {
	my ($self, $patient) = @_;
	confess() unless $patient;
	return $self->{ratio_new_count}->{$patient->name()} if (exists $self->{ratio_new_count}->{$patient->name()});
	my $ratio;
	eval{ $ratio = $self->get_nb_new_count($patient) / $self->get_dp_count($patient); };
	if($@) { confess(); }
	$self->{ratio_new_count}->{$patient->name()} = $ratio;
	return $ratio;
} 

sub get_percent_new_count {
	my ($self, $patient) = @_;
	return $self->get_ratio_new_count($patient) * 100;
}

sub getSvgPlotPath {
	my ($self, $patient) = @_;
	my $path_svg = $patient->getJunctionsAnalysePath().'/../'.$self->annex->{$patient->name()}->{'ensid'}.'/SpliceRes/Figs/';
	my $svg_patient = $path_svg.'/juncPairPos_'.$self->annex->{$patient->name()}->{ensid}.'_'.$self->getChromosome->id().'_'.$patient->name().'_'.$patient->getProject->name().'.svg'; 
	return $svg_patient if (-e $svg_patient);
	return;
}

sub createSashiPlot {
	my ($self, $patient, $locus) = @_;
	my $file = $self->getSashimiPlotPath($patient, $locus);
	return if -e $file;
	my $patient_name = $patient->name();
	my $path_analysis = $patient->getJunctionsAnalysePath();
	my $ensg = $self->annex->{$patient->name()}->{'ensid'};
	my $score = $self->get_score($patient);
	my $project_name = $patient->getProject->name();
	my $cmd = $patient->getProject->buffer->software('ggsashimi');
	my $bam_path = $patient->getJunctionsAnalysePath()."/../resGenes/$ensg/rmdup/$patient_name\_$project_name\_rmdup.bam";
	if (not -e $bam_path) { $bam_path = $patient->getBamFiles->[0]; }
	$cmd .= " -b $bam_path";
	if ($self->getProject->is_human_genome()) { $cmd .= " -c chr".$locus; }
	else { $cmd .= " -c ".$locus; }
	
	$cmd .= " -o ".$file;
#	if ($score and $score >= 100) { $cmd .= " -P /data-isilon/bipd-src/mbras/ggsashimi/ggsashimi-master/colors/red.txt"; }
#	elsif ($score and $score >= 10) { $cmd .= " -P /data-isilon/bipd-src/mbras/ggsashimi/ggsashimi-master/colors/orange.txt"; }
#	elsif ($score and $score >= 1) { $cmd .= " -P /data-isilon/bipd-src/mbras/ggsashimi/ggsashimi-master/colors/green.txt"; }
#	else { $cmd .= " -P /data-isilon/bipd-src/mbras/ggsashimi/ggsashimi-master/colors/black.txt"; }
	$cmd .= " -P /data-isilon/bipd-src/mbras/ggsashimi/ggsashimi-master/colors/black.txt";
	$cmd .= " -C 1";
	$cmd .= " --shrink --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18";
	$cmd .= " -g ".$patient->getProject->get_gtf_genes_annotations_igv();
	$cmd .= " -F svg";
#	warn "\n";
#	warn $cmd;
	`$cmd`;
	return $file;
}

has can_create_sashimi_plots => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> 0,
);

sub createListSashimiPlots {
	my ($self, $patient) = @_;
	$self->can_create_sashimi_plots(1);
	$self->getListSashimiPlotsPathFiles($patient, '1');
}

sub getSashimiPlotPath {
	my ($self, $patient, $locus) = @_;
	my $locus_text;
	if ($locus) { $locus_text = $locus; }
	else { $locus_text = $self->getChromosome->name().':'.($self->start() - 100).'-'.($self->end() + 100); }
	$locus_text =~ s/chr//;
	$locus_text =~ s/:/-/;
	my $path = $patient->getProject->getProjectPath.'/align/sashimi_plots/';
	unless (-d $path) {
		mkdir $path;
		`chmod 775 $path`;
	}

	$path .= $patient->name().'/';
	unless (-d $path) {
		mkdir $path;
		`chmod 775 $path`;
	}
	my $file_name;
	if (exists $self->annex->{$patient->name()}->{'ensid'}) {
		$file_name = 'sashimi_'.$patient->name().'.'.$self->annex->{$patient->name()}->{'ensid'}.'.'.$locus_text.'.svg';
	}
	else {
		$file_name = 'sashimi_'.$patient->name().'.'.$locus_text.'.svg';
	}
	my $outfile = $path.'/'.$file_name;
	if (-e $path.'/../'.$file_name) {
		my $cmd = qq{mv $path/../$file_name $outfile};
		`$cmd`;

	}
	return $outfile;
}

sub getListSashimiPlotsPathFiles {
	my ($self, $patient) = @_;
	my $locus = $self->getChromosome->id().':'.($self->start() - 100).'-'.($self->end() + 100);
	my $sashimi_file = $self->getSashimiPlotPath($patient, $locus);
	$self->createSashiPlot($patient, $locus) if ($self->can_create_sashimi_plots());
	my @lFiles;
	if ($sashimi_file) {
		push(@lFiles, $sashimi_file);
		my $locus_text = $locus;
		$locus_text =~ s/chr//;
		$locus_text =~ s/:/-/;
		my ($chr_id, $start, $end) = split('-', $locus_text);
		my $i = 1;
		while ($i < 3) {
			$start -= (1000*$i);
			$end += (1000*$i);
			my $locus_extended = $chr_id.':'.$start.'-'.$end;
			my $sashimi_plot_file = $self->getSashimiPlotPath($patient, $locus_extended);
			$self->createSashiPlot($patient, $locus_extended) if ($self->can_create_sashimi_plots());
			push(@lFiles, $sashimi_plot_file);
			$i++;
		}
		return \@lFiles;
	}
	return;
}

has dejavu_percent_coordinate_similar => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> 98,
);

sub dejavu_nb_int_this_run_patients {
	my ($self, $patient_view, $min_percent) = @_;
	my $nb = 0;
	foreach my $patient (@{$self->getPatients()}) {
		next if ($patient_view and $patient_view->name() eq $patient->name());
		next if ($min_percent and $self->get_percent_new_count($patient) < $min_percent);
		$nb++;
	}
	return $nb;
}

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
	return;
}

sub get_noise_score {
	my ($self, $patient) = @_;
	my $h_noise = $self->get_hash_noise($patient);
	my $max_noise_in_tr = 0;
	return $max_noise_in_tr unless $h_noise;
	foreach my $enst (keys %{$h_noise}) {
		next if $enst eq 'all';
		foreach my $pos (keys %{$h_noise->{$enst}}) {
			my $this_noise = scalar keys %{$h_noise->{$enst}->{$pos}};
			$max_noise_in_tr = $this_noise if $this_noise > $max_noise_in_tr;
		}
	}
	my $noise_all = scalar keys %{$h_noise->{all}};
	return ($max_noise_in_tr, $noise_all);
}

sub junction_score_penality_ratio {
	my ($self, $patient) = @_;
	my $score_penality = 0;
	my $ratio = $self->get_percent_new_count($patient);
	if ($ratio <= 1) { $score_penality += 8; }
	elsif ($ratio <= 5) { $score_penality += 7; }
	elsif ($ratio <= 10) { $score_penality += 5; }
	elsif ($ratio <= 15) { $score_penality += 2.5; }
	elsif ($ratio <= 20) { $score_penality += 2; }
	elsif ($ratio <= 30) { $score_penality += 1; }
	return $score_penality;
}

sub junction_score_penality_dp {
	my ($self, $patient) = @_;
	return 0;
	
	my $score_penality = 0;
	my $dp = $self->get_dp_count($patient);
	if ($self->get_percent_new_count($patient) == 100) {
		if ($dp <= 5) { $score_penality += 8; }
		elsif ($dp <= 20) { $score_penality += 3; }
		elsif ($dp <= 50) { $score_penality += 1.5; }
		elsif ($dp <= 100) { $score_penality += 0.5; }
	}
	else {
		if ($dp <= 5) { $score_penality += 8; }
		elsif ($dp <= 10) { $score_penality += 3; }
		elsif ($dp <= 20) { $score_penality += 1.5; }
		elsif ($dp <= 30) { $score_penality += 0.5; }
	}
	return $score_penality;
}

sub junction_score_penality_new_junction {
	my ($self, $patient) = @_;
	my $score_penality = 0;
	my $nb_new = $self->get_nb_new_count($patient);
	if ($nb_new == 1) { $score_penality += 8; }
	elsif ($nb_new < 3) { $score_penality += 5; }
	elsif ($nb_new < 5) { $score_penality += 2; }
	elsif ($nb_new < 10) { $score_penality += 1; }
	return $score_penality;
}

sub junction_score_penality_noise {
	my ($self, $patient) = @_;
	my ($max_noise_in_tr, $noise_global) = $self->get_noise_score($patient);
	my ($penality_1, $penality_2, $penality_3, $penality_4, $penality_5) = (8, 6, 5, 3, 2);
	my $length = $self->length();
	if ($length >= 10000) {
		if ($max_noise_in_tr >= 30) { return 3; }
		elsif ($max_noise_in_tr >= 20) { return 2; }
		elsif ($max_noise_in_tr >= 15) { return 1; }
	}
	elsif ($length >= 1000) {
		if ($max_noise_in_tr >= 30) { return 4; }
		elsif ($max_noise_in_tr >= 25) { return 3; }
		elsif ($max_noise_in_tr >= 20) { return 2.5; }
		elsif ($max_noise_in_tr >= 15) { return 1.5; }
		elsif ($max_noise_in_tr >= 10) { return 1; }
	}
	else {
		if ($max_noise_in_tr >= 30) { return 8; }
		elsif ($max_noise_in_tr >= 25) { return 6; }
		elsif ($max_noise_in_tr >= 20) { return 5; }
		elsif ($max_noise_in_tr >= 15) { return 6; }
		elsif ($max_noise_in_tr >= 10) { return 2; }
		elsif (defined($noise_global) and $noise_global >= 30) { return 4; }
		elsif (defined($noise_global) and $noise_global >= 25) { return 3; }
		elsif (defined($noise_global) and $noise_global >= 20) { return 2; }
		elsif (defined($noise_global) and $noise_global >= 15) { return 1.5; }
		elsif (defined($noise_global) and $noise_global >= 10) { return 1; }
	}
	
	return 0;
}

sub junction_score_penality_dejavu {
	my ($self, $patient_name) = @_;
	my $my_ratio = $self->get_percent_new_count($patient_name);
	#return $self->junction_score_penality_dejavu_high_ratio($patient) if ($my_ratio >= 70);
	return $self->junction_score_penality_dejavu_low_ratio($patient_name);
}

sub junction_score_penality_dejavu_low_ratio {
	my ($self, $patient) = @_;
	my $use_percent_dejavu = $self->dejavu_percent_coordinate_similar();
	my $dv_ratio_all = 	$self->project->dejavuJunctionsResume->get_nb_junctions($self->getChromosome->id(), $self->start(),$self->end(), $use_percent_dejavu, $patient->name);
	my $dv_ratio_10 = $self->project->dejavuJunctionsResume->get_nb_junctions_ratio10($self->getChromosome->id(), $self->start(),$self->end(), $use_percent_dejavu, $patient->name);
	my $dv_ratio_20 = $self->project->dejavuJunctionsResume->get_nb_junctions_ratio20($self->getChromosome->id(), $self->start(),$self->end(), $use_percent_dejavu, $patient->name);
	if    ($dv_ratio_20 >= 10 or $dv_ratio_10 >= 20) { return 8; }
	elsif ($dv_ratio_20 >= 5 or $dv_ratio_10 >= 10)  { return 5; }
	elsif ($dv_ratio_20 >= 2 or $dv_ratio_10 >= 5)   { return 3; }
	elsif ($dv_ratio_20 >= 1 or $dv_ratio_10 >= 3)   { return 2; }
	return 0;
}

sub junction_score_penality_dejavu_high_ratio {
	my ($self, $patient) = @_;
	my $dv_ratio_all = $self->dejavu_nb_other_patients($patient);
	my $dv_ratio_70 = $self->dejavu_nb_other_patients_min_ratio_70($patient);
	my $dv_ratio_90 = $self->dejavu_nb_other_patients_min_ratio_90($patient);
	
	if ($dv_ratio_90 >= 20) { return 8; }
	elsif ($dv_ratio_90 >= 10) { return 5; }
	elsif ($dv_ratio_70 >= 20) { return 5; }
	elsif ($dv_ratio_70 >= 15) { return 2; }
	elsif ($dv_ratio_70 >= 5) { return 1; }
	return 0;
}

sub junction_score_penality_dejavu_inthisrun {
	my ($self, $patient) = @_;
	my $score_penality = 0;
	my $dv_run_max_30 = $self->dejavu_nb_int_this_run_patients($patient, 30);
	if ($dv_run_max_30 >= 1) { $score_penality += (1 * ($dv_run_max_30)); }
	my $dv_run_max_20 = $self->dejavu_nb_int_this_run_patients($patient, 20);
	if ($dv_run_max_20 >= 1) { $score_penality += (0.5 * ($dv_run_max_20)); }
	my $dv_run_max_15 = $self->dejavu_nb_int_this_run_patients($patient, 15);
	if ($dv_run_max_15 > 1) { $score_penality += 0.5; }
	my $dv_run_max_10 = $self->dejavu_nb_int_this_run_patients($patient, 10);
	if ($dv_run_max_10 > 1) { $score_penality += 0.5; }
	return $score_penality;
}

sub is_junction_with_known_coordinates {
	my ($self) = @_;
	#$h_exons_introns->{$t_id}->{by_pos}->{$i->start} = $i_id;
	my $h_exons_introns = $self->get_hash_exons_introns();
	my ($start_known, $end_known);
	foreach my $t_id (keys %{$h_exons_introns}) {
		foreach my $pos (keys %{$h_exons_introns->{$t_id}->{by_pos}}) {
			$start_known = 1 if $pos == $self->start();
			$end_known = 1 if $pos == $self->end(); 
		}
	}
	return 1 if $start_known and $end_known;
	return;
}

sub junction_score_penality_known_coordinates {
	my ($self) = @_;
	return 100000 if $self->is_junction_with_known_coordinates();
	return 0;
}

sub junction_score_without_dejavu_global {
	my ($self, $patient) = @_;
	my $score = 10;
	$score = 99999 if ($self->isCanonique());
	$score -= $self->junction_score_penality_ratio($patient);
	$score -= $self->junction_score_penality_dp($patient);
	$score -= $self->junction_score_penality_new_junction($patient);
	$score -= $self->junction_score_penality_noise($patient);
	$score -= $self->junction_score_penality_dejavu_inthisrun($patient);
	return $score;
}

sub junction_score {
	my ($self, $patient, $use_percent_dejavu) = @_;
	$self->dejavu_percent_coordinate_similar($use_percent_dejavu) if $use_percent_dejavu;
	my $gene;
	if (exists $self->annex->{$patient->name()}->{ensid}) {
		foreach my $g (@{$self->getGenes()}) {
			my @ltmp = split('_', $g->id());
			next if ($ltmp[0] ne $self->annex->{$patient->name()}->{ensid});
			$gene = $g;
		}
		return -999 unless $gene;
		confess() unless $gene;
	}
	my $score = $self->junction_score_without_dejavu_global($patient);
	$score -= $self->junction_score_penality_dejavu($patient);
	$score -= $self->junction_score_penality_known_coordinates();
	return $score;
}

1;