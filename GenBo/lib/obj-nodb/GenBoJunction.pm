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
has sj_id => (
	is		=> 'ro',
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
		my $h_exons_introns;
		my (@lExons, @lIntrons);
		my $intspan_junction = $self->getStrictGenomicSpan();
		my $i = 0;
		foreach my $g (@{$self->getGenes()}) {
			my @ltmp = split('_', $g->id());
			foreach my $t (@{$g->getTranscripts()}) {
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
#	default	=> sub {
#		my $self = shift;
#
#		if ($self->is_dragen() or $self->is_star()) {
#			die();
#			my $id = $self->getChromosome->id().'_'.$self->start().'_'.$self->end();
#			return 1 if $self->getChromosome()->get_lmdb_junctions_canoniques('r')->get($id);
#		}
#		else {
#			my $id2 = $self->getChromosome->id().'_'.($self->start() + 1).'_'.($self->end() - 1);
#			return 1 if $self->getChromosome()->get_lmdb_junctions_canoniques('r')->get($id2);
#		}
#		return;
#	},
);

sub getTypeDescription {
	my ($self, $patient) = @_;
	confess() unless ($patient);
	
	return $self->annex->{$patient->name()}->{type}; 
}

sub isRegtools {
	my ($self, $patient) = @_;
	return 1 if $self->getTypeDescription($patient) eq 'NDA';
	return 1 if $self->getTypeDescription($patient) eq 'DA';
	return 1 if $self->getTypeDescription($patient) eq 'N';
	return 1 if $self->getTypeDescription($patient) eq 'D';
	return 1 if $self->getTypeDescription($patient) eq 'A';
	return;
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
	elsif (exists $self->annex->{$patient->name()}->{alt_count}) {
		$count  = $self->annex->{$patient->name()}->{alt_count};
	}
	else {
		$count += $self->annex->{$patient->name()}->{junc_ri_count} if ($self->isRI($patient) and exists $self->annex->{$patient->name()}->{junc_ri_count});
		$count += $self->annex->{$patient->name()}->{junc_se_count} if ($self->isSE($patient) and exists $self->annex->{$patient->name()}->{junc_se_count});
	}
	$self->{nb_new_count}->{$patient->name()} = $count;
	return $count;
}

sub get_canonic_count {
	my ($self, $patient) = @_;
	confess() unless $patient;
	if ($self->isCanonique()) { 
		$self->{nb_new_count}->{$patient->name()} = 0;
		return 0;
	} 
	if ($self->isRI($patient) or $self->isSE($patient)) {
		my $nb;
		$nb = $self->{annex}->{$patient->name}->{canonic_count} if (exists $self->annex->{$patient->name()}->{canonic_count});
		return $nb if $nb and $nb > 0;
		return $self->getCoverageInInterval($patient);
	}
	return $self->{annex}->{$patient->name}->{canonic_count} if (exists $self->annex->{$patient->name()}->{canonic_count});
	return 0;
	confess();
}

sub get_dp_count {
	my ($self, $patient) = @_;
	confess() unless $patient;
	return $self->{dp_count}->{$patient->name()} if (exists $self->{dp_count}->{$patient->name()});
	my $new_count = $self->get_nb_new_count($patient);
	$new_count = 0 if $self->get_nb_new_count($patient) eq '---';
	
	my $normal_count = $self->get_canonic_count($patient);
	$normal_count = 0 if $self->get_canonic_count($patient) eq '---';
	
#	warn $new_count." ".$normal_count." ".$self->id;
	$self->{dp_count}->{$patient->name()} = ($new_count + $normal_count);
	return $self->{dp_count}->{$patient->name()};
}

sub ratio {
	my ($self, $patient) = @_;
	confess() unless $patient;
	return $self->{ratio}->{$patient->name()} if (exists $self->{ratio}->{$patient->name()});
	if ($self->isCanonique()) { 
		$self->{ratio}->{$patient->name()} = 0;
		return 0;
	} 
	if ($self->get_dp_count($patient) == 0){
		warn "\n\n";
		warn $self->id;
		warn $patient->name();
		warn 'NEW: '.$self->get_nb_new_count($patient);
		warn 'DP: '.$self->get_dp_count($patient);
		warn Dumper $self->annex();
		warn "\n\n";
	}
	
	
	my $ratio  = $self->get_nb_new_count($patient) / $self->get_dp_count($patient);
	$self->{ratio}->{$patient->name()} = $ratio;
	return $ratio;
}
sub get_ratio_new_count {
	my ($self, $patient) = @_;
	confess();
	return $self->ratio($patient);
} 

sub get_percent_new_count {
	my ($self, $patient) = @_;
	return $self->ratio($patient) * 100;
}

sub getSvgPlotPath {
	my ($self, $patient) = @_;
	my $path_svg = $patient->getJunctionsAnalysePath().'/../'.$self->annex->{$patient->name()}->{'ensid'}.'/SpliceRes/Figs/';
	my $svg_patient = $path_svg.'/juncPairPos_'.$self->annex->{$patient->name()}->{ensid}.'_'.$self->getChromosome->id().'_'.$patient->name().'_'.$patient->getProject->name().'.svg'; 
	return $svg_patient if (-e $svg_patient);
	return;
}

sub createSashiPlot {
	my ($junction, $patient, $locus, $bam_path) = @_;
	my $file = $junction->getSashimiPlotPath($patient, $locus);
	warn "\n# check $file\n";
	return if -e $file;
	my $nb_new = $junction->get_nb_new_count($patient);
	my $patient_name = $patient->name();
	my $path_analysis = $patient->getJunctionsAnalysePath();
	my $ensg = $junction->annex->{$patient->name()}->{'ensid'};
	my $score = $junction->get_score($patient);
	my $project_name = $patient->getProject->name();
	if (not $bam_path) { confess(); $bam_path = $patient->getBamFiles->[0]; }
	my $cmd = $patient->getProject->buffer->software('ggsashimi');
	$cmd .= " -b $bam_path ";
	if ($junction->getProject->is_human_genome()) { $cmd .= " -c chr".$locus; }
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
	if ($nb_new >= 1000) { $cmd .= " -M 100"; }
	elsif ($nb_new >= 100)  { $cmd .= " -M 20"; }
	elsif ($nb_new >= 50)   { $cmd .= " -M 10"; }
	elsif ($nb_new >= 20)   { $cmd .= " -M 5"; }
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
	my ($self, $patient, $bam_tmp) = @_;
	my $locus = $self->getChromosome->id().':'.($self->start() - 100).'-'.($self->end() + 100);
	my $sashimi_file = $self->getSashimiPlotPath($patient, $locus);
	$self->createSashiPlot($patient, $locus, $bam_tmp) if ($self->can_create_sashimi_plots());
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
			$self->createSashiPlot($patient, $locus_extended, $bam_tmp) if ($self->can_create_sashimi_plots());
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



sub setPatients {
	my $self = shift;
	my $h;
	foreach my $patient (@{$self->getProject->getPatients()}) {
		$h->{$patient->id} = undef if exists $self->annex->{$patient->name()};
	}
	return $h;
}
sub setGenes {
	my $self = shift;
#	return [keys %{$self->vannot_chr->{genes}}];

	my $hGenesid;
	map {$hGenesid->{$_} ++} @{$self->getChromosome()->genesIntervalTree->fetch($self->start-10,$self->start+10)};
	map {$hGenesid->{$_} ++} @{$self->getChromosome()->genesIntervalTree->fetch($self->end-10,$self->end+10)};
	my @t = keys  %$hGenesid;
	$self->vannot({}) unless @t;
	return \@t;
}

sub getGenes {
	my $self = shift;
	my (@lGenes, @lGenes_default);
	foreach my $g (@{$self->getProject()->myflushobjects($self->genes_object(), "genes")}) {
		my $get;
		foreach my $t (@{$g->getTranscripts()}) {
			my $intspan_coding = $t->genomic_span();
			if ($intspan_coding->contains($self->start()) or $intspan_coding->contains($self->start()-1)) { $get = 1; }
			elsif ($intspan_coding->contains($self->end()) or $intspan_coding->contains($self->end()+1))  { $get = 1; }
		}
		if ($get) { push(@lGenes, $g); }
		elsif (scalar(@lGenes) == 0) {
			if ($g->getGenomicSpan->contains($self->start()))  { push(@lGenes_default, $g) ; }
			elsif ($g->getGenomicSpan->contains($self->end())) { push(@lGenes_default, $g); }
		}
	}
	return \@lGenes if (scalar(@lGenes) > 0);
	return \@lGenes_default;
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


sub multiple_align_count {
	my ($self,$patient) = @_;
	return 0 unless exists  $self->annex->{$patient->name}->{junc_multiple_count_sj};
	return $self->annex->{$patient->name}->{junc_multiple_count_sj};
}

sub is_sj {
	my ($self,$patient) = @_;
	return exists $self->annex->{$patient->name}->{is_sj}
}
sub junction_score_penality_new_junction {
	my ($self, $patient) = @_;
	my $score_penality = 0;
	my $nb_new = $self->get_nb_new_count($patient);
	$nb_new -= $self->multiple_align_count($patient);
	
	if ($nb_new <= 2) { $score_penality = 9; }
	elsif ($nb_new <= 3) { $score_penality += 5; }
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
	my $dv_ratio_all = 	$self->dejavu_patients("all",$patient);
	my $dv_ratio_10 = 	$self->dejavu_patients(10,$patient);
	
	my $dv_ratio_20 = $self->dejavu_patients(20,$patient);
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
	my $dv_run_max_30 = $self->in_this_run_patients(30,$patient);
	if ($dv_run_max_30 >= 1) { $score_penality += (1 * ($dv_run_max_30)); }
	my $dv_run_max_20 = $self->in_this_run_patients(20,$patient);
	if ($dv_run_max_20 >= 1) { $score_penality += (0.5 * ($dv_run_max_20)); }
	my $dv_run_max_15 = $self->in_this_run_patients(15,$patient);
	if ($dv_run_max_15 > 1) { $score_penality += 0.5; }
	my $dv_run_max_10 = $self->in_this_run_patients( 10,$patient);
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
	$score = 0 if ($self->isCanonique());
	$score = 0 if ($self->getTypeDescription($patient) eq 'DA');
	$score -= $self->junction_score_penality_ratio($patient);
	$score -= $self->junction_score_penality_dp($patient);
	$score -= $self->junction_score_penality_new_junction($patient);
	$score -= $self->junction_score_penality_noise($patient);
	$score -= $self->junction_score_penality_dejavu_inthisrun($patient);
	if ($self->length < 50) {
		my $short_j_interesting = 1;
		my $h_e_i = $self->get_hash_exons_introns();
		my $h_only_exon_intron;
		foreach my $enst (keys %{$h_e_i}) {
			my @lPos =  keys %{$h_e_i->{$enst}->{'by_pos'}};
			if (scalar @lPos > 1) {
				$h_only_exon_intron = undef;
				last;
			}
			$h_only_exon_intron->{$enst} = $h_e_i->{$enst}->{'by_pos'}->{$lPos[0]};
		}
		if ($h_only_exon_intron) {
			$short_j_interesting = undef;
			foreach my $enst (keys %{$h_only_exon_intron}) {
				my $found;
				my $t = $self->getProject->newTranscript($enst);
				foreach my $exon (@{$t->getExons}) {
					next if $exon->id() ne $t->id.$h_only_exon_intron->{$enst};
					$short_j_interesting = 1 if ($self->start >= $exon->start - 20) and ($self->start <= $exon->start + 20);
					$short_j_interesting = 1 if ($self->start >= $exon->end - 20) and ($self->start <= $exon->end + 20);
					$short_j_interesting = 1 if ($self->end >= $exon->start - 20) and ($self->end <= $exon->start + 20);
					$short_j_interesting = 1 if ($self->end >= $exon->end - 20) and ($self->end <= $exon->end + 20);
				}
			}
		}
		$score = $score - 8  if not $short_j_interesting;
	}
	return $score;
}

sub junction_score {
	my ($self, $patient, $use_percent_dejavu) = @_;
	$self->dejavu_percent_coordinate_similar($use_percent_dejavu) if $use_percent_dejavu;
	my @lGenes = @{$self->getGenes()};
	my $score = $self->junction_score_without_dejavu_global($patient);
	$score -= $self->junction_score_penality_dejavu($patient);
	$score -= $self->junction_score_penality_known_coordinates();
	$score -= 999 if not @lGenes; 
	return $score;
}


#########
# DejaVu
###########

sub dejavu_nb_int_this_run_patients {
	my ($self, $patient_view, $min_percent) = @_;
	confess();
}

sub hash_in_this_run_patients {
	my($self,$ratio) = @_;
	$ratio = 0 unless $ratio;
	return $self->{inthisrun}->{$ratio} if exists $self->{inthisrun}->{ratio};
	$self->{inthisrun}->{$ratio} = undef;
	foreach my $patient (@{$self->getPatients()}) {
		my $fam_name = $patient->getFamily->name();
		$self->{inthisrun}->{$ratio}->{$fam_name}->{$patient->name()} = $ratio  if $self->get_percent_new_count($patient) >= $ratio;
	}
	die($ratio) unless exists $self->{inthisrun}->{$ratio};
	return $self->{inthisrun}->{$ratio};
}

sub in_this_run_patients {
	my($self,$ratio,$patient) = @_;
	$ratio =0 unless $ratio;
	my $hash = $self->hash_in_this_run_patients($ratio);
	if ($patient){
		delete $hash->{$patient->getFamily->name()} if exists $hash->{$patient->getFamily->name()};
	}
	my $nb = scalar (keys %$hash);
	return $nb;	
}
sub in_this_run_ratio {
	my($self,$ratio,$patient) = @_;
	my $a = $self->in_this_run_patients($ratio,$patient);
	warn $a;
	my $nb = scalar(@{$self->project->getPatients});
	return  ($a / $nb);
}

sub dejavu_patients {
	my($self,$ratio,$patient) = @_;
	$ratio = "all" unless $ratio;
	return 0 unless exists $self->dejavu->{$ratio};
	my $nb = 0;
	my $h_details = $self->dejavu_details_by_patnames();
	foreach my $pat_name (keys %{$self->dejavu->{$ratio}}) {
		next if $patient and $patient->name() eq $pat_name;
		my $proj_name = $h_details->{$pat_name};
		next if $proj_name eq $self->getProject->name();
		$nb++;
	}
	return $nb;
}

sub dejavu {
	my($self) =@_;
	return $self->{dejavu} if exists $self->{dejavu};
	$self->{dejavu} = $self->project->dejavuJunctionsResume()->dejavu($self);
	return $self->{dejavu};
}

sub dejavu_details   {
	my($self) =@_;
	my $h = $self->dejavu();
	return {} unless $h;
	my @l_pat = split(';', $h->{details});
	my $res ={};
	foreach my $patinfos (@l_pat) {
		my @lTmp = split('_', $patinfos);
		my $project_name_dv = 'NGS20'.shift(@lTmp).'_'.shift(@lTmp);
		my $patient_name_dv = join('_', @lTmp);
		push(@{$res->{$project_name_dv}},$patient_name_dv);
	}
	return $res;
}

sub dejavu_details_by_patnames   {
	my($self) =@_;
	my $h = $self->dejavu();
	return {} unless $h;
	my @l_pat = split(';', $h->{details});
	my $res ={};
	foreach my $patinfos (@l_pat) {
		my @lTmp = split('_', $patinfos);
		my $project_name_dv = 'NGS20'.shift(@lTmp).'_'.shift(@lTmp);
		my $patient_name_dv = join('_', @lTmp);
		$res->{$patient_name_dv} = $project_name_dv;
	}
	return $res;
}

sub getCoverageInInterval {
	my ($self, $patient) = @_;
	return sprintf("%.2f",$self->get_coverage($patient)->coverage($self->start() +1, $self->end() -1)->{mean});
}

sub getHashSpliceAiNearStartEnd {
	my ($self) = @_;
	my $h;
	return $h if $self->isCanonique();
	my ($max_start_score, $max_start_infos, $h_start_details) = $self->getHashSpliceAiInInterval($self->start() - 10, $self->start() + 10);
	my ($max_end_score, $max_end_infos, $h_end_details) = $self->getHashSpliceAiInInterval($self->end() - 10, $self->end() + 10);
	if ($max_start_score > 0) {
		$h->{start}->{max_score} = $max_start_score;
		$h->{start}->{max_infos} = $max_start_infos;
		$h->{start}->{all} = $h_start_details;
	}
	if ($max_end_score > 0) {
		$h->{end}->{max_score} = $max_end_score;
		$h->{end}->{max_infos} = $max_end_infos;
		$h->{end}->{all} = $h_end_details;
	}
	return $h;
}

sub getHashSpliceAiInInterval {
	my ($self, $start, $end) = @_;
	my @lPositions = ($start..$end);
	my $i = $start;
	my $max_score = 0;
	my $max_score_infos = '';
	
	my $b = new GBuffer;
	my $p = $b->newProject( -name => $self->getProject->name());
	my $chr = $p->getChromosome($self->getChromosome->id());

	my $h;
	foreach my $pos (@lPositions) {
		my $h2 = $chr->get_lmdb_spliceAI()->get($i);
		foreach my $alt (keys %$h2) {
			foreach my $gene_name (keys %{$h2->{$alt}}) {
				my @data = unpack( "W4 C4", $h2->{$alt}->{$gene_name});
				my @type = ( "AG", "AL", "DG", "DL" );
				for ( my $j = 0 ; $j < 4 ; $j++ ) {
					my $score = sprintf("%.2f", $data[$j] / 100);
					next if $score == 0;
					$h->{$pos}->{$alt}->{$gene_name}->{$type[$j]} = $score;
					if ($score > $max_score) {
						$max_score = $score;
						$max_score_infos = $chr->id().':'.$pos.';'.$alt.';'.$type[$j].':'.$score;
					}
				}
			}
		}
	}
	return ($max_score, $max_score_infos, $h);
}


1;
