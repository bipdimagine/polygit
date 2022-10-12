package GenBoJunction;

use strict;
use Moose;
use MooseX::Method::Signatures;
use Data::Dumper;
use Config::Std;
use GenBoCapture;
use Position;
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
		my $ensid = $self->annex->{$lPat[0]->name()}->{ensid};
		my $h_exons_introns;
		my (@lExons, @lIntrons);
		my $intspan_junction = $self->getStrictGenomicSpan();
		my $i = 0;
		foreach my $g (@{$self->getGenes()}) {
			my @ltmp = split('_', $g->id());
			next if ($ltmp[0] ne $ensid);
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

sub isCanonique {
	my ($self, $patient) = @_;
	confess() unless ($patient);
	my $type = lc($self->getTypeDescription($patient)); 
	return 1 if ($type =~ /canonique/);
	return;
}

sub getTypeDescription {
	my ($self, $patient) = @_;
	confess() unless ($patient);
	return $self->annex->{$patient->name()}->{type}; 
}

sub isRI {
	my ($self, $patient) = @_;
	return 1 if $self->getTypeDescription($patient) =~ /RI/; 
	return 1 if $self->annex->{$patient->name}->{type_origin_file} eq 'RI';
	return;
}

sub isSE {
	my ($self, $patient) = @_;
	return 1 if $self->getTypeDescription($patient) =~ /SE/; 
	return 1 if $self->annex->{$patient->name}->{type_origin_file} eq 'SE';
	return;
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
	my $count = 0;
	$count += $self->annex->{$patient->name()}->{junc_ri_count} if ($self->isRI($patient) and exists $self->annex->{$patient->name()}->{junc_ri_count});
	$count += $self->annex->{$patient->name()}->{junc_se_count} if ($self->isSE($patient) and exists $self->annex->{$patient->name()}->{junc_se_count});
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
	my $new_count = $self->get_nb_new_count($patient);
	$new_count = 0 if $self->get_nb_new_count($patient) eq '---';
	my $normal_count = $self->get_nb_normal_count($patient);
	$normal_count = 0 if $self->get_nb_normal_count($patient) eq '---';
	return ($new_count + $normal_count);
}

sub get_ratio_new_count {
	my ($self, $patient) = @_;
	confess() unless $patient;
	my $ratio = $self->get_nb_new_count($patient) / $self->get_dp_count($patient);
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
	$cmd .= " -b $bam_path";
	$cmd .= " -c chr".$locus;
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
	else { $locus_text = $self->getChromosome->id().':'.($self->start() - 100).'-'.($self->end() + 100); }
	$locus_text =~ s/chr//;
	$locus_text =~ s/:/-/;
	my $path = $patient->getProject->getProjectPath.'/align/sashimi_plots/';
	unless (-d $path) {
		mkdir $path;
		`chmod 775 $path`;
	}
	my $outfile = $path.'/sashimi_'.$patient->name().'.'.$self->annex->{$patient->name()}->{'ensid'}.'.'.$locus_text.'.svg';
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

sub get_dejavu_list_similar_junctions {
	my ($self, $identity) = @_;
	$identity = $self->dejavu_percent_coordinate_similar() unless $identity;
	my $chr_id = $self->getChromosome->id();
	my $type = 'all';
	confess() unless $type;
	return $self->getProject->dejavuJunctions->get_junction($chr_id, $self->start(), $self->end(), $type, 'all', $identity);
}

sub get_dejavu_list_similar_junctions_resume {
	my ($self, $identity) = @_;
	$identity = $self->dejavu_percent_coordinate_similar() unless $identity;
	my $chr_id = $self->getChromosome->id();
	my $type = 'all';
	confess() unless $type;
	return $self->getProject->dejavuJunctions->get_junction_resume($chr_id, $self->start(), $self->end(), $type, 'all', $identity);
}

has parse_nb_projects_patients => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return if (not $self->getProject->annotation_genome_version() =~ /HG19/);
		my ($h_proj, $h_pat, $hpatrun);
		my $similar = $self->dejavu_percent_coordinate_similar();
		foreach my $h (@{$self->get_dejavu_list_similar_junctions_resume($similar)}) {
			my @list = split(';', $h->{projects});
			my @list2 = split(';', $h->{patients});
			my $i = 0;
			foreach my $proj (@list) {
				$h_proj->{$proj} = undef;
				foreach my $pat (split(',', $list2[$i])) {
					$h_pat->{$proj.'_'.$pat} = undef;
					#if ($h->{same_as} eq '100%' and 'NGS20'.$proj eq $self->getProject->name()) {
					if ('NGS20'.$proj eq $self->getProject->name()) {
						$hpatrun->{$proj.'_'.$pat} = undef;
					}
				}
				$i++;
			}
		}
		my $h;
		$h->{projects} = scalar(keys %$h_proj);
		$h->{patients} = scalar(keys %$h_pat);
		#$h->{patients_inthisrun} = scalar(keys %$hpatrun);
		$h->{others_patients} = scalar(keys %$h_pat) - scalar(keys %$hpatrun);
		#$h->{others_patients} = scalar(keys %$h_pat);
		return $h;
	},
);

has dejavu_nb_projects => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return if (not $self->getProject->annotation_genome_version() =~ /HG19/);
		return $self->parse_nb_projects_patients->{projects};
	},
);

has dejavu_nb_patients => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return if (not $self->getProject->annotation_genome_version() =~ /HG19/);
		return $self->parse_nb_projects_patients->{patients};
	},
);

has dejavu_nb_others_patients => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return if (not $self->getProject->annotation_genome_version() =~ /HG19/);
		return $self->parse_nb_projects_patients->{others_patients};
	},
);

sub dejavu_nb_int_this_run_patients {
	my ($self, $min_percent) = @_;
	my $nb = scalar(@{$self->getPatients()});
	return $nb unless ($min_percent);
	foreach my $patient (@{$self->getPatients()}) {
		$nb-- if ($self->get_percent_new_count($patient) < $min_percent);
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
	my $h_noise;
	my $h_exons_introns = $self->get_hash_exons_introns();
	return $h_noise unless $h_exons_introns;
	my @lJunctions;
	my $this_chr = $self->getChromosome();
	my $i = $self->vector_id();
	my $min = $self->vector_id() - 15;
	$min = 0 if ($min < 0);
	while ($i >= $min) {
		my $junction2 = $this_chr->getVarObject($i);
		push(@lJunctions, $junction2) if $junction2->get_hash_exons_introns();
		$i--;
	}
	$i = $self->vector_id();
	my $max = $self->vector_id() + 15;
	$max = $this_chr->size_vector() if ($max >= $this_chr->size_vector());
	while ($i < $max) {
		my $junction2 = $this_chr->getVarObject($i);
		push(@lJunctions, $junction2) if $junction2->get_hash_exons_introns();
		$i++;
	}
	foreach my $junction2 (@lJunctions) {
		foreach my $tid (keys %{$junction2->get_hash_exons_introns}) {
			next if $tid eq 'by_pos';
			next unless (exists $h_exons_introns->{$tid});
			my @lPos = (sort keys %{$junction2->get_hash_exons_introns->{$tid}->{by_pos}});
			next if (@lPos == 0);
			foreach my $pos (keys %{$junction2->get_hash_exons_introns->{$tid}->{by_pos}}) {
				if (exists $h_exons_introns->{$tid}->{by_pos}->{$pos}) {
					$h_noise->{$tid}->{$pos}++;
				}
			}
		}
	}
	return $h_noise;
}

sub junction_score {
	my ($self, $patient) = @_;
	my $gene;
	foreach my $g (@{$self->getGenes()}) {
		my @ltmp = split('_', $g->id());
		next if ($ltmp[0] ne $self->annex->{$patient->name()}->{ensid});
		$gene = $g;
	}
	confess() unless $gene;
	my $score = $gene->score();
	
	#ratio
	my $ratio = $self->get_percent_new_count($patient);
	if ($ratio >= 30) { $score += 2; }
	elsif ($ratio >= 20) { $score += 1.5; }
	elsif ($ratio >= 10) { $score += 1; }
	
	#dp
	my $dp = $self->get_dp_count($patient);
	if ($dp >= 75) { $score += 2; }
	elsif ($dp >= 50) { $score += 1.5; }
	elsif ($dp >= 30) { $score += 1; }
	elsif ($dp <= 10) { $score -= 1; }
	elsif ($dp <= 5) { $score -= 3; }
	
	#nb new
	my $nb_new = $self->get_nb_new_count($patient);
	if ($dp >= 15) { $score += 1; }
	elsif ($dp <= 5) { $score -= 1; }
	elsif ($dp <= 2) { $score -= 3; }
	
	#noise
	my $h_noise = $self->get_hash_noise($patient);
	my $max_noise = 0;
	foreach my $enst (keys %{$h_noise}) {
		foreach my $pos (keys %{$h_noise->{$enst}}) {
			my $this_noise = $h_noise->{$enst}->{$pos};
			$max_noise = $this_noise if $this_noise > $max_noise;
		}
	}
	if ($max_noise >= 20) { $score -= 3; }
	elsif ($dp >= 12) { $score -= 2; }
	elsif ($dp >= 5) { $score -= 1; }
	return $score;
}

1;