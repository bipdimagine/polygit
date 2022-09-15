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
	$cmd .= " -b $path_analysis/../$ensg/rmdup/$patient_name\_$project_name\_rmdup.bam";
	$cmd .= " -c chr".$locus;
	$cmd .= " -o ".$file;
	if ($score and $score >= 100) { $cmd .= " -P /data-isilon/bipd-src/mbras/ggsashimi/ggsashimi-master/colors/red.txt"; }
	elsif ($score and $score >= 10) { $cmd .= " -P /data-isilon/bipd-src/mbras/ggsashimi/ggsashimi-master/colors/orange.txt"; }
	elsif ($score and $score >= 1) { $cmd .= " -P /data-isilon/bipd-src/mbras/ggsashimi/ggsashimi-master/colors/green.txt"; }
	else { $cmd .= " -P /data-isilon/bipd-src/mbras/ggsashimi/ggsashimi-master/colors/black.txt"; }
	$cmd .= " -C 1";
	$cmd .= " --shrink --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18";
	$cmd .= " -g ".$patient->getProject->get_gtf_genes_annotations_igv();
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
	my $outfile = $path.'/sashimi_'.$patient->name().'.'.$self->annex->{$patient->name()}->{'ensid'}.'.'.$locus_text.'.pdf';
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
		while ($i < 5) {
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

sub get_dejavu_list_similar_junctions {
	my ($self, $identity) = @_;
	$identity = 98 unless $identity;
	my $chr_id = $self->getChromosome->id();
	my $type = 'all';
	confess() unless $type;
	return $self->getProject->dejavuJunctions->get_junction($chr_id, $self->start(), $self->end(), $type, 'all', $identity);
}

sub get_dejavu_list_similar_junctions_resume {
	my ($self, $identity) = @_;
	$identity = 98 unless $identity;
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
		foreach my $h (@{$self->get_dejavu_list_similar_junctions_resume(98)}) {
			my @list = split(';', $h->{projects});
			my @list2 = split(';', $h->{patients});
			my $i = 0;
			foreach my $proj (@list) {
				$h_proj->{$proj} = undef;
				foreach my $pat (split(',', $list2[$i])) {
					$h_pat->{$proj.'_'.$pat} = undef;
					if ($h->{same_as} eq '100%' and 'NGS20'.$proj eq $self->getProject->name()) {
						$hpatrun->{$proj.'_'.$pat} = undef;
					}
				}
				$i++;
			}
		}
		my $h;
		$h->{projects} = scalar(keys %$h_proj);
		$h->{patients} = scalar(keys %$h_pat);
		$h->{patients_inthisrun} = scalar(keys %$hpatrun);
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

has inthisrun_nb_patients => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return if (not $self->getProject->annotation_genome_version() =~ /HG19/);
		return $self->parse_nb_projects_patients->{patients_inthisrun};
	},
);

1;