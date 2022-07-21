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

has isRI => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> 0,
);

has isSE => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> 0,
);

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

has score => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> 0,
);

has nb_new_count => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> 0,
);

has nb_normal_count => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> 0,
);

has annex => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> 0,
);


#TODO: a faire
has dejaVuInfosForDiag2 => (
	is              => 'rw',
	lazy    => 1,
	default => sub {
		return;
	}
);

#TODO: a faire
sub dejavu_all {
	my $self = shift;
	return '.';
}

#TODO: a faire
sub dejavu_in_this_run {
	my $self = shift;
	return '.';
}

sub getSvgPlotPath {
	my ($self, $patient) = @_;
	my $path_svg = $patient->getJunctionsAnalysePath().'/../'.$self->annex->{'ensid'}.'/SpliceRes/Figs/';
	my $svg_patient = $path_svg.'/juncPairPos_'.$self->annex->{ensid}.'_'.$self->getChromosome->id().'_'.$patient->name().'_'.$patient->getProject->name().'.svg'; 
	return $svg_patient if (-e $svg_patient);
	return;
}

sub createSashiPlot {
	my ($self, $patient, $locus) = @_;
	my $file = $self->getSashimiPlotPath($patient, $locus);
	return if -e $file;
	my $patient_name = $patient->name();
	my $path_analysis = $patient->getJunctionsAnalysePath();
	my $ensg = $self->annex->{'ensid'};
	my $score = $self->score();
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

sub createListSashimiPlots {
	my ($self, $patient) = @_;
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
	my $outfile = $path.'/sashimi_'.$patient->name().'.'.$self->annex->{'ensid'}.'.'.$locus_text.'.pdf';
	return $outfile;
}

sub getListSashimiPlotsPathFiles {
	my ($self, $patient, $to_create) = @_;
	my $locus = $self->getChromosome->id().':'.($self->start() - 100).'-'.($self->end() + 100);
	my $sashimi_file = $self->getSashimiPlotPath($patient, $locus);
	$self->createSashiPlot($patient, $locus) if ($to_create);
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
			$self->createSashiPlot($patient, $locus_extended) if ($to_create);
			push(@lFiles, $sashimi_plot_file);
			$i++;
		}
		return \@lFiles;
	}
	return;
}

1;