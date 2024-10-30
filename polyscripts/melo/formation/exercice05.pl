#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use Data::Dumper;
use Getopt::Long;
use Carp;


use GBuffer;


my $project_name = 'NGS2015_0794';

my $buffer = new GBuffer;
my $project = $buffer->newProjectCache(-name => $project_name);

my $pat = $project->getPatient('GUI_REF');

my $gene = $project->newGene('POLR3A');

my $chr = $gene->getChromosome;
print ref($chr) . "\t" . $chr->name . "\n";

my $vector_gene = $gene->getVariantsVector->Clone();
print ref($vector_gene) . "\n";

my $vector_pat = $pat->getVariantsVector($chr);
print ref($vector_pat) . "\n";

my $vector_variants = $gene->getVariantsVector;
$vector_variants->Intersection($vector_gene, $vector_pat);
print $chr->countThisVariants($vector_variants) . "\n";

my $variants = $chr->getListVarObjects($vector_variants);
foreach my $var (@$variants) {
	print $var->id . "\n";
	my $transcripts = $var->getTranscripts;

	foreach my $transcript (@$transcripts) {
		print "\t" . $transcript->id . "\t";
		print $var->variationTypeInterface($transcript) . "\n";
	}
#	print "Gnomad AC: " . $var->gnomad->{'populations'}->{'all'}->{'AC'} . "\n";
	print "Gnomad AC: " . $var->getGnomadAC . "\tfrequency: " . $var->frequency . "\n";
	
	print "DejaVu: " . keys(%{$var->deja_vu}) . "\n";
}













# Projet: NGS2015_0794
# Patient: GUI_REF
# Gene: POLR3A
# nombre de variants
# leur conséquences par transcrits
# leur frequences gNOMAD AC
# leur DejaVu
# Avec les objets Caches et vecteur pour le nombre de variants et pour les obtenirs.
