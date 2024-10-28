#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
use Data::Dumper;
use Getopt::Long;
use Carp;

use GBuffer;


my $project_name = 'NGS2015_0794';

my $buffer = new GBuffer;
my $project = $buffer->newProjectCache(-name => $project_name);


my $chr = $project->getChromosome('1');
my $pat = $project->getPatient('GUI_YOS');

my $fam = $pat->getFamily();

my $pat2 = $project->getPatient('GUI_REF');

my $gene = $project->newGene('NPHS2');

warn ref($project);
warn ref($chr);
warn ref($pat);
warn ref($gene);

my $vector = $chr->getVariantsVector();
warn ref($vector);
warn $chr->countThisVariants($vector);

my $vector_pat = $pat->getVariantsVector($chr);
warn ref($vector_pat);
warn $chr->countThisVariants($vector_pat);

my $vector_pat_2 = $pat2->getVariantsVector($chr);
warn ref($vector_pat_2);
warn $chr->countThisVariants($vector_pat_2);


print "Common variants pat1 et 2\n";
$vector_pat->Intersection($vector_pat, $vector_pat_2);
warn $chr->countThisVariants($vector_pat);

#print "Common variants pat1 et 2 et gene\n";
#$vector_pat->Intersection($vector_pat, $gene->getVariantsVector());
#warn $chr->countThisVariants($vector_pat);

foreach my $var (@{$chr->getListVarObjects($vector_pat)}) {
	#print $var->id()."\t".$var->variationTypeInterface()."\t".$var->getTransmissionModelType($fam, $pat2)."\n";
}

die;

my $chromosomes = $project->getChromosomes; # liste de GenBoChromosome
foreach my $chr (@$chromosomes) {
	my $chr_name = $chr->name; # name = 1, fasta_name = chr1
	
	
	print "chr " . $chr_name . ":\n";
	
	my $genes = $chr->getGenes; # liste de GenBoGene
	print "\t" . scalar(@$genes) . " genes\n";
	foreach my $gene (@$genes) {
#		print "\t" . $gene->name . "\n";
	}
	print "\n";
	
	my $variants = $chr->getStructuralVariations();
	print "\t" . scalar(@$variants) . " variants:\n";
	foreach my $variant (@$variants){
		
		my $var_genes = $variant->getGenes;
		foreach my $var_gene (@$var_genes) {
			print "\t" . $variant->name . "\t";
			print $var_gene->name . "\t";
			print $var_gene->external_name . "\t";
			print $variant->variationTypeInterface($var_gene) . "\n";
			
			
			warn ref($chr);
			warn ref($var_gene);
			warn ref($variant);
			die;
			
		}
	}
	print "\n";
}











# A partir du projet NGS2016_1231 (attention, bien ce projet sinon cela va être
# méga long et le serveur risque de ne pas aimer)
# Pour tous les patients, faire une boucle sur chaque chromosome,
# donner les gènes et les variants
# variationTypeInterface

