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
my $project = $buffer->newProject(-name => $project_name);


my $chromosomes = $project->getChromosomes; # liste de GenBoChromosome
foreach my $chr (@$chromosomes) {
	my $chr_name = $chr->name; # name = 1, fasta_name = chr1
	
	next if $chr_name ne '21';
	
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
			print $variant->variationTypeInterface($var_gene) . "\t";
			print "Frequency: " . $variant->frequency . "\n";
			
			my $transcripts = $var_gene->getTranscripts;
			foreach my $transcript (@$transcripts) {
				print "\t\t" . $transcript->name . "\t";
				print $variant->variationTypeInterface($transcript) . "\n";
			}
		}
	}
	print "\n";
}











# A partir du projet NGS2016_1231 (attention, bien ce projet sinon cela va être
# méga long et le serveur risque de ne pas aimer)
# Pour tous les patients, faire une boucle sur chaque chromosome,
# donner les gènes et les variants
# variationTypeInterface
# fréquences
# transcripts + variationTypeInterface

