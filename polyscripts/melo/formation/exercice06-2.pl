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




my $pat1 = $project->getPatient('GUI_REF');
my $fam1 = $pat1->getFamily;
my $pat2 = $project->getPatient('LEF_MEL');
my $fam2 = $pat2->getFamily;

my $chromosomes = $project->getChromosomes;
foreach my $chr (@$chromosomes) {
	next unless ($chr->name eq '1');
	print "chr " . $chr->name ."\n";

	my $vec_pat1 = $pat1->getVariantsVector($chr)->Clone;
	my $vec_pat2 = $pat2->getVariantsVector($chr)->Clone;
	
	my $vec_denovo_fam1 = $fam1->getModelVector_fam_denovo($chr)->Clone;
	my $vec_denovo_fam2 = $fam2->getModelVector_fam_denovo($chr)->Clone;

	$vec_pat1->And($vec_pat1, $vec_pat2);
	$vec_denovo_fam1->And($vec_denovo_fam1, $vec_denovo_fam2);
	$vec_pat1->And($vec_pat1, $vec_denovo_fam1);

	print $chr->countThisVariants($vec_pat1) . "\n";
	
	my $var_denovo = $chr->getListVarObjects($vec_pat1);
#	foreach my $var (@$var_denovo) {
#		print "\t" . $var->name . "\n";
#	}
}







# Résultats différents car il y a deux enfants dans la famille GUI, dont 1 sain.












#for my $pat_name ('GUI_REF', 'LEF_MEL') {
#	
#	my $pat = $project->getPatient($pat_name);
#	my $fam = $pat->getFamily;
#	print $fam->name . "\n";
#	
#	my $chromosomes = $project->getChromosomes;
#	foreach my $chr (@$chromosomes) {
#		next unless ($chr->name eq '1');
#		print "chr " . $chr->name ."\n";
#		
#		my $vec_pat = $pat->getVariantsVector($chr)->Clone;
#		my $vec_denovo_fam = $fam->getModelVector_fam_denovo($chr)->Clone;
#		$vec_pat->And($vec_pat, $vec_denovo_fam);
#		
#		my $var_denovo = $chr->getListVarObjects($vec_pat);
#		
##		my $variants = $chr->getListVarObjects($vec_pat);
##		foreach my $var (@$variants) {
##			if ($var->isDenovoTransmission($fam,$pat)) {
###				print "\t" . $var->name . "\n";
##				push(@$var_denovo, $var->name);
##			}
##		}
#		print scalar(@$var_denovo) . "\n";
#	}
#}








# Projet: NGS2015_0794
# Type: famille
# Model: denovo => 2 méthodes: variant->isDenovoTransmission, ou famille->getModelVector_fam_denovo
#
# DB public ≤ 1%, DejaVu ≤ 30
# Impact factor: Medium => Mature miRNA, Splice Acc/Don, Frameshift, Stop-gained,
#	 (start/Stop)-lost, Splice Region, Predicted Splice Region, Missense, No-framshift
# 
# 
# Avec les objets Caches et vecteur pour le nombre de variants et pour les obtenirs.
