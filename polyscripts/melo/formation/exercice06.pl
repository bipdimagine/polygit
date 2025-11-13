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
	$vec_pat1->And($vec_pat1, $vec_pat2);
	
	my $variants = $chr->getListVarObjects($vec_pat1);
	my $var_denovo;
	foreach my $var (@$variants) {
		if ($var->isDenovoTransmission($fam1,$pat1) and $var->isDenovoTransmission($fam2,$pat2)) {
			print "\t" . $var->name . "\n";
			push(@$var_denovo, $var->name);
		}
	}
	print scalar(@$var_denovo) . "\n";
	
}



















	


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
#		my $vec_pat = $pat->getVariantsVector($chr)->Clone;
#		
#		#my $vect_chr = $chr->getVariantsVector->Clone();
#		
#		my $variants = $chr->getListVarObjects($vec_pat);
#		my $var_denovo;
#		foreach my $var (@$variants) {
#			if ($var->isDenovoTransmission($fam,$pat)) {
##				print "\t" . $var->name . "\n";
#				push(@$var_denovo, $var->name);
#			}
#		}
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
