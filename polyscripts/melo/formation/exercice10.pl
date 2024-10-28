#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
use Data::Dumper;
use Getopt::Long;
use Carp;

use GBuffer;

my $buffer = new GBuffer;

my $project = $buffer->newProjectCache(-name=>'NGS2016_1231');
my $patients = $project->getPatients;
#foreach my $pat (@$patients) {print $pat->name ."\n"}
# 107589, 108317, 108004, 107884, 107504, 107628, 107598, 107567,
# 108148, 108325, 107649, 108035, 107878, 107595, 108007, 108372,
# 108143, 107880, 108330, 107624, 107566, 116283, 108179, 107974

#my $pat = $patients->[0];
#my $pat = $project->getPatient('107589');

my $i = 0;
foreach my $pat (@$patients) {
	print "Patient: ". $pat->name ."\n";
	
#	my $variations = $pat->getVariations;
#	print scalar @$variations ." variations\t";
	
	my $genes = $pat->getGenes;
	print scalar @$genes ." genes\n";
	
	foreach my $gene (@$genes) {
#	my $gene = $genes->[0];
		print $gene->name . "\t" . $gene->external_name . "\tchr" . $gene->chromosome->name ."\n";
		
#		print $pat->minimum($gene) . "\t" . $pat->mean($gene) ."\n";
	
		my $chr_name = $gene->chromosome->name;
		
#		print "\t";
		print $pat->minDepth($chr_name, $gene->start, $gene->end) . "\t";
		print $pat->maxDepth($chr_name, $gene->start, $gene->end) . "\t";
		print $pat->meanDepth($chr_name, $gene->start, $gene->end) ."\n";
		
		my $exons = $gene->getExons;
		foreach my $exon (@$exons) {
			print "\t". $exon->name .":\t";
			my $start = $exon->start;
			my $end = $exon->end;
			print $pat->minDepth($chr_name, $start, $end) . "\t";
			print $pat->maxDepth($chr_name, $start, $end) . "\t";
			print $pat->meanDepth($chr_name, $start, $end) ."\n";
		}
#		my $introns = $gene->getIntrons;
#		foreach my $intron (@$introns) {
#			print "\t". $intron->name ."\n";
#			print "\t";
#			my $start = $intron->start;
#			my $end = $intron->end;
#			print $pat->minDepth($chr_name, $start, $end) . "\t";
#			print $pat->maxDepth($chr_name, $start, $end) . "\t";
#			print $pat->meanDepth($chr_name, $start, $end) ."\n";
#		}
		print "\n";
	}
	print "\n\n";
	
	$i++;
	last if $i == 3;
}

print $project->getCapture->analyse ."\n";








# Sur le projet NGS2016_1231, un patient au choix
# pour chaque gène (le faire déjà pour un) ayant au moins 1 variant
# calculer la couverture minimale en reads
# calculer la couverture maximale en reads
# calculer la couverture moyenne en reads
