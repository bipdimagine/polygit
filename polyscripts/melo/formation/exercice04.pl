#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use Data::Dumper;
use Getopt::Long;
use Carp;

use GBuffer;


my $project_name = 'NGS2016_1231';
my $patient_name = '107504';

GetOptions(
			'patient=s' => \$patient_name,
);

#confess("ERROR: -patient mandatory. Die.\n") unless ($patient_name);

my $buffer = new GBuffer;
my $project = $buffer->newProject(-name => $project_name);

my $patients = $project->getPatients;
my $patients_name;
#print scalar(@$patients). " patients:\n";
foreach my $patient (@$patients) {
	push(@$patients_name, $patient->name);
}
#print join(", ", @$patients_name) . "\n";

my $patient = $project->getPatient($patient_name);
print "Patient: " . $patient->name . "\n";


my $variants = $patient->getStructuralVariations;
print scalar @$variants . " variants\t";


my $variations = $patient->getVariations; # liste de GenBoVariation
my $indels = $patient->getIndels; # liste de GenBoInsertion, GenBoDeletion
my $cnvs = $patient->getCnvs; # liste de GenBoCnv

#my @variants = (@$variations, @$indels, @$cnvs);
#print scalar @variants . " variants ";
print " (" . scalar @$variations . " substitutions, ";
print scalar @$indels . " indels, ";
print scalar @$cnvs . " cnv):\n";

my $variants_name;
foreach my $variant (@$variants) {
	push(@$variants_name, $variant->name); # name -> avec des '-', id -> avec des '_'
}
print join("\n", @$variants_name) . "\n";







# A partir du projet NGS2016_1231 (attention, bien ce projet sinon cela va être
# méga long et le serveur risque de ne pas aimer)
# me donner tous les variants avec leur ID du patient de ton choix