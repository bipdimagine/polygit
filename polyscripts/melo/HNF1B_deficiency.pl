#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
#print "$RealBin\n";
use lib "$RealBin";
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
use lib "$RealBin/../../polypipeline/dragen/scripts/";
use lib "$RealBin/../../polypipeline/packages/";
use dragen_util; 
use file_util; 
use Data::Dumper;
use Getopt::Long;
use Carp;
use IO::Prompt;
use Net::SFTP;
use File::Find::Rule ;

use GBuffer;
my $buffer = new GBuffer;

my $project_name = 'NGS2023_7054';
#my $patient_name;
GetOptions(
#	'project=s'		=> \$project_name,		# Projet
#	'patient=s'		=> \$patient_name,
);

my $project = $buffer->newProjectCache(-name=>$project_name) or confess ("can\'t open project cache '$project_name'");
warn $project->name;

my $patient = $project->getPatient('23A2_PATIENT_Het');
my $brother = $project->getPatient('23A3_BROTH_Het');
my $mother = $project->getPatient('23A1_MOTH_Het');
my $father = $project->getPatient('23NA1_FATH_WT');
my $fam = $patient->getFamily;

my $chromosomes = $project->getChromosomes;
foreach my $chr (@$chromosomes) {
	print "chr".$chr->name."\t";
	my $vec_pat = $patient->getVariantsVector($chr)->Clone;
	my $vec_brother = $brother->getVariantsVector($chr)->Clone;
	my $vec_father = $father->getVariantsVector($chr)->Clone;
	my $vec_mother = $mother->getVariantsVector($chr)->Clone;
	
	$vec_pat->And($vec_pat, $vec_brother);
#	$vec_pat->And($vec_pat, $vec_father);
	# todo: exclude mother
	
	my $variants = $chr->getListVarObjects($vec_pat);
	my $filtered;
	foreach my $var (@$variants) {
		if (
#			$var->isCompoundTransmission($fam,$patient,$var->getGene,$father)
#			$var->is_familial_recessif_compound($fam,$patient)
#			$var->is_familial_recessif($fam,$patient)
			$var->isFatherTransmission($fam,$patient)
			and $var->getGnomadAC <= 0.01
			) {
			push(@$filtered, $var);
		}
	}
	print scalar @$filtered."\n";
	
}






















































































