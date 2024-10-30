#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use Data::Dumper;
use Getopt::Long;
use Carp;

use GBuffer;

my $buffer = new GBuffer;

my $query = $buffer->getQuery;

my $projects = $query->listProjects;



my $captures;

my $i = 0;
foreach my $project_name (@$projects) {
	
	my $project = $buffer->newProject(-name=>$project_name);
#	print $project->name ."\n";
	
	my $capture = $project->getCaptures;
	foreach my $capt (@$capture) {
#		print ref($capt) . "\n";
		my $capture_name = $capt->name;
		$captures->{$capture_name}->{projects}++;
		}

	
	my $patients = $project->getPatients;
	foreach my $pat (@$patients) {
		my $capture_name = $pat->getCapture->name;
		$captures->{$capture_name}->{patients}++;
	}
	
	$i++;
	last if $i == 100;
}

print Dumper keys %{$captures};

if (exists $captures->{ciliomev3}) { print "ciliomev3 OK \n"; }
else { print "ciliomev3 NOT \n"; }

if (exists $captures->{ciliomev1545}) { print "ciliomev1545 OK \n"; }
else { print "ciliomev1545 NOT \n"; }

print "\n";



































# Vérifier chaque projet, si l'espèce est bien la même pour chaque patient
# Ceux qui ont un souci, peux tu sortir le numero de projet (NGS20...) et sa description stp ?

# tu peux faire un getCaptures pour chaque patient.
# tu auras un objet GenBoCapture pour chaque capture du patient (capture = liste de genes)
# -> chaque capture possede un nom et donc une methode name
# Et peux tu sortir les projets qui ont des patients avec des captures différentes sur les projets stp ?

# parmis ceux là, combien ont une incoherence exome / genome / capture.
# C'est a dire un exome avec une capture ou un exome avec un genome ?
