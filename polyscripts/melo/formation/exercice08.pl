#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use Data::Dumper;
use Getopt::Long;
use Carp;
use List::MoreUtils qw(uniq);

use GBuffer;

my $buffer = new GBuffer;

my $query = $buffer->getQuery;

my $projects = $query->listProjects;


my @proj_multiple_captures;
my @proj_incoherent;

foreach my $project_name (@$projects) {
	my $project = $buffer->newProject(-name=>$project_name);
#	print $project->name ."\n";
#	print $project->species_id ."\t". $project->getVersion ."\t";
	
	my $patients = $project->getPatients;
	my $captures;
	foreach my $pat (@$patients) {
		my $capture_name = $pat->getCapture->name;
		$capture_name .= '_'.$pat->getCapture->analyse();
		push (@$captures, $capture_name);
#		print $pat->name ."\t". $capture_name ."\n";
	}
	my @capture_uniq = uniq(@$captures);
#	print @capture_uniq;
	if (scalar @capture_uniq  > 1) {
		push (@proj_multiple_captures, $project_name);
		
		print $project_name ."\t";
#		print $project->isGenome ."\t". $project->isExome ."\t";
		print join(', ', @capture_uniq) ."\n";
		
		my ($is_genome, $is_exome, $is_capture);
		foreach my $capture (@capture_uniq) {
			if (lc($capture) =~ /genome/) {$is_genome++}
			elsif (lc($capture) =~ /exome/) {$is_exome++}
			else {$is_capture++}
		}
		if (($is_capture and  $is_genome) or ($is_capture and  $is_exome) or ($is_genome and $is_exome)) {
			push (@proj_incoherent, $project_name);
		}
		
		
	}
}

print "\n";
print scalar(@proj_multiple_captures) . " projects have patients with different captures\n";
#print join("\n", @proj_multiple_captures);

print "\n";
print "Including " . scalar(@proj_incoherent) . " projects with genome/exome/capture incoherence:\n";
print join("\n", @proj_incoherent);

print "\n";



































# Vérifier chaque projet, si l'espèce est bien la même pour chaque patient
# Ceux qui ont un souci, peux tu sortir le numero de projet (NGS20...) et sa description stp ?

# tu peux faire un getCaptures pour chaque patient.
# tu auras un objet GenBoCapture pour chaque capture du patient (capture = liste de genes)
# -> chaque capture possede un nom et donc une methode name
# Et peux tu sortir les projets qui ont des patients avec des captures différentes sur les projets stp ?

# parmis ceux là, combien ont une incoherence exome / genome / capture.
# C'est a dire un exome avec une capture ou un exome avec un genome ?
