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
print scalar(@$projects) . " projects\n";

my ($nb_genome, $nb_exome, $nb_ciliome, $nb_rna) = (0,0,0,0);
my $nb_other;
foreach my $project_name (@$projects) {
#	print $project_name . "\t";
	
	my $project = $buffer->newProject(-name=>$project_name);
#	print $project->description . "\t";

	if ($project->isGenome) {
#		print "Genome";
		$nb_genome++;
	}
	elsif ($project->isExome) {
#		print "Exome";
		$nb_exome++;
	}
#	elsif ($project->isCiliome) {
##		print "Ciliome";
#		$nb_ciliome++;
#	}
	elsif ($project->isRnaSeq) {
#		print "RNAseq";
		$nb_rna++;
	}
	else {
		foreach my $c ( @{ $project->getCaptures } ) {
			my $analyse = $c->analyse;
#			print $analyse;
			$nb_other->{$analyse} ++;
		}
		$nb_other->{'total'} ++;
	}
#	print "\n";
}

print $nb_genome ." genome,\t";
print $nb_exome ." exome,\t";
#print $nb_ciliome ." ciliome,\t";
print $nb_rna ." RNA,\n";

print "Other:\n";

my @nb_other_sorted = sort { $nb_other->{$b} <=> $nb_other->{$a} } keys %$nb_other;
foreach my $cle (@nb_other_sorted) {
    print "\t$cle => $nb_other->{$cle}\n";
}
#print Dumper $nb_other;








# Dans QueryMoose.pm, utiliser list(All)Projects pour avoir tous les projets,
# récuper leur description et donner leur type (genome, capture, exome)

# 7588 projects: 1633 genome,	1633 exome,	410 RNA,	3912 other