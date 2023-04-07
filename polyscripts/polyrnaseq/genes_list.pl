#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
use GBuffer;
use Set::IntSpan::Fast::XS;
use Getopt::Long;
my $fork = 1;

my $project_name;
my $patient_name;
my $dir_out;
GetOptions(
	'project=s'    => \$project_name,
	'patient=s'  => \$patient_name,
	'dir=s'  => \$dir_out,
);
my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name);


my $capture = $project->getCapture();
warn $project->gencode_version();
		my $genes;
foreach my $tr (@{ $capture->transcripts_name() }) {
		#warn $tr;
	#	next if $tr =~ /enh/;
			my $t = $project->newTranscript($tr);
			
			$genes->{$t->getGene->id} = $t->getGene
}

foreach my $g (values(%$genes)){
	
	foreach my $t (@{$g->getTranscripts}){
		my $level =1;
		$level = 2 if exists  $t->{tag}->{CCDS};
		$level = 2 if exists  $t->{tag}->{basic};
		$level = 2 if exists  $t->{tag}->{appris_principal};
		foreach my $exon (@{$t->getExons}){
			
			my $strand = "+";
			$strand = "-" if $g->strand == -1;
			print $g->getChromosome->name."\t".$g->name."\t".$t->name."\t".$t->name.$exon->name."\t".$exon->start."\t".$exon->end."\t".$g->external_name."\t"."$level\n";
			
		}
		
	}
	
}