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
my $cmd;
my ($project_name, $patient_name,$vid);
my $file;

GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'cmd=s'  => \$cmd,
	'variant=s'  => \$vid,
	'patient=s'  => \$patient_name,
);

my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name);
my $patient= $project->getPatient($patient_name);
my $chr = $project->getChromosome(21);
#my $vs  = $patient->getCnvs();
my $vs  = $chr->getBoundaries();
warn scalar(@$vs);
foreach my $v (@$vs){
#	next unless $v->isLarge;
	next if $v->isLargeDeletion;
	next if $v->isLargeDuplication;
	warn $v->id;
	warn $v->name();
	my $gs = $v->getTranscripts();
	#foreach my $g (@$gs){
	#	warn "\t".$g->name."\t".$v->variationType($g);
	#}
}
exit(0);

