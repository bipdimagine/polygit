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
	
);

my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name);

my $variant= $project->_newVariant($vid);
warn $variant->name;
warn Dumper $variant->getChromosome->score_gene_spliceAI( $variant->spliceAI(),"COL4A5");
