#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../polypipeline/packages/";
use lib "$Bin/../../polypipeline/dragen/scripts/";
use dragen_util; 
use Data::Dumper;
use Getopt::Long;
use Carp;
use IO::Prompt;
use Net::SFTP;
use File::Find::Rule ;
use autodie;
use file_util;


use GBuffer;
my $buffer = new GBuffer;

my $project_name; # Test ATACseq : NGS2025_09521
my $patient_name;
my $gene_name;
my $cpu = 20;
GetOptions(
	'project=s'		=> \$project_name,
	'patients=s'	=> \$patient_name,
	'gene_name=s'	=> \$gene_name,
	'cpu=i'			=> \$cpu,
) || die("Error in command line");


my $project = $buffer->newProjectCache( -name => $project_name );
my $patients = $project->getPatients();
#my $patient = $project->getPatient($patient_name);
#my $pname = $patient->name;


my $gene = $project->newGene($gene_name);

warn ref $gene;
warn $gene->external_name;

warn $project->getGenBoId($gene_name);


