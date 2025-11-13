#!/usr/bin/perl
use strict;
<<<<<<< HEAD
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
=======
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
use Data::Dumper;
use Getopt::Long;
use Carp;
use IO::Prompt;
use Net::SFTP;
use File::Find::Rule ;

use GBuffer;
my $buffer = new GBuffer;

my $project_name;
my $patient_name;
GetOptions(
	'project=s'		=> \$project_name,		# Projet
	'patient=s'		=> \$patient_name,
);

$project_name = 'NGS2024_8239' unless ($project_name);
my $project = $buffer->newProject(-name=>$project_name) or confess ("can\'t open project '$project_name'");
#my $project = $buffer->newProjectCache(-name=>$project_name) or confess ("can\'t open project cache '$project_name'");
warn "Project: $project_name";

#my $query = $buffer->getQuery;

my $patients = $project->getPatients unless ($patient_name);
print scalar(@$patients) ." patients:\n";
foreach my $pat (@$patients) {
	my $pn = $pat->name();
 	print "$pn\n";
	my $dir1 = $project->getProjectPath();
	my @files = File::Find::Rule->file()
				->name("$pn"."_aberrations.bed.gz*", "$pn"."_bins.bed.gz*")
				->in($dir1);
	print Dumper @files;
	
}

#my $pat = $project->getPatients()->[0] unless ($patient_name);
#my $pat = $project->getPatient($patient_name) if ($patient_name);
#warn "Patient: ".$pat->name;















>>>>>>> branch 'dev-mperin' of https://github.com/bipdimagine/polygit.git
