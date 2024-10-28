#!/usr/bin/perl
use strict;
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















