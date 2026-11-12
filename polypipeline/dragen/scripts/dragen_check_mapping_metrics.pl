#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use Data::Dumper;
use Getopt::Long;
use Carp;
use IO::Prompt;
use Term::ANSIColor;

use GBuffer;

my $project_name;
my $patient_name;
my $version;
my $max_dup; # 0-100%
my $min_map; # 0-100%
my $verbose;
GetOptions(
	'project=s'			=> \$project_name,
	'patient=s'			=> \$patient_name,
	'version=s'			=> \$version,
	'max_duplicates=i'	=> \$max_dup,
	'min_mapping=i'		=> \$min_map,
	'verbose'			=> \$verbose,
) || confess("Error in command line arguments");

confess("ERROR: -project mandatory.\n") unless ($project_name);
confess("ERROR: -patient mandatory.\n") unless ($patient_name);
confess("ERROR: -max_duplicates must be between 0 and 100.\n") if ($max_dup > 100 or $max_dup < 0);
confess("ERROR: -min_mapping must be between 0 and 100.\n") if ($min_map > 100 or $min_map < 0);

my $buffer = new GBuffer;
my $project = $buffer->newProject(-name => $project_name, -version =>$version);
my $patient = $project->getPatient($patient_name);
unless ($max_dup) {
	$max_dup = 20;
	$max_dup = 30 if ($project->isExome);
	$max_dup = 40 if ($project->isRnaseq);
	$max_dup = 50 if ($patient->getSampleProfile() =~ /epigenomic/);
	$max_dup = 60 if ($project->isDiagnostic);
}
unless ($min_map) {
	$min_map = 95;
}

# Mapping stats
my $dir_dragen_pipeline = $patient->getDragenDirName("pipeline");
my $mapping_metrics_file_1 = $dir_dragen_pipeline."/$patient_name.mapping_metrics.csv";
my $dir_stats_dragen = $patient->project->getAlignmentStatsDir("dragen-align");
my $mapping_metrics_file_2 = $dir_stats_dragen."/$patient_name.mapping_metrics.csv";
my $mapping_metrics_file = $mapping_metrics_file_2 if (-e $mapping_metrics_file_2);
$mapping_metrics_file = $mapping_metrics_file_1 if (-e $mapping_metrics_file_1);
confess("ERROR: No mapping metrics file found: '$mapping_metrics_file_1' or '$mapping_metrics_file_2'") unless (-e $mapping_metrics_file);
open(my $file, '<', $mapping_metrics_file) or confess ("ERROR: Could not open file '$mapping_metrics_file': $!");
my $duplicate;
my $mapping;
while (my $line = readline($file)) {
	next unless ($. == 2 or $. == 8);
	chomp $line;
	# l.2:MAPPING/ALIGNING SUMMARY,,Number of duplicate marked reads,4532637,4.86
	# l.8:MAPPING/ALIGNING SUMMARY,,Mapped reads,46106423,49.46
	if ($line =~ m{^MAPPING/ALIGNING SUMMARY,,Number of duplicate marked reads,\d+,([0-9.]+)$}) {
		$duplicate = $1;
	}
	elsif ($line =~ /^MAPPING\/ALIGNING SUMMARY,,Mapped reads,\d+,([0-9.]+)$/) {
		$mapping = $1;
	}
}
close($file);
confess("Error parsing '$mapping_metrics_file':/nno number of mapped reads and/or duplicate marked reads found.") unless ($duplicate and $mapping);
warn ("$patient_name: mapping = $mapping%, duplicates = $duplicate%") if ($verbose);


# Vérification de l'erreur critique (mapping trop bas) - prioritaire
if ($mapping < $min_map) {
    print "$mapping%\t$duplicate%";
    exit(2);  # Code 2 = Error (rouge)
}

# Vérification du warning (duplicates élevés)
if ($duplicate > $max_dup) {
    print "$mapping%\t$duplicate%";
    exit(1);  # Code 1 = Warning (orange/jaune)
}

# Tout est OK
print "$mapping%\t$duplicate%";
exit(0);  # Code 0 = Success (vert)

