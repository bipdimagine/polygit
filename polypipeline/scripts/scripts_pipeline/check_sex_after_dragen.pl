#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use Data::Dumper;
use Getopt::Long; # librairie pour récuper les options depuis le terminal
use Carp; # librairie pour `confess`

use GBuffer;

my $project_name;
my $patient_name;
GetOptions(
	'project=s' => \$project_name,
	'patient=s' => \$patient_name,
);

confess("ERROR: -project mandatory. Die.\n") unless ($project_name);
confess("ERROR: -patient mandatory. Die.\n") unless ($patient_name);


my $buffer = new GBuffer;

# Vérifie que le project existe
my $project;
eval{$project = $buffer->newProject(-name => $project_name)} //
	confess("ERROR: Could not find project '$project_name'\n");

# Vérifie que le patient existe
my $patient;
eval{$patient = $project->getPatient($patient_name)} //
	confess("ERROR: Could not find patient '$patient_name' in project '$project_name'\n");


my $dragen_pipeline_dir = $project->project_dragen_pipeline_path;

#TODO: a enlever
if ($project_name eq 'NGS2019_2429' and $dragen_pipeline_dir =~ m/\/HG19\/$/) {
	$dragen_pipeline_dir =~ s/HG19/HG38_CNG/;
} 

$dragen_pipeline_dir .= 'pipeline/';

if (not -d $dragen_pipeline_dir . "$patient_name/") {
	warn ("No directory found in dragen pipeline for patient '$patient_name': $dragen_pipeline_dir$patient_name/ does not exists");
	next;
}
my $ploidy_estimation_file = $dragen_pipeline_dir . "$patient_name/$patient_name.ploidy_estimation_metrics.csv";
my $ploidy_estimation_file2 = $project->getAlignmentStatsDir."$patient_name.ploidy_estimation_metrics.csv";
$ploidy_estimation_file = $ploidy_estimation_file2 if (-e $ploidy_estimation_file2 and not -e $ploidy_estimation_file);

my $output_path = $project->project_path;

#TODO: a enlever
if ($project_name eq 'NGS2019_2429' and $output_path =~ m/\/HG19\/$/) {
	$output_path =~ s/HG19/HG38_CNG/;
}

$output_path .= "sex_control/";
if (not -d $output_path) {mkdir $output_path}

open (my $output, '>', $output_path . "$patient_name.sex_control.log");
print {$output} "#Patient\tSex\tSex chromosomes\tMatching\n";
print {$output} $patient_name ."\t";


# Patient sex
my $patient_sex = $patient->sex;
print {$output} $patient_sex ."\t";


# Estimated sex chr
open(my $file, '<', $ploidy_estimation_file) or confess ("ERROR: Could not open file '$ploidy_estimation_file': $!.\n");
my $sex_chr;
while (my $line = readline($file)) {
	chomp $line;
	$sex_chr = substr($line, -2); # sur la dernière ligne
}
close($file);
print {$output} $sex_chr ."\t";


# Comparison
if (($patient_sex == 1 and $sex_chr eq 'XY') or ($patient_sex == 2 and $sex_chr eq 'XX')) {
	print {$output} "1";
}
elsif (not $sex_chr =~ /X[XY]/) {
	warn "Anormal ploidy on sex chromosomes for patient '$patient_name': $sex_chr\n"
}
else {
	print {$output} "0";
	confess "Estimated and patient sex not matching for patient '$patient_name': sex: $patient_sex, sex chr: $sex_chr\n"
}


close ($output);













