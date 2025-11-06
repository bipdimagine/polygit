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

die("ERROR: -project mandatory. Die.\n") unless ($project_name);
die("ERROR: -patient mandatory. Die.\n") unless ($patient_name);


my $buffer = new GBuffer;
my $project = $buffer->newProject(-name => $project_name);
my $patient = $project->getPatient($patient_name);
my $patient_sex = $patient->sex;

# Estimated sex chr
my $dragen_stats = $project->getAlignmentStatsDir("dragen-align");
my $ploidy_estimation_file = $dragen_stats."$patient_name.ploidy_estimation_metrics.csv";
open(my $file, '<', $ploidy_estimation_file) or die ("ERROR: Could not open file '$ploidy_estimation_file': $@.\n");
my $sex_chr;
while (my $line = readline($file)) {
	chomp $line;
	$sex_chr = substr($line, -2); # sur la dernière ligne
}
close($file);


my $output_path = $project->getAlignmentDir("dragen-align")."sex_control/";
mkdir $output_path unless (-d $output_path);
open (my $output, '>', $output_path . "$patient_name.sex_control.log");
print {$output} "#Patient\tSex\tSex chromosomes\tMatching\n";
print {$output} $patient_name ."\t". $patient_sex ."\t". $sex_chr ."\t";

# Comparison
if (($patient_sex == 1 and $sex_chr eq 'XY') or ($patient_sex == 2 and $sex_chr eq 'XX')) {
	print {$output} "1";
	my $tmp = $project->getAlignmentPipelineDir("dragen-align")."sex_control/";
#	system("touch $tmp$patient_name.ok")
}
elsif (not $sex_chr =~ /X[XY]/) {
	warn "Anormal ploidy on sex chromosomes for patient '$patient_name': $sex_chr\n"
}
else {
	print {$output} "0";
	die "Estimated sex and patient sex not matching for patient '$patient_name': sex: $patient_sex, sex chr: $sex_chr\n"
}


close ($output);













