#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use Data::Dumper;
use Getopt::Long; # librairie pour récuper les options depuis le terminal
use Carp; # librairie pour `confess`
use IO::Prompt;
use Term::ANSIColor;
use Text::CSV_XS qw( csv );
use DBI;

use GBuffer;


my $project_name;
my $patient_name;
my $version;
GetOptions(
	'project=s'	=> \$project_name,
	'patient=s'	=> \$patient_name,
	'version=s'	=> \$version,
) || confess("Error in command line arguments");

confess("ERROR: -project mandatory.\n") unless ($project_name);
confess("ERROR: -patient mandatory.\n") unless ($patient_name);


my $buffer = new GBuffer;
my $project = $buffer->newProject(-name => $project_name, -version =>$version);
my $project_id = $project->id;
my $patient = $project->getPatient($patient_name);

# Get sex from db
my $dbh = $buffer->dbh();
my $pat_id = $patient->id;
my $where_cmd_dbh = "patient_id = '$pat_id' and project_id = '$project_id' and name = '$patient_name' ";
my $cmd_dbh_sel = "SELECT name, patient_id, sex FROM PolyprojectNGS.patient WHERE $where_cmd_dbh";
my $sel = $dbh->selectall_hashref($cmd_dbh_sel, 'name');
confess("ERROR selecting patient from database:\n$cmd_dbh_sel:\n".Dumper $sel) unless (keys %$sel == 1);
my $pat_sex = $sel->{$patient_name}->{sex};
warn $patient_name.' sex from db = '.$pat_sex;


# Estimated sex chr
my $computed_sex;
my $sex_chr;
my $average_cov;
my $cov_SRY;
my $dir_dragen_pipeline = $patient->getDragenDirName("pipeline");
my $dir_stats_dragen = $patient->project->getAlignmentStatsDir("dragen-align");

if ($project->isGenome or $project->isExome) {
	my $ploidy_estimation_file = $dir_dragen_pipeline."/$patient_name.ploidy_estimation_metrics.csv";
	my $ploidy_estimation_file_2 = $ploidy_estimation_file =~ s/^$dir_dragen_pipeline/$dir_stats_dragen/r;
	$ploidy_estimation_file = $ploidy_estimation_file_2 unless (-e $ploidy_estimation_file);
#	warn $ploidy_estimation_file;
	confess("ERROR: No ploidy estimation file found: '$ploidy_estimation_file' or '$ploidy_estimation_file_2'") unless (-e $ploidy_estimation_file);
	my $aoa = csv (in => $ploidy_estimation_file);
	$sex_chr = $aoa->[-1]->[-1] if ($aoa->[-1]->[2] eq 'Ploidy estimation');
	warn 'Dragen ploidy estimation = '.$sex_chr;
	confess("Error parsing '$ploidy_estimation_file':\nno ploidy estimation found for sex chromosomes.") unless ($sex_chr);
	
	if ($sex_chr !~ /X[XY]/) {
		confess ("Anormal ploidy on sex chromosomes for patient '$patient_name': $sex_chr\n");
#		die unless (prompt("Continue anyway ?  (y/n)  ", -yes_no));
	}
	else {
		$computed_sex = 1 if ($sex_chr eq 'XY');
		$computed_sex = 2 if ($sex_chr eq 'XX');
	}
}
elsif ($project->isDiagnostic) {
	# /data-pure/sequencing/ngs/NGS2025_09810/HG19_MT/align/dragen-align/stats/DOL_Ele.target_bed_coverage_metrics.csv
	my $target_cov_metrics_file = $dir_dragen_pipeline."/$patient_name.target_bed_coverage_metrics.csv";
	my $target_cov_metrics_file_2 = $target_cov_metrics_file =~ s/^$dir_dragen_pipeline/$dir_stats_dragen/r;
	$target_cov_metrics_file = $target_cov_metrics_file_2 unless (-e $target_cov_metrics_file);
#	warn $target_cov_metrics_file;
	confess("ERROR: No target bed coverage metrics file found: '$target_cov_metrics_file' or '$target_cov_metrics_file_2'") unless (-e $target_cov_metrics_file);
	my $aoa = csv (in => $target_cov_metrics_file);
	$average_cov = $aoa->[2]->[3] if ($aoa->[2]->[2] eq 'Average alignment coverage over target region');
	confess("Error parsing '$target_cov_metrics_file':\nno average alignment coverage over target region found.") unless ($average_cov);
	$cov_SRY = $patient->coverage_SRY;
	$computed_sex = 2 if ($average_cov < 5);
	$computed_sex = 1 if ($cov_SRY > 30);
	$computed_sex = 2 if ($cov_SRY < 0.05 * $average_cov);
	$computed_sex = 1 if ($cov_SRY > 0.1 * $average_cov);
	warn $cov_SRY.'/'.$average_cov.' = '.($cov_SRY/$average_cov*100).'% => '.$computed_sex;
	confess("ERROR could not compute sex : average cov = $average_cov, SRY cov = $cov_SRY") unless ($computed_sex);
}
else {
	confess("Project is not WGS, WES or capture");
}


# Change sex in db if sex unknown (= 0)
if ($pat_sex == 0) {
	$where_cmd_dbh .= "and sex = 0";
	my $cmd_dbh_update = "UPDATE PolyprojectNGS.patient SET sex = $computed_sex where $where_cmd_dbh";
	my $sth = $dbh->do($cmd_dbh_update) or confess("ERROR changing the sex of patient '$patient_name' to $computed_sex: ".$dbh->errstr);
	warn("'$patient_name' sex changed from $pat_sex to $computed_sex");
	
}
# Not matching
elsif ( not $pat_sex == $computed_sex ) {
	my $msg_err = "Estimated sex and patient sex not matching for patient '$patient_name': sex: $pat_sex, estimated sex: $computed_sex\n";
	$msg_err .= "(cov SRY/avg = $cov_SRY/".sprintf("%d", $average_cov)." = ".sprintf("%.1f%%", ($cov_SRY/$average_cov*100)).")\n" if ($project->isDiagnostic);
	$msg_err .= "(Dragen ploidy estimation = $sex_chr)\n" if ($project->isGenome or $project->isExome);
	confess($msg_err);
}


exit(0);


