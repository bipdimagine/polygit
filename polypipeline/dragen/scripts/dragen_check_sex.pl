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
use Text::CSV_XS qw( csv );
use DBI;

use GBuffer;


my $project_name;
my $patient_name;
my $version;
my $verbose;
GetOptions(
	'project=s'	=> \$project_name,
	'patient=s'	=> \$patient_name,
	'version=s'	=> \$version,
	'verbose'	=> \$verbose,
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
warn $patient_name.' sex from db = '.$pat_sex if ($verbose);


# Estimated sex chr
my $computed_sex;
my $sex_chr;
my $average_cov;
my $cov_SRY;
my $cov_chrX;
my $cov_chrY;
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
	warn 'Dragen ploidy estimation = '.$sex_chr if ($verbose);
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
	my $target_cov_metrics_file_1 = $dir_dragen_pipeline."/$patient_name.target_bed_coverage_metrics.csv";
	my $target_cov_metrics_file_2 = $dir_stats_dragen."/$patient_name.target_bed_coverage_metrics.csv";
	my $target_cov_metrics_file = $target_cov_metrics_file_2 if (-e $target_cov_metrics_file_2);
	$target_cov_metrics_file = $target_cov_metrics_file_1 if (-e $target_cov_metrics_file_1);
#	warn $target_cov_metrics_file;
	confess("ERROR: No target bed coverage metrics file found: '$target_cov_metrics_file' or '$target_cov_metrics_file_2'") unless (-e $target_cov_metrics_file);
	my $aoa = csv (in => $target_cov_metrics_file);
	
	$cov_chrY = $aoa->[27]->[3] if ($aoa->[27]->[2] eq 'Average chr Y coverage over target region');
#	confess("Error parsing '$target_cov_metrics_file':\nno average chr Y coverage over target region: $cov_chrY.") if ($cov_chrY eq 'NaN');
#	confess("Error parsing '$target_cov_metrics_file':\nno average chr Y coverage over target region found.") unless ($cov_chrY);
	$cov_chrX = $aoa->[26]->[3] if ($aoa->[26]->[2] eq 'Average chr X coverage over target region');
#	confess("Error parsing '$target_cov_metrics_file':\nno average chr X coverage over target region: $cov_chrX.") if ($cov_chrX eq 'NaN');
#	confess("Error parsing '$target_cov_metrics_file':\nno average chr X coverage over target region found.") unless ($cov_chrX);
	$computed_sex = 2;
	$computed_sex = 1 if ($cov_chrY/$cov_chrX >= 0.2);
	warn 'cov chrY/chrX = '.sprintf("%.2f",$cov_chrY).'/'.sprintf("%.2f",$cov_chrX).' = '.sprintf("%.2f",$cov_chrY/$cov_chrX).' => '.$computed_sex if ($verbose);
	confess("ERROR could not compute sex : Average chr Y coverage = $cov_chrY, Average chr X coverage = $cov_chrX") unless ($computed_sex);
	
	unless ($cov_chrY >= 0 and $cov_chrX > 0) {
		$average_cov = $aoa->[2]->[3] if ($aoa->[2]->[2] eq 'Average alignment coverage over target region');
		confess("Error parsing '$target_cov_metrics_file':\nno average alignment/chr X/chr Y coverage over target region found.") unless ($average_cov);
		$cov_SRY = $patient->coverage_SRY;
		$computed_sex = 2 if ($average_cov < 5);
		$computed_sex = 1 if ($cov_SRY > 30);
		$computed_sex = 2 if ($cov_SRY/$average_cov < 0.05 );
		$computed_sex = 1 if ($cov_SRY/$average_cov > 0.1 );
		warn 'cov SRY/average = '.$cov_SRY.'/'.$average_cov.' = '.sprintf("%.2f",$cov_SRY/$average_cov*100).'% => '.$computed_sex if ($verbose);
		confess("ERROR could not compute sex : average cov = $average_cov, SRY cov = $cov_SRY") unless ($computed_sex);
	}
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
# Anormal ploidy
unless ($computed_sex) {
	my $msg_err = "Anormal ploidy estimation for patient '$patient_name': $sex_chr (sex in db: $pat_sex)";
	confess($msg_err);
} 
# Not matching
elsif ( not $pat_sex == $computed_sex ) {
	my $msg_err = "Estimated sex and patient sex not matching for patient '$patient_name': sex: $pat_sex, estimated sex: $computed_sex ";
	$msg_err .= "(cov SRY/avg = $cov_SRY/".sprintf("%d", $average_cov)." = ".sprintf("%.1f%%", ($cov_SRY/$average_cov*100)).")" if ($project->isDiagnostic and not ($cov_chrY >= 0 and $cov_chrX > 0));
	$msg_err .= "(cov chrY/chrX = ".sprintf("%.d",$cov_chrY).'/'.sprintf("%.d",$cov_chrX)." = ".sprintf("%.2f%%", ($cov_chrY/$cov_chrX)).")" if ($project->isDiagnostic and ($cov_chrY >= 0 and $cov_chrX > 0));
	$msg_err .= "(Dragen ploidy estimation = $sex_chr)" if ($project->isGenome or $project->isExome);
	confess($msg_err);
}


exit(0);


