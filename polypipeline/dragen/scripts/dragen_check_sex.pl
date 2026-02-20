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
warn $patient_name if ($verbose);

# Get sex from db
my $dbh = $buffer->dbh();
my $pat_id = $patient->id;
my $where_cmd_dbh = "patient_id = '$pat_id' and project_id = '$project_id' and name = '$patient_name' ";
my $cmd_dbh_sel = "SELECT name, patient_id, sex FROM PolyprojectNGS.patient WHERE $where_cmd_dbh";
my $sel = $dbh->selectall_hashref($cmd_dbh_sel, 'name');
confess("ERROR selecting patient from database:\n$cmd_dbh_sel:\n".Dumper $sel) unless (keys %$sel == 1);
my $pat_sex = $sel->{$patient_name}->{sex};
warn 'Sex from database = '.$pat_sex if ($verbose);


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
	my $ploidy_estimation_file_1 = $dir_dragen_pipeline."/$patient_name.ploidy_estimation_metrics.csv";
	my $ploidy_estimation_file_2 = $dir_stats_dragen."/$patient_name.ploidy_estimation_metrics.csv";
	my $ploidy_estimation_file = $ploidy_estimation_file_2 if (-e $ploidy_estimation_file_2);
	$ploidy_estimation_file = $ploidy_estimation_file_1 if (-e $ploidy_estimation_file_1);
#	warn $ploidy_estimation_file;
	confess("ERROR: No ploidy estimation file found: '$ploidy_estimation_file_1' or '$ploidy_estimation_file_2'") unless (-e $ploidy_estimation_file);
	my $aoa = csv (in => $ploidy_estimation_file);
	$sex_chr = $aoa->[-1]->[-1] if ($aoa->[-1]->[2] eq 'Ploidy estimation');
	warn 'Dragen ploidy estimation = '.$sex_chr if ($verbose);
	
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
	my $target_cov_metrics_file_1 = $dir_dragen_pipeline."/$patient_name.target_bed_coverage_metrics.csv";
	my $target_cov_metrics_file_2 = $dir_stats_dragen."/$patient_name.target_bed_coverage_metrics.csv";
	my $target_cov_metrics_file = $target_cov_metrics_file_2 if (-e $target_cov_metrics_file_2);
	$target_cov_metrics_file = $target_cov_metrics_file_1 if (-e $target_cov_metrics_file_1);
	warn $target_cov_metrics_file if ($verbose);
	confess("ERROR: No target bed coverage metrics file found: '$target_cov_metrics_file_1' or '$target_cov_metrics_file_2'") unless (-e $target_cov_metrics_file);
	my $aoa = csv (in => $target_cov_metrics_file);
	
	$cov_chrX = $aoa->[26]->[3] if ($aoa->[26]->[2] =~ /(Average|Median) chr X coverage( \(ignore 0x regions\))? over target region/);
	warn 'cov chrX: '.$cov_chrX if ($verbose);
	confess("Error parsing '$target_cov_metrics_file': no average chr X coverage over target region found.") unless ($cov_chrX);
	$cov_chrY = $aoa->[27]->[3] if ($aoa->[27]->[2] =~ /(Average|Median) chr Y coverage( \(ignore 0x regions\))? over target region/);
	warn 'cov chrY: '.$cov_chrY if ($verbose);
	confess("Error parsing '$target_cov_metrics_file': no average/median chr Y coverage over target region found.") unless ($cov_chrY);
	$average_cov = $aoa->[2]->[3] if ($aoa->[2]->[2] eq 'Average alignment coverage over target region');
	confess("Error parsing '$target_cov_metrics_file': no average alignment coverage over target region.") unless ($average_cov > 0);
	warn 'avg cov: '.$average_cov if ($verbose);
	$cov_SRY = $patient->coverage_SRY;
	warn 'cov SRY: '.$cov_SRY if ($verbose);
	
	my $computed_sex_1;
	if ($cov_SRY >= 0) { # if chrY is in capture
#		$computed_sex_1 = 2
		$computed_sex_1 = 1 if ($cov_SRY > 300 and $cov_SRY/$average_cov > 0.8);
		$computed_sex_1 = 2 if ($cov_SRY/$average_cov < 0.05 );
		$computed_sex_1 = 1 if ($cov_SRY/$average_cov > 0.1 );
		$computed_sex_1 = 2 if ($average_cov < 5);
		$computed_sex_1 = 1 if ($cov_SRY > 30);
		warn 'cov SRY/average = '.$cov_SRY.'/'.$average_cov.' = '.sprintf("%.2f",$cov_SRY/$average_cov).' => '.$computed_sex_1 if ($verbose);
		confess("ERROR $patient_name: could not compute sex: SRY cov = $cov_SRY, average cov = $average_cov") unless ($computed_sex_1);
		
	}
	my $computed_sex_2;
	if ($cov_chrX >= 0 and $cov_SRY < 0) {
		$computed_sex_2 = 1 if ($cov_chrX/$average_cov <= 0.7);
		$computed_sex_2 = 2 if ($cov_chrX/$average_cov > 0.7);
		warn 'cov chrX/average = '.$cov_chrX.'/'.$average_cov.' = '.sprintf("%.2f",$cov_chrX/$average_cov).' => '.$computed_sex_2 if ($verbose);
		confess("ERROR $patient_name: could not compute sex: average chrX cov = $cov_chrX, average cov = $average_cov") unless ($computed_sex_2);
		confess("Possible anormal ploïdy for $patient_name: cov chrX / avg = ".sprintf("%.d",$cov_chrX).' / '.sprintf("%.d",$average_cov).' = '.sprintf("%.2f",$cov_chrX / $average_cov)) if ($cov_chrX / $average_cov > 1.3);
	}
	confess("ERROR $patient_name: could not compute sex: no median SRY or chrX cov: $cov_SRY, $cov_chrX") unless ($cov_SRY >= 0 or $cov_chrX > 0);
	confess("$patient_name computed sexes using avg cov SRY/avg vs chrX/avg cov are differents: $computed_sex_1 vs $computed_sex_2.") if ($computed_sex_1 and $computed_sex_2 and $computed_sex_1 != $computed_sex_2);
	$computed_sex = $computed_sex_1 if ($computed_sex_1);
	$computed_sex = $computed_sex_2 if ($computed_sex_2);
	confess("ERROR $patient_name: could not compute sex") unless ($computed_sex);
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
	$msg_err .= "(cov SRY/avg = $cov_SRY/".sprintf("%d", $average_cov)." = ".sprintf("%.2f", ($cov_SRY/$average_cov)).")" if ($project->isDiagnostic and $cov_SRY >= 0);
	$msg_err .= "(cov X/avg = ".sprintf("%.d",$cov_chrX).'/'.sprintf("%d", $average_cov)." = ".sprintf("%.2f", ($cov_chrX/$average_cov)).")" if ($project->isDiagnostic and $cov_SRY < 0 and $cov_chrX > 0);
	$msg_err .= "(Dragen ploidy estimation = $sex_chr)" if ($project->isGenome or $project->isExome);
	confess($msg_err);
}
# Matching (+ verbose)
elsif ($pat_sex == $computed_sex and $verbose) {
	my $msg = "Estimated sex and patient sex matching for patient '$patient_name': sex: $pat_sex, estimated sex: $computed_sex ";
	$msg .= "(cov SRY/avg = $cov_SRY/".sprintf("%d", $average_cov)." = ".sprintf("%.2f", ($cov_SRY/$average_cov)).")" if ($project->isDiagnostic and $cov_SRY >= 0);
	$msg .= "(cov X/avg = ".sprintf("%.d",$cov_chrX).'/'.sprintf("%d", $average_cov)." = ".sprintf("%.2f", ($cov_chrX/$average_cov)).")" if ($project->isDiagnostic and $cov_SRY < 0 and $cov_chrX > 0);
	$msg .= "(Dragen ploidy estimation = $sex_chr)" if ($project->isGenome or $project->isExome);
	confess($msg);
	warn($msg);
}


exit(0);


