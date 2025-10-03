#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";

use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use colored;

use Storable qw(store retrieve freeze);
use Term::ANSIColor;

my $filein;
my $dir;
my $file_bed;
 
my $project_name;
my $patient_name;
my $bam_file;
my $log_file;
my $vcf_final;
my $type;
my $fork;


GetOptions(
	'project=s'   => \$project_name,
	"log=s" =>\$log_file,
	"patient=s" => \$patient_name,
	"fork=s" => \$fork
);
$fork =1 unless $fork;

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatient($patient_name);
my $bam_prod = $patient->getAlignFileName();

my $stat_dir = $project->getMetricsDir();
mkdir $stat_dir unless -e $stat_dir;
my $fileout = $patient->getMetricsFile();
my $metricsHeader = $project->getMetricsHeader();

#my $metricsFile = $patient->getCaptureMetricsFile($metricsHeader);
my $metricsFile = doCaptureMetricsFile($patient, $metricsHeader);

my $bed = $patient->getCaptureFile();
my $picard =  $project->getSoftware('java')." -jar ".$project->getSoftware('picard_path');
my $version = $project->getVersion();
my $public_data = $project->dirGenome();

my $reference_bwa = $project->genomeFasta;#$public_data."/fasta/all.fa";

warn "\n\n";
warn "### PICARD STATS ###\n";
warn "# PROJECT: ".$project->name()."\n";
warn "# PATIENT: ".$patient->name()."\n";
warn "  -> bam_prod: $bam_prod\n";
warn "  -> stat_dir: $stat_dir\n";
warn "  -> fileout: $fileout\n";
warn "  -> metricsHeader: $metricsHeader\n";
warn "  -> bed: $bed\n";
warn "  -> picard: $picard\n";
warn "  -> version: $version\n";
warn "  -> public_data: $public_data\n";
warn "  -> reference_bwa: $reference_bwa\n";
warn "\n\n";

if (-e $bam_prod ){
	my $cmd = qq{$picard CollectHsMetrics BAIT_INTERVALS=$metricsFile VALIDATION_STRINGENCY=LENIENT TARGET_INTERVALS=$metricsFile INPUT=$bam_prod OUTPUT=$fileout REFERENCE_SEQUENCE=$reference_bwa };
	warn 'CMD: '.$cmd;
	system($cmd);
}
else {die("can't find bam file ") ;}

colored::stabilo('green'," ---- END STATS PICARD : ".$patient->name."----");
colored::stabilo('green',"  ----------------------------------------------------");
exit(0);


sub doCaptureMetricsFile {
	my ( $patient, $headerHSmetrics ) = @_;
	my $dir_out = $patient->project->getCallingPipelineDir("picard_stats");
	my @lCaptBedGz = @{$patient->captureFiles()};
	my $gz = $lCaptBedGz[0];
	my $cmd_cp = "cp $gz $dir_out/.";
	warn $cmd_cp;
	`$cmd_cp`;
	my @lTmp_gz = split('/', $gz);
	$gz = $dir_out.'/'.$lTmp_gz[-1];
	my $cmd_zip = "gunzip -f $gz";
	warn $cmd_zip;
	`$cmd_zip`;
	my $bedFileName = $gz;
	$bedFileName =~ s/\.gz//;
	my $metricsFileName = $bedFileName;
	$metricsFileName =~ s/\.bed/\.metrics/;
	my @lTmp = split('/', $metricsFileName);
	my $intervalFileName = $dir_out.'/'.$lTmp[-1] . ".intervals";
	my $cmd = "cat  $bedFileName |awk  -F \"\\t\" 'BEGIN {OFS=\"\\t\"} {print \$1,\$2,\$3,\"+\",(\$3-\$2)+1}' > $metricsFileName";
	`$cmd`;
	my $cmd2 = "cat $headerHSmetrics $metricsFileName > $intervalFileName";
	`$cmd2`;
	return $intervalFileName;
}
	

