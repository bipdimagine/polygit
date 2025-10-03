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
use autodie;
use Parallel::ForkManager;

use GBuffer;
my $buffer = new GBuffer;

my @project_names;
my @patient_names;
my $prop = 30/400;
my $fork = 7;
my $no_exec;
GetOptions(
	'projects=s{1,}'			=> \@project_names,	# NGS2023_7287,NGS2023_7288,NGS2023_6681,NGS2023_6732,NGS2020_3222,NGS2025_08859,NGS2025_08858,NGS2024_7753,NGS2024_7752,NGS2020_3005,NGS2022_5805
	'patients=s{1,}'			=> \@patient_names,	# WIHU23-1711,TIQU23-0650,AFO-ILA,VOSA25-0485,HAMO24-0957,KAM-YAN,WABE22-0467
	'probability|proportion=f'	=> \$prop,
	'fork=i'					=> \$fork,
	'no_exec'					=> \$no_exec,
) || die ("Error in command line arguments\n");

#warn $prop;
@patient_names = split(/,/, join(',',@patient_names));
@project_names = split(/,/, join(',',@project_names));
my $java = $buffer->getSoftware('java');
my $picard = $buffer->getSoftware('picard');

my $jobs;
my $pm = new Parallel::ForkManager($fork);

foreach my $project_name (@project_names) {
	my $pid = $pm->start and next;
	warn $project_name;
	my $project = $buffer->newProject(-name=>$project_name) || die ("Can't open project '$project_name': $!");
	my $patients = $project->get_only_list_patients(join(',',@patient_names));
	die("No patient in project '$project_name'") unless ($patients);
	my $tmp = $project->getAlignmentPipelineDir("DownsampleSam");
	
	
	foreach my $pat (@$patients) {
		my $pname = $pat->name;
		warn $pname;
		my $bam_in = $pat->getBamFile;
		my $bam_out = $bam_in =~ s/$pname\.bam$/$pname\_subsampled_30X\.bam/r;
		$bam_out =~ s/$project_name/NGS2025_090XX/r;
		my $cmd = "$java -jar $picard DownsampleSam -P $prop -I $bam_in -O $bam_out --CREATE_INDEX true --TMP_DIR $tmp";
		warn $cmd;
		system($cmd) unless ($no_exec);
#		$jobs .= $cmd."\n";
	}
	print "\n";
	$pm->finish;
}
$pm->wait_all_children();

#warn $jobs;
#system("echo $jobs L run_cluster.pl -cpu=")