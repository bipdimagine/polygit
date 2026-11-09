#! /usr/bin/env perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use GBuffer; 
use Getopt::Long;
use Carp;
use autodie qw(system open);


my $project_name;
my $patient_name;
my $force;
my $fork = 0;
GetOptions(
	'project=s'	=> \$project_name,
	'patient=s'	=> \$patient_name,
	'fork=i'	=> \$fork,
	'force'		=> \$force,
) || confess("Error in command line arguments\n");

die("Error: -project is mandatory") unless ($project_name);
die("Error: -patient is mandatory") unless ($patient_name);

my $buffer				= GBuffer->new();
my $project				= $buffer->newProject( -name => $project_name );
my $patient				= $project->getPatient($patient_name);
my $sequencing_method	= $patient->getRun->sequencing_method;
warn $sequencing_method;
my $bam_in;
if ($sequencing_method =~ /^kinnex-(\d+)$/) {
	$bam_in = $project->getAlignmentPipelineDir('skera') . $patient->name.'.segmented.bam';
}
else {
	$bam_in = $patient->uBam;
}
warn $bam_in;
confess ("Error: Input bam '$bam_in' not found") unless (-f $bam_in);

my $dir_out			= $project->getAlignmentPipelineDir('lima_isoseq');
my $bam_out_prefix	= $dir_out . $patient->name.'.fl.bam';
my $bam_out			= $dir_out . $patient->name.'.fl.5p--3p.bam';
my $summary_file 	= $dir_out . $patient->name.'.fl.lima.summary';

if (-e  $bam_out){
	if ($force){
		unlink $bam_out;
		unlink $bam_out.'pbi';
	}
	else {
		warn "Already done: $bam_out";
		check_summary($summary_file);
		exit(0);
	}
}

my $kit = $patient->getChemistry;
my $profile 		= $patient->getSampleProfile;
my $capture_name	= $patient->getCapture->name;
my $primers;
if ($profile =~ /^singlecell/ or $capture_name eq 'transcriptome_10X') {
	if ($kit =~ /^SC3P/) {
		$primers = $buffer->config_path('root','public_data')."pacbio/10x_3p_primers.fa";
	}
	elsif ($kit =~ /^SC5P/) {
		$primers = $buffer->config_path('root','public_data')."pacbio/10x_5p_primers.fa";
	}
	else {
		confess("Error determinig the primers fasta for chemistry/kit '$kit'");
	}
}
else {
	$primers = $buffer->config_path('root','public_data')."pacbio/isoseq_primers.fa";
}
confess("Error: '$primers' doesn't exists.") unless (-f $primers);

my $cmd = qq{cd $dir_out && exec_singularity.sh lima-isoseq.sif lima --isoseq --split-named $bam_in $primers $bam_out_prefix -j $fork};

print $cmd."\n";
system($cmd);
die("Error running lima: no output '$bam_out'") unless -e $bam_out;

check_summary($summary_file);

exit(0);



sub check_summary {
	my ($summary_file, $min_passing_percent) = @_;
	$min_passing_percent //= 95; # Seuil par défaut
	confess("Error: no output summary file '$summary_file'") unless (-e $summary_file);

	open(my $fh, '<', $summary_file) or die "Can't open '$summary_file' : $!";
	my ($total_reads, $passing_reads);
	while (my $line = <$fh>) {
		if ($line =~ /^Reads input\s*:\s*(\d+)/) {
			$total_reads = $1;
		}
		elsif ($line =~ /^Reads above all thresholds\s*\(A\)\s*:\s*(\d+)/) {
			$passing_reads = $1;
	}
	}
	close($fh);
	confess ("Error: Can't exctract 'Reads input' or 'Reads above all thresholds' in '$summary_file'") if (!defined $total_reads || !defined $passing_reads);
	confess ("Error: No reads input in '$summary_file'") if ($total_reads == 0);

	my $pass_rate = ($passing_reads / $total_reads) * 100;
	warn(sprintf("Lima summary: %.2f%% passing reads (%d / %d)", $pass_rate, $passing_reads, $total_reads));
	confess (sprintf("Error control: Reads passing (%.2f%%) is below threshold (%d%%). Check primers sequence: '$primers'.", $pass_rate, $min_passing_percent)) if ($pass_rate < $min_passing_percent);
}



