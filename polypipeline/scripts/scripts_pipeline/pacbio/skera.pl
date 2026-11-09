#! /usr/bin/env perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use GBuffer; 
use Getopt::Long;
use Carp;
use autodie qw(system);


my $project_name;
my $patient_name;
my $force;
my $fork = 10;
GetOptions(
	'project=s'	=> \$project_name,
	'patient=s'	=> \$patient_name,
	'fork=i'	=> \$fork,
	'force'		=> \$force,
) || confess("Error in command line arguments\n");

die("Error: -project is mandatory") unless ($project_name);
die("Error: -patient is mandatory") unless ($patient_name);

my $buffer			= GBuffer->new();
my $project			= $buffer->newProject( -name => $project_name );
my $patient			= $project->getPatient($patient_name);
my $ubam			= $patient->uBam;

my $sequencing_method = $patient->getRun->sequencing_method;
warn $sequencing_method;
$sequencing_method =~ /^kinnex-(\d+)$/;
my $concat = $1;
confess("Error: concat can be one of (8, 12, 16), found '$concat'.") unless (grep {$concat} (8, 12, 16) );

my $dir_out			= $project->getAlignmentPipelineDir('skera');
my $bam_out			= $dir_out.$patient->name.'.segmented.bam';

if (-e  $bam_out){
	if ($force){
		unlink $bam_out;
		unlink $bam_out.'pbi';
	}
	else {
		warn "Already done: $bam_out";
		exit(0);
	}
}

my $MAS_adapters	= $buffer->config_path('root','public_data')."pacbio/MAS_adapters/mas${concat}_primers.fasta";

my $cmd = qq{cd $dir_out && exec_singularity.sh smrttools.sif skera split $ubam $MAS_adapters $bam_out -j $fork};
print $cmd."\n";
system($cmd);
die("Error running skera: no output '$bam_out'") unless -e $bam_out;

exit(0);




