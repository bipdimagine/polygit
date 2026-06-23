#!/usr/bin/env perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
#use Set::IntSpan;
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use Carp;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use String::ProgressBar;
use List::Util qw(sum);




 my $project_name;
 my $patient_name;
 #my $low_calling;
 my $method;
 

 
my $fork = 5;
GetOptions(
	'project=s'   => \$project_name,
	"patient=s" => \$patient_name,
	"fork=s" => \$fork,
);
die("miss fork") unless $fork;


my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );


my $singularity = "run_singularity.sh";

my @ds;

my $deeptools = "deeptools.sif ";
my $patient = $project->getPatient($patient_name);
	my $run = $patient->getRun();
my $ref =  $project->genomeFasta();
my $dir = $project->getCoverageDir();
my $outf = $dir."/".$patient->name.".bw";
#mosdepth -t 8 -x -b 1000 -Q 20 "${out_path}/${sample_id}" "${bam_path}";
my $align = $patient->getAlignmentFile();
my $cmd = qq{$singularity $deeptools bamCoverage -b $align -o $outf --binSize 50 -p $fork};
#my $cmd = qq{$singularity $deeptools bamCoverage -b $align -o $outf  -p $fork --binSize 50 --normalizeUsing None --extendReads 0 --minMappingQuality 10};
warn $cmd;
system($cmd);

exit(0);


