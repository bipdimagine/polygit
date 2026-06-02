#!/usr/bin/perl
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
 my $fork;
 my $callable_intspan_file;
 my $patient_name;
 #my $low_calling;
 my $method;
 
 
my $fork = 5;
GetOptions(
	'project=s'   => \$project_name,
	"patient=s" => \$patient_name,
);


my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatient($patient_name);
my $dir_ubam = $patient->getSequencesDirectory();
my $ubam_file = $dir_ubam."/".$patient->name."/".$patient->name."_merged.bam";

my $bam_prod = $patient->getBamFileName("pbmm2");

my $ref               = $project->genomeFasta();

my $cmd = qq{run_singularity.sh pacbio_analysis.sif pbmm2 align $ref $ubam_file $bam_prod --sort };

system("$cmd");
die() unless -e $bam_prod.".bai";



