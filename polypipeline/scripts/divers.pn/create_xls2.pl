#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "/software/polyweb/poly-disk/poly-src/GenBo/lib/obj-nodb/";
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use File::Basename;
use  File::Temp;



my $fork = 1;
my $force;
my ($project_name, $chr_name, $annot_version);
my $patient_name;
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'patient=s'        => \$patient_name,
	'force=s'  => \$force,
);

$patient_name = 'all' unless $patient_name;

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );

my $patients = $project->getPatients();
my $dir = $Bin."/xls.final"."/";
mkdir $dir unless -e $dir;
foreach my $patient (@{$patients}){
	my $fname = $project->name."-".$patient->name;
	next if ($patient_name ne 'all' && $patient->name ne $patient_name);
	warn $fname;
	my $cmd = "/data-xfs/dev/pnitschk/svn-genbo/polymorphism-cgi/polydiag_old/patient_report.pl project=".$project->name
	.qq{ panel= edit_mode=1 never=1 this=6 impact=3 frequence=2 allele_quality=5 report_mode=1 project_summary=1 graphic_coverage=1 sanger_variations=1 validated_variations=1 todo_variations=1 xls=1 project_summary=1 all_coverage=1 table_variations=1 sanger_variations=1 validated_variations=1 all_variations=1 denovo=1 recessive=1 xor=1 both=1 span=20 limit=30 transcripts=all user_name=pnitschk patients=}.
	$patient->name."|  tail -n +4  >$dir".$fname.".xls";
	system($cmd);
}