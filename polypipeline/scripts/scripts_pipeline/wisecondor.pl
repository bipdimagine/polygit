#!/usr/bin/env perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
#use Set::IntSpan;
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use autodie qw(system);


 my $project_name;
 my $fork;
 my $callable_intspan_file;
 my $patient_name;
 #my $low_calling;
 
my $force;
 
 my $fork = 5;
GetOptions(
	'project=s'   => \$project_name,
	"patient=s" => \$patient_name,
	"fork=s" => \$fork,
	"force=s" =>\$force,
);
die("miss fork") unless $fork;

warn "1";
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $wise_sif  = "wisecondor.sif";
my $singularity = "run_singularity.sh";
my $patient = $project->getPatient($patient_name);
 my $bam = $patient->getBamFile(undef,1) ;
 unless (-e $bam){
 	warn "no bam ".$project->name;
 	die();
 }
 my $ref;
 if ($bam =~ /\.cram/ ) {
 	$ref = "--reference ".$project->genomeFasta();
 }
 
 my $chr = $project->getChromosome("1");
 my $fileout = $patient->fileWiseCondor();
 my $fileout_wise = $fileout;
# $fileout_wise = "/data-pure/workspace/tmp/wisecondor/npz/novaseq/".$patient->name;
 $fileout_wise =~ s/\.npz//;
 if (-e $fileout){
 	if ($force){
 		unlink $fileout;
 	}
 	else {
 		exit(0);
 	}
 }
 my $cmd = "exec_singularity.sh wisecondor.sif  wisecondorx convert $bam $fileout_wise $ref";
 system("$cmd");
 
 die() unless -e $fileout;

 