#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../GenBo/lib/obj-nodb/";
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
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"patient=s" => \$patient_name,
);



my $buffer = GBuffer->new();

my $project = $buffer->newProject( -name => $project_name );
my $breakdancer  = $project->getSoftware('breakdancer');
my $bam2cfg  = $project->getSoftware('bam2cfg');

my $tabix  = $project->getSoftware('tabix');
my $ref =  $project->genomeFasta();
my $patient = $project->getPatient($patient_name);
 my $bam = $patient->getBamFile() ;
 my $dirout= $project->getCallingPipelineDir("breakdancer");
 my $dd .="$dirout/".$patient->name();
system("mkdir -p $dirout");

unless (-e "$dd" ){
	system("mkdir -p $dd");
}

 my $cfg = $dd."/".$patient->name.".cfg";
my $cmd = qq{$bam2cfg $bam > $cfg}; 
warn $cmd;
system ($cmd);
my $fork= 15;
my $files;
my @files;
 my $pm = new Parallel::ForkManager($fork);
 foreach my $chr (@{$project->getChromosomes}){
 	warn $chr->name;
 	my $output = $dd."/".$patient->name.".".$chr->name.".out";
 	push (@files,$output);
 	my $pid = $pm->start and next;
 		warn "$breakdancer -o ".$chr->fasta_name." $cfg > $output";
		system ("$breakdancer -o ".$chr->fasta_name." $cfg > $output");
		
		$pm->finish(0);
	}
	$pm->wait_all_children();

 my $final_dir =  $project->getVariationsDir("breakdancer") ;
 my $final_file =  $final_dir."/".$patient->name.".breakdancer.txt" ;
 system("cat ".join(" ",@files).">".$final_file." && bgzip $final_file");
 	
	