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
 my $version;
 my $method;
 my $f1;
 my $f2;
 my $ln;
 my $ls;
 
GetOptions(
	'project=s'   => \$project_name,
	'f1=s' => \$f1,
	'f2=s' => \$f2,
	'ln=s' => \$ln,
	'ls=s' => \$ls,
	"fork=s"  => \$fork,
	"patient=s" => \$patient_name,
	"version=s" => \$version,
);


my $buffer = GBuffer->new();
warn $patient_name;
my $project = $buffer->newProject( -name => $project_name,-version=>$version );
my $patient = $project->getPatient($patient_name);
my $dirout= $project->getCallingPipelineDir("htlv1");
print $patient->name().":\n";
my $fileout = $dirout."/".$patient_name;

#my $bam = $patient->getBamFile();
my $fileout2 =  $project->getVariationsDir("htlv1")."/".$patient_name.".log";
#my $fileout2 =  $project->getVariationsDir("muc1")."/".$patient_name.".aln";
#my $fileout4 =  $project->getVariationsDir("muc1")."/".$patient_name.".log";
#my $fileout3 =  $project->getVariationsDir("muc1")."/".$patient_name.".error";
#my $bam_dir = $patient->getProject->getAlignmentDir("bwa");


 #exit(0) if -e $fileout;
 my $dir_working = "$dirout/$patient_name"  ;
  if (-e $dir_working){
 	system("rm $dirout/*.log*");
 }
 else {
 system ("mkdir $dir_working");
 }
my $singularity = $buffer->software("singularity");
my $image = "/software/distrib/Clonality/SINGULARITY/step4_clonality.sif";
my $script = "/software/distrib/Clonality/cecile_app.sh";
my $cmd = qq{$singularity run --pwd /DATA/htlv1/ -B $dir_working:/TMP -B $dirout:/OUT  $image  $script $f1 $f2 $patient_name $ls $ln  $fork --working_directory  /TMP/};

warn $cmd;

#system($cmd." && touch $dirout/$patient_name.ok");



#system("cp $dirout/$patient_name.log $fileout");

#my $file_aln = $dir_working."/log_".$patient_name.".bam.log.aln";
#my $file_log = $dir_working."/log_".$patient_name.".bam.log";
#system("cp $file_aln $fileout2") if -e $file_aln;
#system("cp $file_log $fileout4") if -e $file_log;
#system("touch $fileout3") unless -e "$dirout/$patient_name.ok";

#die() unless -e "$dirout/$patient_name.ok";
