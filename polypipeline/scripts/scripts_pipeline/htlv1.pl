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
#my $fileout = $dirout."/".$patient_name;

#my $bam = $patient->getBamFile();
my $fileout =  $project->getVariationsDir("htlv1_calling")."/".$patient_name."-clonalityResults.txt";
my $fileout2 =  $project->getVariationsDir("htlv1_calling")."/".$patient_name."-mergedIS.xls";
my $fileout4 =  $project->getVariationsDir("htlv1_calling")."/".$patient_name."-mergedIS.txt";
my $fileout3 =  $project->getVariationsDir("htlv1_calling")."/".$patient_name."-SIMPLIFIED_mergedIS.txt";


 #exit(0) if -e $fileout;
 my $dir_working = "$dirout/$patient_name"  ;
  if (-e $dir_working){
 	system("rm $dirout/*.log*");
 }
 else {
 system ("mkdir $dir_working");
 }
my $singularity = $buffer->software("singularity");
my $image = "/software/distrib/Clonality/SINGULARITY_R4/clonality_R4.sif";
my $script = "/software/distrib/Clonality/fred_app.sh";

my $cmd = qq{$singularity run --pwd /DATA/htlv1/ -B /software:/software -B /data-isilon/:/data-isilon -B $dir_working:/TMP -B $dirout:/OUT  $image  $script $f1 $f2 $patient_name $ls $ln $fork /TMP --working_directory  /TMP/ };

warn $cmd;

system($cmd." && touch $dirout/$patient_name.ok");



#system("cp $dirout/$patient_name.log $fileout");



my $file_result = $dir_working."/".$patient_name.$ln."-clonalityResults.txt";
my $file_merged_xls = $dir_working."/".$patient_name.$ln."-mergedIS.xls";
my $file_merged_txt = $dir_working."/".$patient_name.$ln."-mergedIS.txt";
my $file_simpl = $dir_working."/".$patient_name.$ln."-SIMPLIFIED_mergedIS.txt";
system("cp $file_result $fileout") if -e $file_result;
system("cp $file_merged_xls $fileout2") if -e $file_merged_xls;
system("cp $file_merged_txt $fileout4") if -e $file_merged_txt;
system("cp $file_simpl $fileout3") if -e $file_simpl;
#system("touch $fileout3") unless -e "$dirout/$patient_name.ok";

#die() unless -e "$dirout/$patient_name.ok";
