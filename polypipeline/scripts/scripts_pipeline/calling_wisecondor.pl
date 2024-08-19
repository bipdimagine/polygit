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
 
 my $chr_syno ={
		1=> "chr1",
		2=>"chr2",
		3=>"chr3",
		4=>"chr4",
		5=>"chr5",
		6=>"chr6",
		7=>"chr7",
		8=>"chr8",
		9=>"chr9",
		10=>"chr10",
		11=>"chr11",
		12=>"chr12",
		13=>"chr13",
		14=>"chr14",
		15=>"chr15",
		16=>"chr16",
		17=>"chr17",
		18=>"chr18",
		19=>"chr19",
		20=>"chr20",
		21=>"chr21",
		22=>"chr22",
		X=>"chrX",
		Y=>"chrY",
		MT=>"chrM",
};

 
 my $fork = 5;
GetOptions(
	'project=s'   => \$project_name,
	"patient=s" => \$patient_name,
	"fork=s" => \$fork,
);
die("miss fork") unless $fork;


my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $wise  = $project->getSoftware('wisecondor');
my $patient = $project->getPatient($patient_name);
 my $bam = $patient->getBamFile() ;
 my $npz =  $patient->fileWiseCondor();
 
 my $ref = $project->get_wisecondor_reference;
 #'"/data-beegfs/npz/reference/ref5Kb.npz";
 my $out = $project->getCallingPipelineDir("wiseCondor")."/".$patient->name;
 warn $out;
 my $cmd = "$wise predict  $npz $ref $out --beta 1 --bed";
 system("$cmd");
 
 my $bgzip = $buffer->software("bgzip");
 my $tabix = $buffer->software("tabix");
my $outbed1 = $out."_bins.bed";

my $prod_file = $project->getVariationsDir("wisecondor")."/".$patient->name."_bins.bed.gz";

system("$bgzip $outbed1 && mv $outbed1.gz $prod_file ; $tabix -f -p bed -S 1 $prod_file ");

$outbed1 = $out."_aberrations.bed";
$prod_file = $project->getVariationsDir("wisecondor")."/".$patient->name."_aberrations.bed.gz";
system("$bgzip $outbed1 && mv $outbed1.gz $prod_file ; $tabix -f -p bed -S 1 $prod_file ");

 exit(0);
 