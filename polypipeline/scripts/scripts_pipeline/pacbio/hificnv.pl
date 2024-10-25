#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
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
 my $fork =20;
GetOptions(
	'project=s'   => \$project_name,
	"proc=s"  => \$fork,
	"patient=s" => \$patient_name,
);



my $buffer = GBuffer->new();

my $project = $buffer->newProject( -name => $project_name );

my $hifi  = "hificnv";# $project->getSoftware('hificnv');

my $bam2cfg  = $project->getSoftware('bam2cfg');

my $tabix  = $project->getSoftware('tabix');
my $ref =  $project->genomeFasta();
my $patient = $project->getPatient($patient_name);
 my $bam = $patient->getBamFile() ;
 my $dirout= $project->getCallingPipelineDir("hificnv");
my $name = $patient->name;

my $cmd =qq{$hifi --ref $ref --bam $bam --output-prefix $dirout/$name --threads $fork};
system($cmd);
my ($vcf1) = `ls $dirout/$name*.vcf.gz`;
chomp($vcf1);
my $prod_file = $project->getVariationsDir("hificnv")."/".$name.".vcf.gz";
warn "mv $vcf1 $prod_file && tabix -f -p  vcf $prod_file";
system ("mv $vcf1 $prod_file && tabix -f -p  vcf $prod_file");
exit(0);

