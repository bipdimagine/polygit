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
my $manta  = $project->getSoftware('manta');
my $tabix  = $project->getSoftware('tabix');
my $bgzip  = $project->getSoftware('bgzip');
my $sambamba  = $project->getSoftware('sambamba');
my $samtools  = $project->getSoftware('samtools');
my $lumpy  = $project->getSoftware('lumpy');
my $extractSplitReads_BwaMem = $project->getSoftware("extractSplitReads_BwaMem");

my $ref =  $project->genomeFasta();
my $patient = $project->getPatient($patient_name);
 my $bam = $patient->getBamFile() ;
 my $dirout= $project->getCallingPipelineDir("lumpy");
 my $dd  ="$dirout/".$patient->name();
 if (-e"$dd" ){
	system("rm  $dd/*");
}else {
system("mkdir -p $dd");

}
warn $dd ;
my $discordants = $dd."/".$patient->name().".discordants.bam";
my $f1 = "$dd/sample.discordants.unsorted.bam";
 my $cmd1 = qq{$samtools view -b -F 1294 $bam  > $f1 &&  $sambamba sort -t $fork -o $discordants $f1 ; rm $f1};
 warn $cmd1;
 system($cmd1);
 die("$cmd1") unless -e $discordants;
 my $f2 = "$dd/sample.splitters.unsorted.bam";
 my $splitters = $dd."/".$patient->name().".splitters.bam";
 my $cmd2 =qq{module load python;$samtools view -h $bam | $extractSplitReads_BwaMem -i stdin | samtools view -Sb - > $f2 && $sambamba sort -t $fork -o $splitters $f2;rm $f2;};
 
 warn $cmd2;
 system($cmd2);
 die($cmd2) unless -e $splitters;
 
 my $final_dir =  $project->getVariationsDir("lumpy") ;
 my $out_vcf = $final_dir."/".$patient->name.".vcf";
 my $cmd3 = qq{module load python;$lumpy  -B $bam -S $splitters -D $discordants -o $out_vcf && $bgzip $out_vcf && $tabix -p vcf $out_vcf.gz};   
 warn $cmd3;
 system($cmd3);
 warn $out_vcf;
 
 
 