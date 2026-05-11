#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use FindBin qw($Bin);
#use lib "/software/polyweb/poly-disk/poly-src/polygit/GenBo/lib/obj-nodb/";
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages/";
use  File::Temp;
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use colored;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
 use IPC::Open2;
 use List::MoreUtils qw(natatime);
 
 my $project_name;
my $final_vcf;
my $log_file;
my $patient_name;
my $fork;

GetOptions(
	'project=s'   => \$project_name,
	"log=s" =>\$log_file,
	"vcf=s" => \$final_vcf,
	"patient=s"=>\$patient_name,
	"fork=s" =>\$fork,
	"print=s"  =>\$fork,
);
my $date = `date`;
chomp($date);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatientOrControl($patient_name);
my $bcftools = $buffer->software("bcftools");
my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");
my $dir_gvcf_out= $project->getCallingPipelineDir($patient->name.time.".deepvariant");
my $gatk = $buffer->software("gatk");
my $bed = $dir_gvcf_out."/".$patient->name.".bed";
my $vcf = $patient->getVariationsFileName("mutect-mt");
my $ref               = $project->genomeFasta();
my $bam                 = $patient->getAlignmentFile();

my $vcf_out = $dir_gvcf_out."/".$patient->name.".mt.vcf.gz";

my $dir_gvcf_tmp = $dir_gvcf_out."/tmp.".time;
my $vcf = $patient->getVariationsFileName("mutect-mt");
my $gatksif = "/data-pure/software/SINGULARITY/gatk.sif";

my $singularity = $buffer->software("singularity-run");


my $cmd;

my $chr = $project->getChromosome("MT");
my $cn = $chr->fasta_name;   

$cmd = qq{$singularity $gatksif gatk Mutect2 -R $ref -I $bam -L $cn  --mitochondria-mode -O $vcf_out };
system($cmd);
$cmd = qq{$singularity $gatksif gatk FilterMutectCalls -R $ref -V  $vcf_out -O  $vcf};
system($cmd);
my $vcf1 = $vcf.".old.gz";
warn "mv $vcf $vcf1; mv $vcf.tbi $vcf1.tbi;$bcftools  view -f PASS $vcf1 -Oz -o $vcf;tabix -p -f $vcf";
system ("mv $vcf $vcf1; mv $vcf.tbi $vcf1.tbi;$bcftools  view -f PASS $vcf1 -Oz -o $vcf;tabix -p -f $vcf");

die() unless -e $vcf_out;
exit(0);
