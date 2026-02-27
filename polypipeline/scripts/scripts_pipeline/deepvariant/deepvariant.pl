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
my $bed = $dir_gvcf_out."/".$patient->name.".bed";
unless ($project->isGenome){
my @zbed;
foreach my $chr (@{$project->getChromosomes}) {
	my $intspan = $chr->getIntSpanCaptureForCalling(150);
	push (@zbed,$buffer->intspanToBed($chr,$intspan));
} 
open (BED , ">".$bed);
	print BED join("\n",@zbed);
close(BED);
}
my $ref               = $project->genomeFasta();
my $bam                 = $patient->getBamFile();
my $vcf_out = $dir_gvcf_out."/".$patient->name.".deep.vcf.gz";
my $dir_gvcf_tmp = $dir_gvcf_out."/tmp.".time;
mkdir $dir_gvcf_tmp;
my $cmd;
$fork =20 if $fork >20;
my $deepvariant = $buffer->software("deepvariant-sif");
my $singularity = $buffer->software("singularity-run");
if ($project->isGenome) {
 $cmd = qq{$singularity $deepvariant run_deepvariant  --model_type=WGS --intermediate_results_dir=$dir_gvcf_tmp --ref=$ref --reads=$bam --output_vcf=$vcf_out -num_shards=$fork};
}
else {
 $cmd = qq{$singularity $deepvariant run_deepvariant   --model_type=WES --ref=$ref --reads=$bam --regions=$bed --output_vcf=$vcf_out -num_shards=$fork};
}
system($cmd);
warn $cmd;
warn $vcf_out;
die() unless -e $vcf_out;

my $vcf = $patient->getVariationsFileName("deepvariant");
#if ($project->isExome or $patient->isGenome){
my $cmd2 =qq{bcftools view   -c 1  -e '(QUAL<30 && FORMAT/VAF<0.15 && FORMAT/DP < 10)' $vcf_out -o $vcf -O z && tabix -f -p vcf $vcf};
system($cmd2);
#}
#else {
#	warn "cp $vcf_out $vcf && tabix -f -p vcf $vcf ";
#	system("cp $vcf_out $vcf && tabix -f -p vcf $vcf ");
#}

die() unless -e $vcf.".tbi";
