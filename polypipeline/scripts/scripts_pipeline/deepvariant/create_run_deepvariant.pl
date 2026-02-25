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
my $fork;
my $force;

GetOptions(
	'project=s'   => \$project_name,
	"log=s" =>\$log_file,
	"vcf=s" => \$final_vcf,
	"fork=s" =>\$fork,
	"force=s" =>\$force,
);
my $date = `date`;
chomp($date);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );

my $dir_pipeline = $project->getCallingPipelineDir($project->name.time.".deepvariant");


my $bed = $dir_pipeline."/".$project->name.".bed";
my $cmd_file = $dir_pipeline."/".$project->name.".cmd";
open (CMD , ">".$cmd_file);
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
my $deepvariant = $buffer->software("deepvariant-sif");
my $singularity = $buffer->software("singularity-run");
my $bcftools =  $buffer->software("bcftools-sif");
foreach my $patient (@{$project->getPatients}){
my $bam                 = $patient->getBamFile();

my $dir_gvcf_out= $dir_pipeline."/".$patient->name.time."/";
system("mkdir $dir_gvcf_out && chmod a+rwx $dir_gvcf_out") unless -e $dir_gvcf_out;
my $vcf_out = $dir_gvcf_out."/".$patient->name.".deep.vcf.gz";
my $dir_gvcf_tmp = "/tmp/";
my $cmd;
$fork =20 if $fork >20;


if ($project->isGenome) {
 $cmd = qq{$singularity $deepvariant run_deepvariant  --model_type=WGS --intermediate_results_dir=$dir_gvcf_tmp --ref=$ref --reads=$bam --output_vcf=$vcf_out -num_shards=$fork};
}
else {
 $cmd = qq{$singularity $deepvariant run_deepvariant   --model_type=WES --ref=$ref --reads=$bam --regions=$bed --output_vcf=$vcf_out -num_shards=$fork};
}
my $vcf = $patient->getVariationsFileName("deepvariant");

my $cmd2 =qq{$singularity $bcftools bcftools view   -c 1  -e '(QUAL<30 && FORMAT/VAF<0.15 && FORMAT/DP < 10)' $vcf_out -o $vcf -O z && $singularity $bcftools bcftools tabix -f -p vcf $vcf};
print CMD "$cmd && $cmd2 \n";
}
close (CMD);
print "your file is here : \n".$cmd_file."\n";