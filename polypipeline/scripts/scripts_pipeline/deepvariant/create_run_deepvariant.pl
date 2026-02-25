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
 
my $RED    = "\033[31m";
my $GREEN  = "\033[32m";
my $YELLOW = "\033[33m";
my $BLUE   = "\033[34m";
my $BOLD   = "\033[1m";
my $RESET  = "\033[0m";
 
 my $project_name;
my $final_vcf;
my $log_file;
my $fork;
my $force;
my $type;
GetOptions(
	'project=s'   => \$project_name,
	'type=s'  => \$type,
	"log=s" =>\$log_file,
	"vcf=s" => \$final_vcf,
	"fork=s" =>\$fork,
	"force=s" =>\$force,
);
die("please use -type=deepvariant or -type=deepsomatic or -type=deepvariant,deepsomatic") unless $type;
foreach my $t (split(",",$type)){
die("please use -type=deepvariant or -type=deepsomatic or -type=deepvariant,deepsomatic") if $t ne "deepvariant" and $t ne "deepsomatic";
}

my $date = `date`;
chomp($date);
print "\n" . "=" x 60 . "\n";
	print "---- $RED CREATE  RUN  ".uc($type)." $RESET --  \n ";
print "\n" . "=" x 60 . "\n";	
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

foreach my $t (split(",",$type)){
my $ref               = $project->genomeFasta();
my $deepvariant = $buffer->software("deepvariant-sif");
my $dcmd = "run_deepvariant   --model_type=WES";
if ($type =~ /somatic/){
	$deepvariant = $buffer->software("deepsomatic-sif");
	$dcmd = "run_deepsomatic   --model_type=WES_TUMOR_ONLY";
}
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

my $nsh = "";
$nsh = "-num_shards=$fork" if $fork;
if ($project->isGenome) {
	die();
 $cmd = qq{$singularity $deepvariant  --intermediate_results_dir=$dir_gvcf_tmp --ref=$ref --reads=$bam --output_vcf=$vcf_out $nsh};
}
else {
 $cmd = qq{$singularity $deepvariant $dcmd --ref=$ref --reads=$bam --regions=$bed --output_vcf=$vcf_out $nsh > /dev/null 2>&1};
}
my $vcf = $patient->getVariationsFileName($type);

my $cmd2 =qq{$singularity $bcftools bcftools view   -c 1  -e '(QUAL<30 && FORMAT/VAF<0.15 && FORMAT/DP < 10)' $vcf_out -o $vcf -O z > /dev/null 2>&1 && $singularity $bcftools bcftools tabix -f -p vcf $vcf > /dev/null 2>&1};
print CMD "$cmd && $cmd2 \n";
}
}

close (CMD);
print "$GREEN your file is here : $RESET\n".$cmd_file."\n";
