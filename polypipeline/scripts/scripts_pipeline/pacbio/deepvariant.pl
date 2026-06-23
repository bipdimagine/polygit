#! /usr/bin/env perl
use FindBin qw($Bin);
use strict;
use FindBin qw($Bin);
#use lib "/software/polyweb/poly-disk/poly-src/polygit/GenBo/lib/obj-nodb/";
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages/";
use Getopt::Long;
use GBuffer;
use autodie qw(system);
 
 my $project_name;
my $final_vcf;
my $log_file;
my $patient_name;
my $fork;

GetOptions(
	'project=s'   => \$project_name,
	"patient=s"=>\$patient_name,
	"fork=s" =>\$fork,
);
my $date = `date`;
chomp($date);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $bcftools = $buffer->software("bcftools");
my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");




my $cmd;
my $deepvariant = "deepvariant.sif";#$buffer->software("deeptrio-sif");
my $singularity = "run_singularity.sh";# "/data-bipd/".$buffer->software("singularity-run");
my $ref               = $project->genomeFasta();
my $model = "WGS";
#foreach my $patient (@{$project->getPatients}){
	my $patient = $project->getPatient($patient_name);
	my $dir_out= $project->getAlignmentPipelineDir($patient->name);
	my $bam_out = $dir_out."/".$patient->name.".bam";
	
my $vcf_out = $dir_out."/".$patient->name.".deep.vcf.gz";
my $gvcf_out == $dir_out."/".$patient->name.".deep.gvcf.gz";
exit(0) if -e $vcf_out;

my $gvcf_out =  $patient->getVariationsFileName("deepvariant");#$dir_out."/".$patient->name.".deep.gvcf.gz";
my $dir_gvcf_tmp = $dir_out."tmp.".time;
mkdir $dir_gvcf_tmp;
$fork=64 if $fork>64;

if  ($patient->getRun->isPacBio){
	$model = "PACBIO";
}
elsif  ($patient->getRun->isNanopore){
	$model = "ONT_R104";
}
my $vcf = $patient->getVariationsFileName("deepvariant");
my $gvcf = $patient->gvcfFileName("deepvariant");
 $cmd = qq{ulimit -n 65535 && $singularity $deepvariant run_deepvariant  --model_type=$model --intermediate_results_dir=$dir_gvcf_tmp --ref=$ref --reads=$bam_out --output_vcf=$vcf --output_gvcf=$gvcf --num_shards=$fork};

print $cmd."\n";
system($cmd);

exit(0);
#}

