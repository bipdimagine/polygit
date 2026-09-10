#! /usr/bin/env perl
use FindBin qw($Bin);
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages/";
use Getopt::Long;
use GBuffer;
use autodie qw(system);

my $project_name;
my $final_vcf;
my $log_file;
my $patient_name;
my $fork = 64;
my $force;

GetOptions(
	'project=s'   => \$project_name,
	"patient=s"=>\$patient_name,
	"fork=s" =>\$fork,
	"force=s"  =>\$force,
);
my $date = `date`;
chomp($date);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $bcftools = $buffer->software("bcftools");
my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");


my $deepvariant = "deepvariant.sif";#$buffer->software("deeptrio-sif");
my $singularity = "run_singularity.sh";# "/data-bipd/".$buffer->software("singularity-run");
my $ref = $project->genomeFasta();
my $model = "WGS";
my $patient = $project->getPatient($patient_name);
my $dir_out = $project->getAlignmentPipelineDir($patient->name);
my $bam = $dir_out."/".$patient->name.".bam";
$bam = $patient->getAlignmentFile unless -e $bam;
my $dir_gvcf_tmp = $dir_out."tmp.".time;
mkdir $dir_gvcf_tmp;
$fork=64 if $fork>64;

if ($patient->getRun->isPacBio){
	$model = "PACBIO";
}
elsif ($patient->getRun->isNanopore){
	$model = "ONT_R104";
}
my $vcf = $patient->getVariationsFileName("deepvariant");
my $gvcf = $patient->gvcfFileName("deepvariant");
if (-e $vcf){
	if ($force){
		unlink $vcf;
		unlink ($vcf.".tbi");
	}
	else {
		warn "already done: $vcf";
		exit(0);
	}
}



my $cmd = qq{ulimit -n 65535 && $singularity $deepvariant run_deepvariant  --model_type=$model --intermediate_results_dir=$dir_gvcf_tmp --ref=$ref --reads=$bam --output_vcf=$vcf --output_gvcf=$gvcf --num_shards=$fork};

print $cmd."\n";
system($cmd);

exit(0);
#}

