#! /usr/bin/env perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
#use Set::IntSpan;
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use Carp;
use Storable qw(store retrieve freeze);
use List::Util qw(sum);
use autodie qw(system);



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
my $force;

GetOptions(
	'project=s'   => \$project_name,
	"patient=s" => \$patient_name,
	"fork=s" => \$fork,
	"force=s"  =>\$force,
);
die("miss fork") unless $fork;


my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );




#my $singularity= "singularity run --bind /data-pure:/data-pure --bind /data-isilon:/data-isilon --bind /data-beegfs:/data-beegfs " ;
my $singularity = "run_singularity.sh";# "/data-bipd/".$buffer->software("singularity-run");
my @ds;
my $sawfish = "sawfish.sif ";
my $patient = $project->getPatient($patient_name);
my $run = $patient->getRun();
my $ref =  $project->genomeFasta();
my $dir_out= $project->getAlignmentPipelineDir($patient->name);
my $bam_out = $dir_out."/".$patient->name.".bam";
die() unless -e $bam_out;


my $dir_pipeline = $dir_out."/sawfiwh1".time;
mkdir $dir_pipeline;
my $dir_opt2 = $dir_out."/sawfiwh2".time;
mkdir $dir_opt2;


my $expected = "/data-bipd/data-pure/software/distrib/sawfish/data/expected_cn/expected_cn.hg38.XX.bed";
my $excluded = "/data-bipd/data-pure/software/distrib/sawfish/data/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz";
if ($patient->isMale) {
	 $expected = "/data-bipd/data-pure/software/distrib/sawfish/data/expected_cn/expected_cn.hg38.XY.bed";
}
my $vcf = $dir_out."/".$patient->name.".vcf.gz" ;
 $vcf = $patient->vcfFileName("deepvariant") unless -e $vcf;
my $prod_file = $dir_out."/".$patient->name.".vcf.gz";
if (-e $prod_file){
	if ($force){
		unlink $prod_file;
		unlink $prod_file.".tbi";
	}
	else {
		warn "already done";
		exit(0);
	}
}
my $cmd = qq{$singularity $sawfish sawfish discover --threads $fork --ref $ref  --bam $bam_out  --expected-cn $expected --cnv-excluded-regions $excluded --output-dir $dir_pipeline --maf-sample-name $vcf };

my $cmd2 = qq{$singularity $sawfish sawfish joint-call --threads $fork --sample $dir_pipeline --output-dir $dir_opt2 };

 # --cnv-excluded-regions ${DISTRO_ROOT_DIR}/data/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz \
 # --output-dir HG002_discover_dir
#warn $cmd;
warn "$cmd && $cmd2";
warn "$dir_opt2/genotyped.sv.vcf.gz";
system("$cmd && $cmd2" ) unless -e "$dir_opt2/genotyped.sv.vcf.gz";
#die();
my $prod_file = $dir_out."/".$patient->name.".vcf.gz";
my $cmd3 = "cp $dir_opt2/genotyped.sv.vcf.gz  $prod_file && tabix -p vcf $prod_file";
system($cmd3);
die() unless -e $prod_file.".tbi";

#print $cmd ."\n";
#print $cmd ." && $cmd2 && $cmd3\n";
exit(0);


