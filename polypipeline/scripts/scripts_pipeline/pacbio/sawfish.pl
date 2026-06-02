#!/usr/bin/perl
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




my $singularity= "singularity run --bind /data-pure:/data-pure --bind /data-isilon:/data-isilon --bind /data-beegfs:/data-beegfs " ;
my @ds;
my $sawfish = "/data-pure/software/SINGULARITY/sawfish.sif ";
foreach my $patient (@{$project->getPatients}){
next if ($patient_name && $patient->name ne $patient_name);
my $run = $patient->getRun();
my $ref =  $project->genomeFasta();
my $bam = $patient->getAlignFileName("pbmm2",1) ;
next unless -e $bam;
my $npz =  $patient->fileWiseCondor();
my $dir_pipeline= $project->getCallingPipelineDir($patient->name.".sawfish");
push(@ds,"--sample ".$dir_pipeline);
my $dir_opt = $dir_pipeline."/".$patient->name;
my $dir_opt2 = $dir_pipeline."/".$patient->name.".gentotype";

#mosdepth -t 8 -x -b 1000 -Q 20 "${out_path}/${sample_id}" "${bam_path}";
my $align = $patient->getAlignmentFile();
my $expected = "/data-pure/software/distrib/sawfish/data/expected_cn/expected_cn.hg38.XX.bed";
my $excluded = "/data-pure/software/distrib/sawfish/data/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz";
if ($patient->isMale) {
	 $expected = "/data-pure/software/distrib/sawfish/data/expected_cn/expected_cn.hg38.XY.bed";
}
my $vcf = $patient->vcfFileName("deepvariant");
my $cmd = qq{$singularity $sawfish sawfish discover --threads $fork --ref $ref  --bam $align  --expected-cn $expected --cnv-excluded-regions $excluded --output-dir $dir_pipeline --maf-sample-name $vcf };

my $cmd2 = qq{$singularity $sawfish sawfish joint-call --threads $fork --sample $dir_pipeline --output-dir $dir_opt2 };

 # --cnv-excluded-regions ${DISTRO_ROOT_DIR}/data/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz \
 # --output-dir HG002_discover_dir
#warn $cmd;
system("$cmd && $cmd2" ) unless -e "$dir_opt2/genotyped.sv.vcf.gz";
#die();
my $prod_file = $project->getVariationsDir("sawfish")."/".$patient->name.".vcf.gz";
my $cmd3 = "cp $dir_opt2/genotyped.sv.vcf.gz  $prod_file && tabix -p vcf $prod_file";
system($cmd3);
die() unless -e $prod_file.".tbi";
#print $cmd ."\n";
#print $cmd ." && $cmd2 && $cmd3\n";
}
exit(0);
my $dir_pipeline2= $project->getCallingPipelineDir($project->name.".sawfish.calling");
my $string =  join(" ",@ds);
my $cmd = qq{$singularity $sawfish sawfish joint-call --threads 16 $string  --output-dir $dir_pipeline2 };
warn $cmd;
system($cmd);
warn $cmd;

