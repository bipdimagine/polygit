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
my $pbsv = "/data-pure/software/SINGULARITY/pbsv.sif ";
foreach my $patient (@{$project->getPatients}){
	my $run = $patient->getRun();
my $ref =  $project->genomeFasta();
my $bam = $patient->getBamFile() ;
my $npz =  $patient->fileWiseCondor();
my $dir_pipeline= $project->getCallingPipelineDir($patient->name.".pbsv");
my $dir_opt = $dir_pipeline."/".$patient->name;
my $fo = $dir_pipeline."/".$patient->name.".svsig.gz";

#mosdepth -t 8 -x -b 1000 -Q 20 "${out_path}/${sample_id}" "${bam_path}";
my $align = $patient->getAlignmentFile();
my $expected = "/data-pure/software/distrib/sawfish/data/expected_cn/expected_cn.hg38.XX.bed";
my $excluded = "/data-pure/software/distrib/sawfish/data/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz";
if ($patient->isMale) {
	 $expected = "/data-pure/software/distrib/sawfish/data/expected_cn/expected_cn.hg38.XY.bed";
}

my $cmd = " $singularity $sawfish pbsv discover $align $fo";



