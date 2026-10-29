#! /usr/bin/env perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use GBuffer; 
use Getopt::Long;
use List::Util qw(sum);
use autodie qw(system);


my $project_name;
my $patient_name;
my $method;
my $force;
my $fork = 5;
GetOptions(
	'project=s'   => \$project_name,
	"patient=s" => \$patient_name,
	"fork=s" =>\$fork,
	"force=s" =>\$force,
);


my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatient($patient_name);
my $dir_ubam = $patient->getSequencesDirectory();

my $ubam_file = $patient->uBam;
my $dir_out= $project->getAlignmentPipelineDir($patient->name);
my $bam_out = $dir_out."/".$patient->name.".bam";
#if (-e  $bam_out){
#	if ($force){
#		unlink $bam_out;
#		unlink $bam_out."bai";
#	}
#	else {
#		warn "already done: $bam_out";
#		exit(0);
#	}
#}
my $cram_prod = $patient->getCramFileName("minimap2");

my $ref               = $project->genomeFasta();
my $ref_index = $project->getGenomeIndex("minimap2").'/all.mmi';
my $ref_splice = $project->getGenomeIndex("minimap2").'/splice.bed.gz';
my $cmd = qq{samtools fastq -@ 16 ${ubam_file} | run_singularity.sh minimap2.sif minimap2  -ax lr:hq    -t $fork   $ref  - | samtools sort -@ 16 -O BAM -o $cram_prod && samtools index ${cram_prod} -@ $fork };
	
#-ax splice -t 40 -ub --MD --junc-bed /data-pure/public-data/repository/HG38/annotations/gencode.v49/bed/gencode.v49.annotation.bed /data-isilon/public-data/genome/HG38_CNG/fasta/all.fa /data-pure/workspace/mhamici/Alice_Luis_transcriptome/9098_Ags3762.fastq
if ($patient->isRna){
	 $cmd = qq{samtools fastq -@ 16 ${ubam_file} | run_singularity.sh minimap2.sif minimap2  -ax  splice --secondary=no --MD  -t $fork --junc-bed $ref_splice  $ref  - | samtools sort -@ 16 -O BAM -o $cram_prod && samtools index ${cram_prod} -@ $fork };
}
	
print $cmd."\n";
system("$cmd") ;
die() unless -e $cram_prod.".crai";

exit(0);




