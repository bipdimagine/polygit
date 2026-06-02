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
my $deeptools = "/data-pure/software/SINGULARITY/deeptools.sif ";
my $patient = $project->getPatient($patient_name);
	my $run = $patient->getRun();
my $ref =  $project->genomeFasta();
my $dir = $project->getCoverageDir();
my $outf = $dir."/".$patient->name.".RGPC.bw";
#mosdepth -t 8 -x -b 1000 -Q 20 "${out_path}/${sample_id}" "${bam_path}";
my $align = $patient->getAlignmentFile();

my $vcf = $patient->vcfFileName("deepvariant");
my $cmd = qq{$singularity $deeptools bamCoverage -b $align -o $outf  -p $fork --binSize 1000 --normalizeUsing BPM --extendReads 0 --minMappingQuality 10};
warn $cmd;
system($cmd);

 # --cnv-excluded-regions ${DISTRO_ROOT_DIR}/data/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz \
 # --output-dir HG002_discover_dir
#warn $cmd;
#system($cmd);
exit(0);


