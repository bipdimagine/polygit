#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use Parallel::ForkManager;
use Storable qw(store retrieve freeze thaw);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Cwd 'abs_path';
use Digest::MD5 qw(md5_hex);
use Statistics::Descriptive::Smoother;
use File::Temp;

use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
#use lib "/software/polyweb/poly-disk/poly-src/GenBo/lib/obj-nodb/packages";
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/packages";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use Set::IntSpan::Fast::XS;
use image_coverage;
 use List::Util qw( min sum max);
 # import specific functions
 
use List::MoreUtils qw(firstidx uniq);
#use Cache_Commons;

my $fork = 1;
my $force;
my ($project_name, $patient_name);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'patient=s'    => \$patient_name,
	'force=s'  => \$force,
	
);


unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
#unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }
my $t = `hostname`;
my $nbErrors = 0;
my $buffer = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProject( -name => $project_name);
my $patient= $project->getPatient($patient_name);
my $bam = $patient->getBamFile();
my $res;
# comment on connait le strand ?
my $output = $project->getCallingPipelineDir("hificnv")."/".$patient->name;
my $fasta = $project->genomeFasta();
my $regtools = qq{singularity run --writable-tmpfs -B /data-isilon:/data-isilon -B /data-beegfs:/data-beegfs /software/distrib/regtools/regtools.sif regtools};
my $regtools1 = qq{ $regtools junctions extract -o $output -s FR $bam };
warn $regtools1;
system($regtools1);
die() unless -e  $output;
my $dir = $project->getJunctionsDir("regtools");
my $prod = "$dir/".$patient->name.".tsv";
my $gtf = $project->gtf_file;
my $regtools2 = qq{$regtools junctions annotate -o $prod $output $fasta $gtf};
warn $regtools2;
#echo "singularity run --writable-tmpfs --pwd /data-isilon/bipd-src/mhamici/Regtools/runs/NGS2024_7582 -B /data-isilon:/data-isilon -B /data-beegfs:/data-beegfs /data-isilon/bipd-src/mhamici/Regtools/regtools_v_1_0_0.sif regtools junctions extract -o extracted_ID29_EA1_Q_RNASEQ.bed -s FR /data-isilon/sequencing/ngs/NGS2024_7582/HG19_MT/align/star/ID29_EA1_Q_RNASEQ.bam " | run_cluster.pl -cpu=40
#echo "singularity run --writable-tmpfs --pwd /data-isilon/bipd-src/mhamici/Regtools/runs/NGS2024_7582 -B /data-isilon:/data-isilon -B /data-beegfs:/data-beegfs /data-isilon/bipd-src/mhamici/Regtools/regtools_v_1_0_0.sif regtools junctions annotate -o annotated_ID29_EA1_Q_RNASEQ.tsv /data-isilon/bipd-src/mhamici/Regtools/runs/NGS2024_7582/extracted_ID29_EA1_Q_RNASEQ.bed ../../data/HG19_MT/HG19_MT.fa /data-isilon/public-data/repository/HG19/annotations/gencode.v43/gtf/gencode.v43lift37.annotation.gtf " | run_cluster.pl -cpu=40
system($regtools2);
#tsv ?
