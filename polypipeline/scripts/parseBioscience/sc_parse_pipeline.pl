#!/usr/bin/perl
use Data::Dumper;
use File::Find;
use Getopt::Long;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/";
use lib "$Bin/../../../GenBo/lib/GenBoDB";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../packages";
use GBuffer;

my $analyse;
my $file;
my $projectName;
my $no_exec;
my $fork =20;
my $p;

GetOptions(
	'project=s' => \$projectName,
	'patient=s' => \$p,
	'no_exec=s' => \$no_exec,
	'fork=s' => \$fork,
 );

 
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
my $patient;
my %poolName;

foreach my $p (@{$project->getPatients()}){
	$patient = $p ;
	my $group = $patient->somatic_group();
	$poolName{$group}=1;
}


my $dir_pipeline = $project->getAlignmentPipelineDir("split-pipe");
my ($fastq1,$fastq2) = get_fastq_file($patient,$dir_pipeline);

sub get_fastq_file {
	my ($patient,$dir_pipeline) = @_;
	my $name=$patient->name();
	my $dir_fastq = $patient->getSequencesDirectory() ;
	#$dir_fastq = $patient->getSequencesDirectory() unless $dir_fastq;
	warn $dir_fastq;
	my $group = $patient->somatic_group();
	my @r1 = glob($dir_fastq."/*".$group."*R1*.fastq.gz");
	my @r2 = glob($dir_fastq."/*".$group."*R2*.fastq.gz");
	my $cmd1 = join(" ",@r1);
	warn $cmd1;
	my $cmd2 = join(" ",@r2);
	warn $cmd2;
	my $fastq1 = $group."_R1_L001.fastq.gz";
	my $fastq2 = $group."_R2_L001.fastq.gz";
	system "cat $cmd1 > ".$dir_pipeline."/".$fastq1 unless -e $dir_pipeline."/".$fastq1;
	system "cat $cmd2 > ".$dir_pipeline."/".$fastq2 unless -e $dir_pipeline."/".$fastq2;
	return  ($fastq1,$fastq2);
}

#singularity run -B /data-isilon/sequencing/ngs/NGS2024_7385/:/PROJECT/ -B /data-isilon/public-data/genome/HG38-EBV/split-pipe/:/REF/ -B /data-isilon/data/sequences/ILLUMINA/NOVASEQ
#/IMAGINE/run456_20231121/:/FASTQ/ /software/distrib/ParseBiosciences-Pipeline/ParseBiosciences-Pipeline.1.1.2/split-pipe.sif split-pipe --mode all --genome_dir /REF --output_dir /P
#ROJECT/multi_library1 --sample LH C7-C11 --sample AG_HD3_IFNa C2-C6 --sample Lantaz D9-D12 --sample AG_HD2_IFNa B4-B8 --sample Gauth D5-D8 --sample AG_HD1_NS A1-A5 --sample AG_HD2_
#NS A11-B3 --sample AG_HD1_IFNa A6-A10 --sample AG_HD3_NS B9-C1 --sample VO C12-D4 --fq1 /FASTQ/multi_library1_R1.fastq.gz --fq2 /FASTQ/multi_library1_R2.fastq.gz --chemistry v2


#my @lPat;
my $patients = $project->get_only_list_patients($p);
my $dir = $project->getProjectRootPath();
my $index = $project->getGenomeIndex("split-pipe");
#my $cmd = "module load python/3.8.2 ; split-pipe --mode all  --genome_dir $index  --output_dir $dir  ";
my @cond;
my $cmd;
my $fastq_dir;
my @names;
foreach my $p (@$patients){
	my $bc = $p->barcode();
	$bc =~ s/-/:/;
	$fastq_dir = $dir_pipeline;
	my $xcond = "--sample ".$p->name." ".$bc;
	push(@cond, $xcond);
	push(@names, "/PROJECT/".$p->name());
	
	
	  
    #--fq1 /newvolume/expdata/202206_ex1/s1_S1_R1_001.fastq.gz \
   # --fq2 /newvolume/expdata/202206_ex1/s1_S1_R2_001.fastq.gz \
   
}


$cmd = qq{singularity run -B  $dir:/PROJECT/ -B $index:/REF/ -B $fastq_dir:/FASTQ/ /software/distrib/ParseBiosciences-Pipeline/ParseBiosciences-Pipeline.1.1.2/split-pipe.sif split-pipe --mode all --genome_dir /REF } ;
$cmd .= join(" ",@cond);
$cmd .= " --fq1 /FASTQ/".$fastq1." --fq2 /FASTQ/".$fastq2." --chemistry v2 --output_dir /PROJECT/"    ;
warn $cmd;


my $sublib = join (" ",@names);
my $cmd2 = qq{singularity run -B  $dir:/PROJECT/ -B $index:/REF/ -B $fastq_dir:/FASTQ/ /software/distrib/ParseBiosciences-Pipeline/ParseBiosciences-Pipeline.1.1.2/split-pipe.sif split-pipe --mode comb --sublibraries $sublib  --output_dir /PROJECT/} ;
warn $cmd2;
