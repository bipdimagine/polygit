#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use FindBin qw($Bin);
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
 
 my $project_name;
my $final_vcf;
my $log_file;
my $patient_name;
my $fork;

GetOptions(
	'project=s'   => \$project_name,
	"log=s" =>\$log_file,
	"vcf=s" => \$final_vcf,
	"patient=s"=>\$patient_name,
	"fork=s" =>\$fork,
);
my $date = `date`;
chomp($date);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $gene_name = "RNU7-1";
#my $gene = $project->newGene($gene_name);
my $patient = $project->getPatient($patient_name);
my $bcftools = $buffer->software("bcftools");
my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");
my $dir_gvcf_out= $project->getCallingPipelineDir($patient->name.".".$gene_name.".deepvariant");

my $bed = $dir_gvcf_out."/".$patient->name.".bed";
my $chromosome = $project->getChromosome(12);

open (BED , ">".$bed);
	print BED $chromosome->fasta_name."\t7052979	7053041\n";
	print BED $project->getChromosome(15)->fasta_name()."\t32718089\t32718621\n";
	#print BED $gene->getChromosome->fasta_name."\t".$gene->start."\t".$gene->end."\n";
close(BED);
warn $bed;
my $ref               = $project->genomeFasta();
my $bam                 = $patient->getBamFile();
my $vcf_out = $dir_gvcf_out."/".$patient->name.".deep.vcf.gz";
my $dir_gvcf_tmp = $dir_gvcf_out."/tmp.".time;

mkdir $dir_gvcf_tmp;
my $cmd;

 $cmd = qq{singularity run  --bind /data-isilon:/data-isilon --bind /data-beegfs:/data-beegfs  /software/distrib/deepvariant/deepvariant.sif /opt/deepvariant/bin/run_deepvariant  --num_shards=$fork --model_type=WES --ref=$ref --reads=$bam --regions=$bed --output_vcf=$vcf_out};

system($cmd);
die() unless -e $vcf_out;
my $vcf = $patient->getVariationsFileName("deepvariant-rnu7");

my $cmd2 =qq{bcftools view  $vcf_out -o $vcf -O z && tabix -f -p vcf $vcf};

system($cmd2);

die() unless -e $vcf.".tbi";
warn $vcf;