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
use Text::CSV;


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
my $project_names;
GetOptions(
	'project=s'   => \$project_names,
	"fork=s" => \$fork,
);
die("miss fork") unless $fork;
my $singularity= "singularity run --bind /data-pure:/data-pure --bind /data-isilon:/data-isilon --bind /data-beegfs:/data-beegfs " ;
my @ds;
my $sawfish = "/data-pure/software/SINGULARITY/sawfish.sif ";
my $dir_final = "/data-beegfs/sawfish/";
system ("mkdir  $dir_final");
my $filename = "/data-beegfs/sawfish/all.csv";
	
my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
		
	my $fh;		
	open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
foreach my $project_name (split(",",$project_names)){
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );






my $dir_pipeline2= $project->getCallingPipelineDir($project->name.".sawfish.calling");


my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
my %uniq;
foreach my $patient (@{$project->getPatients}){
	next if exists $uniq{$patient->name};
	$uniq{$patient->name}++;
	my $run = $patient->getRun();
my $ref =  $project->genomeFasta();
my $bam = $patient->getAlignFileName("pbmm2",1) ;
next unless -e $bam;
my $bam = $patient->getAlignmentFile() ;

my $npz =  $patient->fileWiseCondor();
my $dir_pipeline= $project->getCallingPipelineDir($patient->name.".sawfish");
	$csv->print($fh, [$dir_pipeline,$bam]); 

}

}
close($fh);
my $dir_pipeline2= $dir_final;
my $cmd = qq{$singularity $sawfish sawfish joint-call --threads $fork --sample-csv $filename  --output-dir $dir_pipeline2 };
warn $cmd;
system($cmd);
system(qq{bcftools query -l $dir_pipeline2/genotyped.sv.vcf.gz | while read S; do   bcftools view     -s "\$S"   -f PASS  -e 'GT="./."'     -Oz     -o \${S}.vcf.gz    genotyped.sv.vcf.gz;    bcftools index -t \${S}.vcf.gz; done});

warn $cmd;

