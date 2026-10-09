#!/usr/bin/env perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
#use Set::IntSpan;
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use Carp;

use autodie qw(system);


 my $project_name;
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
warn "coucou";

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );


my $singularity = "exec_singularity.sh";

my @ds;

my $deeptools = "deeptools.sif ";
my $patient = $project->getPatient($patient_name);
	my $run = $patient->getRun();
my $ref =  $project->genomeFasta();
my $dir = $project->getCoverageDir();
my $outf = $dir."/".$patient->name.".bw";
#mosdepth -t 8 -x -b 1000 -Q 20 "${out_path}/${sample_id}" "${bam_path}";
my $align = $patient->getAlignmentFile();

my $dir_out= $project->getAlignmentPipelineDir($patient->name);
my $bam_out = $patient->getBamFile();
#my $bam_out = $dir_out."/".$patient->name.".bam";
my $pname = $patient->name;
my $mosdepth_file =  "${dir_out}/${pname}.regions.bed.gz";
my $cmd = qq{$singularity mosdepth.sif mosdepth -t $fork -x -b 1000 -Q 20 $dir_out/${pname} ${bam_out} -f ${ref}};
system($cmd);# unless -e $mosdepth_file;
warn "----------------";
warn $cmd;
warn "------------------";

system("tabix -f -C -p bed $mosdepth_file");
my $ref               = $project->genomeFasta();
my $vcf = $patient->getVariationsFileName("deepvariant");
my $snf  = $dir_out."/".$patient->name.".snf";
my $snf_vcf_gz = $patient->getVariationsFileName("sniffles"); 
my $snf_vcf = $snf_vcf_gz;
$snf_vcf =~ s/\.gz//;

if  (-e $snf) {
	unlink $snf;
	unlink $snf_vcf;
	
}
my $cmd_sniffles = qq{$singularity sniffles.sif sniffles --input ${bam_out} --snf $snf --vcf $snf_vcf --threads $fork };
unless (-e $snf_vcf_gz){
system("$cmd_sniffles") unless -e $snf;
die() unless -e $snf;

system("bgzip -f $snf_vcf && tabix -f -p vcf $snf_vcf_gz ")  ;

}
my $snf_json  = $dir_out."/".$patient->name.".snf.json";
my $cmd_sniffles2 = qq{$singularity snf2json.sif snf2json $snf $snf_json };

system("$cmd_sniffles2")  unless -e $snf_json;

my $cmd_spectre = qq{$singularity spectre.sif spectre CNVCaller --coverage $mosdepth_file --sample-id $pname --output-dir $dir_out --reference $ref  --snfj $snf_json --metadata /data-bipd/data-pure/public-data/genome/HG38_DRAGEN/fasta/all.spectre.mdr --min-cnv-len 30000  --threads $fork};
#my $cmd = qq{$singularity $deeptools bamCoverage -b $align -o $outf  -p $fork --binSize 50 --normalizeUsing None --extendReads 0 --minMappingQuality 10};
warn $cmd_spectre;
my $spectre_vcf_gz = $patient->getVariationsFileName("spectre");

system($cmd_spectre);
my $vf = $dir_out."/".$pname.".vcf.gz";
die($vf) unless (-e $vf);
system("cp $vf $spectre_vcf_gz;tabix -f -p vcf $spectre_vcf_gz");
die() unless -e $spectre_vcf_gz.".tbi";
exit(0);


