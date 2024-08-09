#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/";
use lib "$Bin/../../../GenBo/lib/GenBoDB";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../../../GenBo/lib/kyoto/";
use lib "$Bin/../../../GenBo/lib/GenBoDB/writeDB";
use lib "$Bin/../../packages";
use lib "$Bin";
use epi2me_util; 
use Logfile::Rotate;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moo;
use GBuffer;
use GenBoProject;
use colored; 
use Config::Std;
use Text::Table;
use file_util;
use File::Temp qw/ tempfile tempdir /;; 

use Term::Menus;
use Proc::Simple;
use Storable;
use JSON::XS;
use Net::SSH::Perl; 
use Carp;
 

my $bin_cecile=qq{$Bin/scripts/scripts_db_polypipeline};
my $bin_script_pipeline = qq{$Bin/scripts/scripts_pipeline};


my $projectName;
my $filename;

my $name;
my $patients_name;
my $steps_name;
my $force;
my $type ="";
my $fastq_ext;
my $exclude_patients;
my $max_cpu ;
my $bds;
 my @running_steps;
 my $predef_steps;
 my $nocluster = 0;
my $low_calling;
my $predef_type;
my $define_steps;
my $step;
my $umi;

my $spipeline;

my $limit;
my $version;
my $rna;
my $phased; 
my $neb;
my $pad;
my $fork=20;

GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'command=s'=>\$spipeline,
	'rna=s' =>\$rna,
);
my $pipeline = {};
foreach my $l (split(",",$spipeline)){
	$pipeline->{$l} ++;
}
my $user = system("whoami");
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName , -version =>$version);
my $tm = $project->getPipelineDir();
#my $run_name = $project->getRun->infosRun->{run_name};

#system ("mkdir -p $dir_dragen/".$project->name );

my $patient = $project->getPatient($patients_name);
my $dir_pipeline = $tm;
my $prefix = $patient->name;
my $bam_prod = $patient->getBamFileName("epi2me");
#warn $bam_prod if -e $bam_prod;
#exit(0) if -e $bam_prod;
my $run_name = $patient->project->getRunFromId($patient->getRunId())->infosRun->{run_name};
my $bam_pipeline = $dir_pipeline."/".$prefix."-".$run_name."sorted.aligned.bam";

my $dir_pipeline_log = $patient->getEpi2meDir("pipeline.log");
my $ok_pipeline = $dir_pipeline_log."/".$prefix.".ok.pipeline.".time; 
my $ok_move = $dir_pipeline_log."/".$prefix.".ok.move.".time;
my $log_pipeline = $dir_pipeline_log."/".$prefix.".pipeline.".time.".log"; 
my $log_error_pipeline = $log_pipeline.".err"; 
if ($rna){
	run_pipeline_rna($pipeline);
}
else {
	run_pipeline($pipeline);# unless -e $bam_pipeline;
}
#die() unless  -e $bam_pipeline;
if ($rna == 1) {
	$spipeline.=",sj";
	
}
warn "move";
die "$Bin/epi2me_move.pl -project=$projectName -patient=$patients_name -command=$spipeline -rna=$rna -version=$version && touch $ok_move";
system("$Bin/epi2me_move.pl -project=$projectName -patient=$patients_name -command=$spipeline -rna=$rna -version=$version && touch $ok_move");
exit(0);

################################################
### PIPELINE RNA
################################################

sub run_pipeline_rna {
my ($pipeline) = @_;
my $param_align = "";
my $ref_fasta = $project->getGenomeFasta();
my $epi2me = $buffer->software("epi2me")."/wf-transcriptomes";
my $gtf = $project->gtf_file();
my $cmd_epi2me = qq{$epi2me --ref_genome $ref_fasta -profile singularity --out_dir $dir_pipeline --minimum_mapping_quality 20 --ref_annotation $gtf --threads $fork};
my $runid = $patient->getRun()->id;
my $bam   = $dir_pipeline."/".$patient->name.".bam";
my $fastq1;
my $buffer_ori = GBuffer->new();
my $project_ori = $buffer_ori->newProject( -name => $projectName );
my $patient_ori = $project_ori->getPatient($patients_name);

if (-e $patient->getBamFileName()) {
	my $opt = "--bam ";
	$bam = $patient_ori->getBamFile();
	$param_align = qq{ $opt $bam };
}
else {
	$fastq1 = epi2me_util::get_fastq_file($patient,$dir_pipeline);
	$param_align = " --fastq $fastq1 --prefix ".$patient->name();
}


if (exists $pipeline->{count}){
	my $gtf =  $project->gtf_file();
	die() unless -e $gtf;

}


##
$cmd_epi2me .= $param_align;
warn $cmd_epi2me;
#$patient->update_software_version("epi2me",$cmd_epi2me);
my $exit = system(qq{$Bin/../../../polyscripts/system_utility/run_cluster.pl -cmd=\"$cmd_epi2me\" -cpu=\"20\"}) ;#unless -e $f1;
unlink $fastq1 if  $fastq1 =~ /pipeline/;
die if $exit != 0;

}


#################################################
#### PIPELINE DNA
#################################################
sub run_pipeline {
my ($pipeline) = @_;
warn Dumper $pipeline;
my $param_align;
my $align;
my $ref_fasta = $project->getGenomeFasta();
my $name = $patient->name();
my $epi2me_align = $buffer->software("epi2me")."/wf-alignment";
my $epi2me_calling = $buffer->software("epi2me")."/wf-human-variation";
my $run_name = $patient->project->getRunFromId($patient->getRunId())->infosRun->{run_name};
my $bam_pipeline = $dir_pipeline."/".$patient->name."-".$run_name.".sorted.aligned.bam";
warn "il existe" if -e $bam_pipeline;

my $bam_prod = $patient->getBamFileName();
warn "il existe en prod " if -e $bam_prod;
if (exists $pipeline->{align}){
	unless(-e $bam_prod || -e $bam_pipeline){
		warn "quest ce que tu fous là";
		my $fastq1 = epi2me_util::get_fastq_file($patient,$dir_pipeline);
		confess() unless $fastq1;
		$param_align = "--fastq $fastq1 --references $ref_fasta  ";
		$align=1;
	}
}
if (-e $bam_pipeline){
	warn "pipeline";
	$param_align = "--bam $bam_pipeline";
	$align=0;
}
else {
	warn "prod";
	$param_align = "--bam $bam_prod";
	confess() unless -e $bam_prod;
	$align=0;
}
my $cmd_epi2me;
#my $cmd_epi2me_align = qq{$epi2me_align --out_dir  $dir_pipeline --prefix $name --threads $fork -profile singularity } if (exists $pipeline->{align}) || !$bam_prod;
#my $cmd_epi2me_calling = qq{$epi2me_calling --out_dir  $dir_pipeline --sample_name $name --threads $fork -profile singularity --ref $ref_fasta };

my $param_gvcf ="";
if (exists $pipeline->{gvcf} ){
	$param_gvcf = qq{--GVCF true } ;
}

my $param_vcf ="";
if (exists $pipeline->{vcf} ){
	$param_vcf = qq{--snp --phased} ;
}
my $param_cnv = "";
if (exists $pipeline->{cnv}){
	$param_cnv = qq{ --cnv };

}

my $param_sv = "";
if (exists $pipeline->{sv}){
	$param_sv = qq{ --sv  };
}

my $param_str = "";
if (exists $pipeline->{str}){
	$param_str = qq{ --str };

}
unless(-e $bam_prod || -e $bam_pipeline){
	$cmd_epi2me = qq{$epi2me_align --out_dir  $dir_pipeline --prefix $name --threads $fork -profile singularity };
	$cmd_epi2me .= $param_align." >$log_pipeline 2>$log_error_pipeline  && touch $ok_pipeline "
}
else{
	$cmd_epi2me = qq{$epi2me_calling --out_dir  $dir_pipeline --sample_name $name --threads $fork -profile singularity --ref $ref_fasta };
	$cmd_epi2me .= $param_align." ".$param_gvcf." ".$param_vcf." ".$param_cnv." ".$param_sv." ".$param_str  ." >$log_pipeline 2>$log_error_pipeline  && touch $ok_pipeline ";
}

#$cmd_epi2me_align .= $param_align." >$log_pipeline 2>$log_error_pipeline  && touch $ok_pipeline " if (exists $pipeline->{align}) || !$bam_prod;;
warn qq{$Bin/../../../polyscripts/system_utility/run_cluster.pl -cmd=\"$cmd_epi2me\" -cpu=\"20\"};
#$cmd_epi2me_calling .= $param_align." ".$param_gvcf." ".$param_vcf." ".$param_cnv." ".$param_sv." ".$param_str  ." >$log_pipeline 2>$log_error_pipeline  && touch $ok_pipeline ";
#my $cmd_epi2me = $cmd_epi2me_align;
#$cmd_epi2me = $cmd_epi2me_calling if -e $bam_prod || -e $bam_pipeline;

#$patient->update_software_version("dragen",$cmd_dragen);
my $exit = system(qq{$Bin/../../../polyscripts/system_utility/run_cluster.pl -cmd=\"$cmd_epi2me\" -cpu=\"20\"}) ;
#unlink $fastq1 if  $fastq1 =~ /pipeline/;
#my $exit = system(qq{$Bin/../run_dragen.pl -cmd=\"$cmd_dragen\"}) ;#unless -e $f1;
#die() unless -e $ok_pipeline;
die if $exit != 0;
}




1;


