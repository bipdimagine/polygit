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
use dragen_util; 
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
my $cram;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'umi=s' => \$umi,
	'command=s'=>\$spipeline,
	'version=s' =>\$version,
	'rna=s' =>\$rna,
	'phased=s' => \$phased,
	'neb=s' => \$neb,
	'padding=s' => \$pad,
	'cram=s' => \$cram,
	#'low_calling=s' => \$low_calling,
);
my $pipeline = {};

my $dragen = "dragen";
foreach my $l (split(",",$spipeline)){
	$pipeline->{$l} ++;
}
my $user = system("whoami");
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName , -version =>$version);

 my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
 my $ssh = Net::SSH::Perl->new($buffer->config->{dragen}->{ip});
 $ssh->login("$username");



my $tm = "/staging/tmp/";

if ($project->isGenome){
	warn "coucou";
}
#system ("mkdir -p $dir_dragen/".$project->name );

my $patient = $project->getPatient($patients_name);
my $dir_pipeline = $patient->getDragenDir("pipeline");
my $prefix = $patient->name;
my $bam_prod = $patient->gvcfFileName("dragen-calling");
#warn $bam_prod if -e $bam_prod;
#exit(0) if -e $bam_prod;

my $bam_pipeline = $dir_pipeline."/".$prefix.".bam";

my $dir_pipeline_log = $patient->getDragenDir("pipeline.log");
my $ok_pipeline = $dir_pipeline_log."/".$prefix.".ok.pipeline.".time; 
my $ok_move = $dir_pipeline_log."/".$prefix.".ok.move.".time;
my $log_pipeline = $dir_pipeline_log."/".$prefix.".pipeline.".time.".log"; 
my $pangenome ="";
 #$pangenome = "--validate-pangenome-reference=false" 
 
 $dragen = "/opt/dragen/4.3.16/bin/dragen" if $project->genome_version =~ /HG19/;
 
 my($version_str, $stderr, $exit2) = $ssh->cmd("$dragen --version");
my @a = split("\n",$version_str);
my @b = split(" ",$a[0]);

my $dragen_version =$b[-1];

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

warn "$Bin/dragen_move.pl -project=$projectName -patient=$patients_name -command=$spipeline -rna=$rna -version=$version -dragen_version=$dragen_version && touch $ok_move";
system("$Bin/dragen_move.pl -project=$projectName -patient=$patients_name -command=$spipeline -rna=$rna -cram=$cram -version=$version -dragen_version=$dragen_version && touch $ok_move");
exit(0);

################################################
### PIPELINE RNA
################################################

sub run_pipeline_rna {
my ($pipeline) = @_;
my $param_align = "";
my $ref_dragen = $project->getGenomeIndex("dragen");
my $param_umi = "";
my $tmp = "/staging/tmp";
my $cmd_dragen = qq{$dragen -f -r $ref_dragen --output-directory $dir_pipeline --intermediate-results-dir $tmp --output-file-prefix $prefix };
my $runid = $patient->getRun()->id;
my $bam   = $dir_pipeline."/".$patient->name.".bam";
my ($fastq1,$fastq2);
my $align = "true";
if (-e $bam){
	$align = "false";
	warn "## $align";
	return;
}

my $buffer_ori = GBuffer->new();
my $project_ori = $buffer_ori->newProject( -name => $projectName );
my $patient_ori = $project_ori->getPatient($patients_name);

if ($version && exists $pipeline->{align} ){
	$patient_ori->{alignmentMethods} =['dragen-align','hisat2'];
	my $bamfile = $patient_ori->getBamFile();
	$param_align = "-b $bamfile --enable-map-align-output true  --enable-rna=true ";
	#$param_align .= "--output-format CRAM " if $version =~/HG38/;
	
}	

elsif (-e $patient->getBamFileName()) {
	my $opt = "--bam-input";
	$bam = $patient_ori->getBamFile();
	#$param_align = qq{ $opt $bam --enable-map-align false --enable-map-align-output false };
	$param_align = qq{ $opt $bam };
}
else {
	
	if($neb){
		my $dir_fastq = $project->getAlignmentPipelineDir("dragen-align");
		#my $dir_fastq = $dir_pipeline."/";
		warn $dir_fastq;
		warn "###NEB!!!!!";
		($fastq1,$fastq2) = dragen_util::get_fastq_file($patient,$dir_pipeline,$dir_fastq) ;
	}
	else{
		 ($fastq1,$fastq2) = dragen_util::get_fastq_file($patient,$dir_pipeline);
	}


 $param_align = " -1 $fastq1 -2 $fastq2 --RGID $runid  --RGSM $prefix --enable-map-align-output $align --enable-rna=true ";


}

if ($umi){
	$param_umi = "--umi-enable true ";
}
elsif(-e $patient->getBamFileName()){
	$param_align .= "--enable-duplicate-marking false ";
}
else{
	$param_align .= "--enable-duplicate-marking true ";
}

if (exists $pipeline->{count}) {
	my $gtf =  $project->gtf_file();
	die() unless -e $gtf;

	$param_align .= "-a $gtf --enable-rna-quantification true";
#PROD
#	$gtf = qq{/data-isilon/public-data/repository/HG19/annotations/gencode.v43/gtf/gencode.v43lift37.annotation.gtf};
#	$param_align .= "-a $gtf --enable-rna-quantification true  --rna-ann-sj-min-len 4";
#
#
}
#
my $param_calling ="";
if (exists $pipeline->{vcf} ){
	$param_calling = qq{--enable-variant-caller true  } ;
	
}


##
$cmd_dragen .= $param_umi." ".$param_align." ".$param_calling;
$patient->update_software_version("dragen",$cmd_dragen);
warn qq{$Bin/../run_dragen.pl -cmd=\"$cmd_dragen\"};
my $exit = system(qq{$Bin/../run_dragen.pl -cmd=\"$cmd_dragen\"}) ;#unless -e $f1;
unlink $fastq1 if  $fastq1 =~ /pipeline/;
unlink $fastq2 if $fastq2 =~ /pipeline/;;
die if $exit != 0;

}


################################################
### PIPELINE DNA
################################################
sub run_pipeline {
my ($pipeline) = @_;

my $param_align;
my $ref_dragen = $project->getGenomeIndex("dragen");
my $param_umi = "";
my $param_phased = ""; 
my ($fastq1,$fastq2);
if (exists $pipeline->{align}){
	warn "1";
	 ($fastq1,$fastq2) = dragen_util::get_fastq_file($patient,$dir_pipeline);
	 warn "2";
	 
}
if (!($fastq1)) {
	my $bamfile = $patient->getBamFileName("bwa");
	$param_align = "-b $bamfile --enable-map-align-output true --enable-duplicate-marking true ";
}
if ($version && exists $pipeline->{align} && !($fastq1)){
	my $buffer_ori = GBuffer->new();
	my $project_ori = $buffer_ori->newProject( -name => $projectName );
	my $patient_ori = $project_ori->getPatient($patients_name);
	
	$patient_ori->{alignmentMethods} =['dragen-align','bwa'];
	my $bamfile = $patient_ori->getBamFile();
	$param_align = "-b $bamfile --enable-map-align-output true --enable-duplicate-marking true ";# unless $cram;
	$param_align .= " --output-format CRAM " if $cram;
	
	if ($umi){
		confess();
	}
}	
elsif (exists $pipeline->{align}){
my ($fastq1,$fastq2) = dragen_util::get_fastq_file($patient,$dir_pipeline);

#	confess() unless $fastq1;
	if ($fastq1) {
		my $runid = $patient->getRun()->id;
		$param_align = " -1 $fastq1 -2 $fastq2 --RGID $runid  --RGSM $prefix --enable-map-align-output true ";
	
		if ($umi){
			$param_align .= qq{ --umi-enable true   --umi-library-type random-simplex  --umi-min-supporting-reads 1 --vc-enable-umi-germline true};
		}
		elsif($patient->getCapture->isPcr()){
			$param_align .= qq{ --enable-duplicate-marking false} ;
		}
		else {
			$param_align .= qq{ --enable-duplicate-marking true };
		}
	}
	else {
		my $bamfile = $patient->getBamFileName("bwa");
		confess("No fastq or (bwa) bam found") unless -e $bamfile;
		$param_align = "-b $bamfile --enable-map-align-output true --enable-duplicate-marking true ";
		 $param_align .= " --output-format CRAM " if $cram;
	}
	$param_align .= " --output-format CRAM " if $cram;
	
}
else {
	my $bam = $patient->getBamFileName();
	my $opt = "--bam-input";
	$bam = $patient->getCramFileName("dragen-align") if $cram;
	warn $bam;
	unless (-e $bam){
		$bam = $patient->getDragenDir("pipeline")."/".$patient->name.".bam";
		$bam = $patient->getDragenDir("pipeline")."/".$patient->name.".cram" if $cram;
		confess() unless -e $bam;
	}

	$opt = "--cram-input" if $bam =~ /cram/;
	$param_align = qq{ $opt $bam --enable-map-align false --enable-map-align-output false };
}
my $param_gvcf = "";
my $tmp = "/staging/tmp";

my $cmd_dragen = qq{$dragen -f -r $ref_dragen --output-directory $dir_pipeline --intermediate-results-dir $tmp --output-file-prefix $prefix };
if (exists $pipeline->{gvcf}){
	
	$param_gvcf = qq{--vc-emit-ref-confidence GVCF } ;
	
}
my $param_bed ="";

unless ($version){
	unless ($project->isGenome ) {
		my $capture_file  = $patient->getCapture->gzFileName();
		confess() unless -e $capture_file;
		if ($project->isExome){
			my $capture_file2 = "$dir_pipeline/".$patient->name.".bed";
			warn "zcat $capture_file /data-isilon/public-data/capture/HG38/exons.HG38.bed.gz | sort -k1,1V -k2,2n | bedtools merge >$capture_file2";
			system ("zcat $capture_file /data-isilon/public-data/capture/HG38/exons.HG38.bed.gz | sort -k1,1V -k2,2n | bedtools merge >$capture_file2" );
			$capture_file = $capture_file2;
		}
		#test prenantome 
		$param_bed .= qq{  --vc-target-bed $capture_file};
		if (defined($pad)){
			$param_bed .= qq{ --vc-target-bed-padding $pad };
		}
		else{
			$param_bed .= qq{  --vc-target-bed-padding 150  };
		}
		
#		if($patient->getCapture->isPcr()){
#			$param_bed .= qq{--enable-dna-amplicon true --amplicon-target-bed $capture_file} ;
#			
#		}
	}
}
my $param_vcf ="";
if (exists $pipeline->{vcf} && exists $pipeline->{gvcf}){
	
	$param_vcf = qq{--vc-enable-vcf-output true } ;
	
}
my $param_calling ="";
if (exists $pipeline->{vcf} or exists $pipeline->{gvcf}){
	$param_calling = qq{--enable-variant-caller true } ;
	
}

my $param_cnv = "";
if (exists $pipeline->{cnv}){
	$param_cnv = qq{ --enable-cnv true --cnv-enable-self-normalization true };

}

my $param_sv = "";
if (exists $pipeline->{sv}){
	$param_sv = qq{   --enable-sv true };
	unless ($project->isGenome) {
	 $param_sv .= " --sv-exome true ";
	}
}

my $param_str = "";
if (exists $pipeline->{str}){
	$param_str = qq{ --repeat-genotype-enable true };

}
$param_phased = "--vc-combine-phased-variants-distance ".$phased if $phased;


$cmd_dragen .= $param_umi." ".$param_align." ".$param_calling." ".$param_gvcf." ".$param_vcf." ".$param_cnv." ".$param_bed." ".$param_sv." ".$param_phased." ".$param_str." $pangenome >$log_pipeline 2>$log_error_pipeline  && touch $ok_pipeline ";
warn qq{$Bin/../run_dragen.pl -cmd=\"$cmd_dragen\"};


$patient->update_software_version("$dragen",$cmd_dragen,$dragen_version);

my $exit = system(qq{$Bin/../run_dragen.pl -cmd=\"$cmd_dragen\"}) ;#unless -e $f1;
die() unless -e $ok_pipeline;

#system("ssh pnitschk\@10.200.27.109 ". $cmd." >$dir_pipeline/dragen.stdout 2>$dir_pipeline/dragen.stderr");
#system("ssh pnitschk\@10.200.27.109 rm $fastq1 $fastq2");
#my ($out,$err,$exit) = $ssh->cmd("$cmd") ;#unless -e $bam_pipeline;
die if $exit != 0;
}




1;


