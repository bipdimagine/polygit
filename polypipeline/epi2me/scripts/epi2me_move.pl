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
use Getopt::Long;
use Data::Dumper;
use GBuffer;
use GenBoProject;
use colored; 
use Text::Table;
use file_util;
use File::Temp qw/ tempfile tempdir /;; 
use Logfile::Rotate; 
use File::Basename;
use Net::SSH::Perl;
 


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

my $spipeline;
my $rna;
my $limit;
my $version;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'step=s'=> \$step,
	'type=s' => \$type,
	'command=s'=>\$spipeline,
	'version=s' =>\$version,
	'rna=s' =>\$rna,
	#'low_calling=s' => \$low_calling,
);
#my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
# my $ssh = Net::SSH::Perl->new("10.1.2.9");
#$ssh->login("$username");


my $pipeline;
foreach my $l (split(",",$spipeline)){
	$pipeline->{$l} ++;
}


#my $user = system("whoami");
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName , -version =>$version);
my $tm = $project->getPipelineDir();
my $run_name = $project->getRun->infosRun->{run_name};

my $patient = $project->getPatient($patients_name);
my $dir_pipeline = $tm;
my $cmd_dir = qq{test -d $dir_pipeline || mkdir -p $dir_pipeline};
my ($out, $err, $exit) = system($cmd_dir);
#my $dir_pipeline = $patient->getDragenDirName("pipeline");
my $prefix = $patient->name;
my $bam_prod = $patient->getBamFileName("epi2me");

if (exists $pipeline->{align}){
	

	my $bam_pipeline = $dir_pipeline."/".$prefix."-".$run_name.".sorted.aligned.bam";
	$bam_pipeline = $dir_pipeline."/".$prefix."_reads_aln_sorted.bam" if $rna;
#	$bam_pipeline =~ s/_reads_aln_sorted// if $rna;
	
	($out, $err, $exit)=  system("test -f $bam_pipeline");
	move_bam($bam_pipeline,$patient);# if ($version );
	my $dir_out = $patient->project->getAlignmentStatsDir("epi2me");
	move_stats($dir_pipeline,$dir_out,$patient);
}
if (exists $pipeline->{gvcf}){
	my $gvcf_pipeline = "$dir_pipeline/".$prefix.".wf.gvcf.gz";
	($out, $err, $exit)=  system("test -f $gvcf_pipeline");
	move_gvcf($gvcf_pipeline,$patient);
}
if (exists $pipeline->{vcf}){
	my $vcf_pipeline = "$dir_pipeline/".$prefix.".wf_snp.vcf.gz";
	($out, $err, $exit)=  system("test -f $vcf_pipeline");
	move_vcf($vcf_pipeline,$patient);
}
if (exists $pipeline->{cnv}){
	my $target_pipeline  ="$dir_pipeline/".$prefix.".wf_cnv.vcf.gz";
	($out, $err, $exit)=  system("test -f $target_pipeline");
	move_cnv($target_pipeline,$patient);
#	if($project->isGenome){
#		my $cnv_file  = "$dir_pipeline/".$prefix.".cnv.vcf.gz";
#		move_cnv($cnv_file,$patient)
#	}
}
if (exists $pipeline->{sv}){
	my $sv_file = $dir_pipeline."/".$prefix.".sv.vcf.gz";
	move_sv($sv_file,$patient);
}
if (exists $pipeline->{str}){
	my $str_file = $dir_pipeline."/".$prefix.".wf_str.vcf.gz";
	move_str($str_file,$patient);
}


if ( $rna ==1){
	my $sv_file = $dir_pipeline."/".$prefix.".SJ.out.tab";
	move_sj($sv_file,$patient);
}


exit(0);



sub move_stats {
	my ($dir1,$dir2,$patient) = @_;
	system("rsync -rav $dir1/$prefix-wf-*-report.html $dir2/ ");
	
}


sub move_bam {
	my ($bam,$patient) = @_;
	my $prod = $patient->getBamFileName("epi2me");
	$prod = $patient->getCramFileName("epi2me") if $bam =~ /cram/;
	system("rsync -rav --remove-source-files $bam $prod ");
	system("rsync -rav  $bam.bai $prod.bai ");
	system("rsync -rav  $bam.cai $prod.cai ") if $bam =~ /cram/;
	
}

sub move_gvcf {
	my ($gvcf,$patient) = @_;
	my $prod = $patient->gvcfFileName("epi2me-calling");
	backup($prod) if -e $prod;
	system("rsync -rav  $gvcf $prod");
	system("rsync -rav  $gvcf.tbi $prod.tbi");
}
sub move_vcf {
	my ($vcf,$patient) = @_;
	my $prod = $patient->getVariationsFileName("epi2me-calling");
	my $prod = $patient->vcfFileName("epi2me-calling");

	backup($prod) if -e $prod;
	system("rsync -rav  $vcf $prod");
	system("rsync -rav  $vcf.tbi $prod.tbi");
}
sub move_count {
	my ($t1,$t2,$patient) = @_;
	my $dir = $patient->project->getTargetCountDir();
	system("rsync -rav  $t1 $dir/");
	system("rsync -rav  $t2 $dir/");
}

sub move_cnv {
	my ($t1,$patient) = @_;
	my $dir = $patient->project->getVariationsDir("epi2me-cnv");
	system("rsync -rav  $t1 $dir/");
	system("rsync -rav  $t1.tbi $dir/");
}

sub move_sv {
	my ($t1,$patient) = @_;
	my $dir = $project->getVariationsDir("epi2me-sv");
	warn "rsync -rav  $t1 $dir/";
	system("rsync -rav  $t1 $dir/");
	system("rsync -rav  $t1.tbi $dir/");
	
}

sub move_sv {
	my ($t1,$patient) = @_;
	my $dir = $project->getVariationsDir("epi2me-str");
	warn "rsync -rav  $t1 $dir/";
	system("rsync -rav  $t1 $dir/");
	system("rsync -rav  $t1.tbi $dir/");
	
}

sub move_sj {
	my ($t1,$patient) = @_;
	my $dir = $project->getJunctionsDir("epi2me-sj");
	warn $dir;
	system("rsync -rav  $t1 $dir/");
	#system("rsync -rav  $t1.tbi $dir/");
	
}
sub backup {
	my ($final_gz) = @_;
	my $dir =  dirname($final_gz);
	my $dirb = $dir."/backup";
	unless (-e $dirb){
		mkdir $dirb unless -e $dirb;
		system("chmod a+w $dirb");
	}
	if (-e $final_gz){
	my $log = new Logfile::Rotate( File => $final_gz,	
	 								Count => 5,
	 								Gzip  => 'no',
	 								 Flock  => 'no',
	 								 Dir => $dirb
	 								);
		$log->rotate();	 								
		}
		unlink $final_gz;
		unlink $final_gz.".tbi" if -e $final_gz.".tbi";
}
