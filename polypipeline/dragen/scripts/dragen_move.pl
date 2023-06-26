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

my $limit;
my $version;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'step=s'=> \$step,
	'type=s' => \$type,
	'command=s'=>\$spipeline,
	'version=s' =>\$version,
	#'low_calling=s' => \$low_calling,
);
my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
 my $ssh = Net::SSH::Perl->new("10.1.2.9");
$ssh->login("$username");


my $pipeline;
foreach my $l (split(",",$spipeline)){
	$pipeline->{$l} ++;
}


#my $user = system("whoami");
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName , -version =>$version);
my $tm = "/staging/tmp/";

#system ("mkdir -p $dir_dragen/".$project->name );

my $patient = $project->getPatient($patients_name);
my $dir_pipeline = $patient->getDragenDirName("pipeline");
my $cmd_dir = qq{test -d $dir_pipeline || mkdir -p $dir_pipeline};
my ($out, $err, $exit) = $ssh->cmd($cmd_dir);
#my $dir_pipeline = $patient->getDragenDirName("pipeline");
my $prefix = $patient->name;
my $bam_prod = $patient->getBamFileName("dragen-align");

my $url = qq{$username\@10.200.27.109:};
$url ="";
#exit(0) if -e $bam_prod;
#warn "coucou";
if (exists $pipeline->{align}){
	my $bam_pipeline = $dir_pipeline."/".$prefix.".bam";
	if ($version && !(-e $bam_pipeline)){
	$bam_pipeline =~ s/bam/cram/;
	$bam_prod =~ s/bam/cram/;
	}
#warn $bam_pipeline;
	($out, $err, $exit)=  $ssh->cmd("test -f $bam_pipeline");
	move_bam($bam_pipeline,$patient);
}

if (exists $pipeline->{gvcf}){
	my $gvcf_pipeline = "$dir_pipeline/".$prefix.".hard-filtered.gvcf.gz";
	($out, $err, $exit)=  $ssh->cmd("test -f $gvcf_pipeline");
	move_gvcf($gvcf_pipeline,$patient);
}
if (exists $pipeline->{vcf}){
	my $vcf_pipeline = "$dir_pipeline/".$prefix.".hard-filtered.vcf.gz";
	($out, $err, $exit)=  $ssh->cmd("test -f $vcf_pipeline");
	move_vcf($vcf_pipeline,$patient);
}
if (exists $pipeline->{cnv}){
	my $target_pipeline  ="$dir_pipeline/".$prefix.".target.counts.gz";
	my $target_pipeline_gc  = "$dir_pipeline/".$prefix.".target.counts.gc-corrected.gz";
	($out, $err, $exit)=  $ssh->cmd("test -f $target_pipeline_gc");
	move_count($target_pipeline,$target_pipeline_gc,$patient);
	if($project->isGenome){
		my $cnv_file  = "$dir_pipeline/".$prefix.".cnv.vcf.gz";
		move_cnv($cnv_file,$patient)
	}
}
if (exists $pipeline->{sv}){

	my $sv_file = $dir_pipeline."/".$prefix.".sv.vcf.gz";
	move_sv($sv_file,$patient);
}

if (exists $pipeline->{dragen_count}){
	warn "toto";

	my $count_file = $dir_pipeline."/".$prefix.".quant.sf";
	my $splice_file = $dir_pipeline."/".$prefix.".SJ.out.tab";
	warn $splice_file;
	move_rna_count($count_file,$splice_file,$patient);
}







#die($gvcf_pipeline ." probleme no gvcf") unless  $exit ==0;


#die($target_pipeline_gc ." probleme no target gc") unless  $exit ==0;
#die() unless -e $target_pipeline;







exit(0);






sub move_bam {
	my ($bam,$patient) = @_;
	my $prod = $patient->getBamFileName("dragen-align");
	 $prod = $patient->getCramFileName("dragen-align") if $bam =~ /cram/;
	system("rsync -rav  $url"."$bam $prod ");
	system("rsync -rav  $url"."$bam.bai $prod.bai ");
	system("rsync -rav  $url"."$bam.cai $prod.cai ");
	
}

sub move_gvcf {
	my ($gvcf,$patient) = @_;
	my $prod = $patient->gvcfFileName("dragen-calling");
	backup($prod) if -e $prod;
	system("rsync -rav  $url"."$gvcf $prod");
	system("rsync -rav  $url"."$gvcf.tbi $prod.tbi");
}
sub move_vcf {
	my ($vcf,$patient) = @_;
	my $prod = $patient->getVariationsFileName("dragen-calling");
	my $prod = $patient->vcfFileName("dragen-calling");

	backup($prod) if -e $prod;
	system("rsync -rav  $url"."$vcf $prod");
	system("rsync -rav  $url"."$vcf.tbi $prod.tbi");
}
sub move_count {
	my ($t1,$t2,$patient) = @_;
	my $dir = $patient->project->getTargetCountDir();
	system("rsync -rav  $url"."$t1 $dir/");
	system("rsync -rav  $url"."$t2 $dir/");
}


sub move_cnv {
	my ($t1,$patient) = @_;
	my $dir = $patient->project->getVariationsDir("dragen-cnv");
	system("rsync -rav  $url"."$t1 $dir/");
	system("rsync -rav  $url"."$t1.tbi $dir/");
}

sub move_sv {
	my ($t1,$patient) = @_;
	my $dir = $project->getVariationsDir("dragen-sv");
	system("rsync -rav  $url"."$t1 $dir/");
	system("rsync -rav  $url"."$t1.tbi $dir/");
	
}

sub move_rna_count {
	my ($t1,$t2,$patient) = @_;
	my $dir = $project->getCountingDir("dragen-count");
	system("rsync -rav  $url"."$t1 $dir/");
	system("rsync -rav  $url"."$t2 $dir/");
	
}

sub move_rna_metrics {
	my ($t1,$patient) = @_;
	my $dir = $project->getCountingDir("dragen-count");
	system("rsync -rav  $url"."$t1 $dir/");
	
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

