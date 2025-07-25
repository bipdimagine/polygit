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
 
my $ip_dragen = "10.1.2.10";

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
my $cram;
my $dragen_version;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'step=s'=> \$step,
	'type=s' => \$type,
	'command=s'=>\$spipeline,
	'version=s' =>\$version,
	'rna=s' =>\$rna,
	'cram=s' =>\$cram,
	'dragen_version' =>\$dragen_version,
	#'low_calling=s' => \$low_calling,
);
my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
 my $ssh = Net::SSH::Perl->new($ip_dragen);
$ssh->login("$username");


my $pipeline;
foreach my $l (split(",",$spipeline)){
	$pipeline->{$l} ++;
}


#my $user = system("whoami");
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName , -version =>$version);
$version = $project->genome_version unless $version;

my $tm = "/staging/tmp/";

#system ("mkdir -p $dir_dragen/".$project->name );

my $patient = $project->getPatient($patients_name);
my $dir_pipeline = $patient->getDragenDirName("pipeline");
my $cmd_dir = qq{test -d $dir_pipeline || mkdir -p $dir_pipeline};
my ($out, $err, $exit) = $ssh->cmd($cmd_dir);
#my $dir_pipeline = $patient->getDragenDirName("pipeline");
my $prefix = $patient->name;
my $bam_prod = $patient->getBamFileName("dragen-align");

my $url = qq{$username\@}.qq{$ip_dragen:};
$url ="";
#exit(0) if -e $bam_prod;
#warn "coucou";
if (exists $pipeline->{align}){
	
	my $bam_pipeline = $dir_pipeline."/".$prefix.".bam";
	if ($cram){
	$bam_pipeline =~ s/bam/cram/;
	$bam_prod =~ s/bam/cram/;
	
	}
#warn $bam_pipeline;
	($out, $err, $exit)=  $ssh->cmd("test -f $bam_pipeline");
	move_bam($bam_pipeline,$patient);# if ($version );
	my $dir_out = $patient->project->getAlignmentStatsDir("dragen-align");
	move_stats($dir_pipeline,$dir_out,$patient);
}
if (exists $pipeline->{gvcf}){
	my $gvcf_pipeline = "$dir_pipeline/".$prefix.".hard-filtered.gvcf.gz";
	($out, $err, $exit)=  $ssh->cmd("test -f $gvcf_pipeline");
	move_gvcf($gvcf_pipeline,$patient);
}
if (exists $pipeline->{vcf}){
	my $vcf_pipeline = "$dir_pipeline/".$prefix.".hard-filtered.vcf.gz";
	($out, $err, $exit)=  $ssh->cmd("test -f $vcf_pipeline");
	my $bcftools = $buffer->software("bcftools");
	my $tabix = $buffer->software("tabix");
	my $vcf2 = "$dir_pipeline/".$prefix.".soft-filtered.vcf.gz";
	system(qq{$bcftools view   -c 1  -i ' (QUAL<10 && (FORMAT/AF<0.3 || FORMAT/DP<10)) || (QUAL<20 && (FORMAT/AF<0.25 || FORMAT/DP<5)) || (QUAL<50 && FORMAT/DP<3 )' $vcf_pipeline -o $vcf2 -O z && $tabix -f -p vcf $vcf2});
	 ($out, $err, $exit)=  $ssh->cmd("test -f $vcf2.tbi");
	move_vcf($vcf2,$patient);
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
my $ploidy_file = $dir_pipeline."/".$prefix."ploidy.vcf.gz";
if (-e $ploidy_file){
	move_ploidy($ploidy_file,$patient);
}

if ( $rna ==1){
	my $sv_file = $dir_pipeline."/".$prefix.".SJ.out.tab";
	move_sj($sv_file,$patient);
}






#die($gvcf_pipeline ." probleme no gvcf") unless  $exit ==0;


#die($target_pipeline_gc ." probleme no target gc") unless  $exit ==0;
#die() unless -e $target_pipeline;






exit(0);



sub move_stats {
	my ($dir1,$dir2,$patient) = @_;
	system("rsync -rav $dir1/".$patient->name."*csv $dir2/ ");
	
}

sub move_file {
	my ($file,$patient,$method,$type,$idxtype,$backup) = @_;	
	
	my $physical_name = $patient->getPhysicalFileName($method,$version,$type);
	if($backup && -e $physical_name ) {
		backup($physical_name);
	}
	#--remove-source-files
	warn "rsync -rav  $url"."$file $physical_name ";
	system("rsync -rav  $url"."$file $physical_name ");
	system ("rsync -rav  $url"."$file.".$idxtype." $physical_name.".$idxtype);
	return $physical_name;
}

sub move_cram {
	my ($bam,$patient,$type,$idxtype) = @_;
	my $physical_name = move_file($bam,$patient,"dragen-align","cram","crai");
	my $prod = $patient->getCramFileName("dragen-align",$version);
	system("ln -sf  $physical_name $prod");
	system("ln -sf  $physical_name.crai $prod.crai");
	 my $filename = $prod;
	 $filename =~ s/cram/idxstats/;
	 my $samtools = $buffer->software("samtools");
	 system("$samtools idxstats $physical_name  -\@ 20  > $filename");
}


sub move_bam {
	my ($bam,$patient) = @_;
	my $prod = $patient->getBamFileName("dragen-align");
	if( $cram){
		move_cram($bam,$patient);
		return;
	}
	my $physical_name = move_file($bam,$patient,"dragen-align","bam","bai");
	system("ln -sf  $physical_name $prod");
	system("ln -sf  $physical_name.bai $prod.bai");
}


sub move_gvcf {
	my ($gvcf,$patient) = @_;
	my $physical_name = move_file($gvcf,$patient,"dragen-calling","gvcf.gz","tbi","backup");
	my $prod = $patient->gvcfFileName("dragen-calling");
	
	system("ln -sf  $physical_name $prod");
	tabix($prod,"vcf");
}

sub move_ploidy {
	my ($vcf,$patient) = @_;
	my $physical_name = move_file($vcf,$patient,"dragen-ploidy","vcf.gz","tbi");
	
	my $prod = $patient->vcfFileName("dragen-ploidy");
	system("ln -sf  $physical_name $prod");
	tabix($prod,"vcf");
}

sub move_vcf {
	my ($vcf,$patient) = @_;
	my $physical_name = move_file($vcf,$patient,"dragen-calling","vcf.gz","tbi","backup");
	
	my $prod = $patient->getVariationsFileName("dragen-calling");
	my $prod = $patient->vcfFileName("dragen-calling");
	system("ln -sf  $physical_name $prod");
	tabix($prod,"vcf");
}
sub move_count {
	my ($t1,$t2,$patient) = @_;
	my $dir = $patient->project->getTargetCountDir();
	system("rsync -rav  $url"."$t1 $dir/");
	system("rsync -rav  $url"."$t2 $dir/");
}


sub move_cnv {
	my ($t1,$patient) = @_;
	my $physical_name = move_file($t1,$patient,"dragen-cnv","vcf.gz","tbi","backup");
	my $dir = $patient->project->getVariationsDir("dragen-cnv");
	my $prod = "$dir/".$prefix.".cnv.vcf.gz";
	system("ln -fs  $physical_name $prod");
	tabix($prod,"vcf");
	
}
sub tabix {
	my ($file,$type) = @_;
	my $tabix = $buffer->software("tabix");
	warn "$tabix -f -p $type $file";
	system("$tabix -f -p $type $file");
}
sub move_sv {
	my ($t1,$patient) = @_;
	my $dir = $project->getVariationsDir("dragen-sv");
	my $physical_name = move_file($t1,$patient,"dragen-sv","vcf.gz","tbi","backup");
	my $prod = $dir."/".$prefix.".sv.vcf.gz";
	system("ln -sf  $physical_name $prod");
	tabix($prod,"vcf");
	
}
sub move_sj {
	my ($t1,$patient) = @_;
	my $dir = $project->getJunctionsDir("dragen-sj");
	system("rsync -rav  $url"."$t1 $dir/");
	#system("rsync -rav  $url"."$t1.tbi $dir/");
	
}
sub backup {
	my ($final_gz) = @_;
	my $dir =  dirname($final_gz);
	my $dirb = $dir."/.backup";
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
