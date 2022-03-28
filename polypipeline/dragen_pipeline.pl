#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo/lib/";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb/";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../GenBo/lib/kyoto/";
use lib "$Bin/../GenBo/lib/GenBoDB/writeDB";
use lib "$Bin/packages";
use Logfile::Rotate;
 use Cwd;
use PBS::Client;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
use MooseX::Method::Signatures;
use bds_steps;   
use file_util;
use Class::Inspector;
use Digest::MD5::File ;
use GBuffer;
use GenBoProject;
use colored; 
use Config::Std;
use Text::Table;

use File::Temp qw/ tempfile tempdir /;; 

use Term::Menus;
 use Proc::Simple;
 use Storable;
use JSON::XS;
 
 

 
  
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
 my $myproc;
my $low_calling;
my $predef_type;
my $define_steps;
my $step;
my $dir_dragen = "/data-dragen/run";


my $limit;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'step=s'=> \$step,
	'type=s' => \$type,
	#'low_calling=s' => \$low_calling,
);


my $user = system("whoami");

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
system ("mkdir -p $dir_dragen/".$project->name );
my $dir2 = "$dir_dragen/".$project->name;
$patients_name = $project->get_only_list_patients($patients_name);
$patients_name = "all" unless $patients_name;

foreach my $patient (@{$patients_name}) {
	my $name=$patient->name();
	warn $patient->name();
	my $prod_bam = 	$patient->getBamFileName("dragen-align");
	my $prod_vcf = $patient->getVariationsFileName("dragen-calling");
	my $prod_vcf_index = $prod_vcf.".tbi";
	my $prod_index = $prod_bam.".bai";
	next() if -e $prod_index && $step eq "align";
	next() if -e $prod_vcf_index && $step eq "calling";


	my $dir3 = $dir2."/".$patient->name();
	warn $patient->name();
	system("mkdir $dir3");
	my $files_pe1 = file_util::find_file_pe($patient,"");
	my $cmd;
	my @r1;
	my @r2;
	foreach my $cp (@$files_pe1) {
		my $file1 = $patient->getSequencesDirectory()."/".$cp->{R1};
		my $file2 = $patient->getSequencesDirectory()."/".$cp->{R2};
		push(@r1,$file1);
		push(@r2,$file2);
	}
	my $cmd1 = join(" ",@r1);
	my $cmd2 = join(" ",@r2);
	my $fastq1 = $dir3."/".$patient->name.".R1.fastq.gz";
	my $fastq2 = $dir3."/".$patient->name.".R2.fastq.gz";
	if ($step eq "align"){
		system "cat $cmd1 > $fastq1" unless -e $fastq1;
		system "cat $cmd2 > $fastq2" unless -e $fastq2;
	}
	my $pname = $patient->name;
	warn "cat $cmd2 > $dir3/".$patient->name.".R2.fastq.gz";
	my $fasta;
	$fasta = qq{/data-dragen/human/reference/HG19_BIPD/HG19_BIPD.fa.k_21.f_16.m_149} if $type eq "exome";
	$fasta = qq{/data-dragen/human/reference/HG19_CNG/fasta/HG19_CNG.fa.k_21.f_16.m_149} if $type eq "genome";
	#my $fasta = qq{/data-dragen/human/reference/HG38/fasta/HG38.fa.k_21.f_16.m_149/};
#	my $fasta = qq{/data-dragen/mouse/reference};
	my $bed1 = $patient->getCaptureFile();
	my $bed2 = "$dir3/".$patient->name.".bed";
	system("zcat $bed1 > $bed2") if $step eq "align";
	warn $bed2;
	my $dragen_bam = 	 $dir3."/$pname.bam";
	my $cmd11;
	$cmd11 =qq{ dragen -f  -r $fasta  -1 $fastq1  -2 $fastq2  --vc-emit-ref-confidence GVCF --enable-variant-caller true  --RGID Illumina_RGID  --RGSM $pname  --output-directory $dir3 --output-file-prefix $pname  --enable-duplicate-marking true  --enable-map-align-output true --vc-target-bed $bed2 --vc-target-bed-padding 250 } if $type eq "exome";
	$cmd11 =qq{ dragen -f  -r $fasta  -1 $fastq1  -2 $fastq2  --vc-emit-ref-confidence GVCF --enable-variant-caller true  --RGID Illumina_RGID  --RGSM $pname  --output-directory $dir3 --output-file-prefix $pname  --enable-duplicate-marking true  --enable-map-align-output true  } if $type eq "genome";
	
	my $t = time;
	my $dev_index = $dragen_bam.".bai";
	if ($step eq "align"){
	system ("ssh ".$user."\@10.200.27.109 ".$cmd11) unless -e $dev_index ;
	}
	my $dragen_vcf = $dir3."/$pname.vcf.gz";
	my $t1 = abs(time -$t);
	warn "--> $t1";

	if ($step eq "align"){
		my $prod_bam = 	$patient->getBamFileName("dragen-align");
		warn "$dragen_bam $prod_bam";
		system("mv $dragen_bam $prod_bam");
		system("mv $dragen_bam.bai $prod_bam.bai");
		my $dragen_gvcf = 	 $dir3."/$pname.hard-filtered.gvcf.gz";
		my $prod_gvcf = 	$patient->gvcfFileName("dragen-calling");
		warn "$dragen_gvcf $prod_gvcf";
		system("mv $dragen_gvcf $prod_gvcf");
		system("mv $dragen_gvcf.tbi $prod_gvcf.tbi");
	}
	
	if ($step eq "calling"){
		my $cmd22;
		$cmd22 = qq{dragen -f -r $fasta -b $prod_bam --enable-variant-caller true  --output-file-prefix $pname --output-directory $dir3  --enable-map-align false --enable-sort false --vc-target-bed $bed2 --vc-target-bed-padding 250} if $type eq "exome";
		$cmd22 = qq{dragen -f -r $fasta -b $prod_bam --enable-variant-caller true  --output-file-prefix $pname --output-directory $dir3  --enable-map-align false --enable-sort false } if $type eq "genome";

		warn $cmd22;
		$t = time;
		system ("ssh ".$user."\@10.200.27.109 ".$cmd22);
		my $t2 = abs(time -$t);
		die($dragen_vcf) unless -e $dragen_vcf;
		my $dragen_vcf = $dir3."/$pname.vcf.gz";
		die($dragen_vcf) unless -e $dragen_vcf;
		my $prod_vcf = $patient->getVariationsFileName("dragen-calling");
		warn "$dragen_vcf $prod_vcf";
		system("mv $dragen_vcf $prod_vcf");
		system("mv $dragen_vcf.tbi $prod_vcf.tbi");
	}

	exit(0);
}







