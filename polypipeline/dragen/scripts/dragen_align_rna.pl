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
use Logfile::Rotate;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
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

 
 my $projectName;
 my $patient_name;
 
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patient_name,
	#'low_calling=s' => \$low_calling,
);


my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
my $patient = $project->getPatient($patient_name);
my $bam_prod = $patient->getBamFileName("dragen-align");
warn $bam_prod if -e $bam_prod;
exit(0) if -e $bam_prod;
my $dir_pipeline = $patient->getDragenDir("pipeline");

my ($fastq1,$fastq2) = get_fastq_file($patient,$dir_pipeline);
my $ref_dragen = $project->getGenomeIndex("dragen");
my $gtf = "/data-isilon/public-data/repository/HG19/gtf/gencode.v40/gencode.gtf.gz";
my $runid = $patient->getRun()->id;
my $cmd = qq{dragen -f \
-r $ref_dragen
-1 $fastq1 \
-2 $fastq2 \
-a $gtf \
--enable-map-align true \
--enable-sort=true \
--enable-bam-indexing true \
--enable-map-align-output true \
--output-format=BAM \
--RGID=$runid \
--RGSM=$patient_name \
--config-file /opt/edico/config/dragen-user-defaults.cfg \
--enable-rna=true \
--output-directory $dir_pipeline \
--output-file-prefix $patient_name};


sub get_fastq_file {
	my ($patient,$dir_pipeline) = @_;
	
	my $name=$patient->name();
	
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
	my $fastq1 = $dir_pipeline."/".$patient->name.".R1.fastq.gz";
	my $fastq2 = $dir_pipeline."/".$patient->name.".R2.fastq.gz";
	#if ($step eq "align"){
		system "cat $cmd1 > $fastq1";# unless -e $fastq1;
		system "cat $cmd2 > $fastq2";# unless -e $fastq2;
	#}
	return  ($fastq1,$fastq2);
}

