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

 
 my $ssh = Net::SSH::Perl->new("10.200.27.109");
$ssh->login("pnitschk");

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
my ($fastq1,$fastq2) = get_fastq_file($patient,$dir_pipeline);
exit(0);





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


