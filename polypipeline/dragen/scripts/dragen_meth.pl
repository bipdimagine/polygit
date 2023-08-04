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
use Logfile::Rotate;
use dragen_util; 
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

 
my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
my $url = qq{$username\@10.200.27.109};

my $projectName;
my $patients_name;
my $protocol;
my $umi;

GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'protocol=s' => \$protocol,
	'umi=s' => \$umi,
);

die("no protocol defined => pbat, directionnal or non_directionnal") unless $protocol;

my $user = system("whoami");
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
my $patients = $project->get_only_list_patients($patients_name);
#my $patients = $project->getPatients();
my $tmp = "/staging/tmp/";


my $ref_dragen = $project->getGenomeIndex("dragen-meth");


my @cmds;
foreach my $p (@$patients) {
	my $bam_prod = $p->getBamFileName("dragen-meth");
	my $name = $p->name();
	next() if -e $bam_prod;
	my $exit =0; 
	my $directional="";
	my $runid = $p->getRun()->id;
	my $dir_pipeline = $p->getDragenDir("pipeline");
	my ($fastq1,$fastq2) = dragen_util::get_fastq_file($p,$dir_pipeline);
	my $cmd = "dragen -f --enable-sort true --enable-bam-indexing true --enable-duplicate-marking false --enable-methylation-calling true --methylation-protocol $protocol --methylation-generate-mbias-report true --methylation-generate-cytosine-report true --intermediate-results-dir $tmp --ref-dir $ref_dragen --RGID=$runid --RGSM=$name --RGPL illumina --RGPU 1 -1 $fastq1 -2 $fastq2 --output-directory $dir_pipeline --output-file-prefix $name ";
	if ($umi){
		$cmd .= qq{ --umi-enable true   --umi-library-type random-simplex  --umi-min-supporting-reads 1 };
	}
	else {
		$cmd .= qq{ --enable-duplicate-marking true };
	 }
	#--methylation-protocol $directional
	#my $gvcf_pipeline = "$dir_pipeline/".$name.".hard-filtered.gvcf.gz";
#	warn $gvcf_pipeline;
	$exit = system(qq{$Bin/../run_dragen.pl -cmd=\"$cmd\"}) unless -e $bam_prod;
	warn qq{$Bin/../run_dragen.pl -cmd=\"$cmd"};
	die($cmd) unless $exit == 0;
	#	move_gvcf($gvcf_pipeline,$group);
	#
	
}
exit(0);

