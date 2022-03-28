#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";

use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use colored;

use Storable qw(store retrieve freeze);
use Term::ANSIColor;

my $filein;
my $dir;
my $file_bed;
 
my $project_name;
my $patient_name;
my $bam_file;
my $log_file;
my $vcf_final;
my $type;
my $fork;


GetOptions(
	'project=s'   => \$project_name,
	'bam=s' => \$bam_file,
	"log=s" =>\$log_file,
	"patient=s" => \$patient_name,
	"fork=s" => \$fork
);
$fork =1 unless $fork;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $ymd = sprintf("%04d_%02d_%02d_%02d_%02d",$year+1900,$mon+1,$mday,$hour,$min);




my $buffer = GBuffer->new();
 my $dir_trash =  $buffer->config->{project_data}->{root}."/TRASH/";
 system("mkdir -p $dir_trash && chmod a+rwx $dir_trash") unless -e $dir_trash;
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatientOrControl($patient_name);
my$m = $patient->alignmentMethod();
my $dir_prod = $project->getAlignmentDir($m);

my $bam_prod = $patient->getBamFileName();
my $bai_prod = $bam_prod.".bai";
my $bai_prod_bad = $dir_prod."/".$patient->name().".bai";
my $samtools = $buffer->getSoftware("samtools");
my $sambamba = $buffer->getSoftware("sambamba");

if (-e $bam_prod ){
	my $trash_file = $dir_trash."/".$patient->name().".$ymd.bam";
	system("chmod a+w $bam_prod ");
	system("chmod a+w $bai_prod ");
	system("rm $bai_prod ");
	system  ("chmod a+w $bam_prod && mv $bam_prod $trash_file");
	unlink $bai_prod if -e $bai_prod;
	unlink $bai_prod_bad if -e $bai_prod_bad;
}
die("can't move bam file probably protect") if -e $bam_prod ;

system("mv $bam_file $bam_prod  ");

my $bai = $bam_file.".bai";

if (-e $bai){
	system("mv $bai $bai_prod");
}
else {
	$bai = $bam_file;
	$bai =~s/\.bam/\.bai/;
	if (-e $bai){
		system("mv $bai $bai_prod");
	}
	else {
		system("$sambamba index -t $fork $bam_prod");
	}
}

system("chmod a-w $bam_prod && touch $bai_prod");


colored::stabilo('green'," ---- END MV BAM : ".$patient->name."----");
colored::stabilo('green',"  ----------------------------------------------------");
exit(0);
