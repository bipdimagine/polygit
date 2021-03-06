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

my $projectName;
my $patients_name;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
);


my $user = system("whoami");
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
my $patients = $project->get_only_list_patients($patients_name);
my $tm = "/staging/tmp/";

my @gvcfs;
my $dir_pipeline =  $project->getCallingPipelineDir("dragen-calling");#$project->project_dragen_pipeline_path()."/pipeline/";
my $file = $dir_pipeline."/".$project->name.'.'.time;
my $file2 = $dir_pipeline."/".$project->name.'.pedigree.'.time;
open(FILE,">$file");
open(PED,">$file2");
foreach my $patient (sort{$a->getFamily->name cmp $b->getFamily->name} @$patients){
	die unless -e $patient->getGvcfFile("dragen-calling");
	print FILE $patient->getGvcfFile("dragen-calling")."\n";
	print PED $patient->pedigreeLine."\n";
}
close(FILE);
close(PED);
my $ref_dragen = $project->getGenomeIndex("dragen");
my $cmd = "dragen -f -r $ref_dragen --output-directory $dir_pipeline --output-file-prefix $projectName --enable-joint-genotyping true --variant-list $file";
my $filein = $dir_pipeline."/$projectName.hard-filtered.vcf.gz";
warn $cmd;
warn $filein;
my($stdout, $stderr, $exit) = $ssh->cmd($cmd ."&& touch $dir_pipeline/dragen-genotype.ok" ) unless -e $filein;
warn "end";
#warn $stderr;

die() unless -e "$dir_pipeline/dragen-genotype.ok";
my $bin_dev = "$Bin/../../scripts/scripts_pipeline/";
my $cmd1 = "perl $bin_dev/correct_gatk.pl -vcf=$filein -project=$projectName -dir=$dir_pipeline  -patient=$patients_name -method=dragen-calling";

#system($cmd1);
my $cmd2 = "perl $bin_dev/move_vcf.pl  -project=$projectName -vcf_dir=$dir_pipeline -patient=$patients_name -method_calling=dragen-calling";
warn $cmd2;
#system($cmd2);

exit(0);









