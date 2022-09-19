#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
#$Bin="/data-isilon/bipd-src/pnitschk/git/repository/polypipeline/dragen/";
use lib "/software/polyweb/poly-disk/poly-src/polygit/GenBo//lib/";
use lib "/software/polyweb/poly-disk/poly-src/polygit/GenBo//lib/obj-nodb/";
use lib "/software/polyweb/poly-disk/poly-src/polygit/GenBo//lib/obj-nodb/packages/";
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
#use file_util;
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

GetOptions(
	'project=s' => \$projectName,
#	'family=s' => \$familyName,
	'patients=s' => \$patients_name,
);


my $user = system("whoami");
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
#my $patients = $project->get_only_list_patients($patients_name);
my $patients = $project->getPatients();
my $tm = "/staging/tmp/";

my @gvcfs;
my $dir_pipeline =  $project->getCallingPipelineDir("dragen-calling");#$project->project_dragen_pipeline_path()."/pipeline/";
my $file = $dir_pipeline."/".$project->name.'.'.time;
my $file2 = $dir_pipeline."/".$project->name.'.pedigree.'.time;
open(FILE,">$file");
open(PED,">$file2");

#foreach my $patient (sort{$a->getFamily->name cmp $b->getFamily->name} @$patients){
##	warn $patient->getGvcfFile("haplotypecaller");
#	die unless -e $patient->getGvcfFile("haplotypecaller4");
#	print FILE $patient->getGvcfFile("haplotypecaller4")."\n";
#	#print PED $patient->pedigreeLine."\n";
#}
#close(FILE);
#close(PED);
my $ref_dragen = $project->getGenomeIndex("dragen");

my @lgroups = @{$project->getSomaticGroups()};
my @cmds;
foreach my $group (@lgroups) {
	my $group_name = $group->name();
	foreach my $ill (@{$group->getPatientsSomatic()}){
		my $ill_bam =$ill->getBamFile();
		foreach my $healthy (@{$group->getPatientsGerminal()}){
			my $healthy_bam = $healthy->getBamFiles->[0];
			my $cmd = "dragen -f -r $ref_dragen --enable-map-align false --output-directory $dir_pipeline --output-file-prefix $group_name --enable-variant-caller true --tumor-bam-input $ill_bam -b $healthy_bam ";
		#	push(@cmds,$cmd);
			my $gvcf_pipeline = "$dir_pipeline/".$group_name.".hard-filtered.gvcf.gz";
			$exit = system(qq{$Bin/../run_dragen.pl -cmd=\"$cmd\"})unless -e $gvcf_pipeline;
			warn qq{$Bin/../run_dragen.pl -cmd=\"$cmd"};
			die($cmd) unless $exit == 0;
			move_gvcf($gvcf_pipeline,$group);
			exit(0);
		}
	
	}
}
#my %pat_dejavu;
#foreach my $patient (sort{$a->getFamily->name cmp $b->getFamily->name} @$patients){
#	next() if $pat_dejavu{$patient->name()} == 1;
#	my $group = $patient->getSomaticGroup();
#	my @pat_healthy = grep{$_->isHealthy || $_->getSomaticGroup() eq $group } @$patients;
#	my @pat_ill = grep{$_->isHealthy || $_->getSomaticGroup() eq $group } @$patients;
#	$pat_dejavu{$patient->name()} == 1;
#}

#foreach my $gr (keys{$project->getGroups()}){

#	foreach my $gr_pat (@{$gr->getPatients()}){
#		my @tum = grep { $_->is}
#		foreach my $tum ()
#	}
#}

#my $cmd = "dragen -f -r $ref_dragen --output-directory $dir_pipeline --output-file-prefix $projectName --enable-joint-genotyping true --variant-list $file";
#my $filein = $dir_pipeline."/$projectName.hard-filtered.vcf.gz";
#warn $filein;
foreach my $c (@cmds){
#	print $c."\n";
#	my $exit =0; 
##	$exit = system(qq{/data-isilon/bipd-src/pnitschk/git/repository/polypipeline/dragen/run_dragen.pl -cmd=\"$c\"}) ;
#	$exit = system(qq{$Bin/../run_dragen.pl -cmd=\"$cmd\"}) unless -e $gvcf_pipeline;
#	warn qq{$Bin/../run_dragen.pl -cmd=\"$cmd"};
#	warn "end";
#	die() if $exit != 0;
}

#my $bin_dev = "$Bin/../scripts/scripts_pipeline/";
#my $cmd1 = "perl $bin_dev/correct_gatk.pl -vcf=$filein -project=$projectName -dir=$dir_pipeline  -patient=$patients_name -method=dragen-calling";
#system($cmd1);
#my $cmd2 = "perl $bin_dev/move_vcf.pl  -project=$projectName -vcf_dir=$dir_pipeline -patient=$patients_name -method_calling=dragen-calling";
#warn $cmd2;
#system($cmd2);


exit(0);


