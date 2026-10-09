#! /usr/bin/env perl
use strict;
use FindBin qw($Bin);
#use lib "/software/polyweb/poly-disk/poly-src/polygit/GenBo/lib/obj-nodb/";
use lib "$Bin/../GenBo/lib/obj-nodb/";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use Data::Dumper;
use Getopt::Long;
use GBuffer;
use autodie qw(system);
 use lib "$Bin/../polyutil";
 use slurm;
 my $arg_project_name;

my $date = `date`;
chomp($date);
my $fork = 64;
$fork =256 unless $fork;


my $force;
GetOptions(
	'project=s' => \$arg_project_name,
	'force=s' => \$force,
) or confess ("Error in command line arguments");

my $slurm = slurm->new();
my $jobs;
my $dir_pipeline_script = qq{$Bin/scripts/scripts_pipeline/pacbio/};
my $fork = 64;

my $nb=0;
my @steps = ("pbmm2","deepvariant","binary_depth","sawfish","wisecondor","spectre","hificnv","calling_wisecondor");

my @calling =("pbmm2","sawfish","wisecondor","spectre","hificnv","calling_wisecondor");
my @project_steps =("deepvariant_denovo");
my $scripts = {
		pbmm2 => {cmd=>"$dir_pipeline_script/pbmm2.pl"},
		deepvariant => {cmd=>"$dir_pipeline_script/deepvariant.pl",previous=>"pbmm2"},
		binary_depth => {cmd=>"$dir_pipeline_script/../coverage_genome.pl",previous=>"pbmm2"},
		sawfish => {cmd=>"$dir_pipeline_script/sawfish.pl",previous=>"pbmm2"},
		wisecondor =>  {cmd=>"$dir_pipeline_script/wisecondor.pl",previous=>"pbmm2"},
		calling_wisecondor =>  {cmd=>"$dir_pipeline_script/calling_wisecondor.pl",previous=>"wisecondor"},
		spectre =>  {cmd=>"$dir_pipeline_script/spectre.pl",previous=>"pbmm2"},
		hificnv =>  {cmd=>"$dir_pipeline_script/hificnv.pl",previous=>"pbmm2"},
		deepvariant_denovo => {cmd=>"$dir_pipeline_script/../deepvariant/deepvariant_denovo.pl",previous=>"deepvariant"},
};

foreach my $project_name (split(",",$arg_project_name)) {
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $groups;
my $stforce;

foreach my $c (@calling){
#	system("add_calling_methods.sh -project=".$project->name." -method=".$c);
}

my $per_fam;
my $deep;
my $hids;
 foreach my $patient (@{$project->getPatients}) {
	
	my $patient_name = $patient->name();
	foreach my $step (@steps){
		unless (-e $patient->getFileName($step) or $force) {
			my $cmd = $scripts->{$step}->{cmd}.qq{ -project=$project_name -patient=$patient_name -fork=$fork $stforce};#."&&  $Bin/bam2cram.pl -project=$project_name -patient=$patient_name -fork=$fork $stforce";
			my $hh ={cmd=>$cmd,name=>"${step}#".$project->name,type=>$patient->name,cpu=>128};
			if (exists $scripts->{$step}->{previous}){
				my $previous = $scripts->{$step}->{previous};
				if (exists $hids->{$previous}){
					$hh->{previous}=$hids->{$previous}->{$patient->id};
				}
			}
			
			my $id = $slurm->add_job($hh);
			push(@{$hids->{$step}->{$patient->id}}, $id);
		}
		} 
	}
	foreach my $step (@project_steps){
		unless (-e $project->getFileName($step) or $force ){
			if ($force) {
				unlink $project->getFileName($step);
			}
			my $cmd = $scripts->{$step}->{cmd}.qq{ -project=$project_name -fork=$fork $stforce};
			my $hh ={cmd=>$cmd,name=>"${step}#".$project->name,type=>$project->name,cpu=>128};
			if (exists $scripts->{$step}->{previous}) {
				my $previous = $scripts->{$step}->{previous};
				if (exists $hids->{$previous}) {
					foreach my $a  (values %{$hids->{$previous}}){
						#$hh->{previous}=$hids->{$previous}->{$patient->id};
						push(@{$hh->{previous}}, @$a);
					}
				}
			}
			my $id = $slurm->add_job($hh);
			$hids->{$step}->{$project->id} = $id;
		}
	}
}
$slurm->print_jobs();
$slurm->run_slurm;
