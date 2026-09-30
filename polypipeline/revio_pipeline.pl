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
foreach my $project_name (split(",",$arg_project_name)) {
	warn $project_name;
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $groups;
my $stforce;

my @steps = ("pbmm2","deepvariant","binary_depth","sawfish","wisecondor","spectre","hificnv","calling_wisecondor");




 foreach my $patient (@{$project->getPatients}) {
	#binary_depth,coverage,melt,deepvariant
	my $patient_name = $patient->name;
	
	my $previous_dude;
	#Coverage 
		
	my $align = qq{$dir_pipeline_script/pbmm2.pl -project=$project_name -patient=$patient_name -fork=$fork $stforce}."&&  $Bin/bam2cram.pl -project=$project_name -patient=$patient_name -fork=$fork $stforce";
	my $id1;
	
	
	unless (-e $patient->getFileName("pbmm2") or $force ){
			
			$id1 = $slurm->add_job({cmd=>$align,name=>"pbmm2!".$project->name,type=>$patient->name,cpu=>128});
	}
	unless (-e $patient->getFileName("deepvariant")){
			my $cmd_deepvariant = qq{$dir_pipeline_script/deepvariant.pl -project=$project_name -patient=$patient_name -fork=$fork $stforce};
			 $slurm->add_job({cmd=>$cmd_deepvariant,name=>"deep!".$project->name,type=>$patient->name,cpu=>$fork,previous=>[$id1]});
	}
	unless (-e $patient->getFileName("binary_depth")){
			my $cmd_coverage = qq{perl $dir_pipeline_script/../coverage_genome.pl -project=$project_name -patient=$patient_name -fork=$fork -$stforce};
			$slurm->add_job({cmd=>$cmd_coverage,name=>"binary-depth!".$project->name,type=>$patient->name,cpu=>$fork,previous=>[$id1]});
			
	}
	unless (-e $patient->getFileName("sawfish")){
			my $cmd_sawfish = qq{$dir_pipeline_script/sawfish.pl -project=$project_name -patient=$patient_name -fork=$fork $stforce};
			$slurm->add_job({cmd=>$cmd_sawfish,name=>"sawfish!".$project->name,type=>$patient->name,cpu=>$fork,previous=>[$id1]});
	}
	my $wise_id;
	unless (-e $patient->getFileName("wisecondor")){
		
			my $cmd_wsiecondor = qq{$dir_pipeline_script/wisecondor.pl -project=$project_name -patient=$patient_name -fork=$fork $stforce};
			$wise_id = $slurm->add_job({cmd=>$cmd_wsiecondor,name=>"wise1!".$project->name,type=>$patient->name,cpu=>$fork,previous=>[$id1]});
	}
	unless (-e $patient->getFileName("spectre")){
			my $cmd_wsiecondor = qq{$dir_pipeline_script/spectre.pl -project=$project_name -patient=$patient_name -fork=$fork $stforce};
			$wise_id = $slurm->add_job({cmd=>$cmd_wsiecondor,name=>"spectre!".$project->name,type=>$patient->name,cpu=>$fork,previous=>[$id1]});
	}
	unless (-e $patient->getFileName("hificnv")){
			my $cmd_wsiecondor = qq{$dir_pipeline_script/hificnv.pl -project=$project_name -patient=$patient_name -fork=$fork $stforce};
			$wise_id = $slurm->add_job({cmd=>$cmd_wsiecondor,name=>"hificnv!".$project->name,type=>$patient->name,cpu=>$fork,previous=>[$id1]});
	}
	unless (-e $patient->getFileName("calling_wisecondor")){
		my $cmd_wsiecondor2 = qq{$dir_pipeline_script/calling_wisecondor.pl -project=$project_name -patient=$patient_name -fork=$fork $stforce};
		$slurm->add_job({cmd=>$cmd_wsiecondor2,name=>"wiseC!".$project->name,type=>$patient->name,cpu=>$fork,previous=>[$wise_id]});
	}

 }
}

$slurm->print_jobs();
die();
$slurm->run_slurm;
