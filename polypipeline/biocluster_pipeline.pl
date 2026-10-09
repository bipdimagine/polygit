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
 use YAML qw(LoadFile);
 use Term::ANSIColor qw(color colored);
 
 my $arg_project_name;

my $date = `date`;
chomp($date);
my $fork = 64;
$fork =256 unless $fork;
my $pipeline_name;

my $force;
GetOptions(
	'project=s' => \$arg_project_name,
	'pipeline=s' => \$pipeline_name,
	'force=s' => \$force,
) or confess ("Error in command line arguments");

my $slurm = slurm->new();
my $jobs;
my $dir_pipeline_script = qq{$Bin/scripts/scripts_pipeline/};


my $fork = 64;

my $nb=0;

my $jj;
foreach my $project_name (split(",",$arg_project_name)) {
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $groups;
my $stforce;

my $per_fam;
my $deep;
my $hids;
my $yaml;
if ($project->isGenome){
 $yaml = "$Bin/scripts/config_pipeline/${pipeline_name}/genome.yaml";
die("$yaml not found") unless -e $yaml;
}
elsif ($project->isExome){
	$yaml = "$Bin/scripts/config_pipeline/${pipeline_name}/exome.yaml";
}
else {
	die("only genome for now " );
}


my $config = LoadFile($yaml);
my @steps = @{$config->{steps}};
my @project_steps = @{$config->{project_steps}};
my $scripts = $config->{scripts};
warn Dumper $scripts;
my $run; 
 foreach my $patient (@{$project->getPatients}) {
	
	my $patient_name = $patient->name();
	foreach my $step (@steps){
		if ($force && -e $patient->getFileName($step)){
			unlink  $patient->getFileName($step);
		}
		unless (-e $patient->getFileName($step) or $force) {
			warn Dumper $scripts->{$step};
			my $cmd = $dir_pipeline_script.$scripts->{$step}->{script}.qq{ -project=$project_name -patient=$patient_name -fork=$fork};#."&&  $Bin/bam2cram.pl -project=$project_name -patient=$patient_name -fork=$fork $stforce";
			my $hh ={cmd=>$cmd,name=>"${step}#".$project->name,type=>$patient->name,cpu=>128};
			if (exists $scripts->{$step}->{previous}){
				my $previous_string = $scripts->{$step}->{previous};
				foreach my $previous (split(",",$previous_string)) {
					if (exists $hids->{$previous}){
						push(@{$hh->{previous}}, $hids->{$previous}->{$patient->id});
					}
				}
			}
			else {
				print colored("\n => skip $step for $patient_name \n",'magenta');
			}
			
			my $id = $slurm->add_job($hh);
			$hids->{$step}->{$patient->id} = $id;
			#push(@{$hids->{$step}->{$patient->id}}, $id);
		}
		} 
	}
	foreach my $step (@project_steps){
		unless (-e $project->getFileName($step) or $force ){
			if ($force) {
				unlink $project->getFileName($step);
			}
			my $cmd = $dir_pipeline_script.$scripts->{$step}->{script}.qq{ -project=$project_name -fork=$fork $stforce};
			my $hh ={cmd=>$cmd,name=>"${step}#".$project->name,type=>$project->name,cpu=>128};
			if (exists $scripts->{$step}->{previous}) {
				my $previous = $scripts->{$step}->{previous};
				if (exists $hids->{$previous}) {
					foreach my $a  (values %{$hids->{$previous}}){
						#$hh->{previous}=$hids->{$previous}->{$patient->id};
						push(@{$hh->{previous}}, $a);
					}
				}
			}
			my $id = $slurm->add_job($hh);
			$hids->{$step}->{$project->id} = $id;
		}
	}
}
warn $slurm->{jobs};
unless ($slurm->{jobs}){
	print colored("\n *** Nothing to do... your carbon footprint thanks you ... ***\n",'green');
	exit(0);
}
$slurm->print_jobs();
$slurm->run_slurm;
