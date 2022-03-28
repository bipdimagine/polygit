#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo/lib/obj-nodb/";
use lib "$Bin/packages_old";
use GBuffer;
#use GenBoProject;
#use PBS::Client;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
use MooseX::Method::Signatures;
use bds_cache_steps;    
#use file_util;
use Class::Inspector;
#use check_utils;
use Text::Table;
use colored; 
use Term::Menus;
 
my ($projectName, $filename, $name, $patients_name, $steps_name, $force, $type, $fastq_ext, $somatic, $method, $no_cluster, $stdout, $help);
my $fork = 1;
my $yes = 0;
my $nocluster = 0;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'steps=s' => \$steps_name,
	'force=s' => \$force,
	'yes=s' => \$yes,
	'nocluster=s' => \$nocluster,
	'h!' => \$help,
);

#$steps_name = "all" unless $steps_name;
my $report;
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );


unless ($steps_name){
	if ( $project->isDiagnostic()){
		$steps_name = "diag";
	}
	else {
		$steps_name = "query";
	}
}

my $predef_steps;
$predef_steps->{query} = ["global_infos", "store_ids", "store_annotations", "check_store_annotations", "strict_denovo", "check_strict_denovo", "loh", "check_loh","coverage","quality_check"];
$predef_steps->{genome} = ["global_infos", "store_ids", "store_annotations", "check_store_annotations", "strict_denovo", "check_strict_denovo", "loh", "check_loh","quality_check"];
$predef_steps->{query_somatic} = ["global_infos", "store_ids", "store_annotations", "check_store_annotations", "strict_denovo", "check_strict_denovo", "loh", "check_loh","coverage","quality_check"];
$predef_steps->{query_fam} = ["global_infos", "store_ids", "store_annotations", "check_store_annotations", "strict_denovo", "check_strict_denovo", "loh", "check_loh","coverage","quality_check"];
$predef_steps->{diag} = ["global_infos", "store_ids", "polydiag", "store_annotations", "check_store_annotations", "strict_denovo", "check_strict_denovo", "loh", "check_loh","coverage","quality_check"];
$predef_steps->{coverage} = [ "coverage"];
$predef_steps->{polydiag} = [ "polydiag","quality_check"];
$predef_steps->{dejavu} = [ "dejavu" ];
$predef_steps->{quality_check} = [ "quality_check" ];
#$predef_steps->{local_config} = [ "local_config" ];
#$predef_steps->{diag} = [ "coverage"];
#$predef_steps->{all} = ["polydiag"];

my $pipeline = bds_cache_steps->new( project=>$project, nocluster=>$nocluster );

my $steps = {
				"store_ids"=>  sub {$pipeline->store_ids(@_)},
				"store_annotations"=> sub {$pipeline->store_annotations(@_)},
				"strict_denovo"=> sub {$pipeline->strict_denovo(@_)},
				"loh"=> sub {$pipeline->loh(@_)},
				"global_infos"=> sub {$pipeline->global_infos(@_)},
				"coverage" =>sub {$pipeline->coverage(@_)},
				"local_config" => sub {$pipeline->local_config(@_)},
				"dejavu" => sub {$pipeline->dejavu(@_)},
				"check_store_ids"=>  sub {$pipeline->check_store_ids(@_)},
				"check_store_annotations"=> sub {$pipeline->check_store_annotations(@_)},
				"check_strict_denovo"=> sub {$pipeline->check_strict_denovo(@_)},
				"check_loh"=> sub {$pipeline->check_loh(@_)},
				"check_global_infos"=> sub {$pipeline->check_loh(@_)},
				"check_coverage" =>sub {$pipeline->check_coverage(@_)},
				"polydiag" =>sub {$pipeline->polydiag(@_)},
				"quality_check" =>sub {$pipeline->quality_check(@_)},
			};

if ($help) { usage(); }
unless ($steps_name) { usage(); }

$pipeline->yes($yes);
$pipeline->unforce(0) if $force;
$pipeline->add_sample(patient=>$project);



my @running_steps;
if (exists $predef_steps->{$steps_name}) {
	@running_steps = @{$predef_steps->{$steps_name}};
}
else {
	confess($steps_name) unless exists $steps->{$steps_name}  ; 
	@running_steps = ($steps_name);# if exists $steps->{$step}
	
}
foreach my $s (@running_steps){
		unless  (exists $steps->{$s} ){
			warn $s;
			usage();
		}
	}

prepare_calling_jobs();
$pipeline->print_all_steps();
unless ($yes){
	my $choice = prompt("run this/these step(s)   (y/n) ? ");
	die() if ($choice ne "y"); 
}
$pipeline->clean();

if ($pipeline->nocluster ne 2){
$SIG{'INT'} = sub {
	 if (defined $pipeline->daemon){
		$pipeline->daemon->Kill_Daemon($pipeline->pid,15) if $pipeline->daemon->Status($pipeline->pid);
		sleep(2);
	 	system("$Bin//scripts/scripts_pipeline/parse_qstat.pl ".$pipeline->timestamp());
	 	sleep(2);
	}
	elsif (defined $pipeline->process()) {
		warn "process";
	 	$pipeline->process->kill();
	 	sleep(5);
	 	system("$Bin//scripts/scripts_pipeline/parse_qstat.pl ".$pipeline->timestamp());
	 };
	 $pipeline->clean_error;
	 exit(1);
};
}
$pipeline->launch_bds_daemon();
exit(0);



##### METHODS ######



sub prepare_calling_jobs {
	my $next_file = "";
	foreach my $step (@running_steps){
		warn $step;
		($next_file) = $steps->{$step}->((filein=>$next_file));
	}
}

sub usage {
	print colored ['red'],  "\n======================= USAGE ========================\n";
	print colored ['blue'], $0." -project=project_name -steps=(check or polyquery or polydiag or all) -force=(1 : force restart step) -cpu_max = (available CPU number if run without PBS protocol) \n";
	print colored ['green'], "type : polyquery ===> steps : " .join(",",@{$predef_steps->{cache_polyquery}})."\n";
	print colored ['green'], "type : polydiag ===> steps : " .join(",",@{$predef_steps->{cache_polydiag}})."\n";
	print colored ['green'], "type : check_polyquery ===> steps : " .join(",",@{$predef_steps->{check_cache_polyquery}})."\n";
	print colored ['blue'], "================== steps list ================\n";
	print colored ['green'], join(", ",keys %$steps)."\n";
	print colored ['blue'], "===================================================\n";
	if (defined $steps_name && $steps_name ne "all") {
		my @list_steps = split(",",$steps_name);
		foreach my $s (@list_steps){
			unless (grep{/$s/} keys %$steps ){
				print colored ['red'], $s." is not a valid step.\n";
			}
		}
	}
	print colored ['red'],"=================================================\n";
	die();
}