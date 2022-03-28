#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo/lib/obj-nodb/";
use lib "$Bin/packages";
use GBuffer;
use GenBoProject;
use PBS::Client;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
use MooseX::Method::Signatures;
 
use calling_steps;    
use file_util;
use Class::Inspector;
use check_utils;
use Text::Table;
use colored; 

my $projectName;
my $filename;

my $name;
my $patients_name;
my $steps_name;
my $force;
my $type;
my $fastq_ext;
my $somatic;
my $method;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'steps=s' => \$steps_name,
	'force=s' => \$force,
	'type=s' => \$type,

	"somatic" => \$somatic,
	"method=s" => \$method,
);



die("and the method ?") if $method ne "unifiedgenotyper" && $method ne "haplotypecaller_gvcf" && $method ne "haplotypecaller_vcf" && $method ne "haplotypecaller";;
unless ($projectName) { usage(); }
usage() unless $patients_name;



my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName);
my $jobs;
my $dir = $project->getPipelineDir();
my $file = $dir."/".$project->name().".log";

my ($m1,$mtype) = split("_",$method);

my $pipeline = calling_steps->new(log_file=>$file, project=>$project,method_calling=>$m1 );
$pipeline->argument_patient($patients_name);
$pipeline->project($project);

$pipeline->unforce(0) if $force;
$pipeline->somatic(1) if $somatic;

my $steps = {
	"unified_genotyper"				=> sub{$pipeline->calling_unifiedgenotyper(@_)},
	"snp_recalibrator"				=> sub{$pipeline->snp_recalibrator(@_)},
	"snp_apply_recalibration"		=> sub{$pipeline->snp_apply_recalibration(@_)},
	"indels_recalibrator"			=> sub{$pipeline->indels_recalibrator(@_)},
	"indels_apply_recalibration"	=> sub{$pipeline->indels_apply_recalibration(@_)},
	"mendelian_verification"		=> sub{$pipeline->mendelian_verification(@_)},
	"sex_verification"				=> sub{$pipeline->sex_verification(@_)},
	"coverage_verification"			=> sub{$pipeline->coverage_verification(@_)},
	"split_vcf"						=> sub{$pipeline->split_vcf(@_)},
	"correct_vcf"					=> sub{$pipeline->correct_vcf(@_)},
	"check_dbsnp"					=> sub{$pipeline->check_dbsnp(@_)},
	"move_vcf"						=> sub{$pipeline->move_vcf(@_)},
	"plink"							=> sub{$pipeline->plink(@_)},
	"quality_check"					=> sub{$pipeline->quality_check(@_)},
	"haplotypecaller"						=>sub {$pipeline->calling_haplotypecaller_gvcf(@_)},
	#"genotype_gvcf"						=>sub {$pipeline->genotype_gvcf(@_)},
	#"join_gvcf"						=>sub {$pipeline->join_gvcf(@_)},
	"move_snp"						=>sub {$pipeline->move_snp_vcf(@_)},
			    
};

unless ($steps_name) {
	usage();
}

if (defined $steps_name && $steps_name ne "all") {
	my @list_steps = split(",",$steps_name);
	foreach my $s (@list_steps){
		unless (grep{/$s/} keys %$steps ){
			usage();
		}
	}
}


 
#my @running_steps =("unified_geontyper","snp_recalibrator","snp_apply_recalibration",
#"indels_recalibrator","indels_apply_recalibration", "split_vcf",
#"check_dbsnp","mendelian_verification","sex_verification");
my @running_steps;
 @running_steps =("unified_genotyper", "correct_vcf", "quality_check" ) if $method eq "unifiedgenotyper"; #plink
 if ($m1 eq  "haplotypecaller_gvcf"){
 @running_steps =("haplotypecaller" ) if $mtype eq "gvcf"; #plink
 @running_steps =("join_gvcf","genotype_gvcf","correct_vcf", "quality_check" ) if $mtype eq "vcf";
 }
 
if ($steps_name ne "all") {
#	if ($steps_name =~ /move_vcf/) {
#		my $jsonFile = $project->getProjectPath().'../'.$projectName.'.resume.json';
#		my ($nbWarn, $nbErr) = check_utils::getTableFromJsonFile($jsonFile);
#		checkErrors($nbWarn, $nbErr);
#	}
	@running_steps = split(",", $steps_name);
}

my $job = prepare_calling_jobs();



$pipeline->print_all_steps;
my $choice = prompt("run this step (y/n)? ");
die() if ($choice ne "y"); 

my $client = PBS::Client->new();
$client->qsub($job);

print "\n";
print colored ['black ON_BRIGHT_GREEN'],"the log file is here : $file";
print color 'reset';
print "\n";
exit(0);


##### METHODS ######


sub checkErrors {
	my ($nbWarn, $nbErr) = @_;
	if (($nbWarn > 0) or ($nbErr > 0)) {
		print "\n";
		if ($nbWarn > 0) {
			print  colored ['black ON_BRIGHT_YELLOW'], " -> $nbWarn WARNING found !";
			print  color 'reset';
			print "\n";
		}
		if ($nbErr > 0) {
			print  colored ['black ON_BRIGHT_RED'], " -> $nbErr WARNING found !";
			print  color 'reset';
			print "\n";
		}
		print  color 'reset';
		print "\n";
		my $choice = prompt("Do you want to run this analysis and avoid the warning(s) ? (y/n)  ");
		die("\nExit...\n\n") if ($choice ne "y");
		print  "\n\n";
	}
}

sub prepare_calling_jobs {
	my $current_job = $pipeline->start_job();
	my $next_file="";

	foreach my $step (@running_steps){
		warn $step;
		($current_job,$next_file) = $steps->{$step}->((filein=>$next_file,previous=>$current_job));
	}
	warn $current_job;
	#die();
	my $final_job = $pipeline->end_job();
	$final_job->prev({ ok => $current_job} );
	return $final_job;
}

sub usage {
print colored ['red'],  "\n======================= USAGE ========================\n";
print colored ['blue'], $0." -project=project_name -steps=(step or all) -patients=(patient or all)  -force=(1 restart step) -method=haplotypecaller opr unifiedgenotyper\n\n";
print colored ['blue'],  "============= steps lists ================\n";
print colored ['green'],join(", ",keys %$steps)."\n";
print colored ['blue'],"=========================================\n";
if (defined $steps_name && $steps_name ne "all") {
	my @list_steps = split(",",$steps_name);
	foreach my $s (@list_steps){
		unless (grep{/$s/} keys %$steps ){
			print colored ['red'], $s." is not a valid step.\n";
		}
	}
}
print colored ['blue'],"=========================================\n";
	die();
}

