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
 
use bds_calling_steps;    
use file_util;
use Class::Inspector;
use check_utils;
use Text::Table;
use colored; 
 use Term::Menus;
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
my $no_cluster = 0;
my $stdout;
my $fork = 0;
my $yes = 0;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'steps=s' => \$steps_name,
	'force=s' => \$force,
	'type=s' => \$type,
	'yes=s' => \$yes,
	"method=s" => \$method,
	"nocluster=s" => \$no_cluster,
	"fork=s" => \$fork,
);


#unless (-e "/usr/bin/qstat"){
#	$no_cluster = 1;
#}


unless ($projectName) { usage(); }
usage() unless $patients_name;



my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName);
my $jobs;
my $dir = $project->getPipelineDir();
my $file = $dir."/".$project->name().".log";


my @running_steps;
@running_steps =("genotype_gvcf4","correct_vcf","move_vcf","dude");# if $method eq "haplotypecaller";
 

my $pipeline = bds_calling_steps->new( project=>$project,nocluster=>$no_cluster);

#my $pipeline = bds_steps->new( project=>$project,argument_patient=>$patients_name,nocluster=>$nocluster);

if ($patients_name eq "novcf"){
	my $buffer1 = GBuffer->new();
	my $project1 = $buffer->newProject( -name => $projectName);
	my @p;
	foreach my $p (@{$project->getPatients()}){
		#warn 
		#next if  $p->getVariationsFile ne "";
		next if -e  $p->getVariationsFile;
		push(@p,$p->name);
	
	}
	$patients_name = join(",",@p);
	warn $patients_name;
	warn scalar(@p);
}
$pipeline->argument_patient($patients_name);

if ($fork > 0){
	$pipeline->nproc();
	$pipeline->nproc($fork);
}



my $steps = {
	"correct_vcf"					=> sub{$pipeline->correct_vcf(@_)},
	"correct_vcf4"					=> sub{$pipeline->correct_vcf4(@_)},
	"move_vcf"						=> sub{$pipeline->move_vcf(@_)},
	"move_vcf4"						=> sub{$pipeline->move_vcf4(@_)},
	"quality_check"					=> sub{$pipeline->quality_check(@_)},
	"genotype_gvcf4"						=>sub {$pipeline->genotype_gvcf4(@_)},
	"genotype_gvcf"						=>sub {$pipeline->genotype_gvcf(@_)},
	"lifinder"						=>sub {$pipeline->calling_larges_indels(@_)},
	"featureCounts"						=>sub {$pipeline->count_featureCounts(@_)},
	"xhmm"						=>  sub {$pipeline->xhmm(@_)},
	"cache"						=>  sub {$pipeline->cache(@_)},
	"dude"						=>  sub {$pipeline->dude(@_)},
	"cache_test"						=>  sub {$pipeline->cache_test(@_)},
	"mutect2"						=>  sub {$pipeline->mutect2(@_)},
	"lofreq"						=>  sub {$pipeline->lofreq(@_)},
	
};

if (defined $steps_name && $steps_name ne "all") {
	my @list_steps = split(",",$steps_name);
	foreach my $s (@list_steps){
		warn $s;
		unless (grep{/$s/} keys %$steps ){
			usage();
		}
	}
}

 
 my $steps_defined ={
 		"all" =>"genotype_gvcf4,correct_vcf4,move_vcf4,dude",
 		"just_calling4" =>"genotype_gvcf4,correct_vcf4,move_vcf4",
 		"all_V3" =>"genotype_gvcf,correct_vcf,move_vcf,dude",
 		"test" =>"genotype_gvcf4,correct_vcf,move_vcf",
 		"without_move" =>"genotype_gvcf4,correct_vcf4,dude",
 		"move" =>"move_vcf4",
 		"featureCounts" => "featureCounts",
		"exome_umi" => "mutect2,move_vcf",
 		"dude" => "dude",

 } ;
  
unless ($steps_name){
my @list=('all','just_calling4','all_V3','test','without_move','move','featureCounts',"dude","exome_umi");
   my $banner="  Please Pick a pipeline type for $projectName :";
   	print colored ['black ON_BRIGHT_GREEN'];
   my $pick=&pick(\@list,$banner);
   $steps_name = $steps_defined->{$pick};
   print "SELECTION = $steps_name\n";
 	
   }
     
 
if ($steps_name ne "all") {
	@running_steps = split(",", $steps_name);
}

my $haplo = "haplotypecaller";
my $test_haplo = grep {$_ eq "genotype_gvcf4"} @running_steps;
$haplo = "haplotypecaller4" if $test_haplo;
$pipeline->method_calling($haplo);

$pipeline->yes($yes);
$pipeline->unforce(0) if $force;
$pipeline->add_sample(patient=>$project);

 prepare_calling_jobs();
 $pipeline->print_all_steps();
 unless ($yes){
	my $choice = prompt("run this/these step(s)   (y/n) ? ");
	die() if ($choice ne "y"); 
 }
 
$pipeline->clean();



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
$pipeline->launch_bds_daemon();

exit(0);



sub prepare_calling_jobs {
	my $next_file="";
	foreach my $step (@running_steps){
		($next_file) = $steps->{$step}->((filein=>$next_file));
	}
	#die();
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

