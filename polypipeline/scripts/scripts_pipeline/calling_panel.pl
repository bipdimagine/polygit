#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use colored;
#use Set::IntSpan;
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use Carp;
use calling_target;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use String::ProgressBar;
use List::Util qw(sum);



 my $project_name;
 my $fork;
 my $fileout;
 my $callable_intspan_file;
 my $patient_name;
 #my $low_calling;
 my $method;


GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"patient=s" => \$patient_name,
	"fileout=s"=>\$fileout,
	"filein=s"=>\$callable_intspan_file,
#	"low_calling=s"=>\$low_calling,
	"method=s"=>\$method,
);

my $verbose;
$fork =1 unless $fork;
my $debug;
#$debug = "-l off"; 
my $date = `date`;


#exit(1) if ($method ne "samtools");
#exit(1) if ($method = "freebayes");
my $buffer = GBuffer->new();

my $project = $buffer->newProject( -name => $project_name );
my $gatk  = $project->getSoftware('gatk');
my $java = $project->getSoftware('java');

	my $ref =  $project->genomeFasta();

my $patient = $project->getPatientOrControl($patient_name);

 my $bam = $patient->getBamFile() ;
my $intspan;

#if  (-e $callable_intspan_file){
#	 $intspan = retrieve($callable_intspan_file);
#	 
#	 #warn  $intspan->{chr5}->as_string;
#	 #die();
#}
#else {
	#die();
	foreach my $chr (@{$project->getChromosomes}) {
	 $intspan->{$chr->ucsc_name} = $chr->getIntSpanCaptureForCalling(300);
}
#	die("error callable region - file .freeze is absent");

#}

# die() unless -e $callable_intspan_file;
 
 #my $intspan = retrieve($callable_intspan_file);

 my $methods = {
			"unifiedgenotyper"				=> { "method" => sub{calling_target::unifiedGenotyper(@_)}, priority =>1}, 
		#	"cp"				=> { "method" => sub{calling_target::cp(@_)}, priority =>1}, 
			"freebayes"				=> { "method" => sub{calling_target::freebayes(@_)}, priority =>3}, 
			"samtools"		=> { "method" => sub{calling_target::samtools(@_)},priority=>2},
			"platypus"		=> { "method" => sub{calling_target::platypus(@_)},priority=>2},
			"duplicate_region_calling"		=> { "method" => sub{calling_target::duplicate_region_calling(@_)},priority=>2},
			"hpduplicate_region_calling"		=> { "method" => sub{calling_target::duplicate_region_calling2(@_)},priority=>2},
			"haplotypecaller"		=> { "method" => sub{calling_target::haplotypecaller(@_)},priority=>1},
			"haplotypecaller4"		=> { "method" => sub{calling_target::haplotypecaller4(@_)},priority=>1},
			"p1_freebayes" => { "method" => sub{calling_target::p1_freebayes(@_)}, priority =>3}, 
			"eif6_freebayes" => { "method" => sub{calling_target::eif6_freebayes(@_)}, priority =>3}, 
			"mutect2" => { "method" => sub{calling_target::mutect2(@_)}, priority =>3}, 
			"lofreq" => { "method" => sub{calling_target::lofreq(@_)}, priority =>3}, 
};
die($method) unless exists  $methods->{$method};
warn $method;
my $vcf =  $methods->{$method}->{method}($project,$patient,$intspan,$fork);
warn $vcf;
die("problem ") unless -e $vcf;
my $bgzip = $buffer->software("bgzip");
system("$bgzip $vcf");

warn "ok";
$vcf = $vcf.".gz";
warn $vcf;
system("rsync -v $vcf $fileout && rm $vcf");

my $tabix =  $buffer->software("tabix");

system("$tabix -p vcf $fileout");

die() unless -e $fileout.".tbi";

exit(0);