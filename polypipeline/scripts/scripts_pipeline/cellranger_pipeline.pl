#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use colored;
use file_util;
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
use Proc::Simple;


 my $project_name;
 my $fork;
 my $list_patients;
 my $bcl_dir;
 my $sampleSheet;
 my $step;
 my @steps;


GetOptions(
	'project=s'   => \$project_name,
#	"fork=s"  => \$fork,
	"patients=s"=>\$list_patients,
	"bcl_dir=s" => \$bcl_dir,
	"sample_sheet=s" => \$sampleSheet,
	"step=s" => \$step,
);

my $verbose;
$fork =1 unless $fork;
my $debug;
#$debug = "-l off"; 
my $date = `date`;


my $buffer = GBuffer->new();

my $project = $buffer->newProject( -name => $project_name );
my $cell  = $project->getSoftware("cellranger");
my $ref_root = $buffer->config->{public_data}->{$project->getVersion()};
my $index = $ref_root."/".$buffer->index("cellranger");
my @seqDir=@{$project->getSequencesDirectories};
confess ("problem sequence dir") if scalar(@seqDir)>1;
warn  $seqDir[0];
my $patients = $project->getPatients();
#my $patients = file_util::return_patients( $project, $patients_name );
#my $patients =  $project->get_list_patients($list_patients);
#confess("no patients") unless scalar(@$patients);
#$step = "mkfastq,count" if $step eq "all";
#@steps= join(",",$step); 


#foreach my $s (split(",",@steps)){
#	warn $s;
#	`$s`;
#}
mkfastq($project_name,$bcl_dir,$seqDir[0],$sampleSheet);
count($patients,$seqDir[0],$index);


sub mkfastq {
	my ($project_name,$bcl_dir,$seqDir,$sampleSheet) = @_;
	my $id = $project_name."_FCB";
	my $align_command = " $cell  mkfastq --run=$bcl_dir  --id=$id  --output-dir=$seqDir --delete-undetermined --ignore-dual-index --csv=$sampleSheet  --jobmode=slurm";
	warn $align_command;
	system( $align_command);
}




sub count{
	my ($patients,$seqDir,$index) = @_;
	foreach my $p (@$patients){
		my $patientName = $p->name; 
		my $myproc = Proc::Simple->new();
		my $count_command = "$cell count --id=$patientName --transcriptome=$index --fastqs=$seqDir --sample=$patientName --jobmode=slurm";
		print $count_command."\n";
		$myproc->start($count_command); 
	}

}


# my $methods = {
#			"unifiedgenotyper"				=> { "method" => sub{calling_target::unifiedGenotyper(@_)}, priority =>1}, 
#			"freebayes"				=> { "method" => sub{calling_target::freebayes(@_)}, priority =>3}, 
#			"samtools"		=> { "method" => sub{calling_target::samtools(@_)},priority=>2},
#			"platypus"		=> { "method" => sub{calling_target::platypus(@_)},priority=>2},
#			"duplicate_region_calling"		=> { "method" => sub{calling_target::duplicate_region_calling(@_)},priority=>2},
#};
#die($method) unless exists  $methods->{$method};
#my $vcf =  $methods->{$method}->{method}($project,$patient,$intspan,$low_calling,$fork);
#warn $vcf;
#die("problem ") unless -e $vcf;
#my $bgzip = $buffer->software("bgzip");
#system("$bgzip $vcf");
#
#warn "ok";
#$vcf = $vcf.".gz";
#warn $vcf;
#system("rsync -v $vcf $fileout && rm $vcf");
#
#my $tabix =  $buffer->software("tabix");
#
#system("$tabix -p vcf $fileout");
#
#die() unless -e $fileout.".tbi";
#
#exit(0);
