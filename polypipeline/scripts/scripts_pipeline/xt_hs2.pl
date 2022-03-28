#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use File::Basename;
use  File::Temp;
use Sys::Hostname;

my $filein;
my $dir;
my $name;
my $project_name;
my $prog;
my $file1;
my $file2;
my $lane ="";
my $mode ;
my $verbose;
my $fileout;
my $fork = 1;
use List::MoreUtils qw(firstidx);

my $dir_out;
my $log_file;
GetOptions(
	"log=s" =>\$log_file,
	"project=s" =>\$project_name,
	"fork=s" =>\$fork,
	"fileout=s"  =>\$fileout,
	"file1=s"  =>\$file1,
	"file2=s"  =>\$file2,
	"method=s" => \$prog,
	"name=s" => \$name,
	"lane=s" => \$lane,
	"mode=s" => \$mode,
	"verbose=s" => \$mode,
);
 
if ($log_file){
	open (STDOUT,">>".$log_file);
}
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatient($name);
my ($dirin) = $patient->getSequencesDirectory();
my $agent_trimmer = "/software/distrib/AGeNT/AGeNT-latest/lib/trimmer-2.0.3.jar";
my $output_dir = $project->getAlignmentPipelineDir($prog);
system("mkdir -p $output_dir;chmod a+rwx $output_dir") unless -e $output_dir;
my $java = $buffer->software('java');
$java ="java" unless -e $java;
	
my $trimmer_out=$output_dir."/done.txt";
my $cmd = qq{$java -jar $agent_trimmer -fq1 $file1 -fq2 $file2  -v2  -out_loc $output_dir};

warn $cmd;
system( $cmd);


my @res_files = glob("$dirin/".$name."*.fastq*_*.gz");
#die (join(",'",@res_files)) if scalar(@res_files) != 3;
$trimmer_out = $res_files[0];
my $newFile;
my %fileout;
$fileout{R1} =$output_dir."/".$name."_R1.fastq.gz";
$fileout{R2} =$output_dir."/".$name."_R2.fastq.gz";
$fileout{MBC} =$output_dir."/".$name."_RN.txt.gz";
$fileout{STATS} =$output_dir."/".$name."_STATS.txt";

	foreach my $res (@res_files){
		warn $res;
		if( $res =~ m/.*_STATS*.properties/) {
			warn $res;
			$newFile = $fileout{properties};
		}
		if( $res =~ m/.*_MBC_0.txt.gz/) {
			$newFile = $fileout{MBC};
		}
		elsif($res =~ m/_R1_.*_Cut_0.fastq.gz/){
			$newFile = $fileout{R1};	
		}
		elsif($res =~ m/_R2_.*_Cut_0.fastq.gz/){
			$newFile = $fileout{R2};	
		}
		else{
			die("problem with file $res");
		}
		my $cmd2 = " test -e $res  || exit 1 &&  mv $res $newFile  || exit 1";#
		system $cmd2;
	}
	
#
#method mv_agent_files (Str: $filein!) {
#	my $method = $self->patient()->alignmentMethod();
#	my ($dirin) = $self->patient()->getSequencesDirectory();
#	warn $dirin;
#	my $dirout = $self->project->getAlignmentPipelineDir($method);
#	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
#	my $name = $self->patient()->name();
#	my $trimmer_out=$dirout."/".$name."_RN.txt.gz";
#	return $trimmer_out if -e $trimmer_out;
#	my $fileout;
#	my @res_files = glob("$dirin/".$name."*fastq*_*");
#	warn Dumper(@res_files);
#	my $ppn = $self->nproc;# if $self->nocluster;
#	my $project_name =  $self->project->name();
#	#die (join(",'",@res_files)) if scalar(@res_files) != 3;
#	my %fileout;
#	$fileout{R1} =$dirout."/".$name."_R1.fastq.gz";
#	$fileout{R2} =$dirout."/".$name."_R2.fastq.gz";
#	$fileout{MBC} =$dirout."/".$name."_RN.txt.gz";
#	my $newFile;
#	foreach my $res (@res_files){
#		if( $res =~ m/.*_MBC_0.txt.gz/) {
#			$newFile = $fileout{MBC};
#			$fileout = $fileout{MBC};
#		}
#		elsif($res =~ m/_R1_.*_Cut_0.fastq.gz/){
#			$newFile = $fileout{R1};	
#		}
#		elsif($res =~ m/_R2_.*_Cut_0.fastq.gz/){
#			$newFile = $fileout{R2};	
#		}
#		else{
#			die("problem with file $res");
#		}
#		my $cmd = qq{ test -e $res  || exit 1 &&  mv $res $newFile || exit 1"; };
#		my $type = "mv_agent_file";
#		my $stepname = $name."@".$type;
#	
#		my $job_bds = job_bds_tracking->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>$res,fileout=>$fileout{MBC} ,type=>$type,dir_bds=>$self->dir_bds,sample_name=>$name,project_name=>$project_name,software=>$method);
#		$self->current_sample->add_job(job=>$job_bds);
#		if ($self->unforce() && -e $trimmer_out){
#	 		$job_bds->skip();
#		}
#		
#	}
#	return ($fileout);
#}