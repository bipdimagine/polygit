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
my $tso_adaptors;
my $illumina_adaptors;
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
	"tso=s" =>\$tso_adaptors,
	"illumina=s" =>\$illumina_adaptors,
	"fork=s" =>\$fork,
#	"fileout=s"  =>\$fileout,
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
my $flexbar = $buffer->getSoftware('flexbar');
my $output_dir = $project->getAlignmentPipelineDir($prog);
my $output_log =  $project->getAlignmentPipelineDir($prog)."/".$name.".flexbarOut.log";
#my $output_log= $output_dir."/".$name.".flexbarOut.log";
my $target = $output_dir."/".$name."_F".$lane;
system("mkdir -p $output_dir;chmod a+rwx $output_dir") unless -e $output_dir;
#my $java = $buffer->software('java');
#$java ="java" unless -e $java;

my $cmd = qq{$flexbar --reads $file1  --reads2 $file2 --stdout-reads --adapters $tso_adaptors --adapter-trim-end LEFT --adapter-revcomp ON --adapter-revcomp-end RIGHT  --htrim-left GT --htrim-right CA  --htrim-min-length 3 --htrim-max-length 5 --htrim-max-first --htrim-adapter  --min-read-length 2  --threads $fork |  $flexbar --reads -  --interleaved --target $target --adapters $illumina_adaptors --zip-output GZ --output-log $output_log  --adapter-trim-end RIGHT --min-read-length 2 --threads $fork};
warn $cmd;
system($cmd);

exit(0);


#flexbar -r reads.fastq --pre-trim-left 5
#
#
#my @res_files1 = glob("$output_dir/".$name."*_1.fastq.gz");
#my @res_files2 = glob("$output_dir/".$name."*_2.fastq.gz");
#my $all_files1 = join(" ",@res_files1);
#my $all_files2 = join(" ",@res_files2) ;
#my $catfile1 = $output_dir."/".$name."_R1.fastq.gz";
#my $catfile2 = $output_dir."/".$name."_R2.fastq.gz";
#
#my $cmd2 = qq{cat $all_files1 > $catfile1 };
#my $cmd3 = qq{cat $all_files2 > $catfile2 };
#system ($cmd2);
#system ($cmd3);

#my $newFile;
#my %fileout;
#$fileout{R1} =$output_dir."/".$name."_R1.fastq.gz";
#$fileout{R2} =$output_dir."/".$name."_R2.fastq.gz";
#
#
#	foreach my $res (@res_files){
#		if($res =~ m/$name."_1.fastq.gz"/){
#			$cmd = qq{cat $res > $catfile1};
#			system ($cmd);
#			$newFile = $fileout{R1};	
#		}
#		elsif($res =~ m/$name."_2.fastq.gz"/){
#			$cmd = qq{cat $res > $catfile1};
#			$newFile = $fileout{R2};	
#		}
#		else{
#			die("problem with file $res");
#		}
#		my $cmd2 = " test -e $res  || exit 1 &&  mv $res $newFile  || exit 1";#
#		system $cmd2;
#	}
#	