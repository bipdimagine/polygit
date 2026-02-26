#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use Set::IntSpan::Island;
use Set::IntSpan::Fast::XS;
use Array::IntSpan;
use lib $Bin;
 use GenBoBinaryFile;
use GBuffer;
#use Tie::IntegerArray;
use IPC::Open2;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Scalar::Util qw(looks_like_number);
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use colored;
use Bio::DB::Sam;
use Parallel::ForkManager;
use String::ProgressBar;
use Set::IntSpan::Fast::XS;
#use Bio::DB::HTS;
 use JSON::XS;
use List::Util  qw(sum);
use IO::Handle;
use Fcntl 'SEEK_SET'; 
use IO::File ;
use Array::Diff;
#use File::Binary qw($BIG_ENDIAN $LITTLE_ENDIAN $NATIVE_ENDIAN);
use List::MoreUtils qw{ natatime };
use JSON::XS;
#use DB_File ;

my $filein;
my $dir;
my $file_bed;
my $name;
my $fork = 1;
my $project_name;
my $patient_name;
my $verbose;
my $use_samtools;
my $log_file;
my $version;
my $json;
GetOptions(
	"fork=s"   => \$fork,
	"project=s" =>\$project_name,
	"patient=s" =>\$patient_name,
	"verbose=i" =>\$verbose,
	"version=s" =>\$version,
	"json=s" =>\$json,
);

my $buffer = new GBuffer;
my $project = $buffer->newProject ( -name => $project_name , -version =>$version);
my $pm   = new Parallel::ForkManager($fork);
my $tabix = $buffer->software("tabix");
my $bgzip = $buffer->software("bgzip");
my $res;
my $patient = $project->getPatient($patient_name);
my $capture = $patient->getCapture();
my $bed =  $capture->gzFileName();
die() unless -e $bed;
my $bam = $patient->getBamFile();
warn $bed;
my $capture = "";
unless ($project->isGenome){
	$capture = "-b $bed ";
}
open SAM , "samtools depth $capture -a $bam -@ $fork| cut -f 3 |";

	my $s5;
	my $s30;
	my $nb;
	my $s15;
	my $s100;
	my $s1;
	my $s20 =0;
	my $sum = 0;
while (<SAM>){
	chomp();
	my $a = $_;
	$nb ++;
		$s1 ++ if $a >= 1;
		$sum += $a;
		$s5 ++ if $a >= 5;
					$s15 ++ if $a >= 15;
					$s20 ++ if $a >= 20;
					$s30 ++ if $a >= 30;
					$s100 ++ if $a >= 100;
} 
close(SAM);
my $hjson;
my $pid = $patient->id;
	my $coverage_file;
	$coverage_file = $patient->getCoverageFile();
	my $bed_coverage = $patient->getCoverageFile();
	$bed_coverage =~ s/\.gz//;
	open(BED,">$bed_coverage");
	warn $coverage_file;
	warn $name."\n";
	my $z =  $nb;#/$res->{$name}->{nb}));
	 $hjson->{$pid}->{"nb"} = $z; 
	print BED "mean_all\t1\t".$z."\n";
	 $z= ($s5/$nb);
	 $hjson->{$pid}->{"5x"} = int($z*1000)/10;  
	print BED "mean_all\t5\t".$z."\n";
	 $z= ($s15/$nb);
	 $hjson->{$pid}->{"15x"} = int($z*1000)/10; ; 
	print BED "mean_all\t15\t".$z."\n";
	$z= ($s20/$nb);
	 $hjson->{$pid}->{"20x"} = int($z*1000)/10; ; 
	print BED "mean_all\t20\t".$z."\n";
	 $z= ($s30/$nb);
	 $hjson->{$pid}->{"30x"} = int($z*1000)/10; ; 
	print BED "mean_all\t30\t".$z."\n";
	 $z= ($sum/$nb);
	 $hjson->{$pid}->{"mean"} = int($z*10)/10; 
	print BED "mean_all\t99\t".$z."\n";
	 $z= ($s100/$nb);
	 $hjson->{$pid}->{"100x"} = int($z*1000)/10; 
	print BED "mean_all\t100\t".$z."\n";
	close BED;
system("$bgzip -f $bed_coverage; $tabix -b 2 -e 2 -s 1 -f $coverage_file");


#my $hjson; 
#
#foreach my $patient (@{$patients}){
#	my $name = $patient->name;
#	my $pid = $patient->id; 
#	
#	my $coverage_file;
#	$coverage_file = $patient->getCoverageFile();
#	my $bed_coverage = $patient->getCoverageFile();
#	$bed_coverage =~ s/\.gz//;
#	open(BED,">$bed_coverage");
#	warn $coverage_file;
#	#die if -e $coverage_file;
#	#die();
#	warn $name."\n";
#	my $z =  $res->{$name}->{nb};#/$res->{$name}->{nb}));
#	 $hjson->{$pid}->{"nb"} = $z; 
#	print BED "mean_all\t1\t".$z."\n";
#	 $z= (($res->{$name}->{s5}/$res->{$name}->{nb}));
#	 $hjson->{$pid}->{"5x"} = int($z*1000)/10;  
#	print BED "mean_all\t5\t".$z."\n";
#	$z = (($res->{$name}->{s15}/$res->{$name}->{nb}));
#	 $hjson->{$pid}->{"15x"} = int($z*1000)/10; ; 
#	print BED "mean_all\t15\t".$z."\n";
#	$z = (($res->{$name}->{s20}/$res->{$name}->{nb}));
#	 $hjson->{$pid}->{"20x"} = int($z*1000)/10; ; 
#	print BED "mean_all\t20\t".$z."\n";
#	$z =  (($res->{$name}->{s30}/$res->{$name}->{nb}));
#	 $hjson->{$pid}->{"30x"} = int($z*1000)/10; ; 
#	print BED "mean_all\t30\t".$z."\n";
#	$z =  ($res->{$name}->{sum}/$res->{$name}->{nb});
#	 $hjson->{$pid}->{"mean"} = int($z*10)/10; 
#	print BED "mean_all\t99\t".$z."\n";
#	$z =  (($res->{$name}->{s100}/$res->{$name}->{nb}));
#	 $hjson->{$pid}->{"100x"} = int($z*1000)/10; 
#	print BED "mean_all\t100\t".$z."\n";
#	close BED;
#	system("$bgzip -f $bed_coverage; $tabix -b 2 -e 2 -s 1 -f $coverage_file");
#}
#
#if ($json) {
#	print encode_json  $hjson;
#}
exit(0);


