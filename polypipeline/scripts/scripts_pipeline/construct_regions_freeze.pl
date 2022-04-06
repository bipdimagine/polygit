#!/usr/bin/perl
use FindBin qw($Bin $RealBin);
use strict;

use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/"; 
use GBuffer ;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Bio::DB::Sam;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use List::Util qw(sum);
use File::Temp;
 use Time::Elapsed qw( elapsed );
 use Time::ETA;
 
my $filein;
my $dir;
my $file_bed;
 
my $dp_limit = 5;
my $al = 3;
my $project_name;
my $fork;
my $patient_name;
$| =1;
my $log_file;
my $vcf_final;
my $region;
my $start;
my $end;
my $chr;
my $region_name;
my $fork =1;
my $window_length;
my $file;
my $version;
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"log=s" =>\$log_file,
	"patients=s" => \$patient_name,
	"file=s" => \$file,
	#"chr=s"=>\$chr,
	#"start=s"=>\$start,
	#"end=s"=>\$end,
	#"region=s"=>\$region_name,
	#"fork=s" =>\$fork,
	"window=s" => \$window_length,
	"version=s" => \$version,
);
my $date = `date`;
chomp($date);

if ($log_file){
	open (STDOUT,">>".$log_file);
}

if ($region_name){
	my $reste;
	 ($chr,$reste) = split(":",$region_name);
	($start,$end) =  split("-",$reste);
}

my $other_project = [];
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $java = $project->getSoftware('java');
my $samtools = $project->getSoftware('samtools');
my $javac = $project->getSoftware('java');
$javac = "java" unless -e $javac;
my $gatk  = $project->getSoftware('gatk');
my $project = $buffer->newProject( -name => $project_name, -version=>$version );

my $chrs = $project->getChromosomes();
my $reference = $project->getGenomeFasta;
my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");
	my $ref =  $project->genomeFasta();
	 $project->get_only_list_patients($patient_name);

my $patient = $project->getPatientOrControl($patient_name);




#my $fork = 4;
my $dir_out_gvcf= $project->getCallingPipelineDir("gvcf");
#my $dir_out_gvcf = $dir_out."/gvcf";
$dir_out_gvcf .= "/$patient_name/";
system ("mkdir -p $dir_out_gvcf") unless -e $dir_out_gvcf;
unlink $file if -e $file;
my $region = "$chr:$start-$end";
my $output_region = "$chr.$start.$end";
my $bam =  $patient->getBamFile();

	my $recal = $patient->getRecalFile();
	my $recal_string="";
	$recal_string = "--BQSR $recal" if -e $recal;
	my $nct ="-nct 2";
	#$nct="";
#	$nct ="-nct $fork" if $fork >1;
	
	my @cmds;
	my $real_fork = $fork;# int($fork /2);
my $pm = new Parallel::ForkManager($real_fork);

my $all_windows;
my $total = 0;
foreach my $chr (@{$project->getChromosomes()}){
	#next if $chr->name eq "MT";
		 my $windows = $chr->getWindowCaptureForCalling(250,$window_length);
		 warn $chr->name();
			#next if $chr->name eq "MT";
		#	next unless $chr->name eq "18";
		 
		#  my $intspan0 = $chr->getIntSpanCaptureForCalling(0);
		 $all_windows->{$chr->name} = [];
		  foreach my $window (@$windows){
		  	#		warn Dumper $window;
		  		my $outfile = $patient->getWindowGvcf($window);
		  		my $debug;
		  		$debug =1 if $outfile eq "/data-beegfs/sequencing/pipeline/NGS2016_1283/HG19/calling/gvcf//B00GHV2/B00GHV2.chr1.6000001.7000001.g.vcf";
		  	#	next if -e $outfile;
		  		$window->{outfile} = $outfile;
		  			$window->{chr_ucsc} = $chr->fasta_name;
		  
		  			$total ++;
		  			
		  		push(@{$all_windows->{$chr->name}},$window);
		  		
		  		
		  }
}
 $file = $dir_out_gvcf."/$patient_name.regions.$window_length.freeze" unless $file;
 store $all_windows,$file;



