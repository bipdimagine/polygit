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
 my $bam;
 my $patient_name;
 



GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"patient=s" => \$patient_name,
	"fileout=s"=>\$fileout,
	"filein=s"=>\$bam,
);

my $verbose;
$fork =1 unless $fork;
my $debug;
#$debug = "-l off"; 
my $date = `date`;




my $buffer = GBuffer->new();

my $project = $buffer->newProject( -name => $project_name );
my $gatk  = $project->getSoftware('gatk');
my $java = $project->getSoftware('java');

	my $ref =  $project->genomeFasta();

my $patient = $project->getPatientOrControl($patient_name);

 $bam = $patient->getBamFile() unless $bam;


my $dir_out=$project->getCallingPipelineDir("callable");





#
 my $pm1 = new Parallel::ForkManager($fork);
 
 my $is_callable;
 
 my %all_intspan;
$pm1->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    	 # $pr->update($c++);
	
 		#$pr->write() if $verbose;
 		my $chr = $data->{chr};
 		my $ucsc = $data->{ucsc};
 			$all_intspan{$ucsc} = Set::IntSpan::Fast->new() unless exists 	$all_intspan{$ucsc};
 			foreach my $region_callable (@{$data->{array}}){
 				my $pos;
 				$pos->{chr} = $chr;
 		 	($pos->{start},$pos->{end}) = split("-",$region_callable);
 		 	$all_intspan{$ucsc}->add_range($pos->{start},$pos->{end});
 		 	$is_callable = 1;
 			}
    		
    }
  );
  

 foreach my $chr (@{$project->getChromosomes()}){
 	#next if $chr->name ne "1";
 my $set = Set::IntSpan::Fast->new();
 my $intspan = $chr->getIntSpanCaptureForCalling(250);
 my @beds = $buffer->intspanToBed($chr,$intspan);
 next unless scalar @beds;

 my $bed1 =  calling_target::getTmpFile($dir_out,$chr->name,"bed");
 open(BED,">$bed1");
 print BED join("\n",@beds);
 close BED;
 my $summary =  calling_target::getTmpFile($dir_out,$chr->name,".sum.txt");
 my $bed =  calling_target::getTmpFile($dir_out,$chr->name,"bed");
 my $pid = $pm1->start and next;
 my $cmd = qq{$java -jar $gatk -T CallableLoci -R $ref -I $bam -o $bed -l off -summary $summary  --minDepth 10 --minDepthForLowMAPQ 15 -frlmq 0.6  -L $bed1};
 system($cmd);
 
  my $callable_intspan = Set::IntSpan::Fast::XS->new();
 	open (GATK, " cat $bed  | ");
 	 	while (my $line  = <GATK>){
 		chomp($line);
 		if  ($line =~/CALLABLE|POOR_MAPPING_QUALITY/){
 			my ($chr,$start,$end,$type) = split(" ",$line);
 			$callable_intspan->add_range($start, $end);
 		}
 	}
 	close GATK; 
 	my $intspan2 = Set::IntSpan::Fast::XS->new();
 	#warn $callable_intspan->as_string();
 	foreach my $capture (@{$project->getCaptures}){
		$intspan2 = $intspan2->union($capture->getIntSpanForChromosome($chr,50));
	}
 	#my $intspan2 = $chr->getIntSpanCaptureForCalling(0);
 	#warn $intspan2->as_string();
 	#die();
 	my @regions_callable = split(",",$callable_intspan->as_string());
 	
 	
 	my $diff = $intspan2->diff($callable_intspan);
 	my $p = int(scalar($diff->as_array)/scalar($intspan2->as_array)*100);

 	my	$color = 'black ON_BRIGHT_GREEN';
 	my $text = "seems to be correct ";
 	
 	if ($p>10){
 		$color = 'black on_magenta';
 		$text= "COVERAGE WARNING PROBLEM :!!!! ";
 	}
 	if (scalar($intspan2->as_array) == 0 ){
 		$color = 'black ON_BRIGHT_RED';
 		$text= "COVERAGE PROBLEM : NO COVERAGE  => chromosome ".$chr->name;
 	}
	
	print  colored [$color],$patient->name."::".$chr->name." coverage ".$p."% ".$text ;
	
	print  color 'reset';
	print  "\n";
 	my $data;
 	$data->{chr} =  $chr->name;
 	$data->{ucsc} =  $chr->ucsc_name;
 	$data->{array} = \@regions_callable;
 	$pm1->finish(0,$data);
 	
 unlink $bed;
 unlink $summary;
 unlink $bed1;
 
 }
 $pm1->wait_all_children();
 unless ($is_callable){
	my	$color = 'black ON_BRIGHT_RED';
 		my $text= "NO CALLABLE REGION FOR SAMPLE  ".$patient->name;
 		print  $text ;
 		die();
}

 store  \%all_intspan, $fileout;
 my $dir_prod= $project->getCoverageCallable();
 my $fileprod = $dir_prod."/".$patient->name.".callable.bed";
 warn $fileprod;
 open(REG,">$fileprod");
 
 foreach  my $chr (@{$project->getChromosomes()}){
 	my $intspan = $all_intspan{$chr->ucsc_name};
 	next unless $intspan;
 	 my $iter2 = $intspan->iterate_runs();
       while (my ( $from, $to ) = $iter2->()) {      
       	print REG $chr->ucsc_name."\t".$from."\t".$to."\n";
       }
 	#  push(@$regions,map{$chr->ucsc_name.":".$_}split(",",$intspan->as_string));
 }
 close REG;
 exit(0);
 
 
 