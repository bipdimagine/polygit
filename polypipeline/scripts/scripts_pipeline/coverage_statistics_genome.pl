#!/usr/bin/env perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
 use GenBoBinaryFile;
use GBuffer;
#use Tie::IntegerArray;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Scalar::Util qw(looks_like_number);
use Storable qw(store retrieve freeze);
use colored;
use Parallel::ForkManager;
use List::Util  qw(sum);
use JSON::XS;

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

my $patients;
if ($patient_name) { push(@$patients, $project->getPatient($patient_name)); }
else { $patients = $project->getPatients(); }

if ($project->isGenome && $json) {
	my $hjson;
	foreach my $patient (@{$patients}){
		my $pid = $patient->id;
			my $all_sum;
		my $coverage_file;
		$coverage_file = $patient->getCoverageFile();
		warn $coverage_file;
		my $h;
		confess($coverage_file) unless -e $coverage_file;
		my $tabix = $patient->tabix_coverage;
			my $res   =  $tabix->query_full( "mean_all") ;
			while ( my $line = $res->next ) {
				my ( $a, $b, $c ) = split( " ", $line );
				
				$h->{"s".$b} = $c;
			}
			$h->{"nb"} = delete $h->{s1};
			$h->{"mean"} = delete $h->{s99};
	 my $z= (($h->{s5}/$h->{nb}));
	 $hjson->{$pid}->{"5x"} = int($h->{s5}*1000)/10;  
	$z = (($h->{s15}/$h->{nb}));
	 $hjson->{$pid}->{"15x"} = int($h->{s15}*1000)/10; ; 
	$z = (($h->{s20}/$h->{nb}));
	 $hjson->{$pid}->{"20x"} = int($h->{s20}*1000)/10; ; 
	$z =  (($h->{s30}/$h->{nb}));
	 $hjson->{$pid}->{"30x"} = int($h->{s30}*1000)/10; ; 
	$z =  ($h->{sum}/$h->{nb});
	 $hjson->{$pid}->{"mean"} = int($h->{"mean"}*10)/10; 
	 $z =  (($h->{s100}/$h->{nb}));
	 $hjson->{$pid}->{"100x"} = int($h->{"s100"}*1000)/10; 
			
	
	}
		print encode_json  $hjson;
		exit(0);
}

if ($project->isGenome){
	 
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $h ) = @_;

			unless ( defined($h) or $exit_code > 0 ) {
				print
				  qq|No message received from child process $exit_code $pid!\n|;
				die();
				return;
			}
		 my $patient = $h->{patient};
		 $res->{$patient}->{s5} += $h->{s5};
		 $res->{$patient}->{s1} += $h->{s1};
		 $res->{$patient}->{s30} += $h->{s30};
		 $res->{$patient}->{s15} += $h->{s15};
		 $res->{$patient}->{s20} += $h->{s20};
		 $res->{$patient}->{s100} += $h->{s100};
		 $res->{$patient}->{sum} += $h->{sum};
		 $res->{$patient}->{nb} += $h->{nb};
		}
	);
	
	foreach my $patient (@{$patients}){
	foreach my $chr (@{$project->getChromosomes}){
			my $intervals = $buffer->divide_by_chunks(1,$chr->length,50_000_000);
			#warn $from." ".$to;
			#warn Dumper $intervals;
			#die();
			foreach my $interval (@$intervals){
				my $pid = $pm->start and next;
				my $h;
				my $s5;
				my $s30;
				my $nb;
				my $s15;
				my $s100;
				my $s1;
				my $s20 =0;
				my $array = $patient->depth($chr->name,$interval->[0],$interval->[1]);
				my $sum = sum @$array;
				#$all_sum = $sum;
				$nb = scalar(@$array);
				foreach my $a (@$array){
					$s1 ++ if $a >= 1;
					$s5 ++ if $a >= 5;
					$s15 ++ if $a >= 15;
					$s20 ++ if $a >= 20;
					$s30 ++ if $a >= 30;
					$s100 ++ if $a >= 100;
				}
		#warn $chr->name." ".$all_sum/$nb." ".(($s5/$nb)*100)." ".(($s30/$nb)*100);
			 $h = {s5=>$s5,s15=>$s15,s30=>$s30,s100=>$s100,s20=>$s20,patient=>$patient->name,nb=>$nb,sum=>$sum} ;
		 	$pm->finish( 0,$h );
		 }
		
		}
	}
$pm->wait_all_children();
}
else {
$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $h ) = @_;

			unless ( defined($h) or $exit_code > 0 ) {
				print
				  qq|No message received from child process $exit_code $pid!\n|;
				die();
				return;
			}
		 my $patient = $h->{patient};
		 $res->{$patient}->{s5} += $h->{s5};
		 $res->{$patient}->{s1} += $h->{s1};
		 $res->{$patient}->{s30} += $h->{s30};
		 $res->{$patient}->{s15} += $h->{s15};
		  $res->{$patient}->{s20} += $h->{s20};
		 $res->{$patient}->{s100} += $h->{s100};
		 $res->{$patient}->{sum} += $h->{sum};
		 $res->{$patient}->{nb} += $h->{nb};
		 
#		 warn  $patient." ".$res->{$patient}->{sum};
		}
	);
foreach my $patient (@{$patients}){
	my $all_sum;
	
	foreach my $chr (@{$project->getChromosomes}){
			#next if $chr->name ne "1";
			my $intspan = $patient->getCapture->genomic_span($chr);
			my @t = split(",",$intspan->as_string);
			next if $intspan->is_empty();
			my $pid = $pm->start and next;
			my $array = $patient->depthIntspan($chr->name,$intspan);
			#warn $from." ".$to;
			#warn Dumper $intervals;
			#die();
		#	foreach my $interval (@$intervals){
				
				my $s5;
				my $s30;
				my $nb;
				my $s15;
				my $s100;
				my $s1;
				my $s20;
				my $sum = sum @$array;
				$all_sum = $sum;
				$nb = scalar(@$array);
				foreach my $a (@$array){
					$s1 ++ if $a >= 1;
					$s5 ++ if $a >= 5;
					$s20 ++ if $a >= 20;
					$s15 ++ if $a >= 15;
					$s30 ++ if $a >= 30;
					$s100 ++ if $a >= 100;
				}
		#	}
		#warn $chr->name." ".$all_sum/$nb." ".(($s5/$nb)*100)." ".(($s30/$nb)*100);
		$pm->finish( 0, {s20=>$s20,s5=>$s5,s15=>$s15,s30=>$s30,s100=>$s100,patient=>$patient->name,nb=>$nb,sum=>$sum} );
		}
	}
$pm->wait_all_children();
}


my $hjson; 

foreach my $patient (@{$patients}){
	my $name = $patient->name;
	my $pid = $patient->id; 
	
	my $coverage_file;
	$coverage_file = $patient->getCoverageFile();
	my $bed_coverage = $patient->getCoverageFile();
	$bed_coverage =~ s/\.gz//;
	open(BED,">$bed_coverage");
	warn $coverage_file;
	#die if -e $coverage_file;
	#die();
	warn $name."\n";
	my $z =  $res->{$name}->{nb};#/$res->{$name}->{nb}));
	 $hjson->{$pid}->{"nb"} = $z; 
	print BED "mean_all\t1\t".$z."\n";
	 $z= (($res->{$name}->{s5}/$res->{$name}->{nb}));
	 $hjson->{$pid}->{"5x"} = int($z*1000)/10;  
	print BED "mean_all\t5\t".$z."\n";
	$z = (($res->{$name}->{s15}/$res->{$name}->{nb}));
	 $hjson->{$pid}->{"15x"} = int($z*1000)/10; ; 
	print BED "mean_all\t15\t".$z."\n";
	$z = (($res->{$name}->{s20}/$res->{$name}->{nb}));
	 $hjson->{$pid}->{"20x"} = int($z*1000)/10; ; 
	print BED "mean_all\t20\t".$z."\n";
	$z =  (($res->{$name}->{s30}/$res->{$name}->{nb}));
	 $hjson->{$pid}->{"30x"} = int($z*1000)/10; ; 
	print BED "mean_all\t30\t".$z."\n";
	$z =  ($res->{$name}->{sum}/$res->{$name}->{nb});
	 $hjson->{$pid}->{"mean"} = int($z*10)/10; 
	print BED "mean_all\t99\t".$z."\n";
	$z =  (($res->{$name}->{s100}/$res->{$name}->{nb}));
	 $hjson->{$pid}->{"100x"} = int($z*1000)/10; 
	print BED "mean_all\t100\t".$z."\n";
	close BED;
	system("$bgzip -f $bed_coverage; $tabix -b 2 -e 2 -s 1 -f $coverage_file");
}

if ($json) {
	print encode_json  $hjson;
}
exit(0);

sub tabix {
	my ($paramparam) 
	
}
