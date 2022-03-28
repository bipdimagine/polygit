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
GetOptions(
	"fork=s"   => \$fork,
	"project=s" =>\$project_name,
	"patient=s" =>\$patient_name,
	"verbose=i" =>\$verbose,
);

my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => $project_name );
my $pm   = new Parallel::ForkManager($fork);
my $tabix = $buffer->software("tabix");
my $bgzip = $buffer->software("bgzip");
my $res;
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
		 $res->{$patient}->{s100} += $h->{s100};
		 $res->{$patient}->{sum} += $h->{sum};
		 $res->{$patient}->{nb} += $h->{nb};
		 warn  $patient." ".$res->{$patient}->{sum};
		}
	);
$patient_name="all" unless $patient_name;	
my $patients = $project->get_list__controls_patients($patient_name);
foreach my $patient (@{$patients}){
	my $all_sum;
	
	foreach my $chr (@{$project->getChromosomes}){
			my $intervals = $buffer->divide_by_chunks(1,$chr->length,50_000_000);
			#warn $from." ".$to;
			#warn Dumper $intervals;
			#die();
			foreach my $interval (@$intervals){
				my $pid = $pm->start and next;
				my $s5;
				my $s30;
				my $nb;
				my $s15;
				my $s100;
				my $s1;
				my $array = $patient->depth($chr->name,$interval->[0],$interval->[1]);
				my $sum = sum @$array;
				$all_sum = $sum;
				$nb = scalar(@$array);
				foreach my $a (@$array){
					$s1 ++ if $a >= 1;
					$s5 ++ if $a >= 5;
					$s15 ++ if $a >= 15;
					$s30 ++ if $a >= 30;
					$s100 ++ if $a >= 100;
				}
		#warn $chr->name." ".$all_sum/$nb." ".(($s5/$nb)*100)." ".(($s30/$nb)*100);
		$pm->finish( 0, {s5=>$s5,s15=>$s15,s30=>$s30,s100=>$s100,patient=>$patient->name,nb=>$nb,sum=>$sum} );
		}
	}
}
$pm->wait_all_children();

foreach my $patient (@{$patients}){
	my $name = $patient->name;
	my $coverage_file;
	$coverage_file = $patient->getCoverageFile();
	my $bed_coverage = $patient->getCoverageFile();
	$bed_coverage =~ s/\.gz//;
	open(BED,">$bed_coverage");
	warn $coverage_file;
	#die if -e $coverage_file;
	#die();
	print $name."\n";
	my $z =  $res->{$name}->{nb};#/$res->{$name}->{nb}));
	print BED "mean_all\t1\t".$z."\n";
	 $z= (($res->{$name}->{s5}/$res->{$name}->{nb}));
	print BED "mean_all\t5\t".$z."\n";
	$z = (($res->{$name}->{s15}/$res->{$name}->{nb}));
	print BED "mean_all\t15\t".$z."\n";
	$z =  (($res->{$name}->{s30}/$res->{$name}->{nb}));
	print BED "mean_all\t30\t".$z."\n";
	$z =  ($res->{$name}->{sum}/$res->{$name}->{nb});
	print BED "mean_all\t99\t".$z."\n";
	$z =  (($res->{$name}->{s100}/$res->{$name}->{nb}));
	print BED "mean_all\t100\t".$z."\n";

	close BED;
	system("$bgzip -f $bed_coverage; $tabix -b 2 -e 2 -s 1 -f $coverage_file");
}