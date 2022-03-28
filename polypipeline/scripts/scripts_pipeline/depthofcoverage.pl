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
use Storable qw(store retrieve freeze);
use List::MoreUtils qw(natatime);

my $project_name;
my $patient_name;
my $fileout;
my $filein;

GetOptions(
	'project=s'   => \$project_name,
	"patients=s" => \$patient_name,
	"fileout=s" => \$fileout,
	"filein=s" => \$filein,
);



my $other_project = [];
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatient($patient_name);
my $bam = $patient->getBamFile();

my $gatk  = $project->getSoftware('gatk');
my $reference = $project->getGenomeFasta;
my $javac = $project->getSoftware('java');
$javac = "java" unless -e $javac;

my $capture = $patient->getCapture();
my @beds;
my $dir_out_gvcf= $project->getCallingPipelineDir("depthofcoverage");
warn $dir_out_gvcf;
my $pm = new Parallel::ForkManager(16);
my $first_line;
my $all_data;
$pm->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    		my $chr_name = $data->{chromosome};
    		my $file = $data->{file}.".sample_interval_summary";
    		my $start = $data->{start};
    		if (-e $file){
    			my @d = `cat $file`;
    				unless ($first_line){
    				$first_line = shift(@d);
    			}
    			else {
    				shift(@d);
    			}
    			
    			$all_data->{$chr_name}->{$start} =\@d;
    			system("rm ".$data->{file}."*");
    			
    		}
    }
  );
my $nb;


foreach my $chr (@{$project->getChromosomes}){
		warn "start ".$chr->name;
		my $intspan = $capture->getIntSpanForChromosome($chr,50);
		my $bed_file = $dir_out_gvcf."/".$patient->name.".".$chr->name.".bed";
			my $chr_name = $chr->name;
		#my $chr_ucsc_name = $chr->ucsc_name;
#		open(BED,">".$bed_file);
#
#		print BED join("\n",$buffer->intspanToBed($chr,$intspan));
#		close BED;
	#last if $chr->name eq "3";
		my @beds = $buffer->intspanToBed($chr,$intspan);
		 my $iter = natatime 1000, @beds;
		 while( my @tmp = $iter->() ){
		 	my $pid = $pm->start and next;
		 	my ($chr_ucsc_name,$start,$x) = split(" ",$tmp[0]);
		 	my ($chr_ucsc_name2,$start2,$end) = split(" ",$tmp[-1]);
			my $bed_file = $dir_out_gvcf."/".$patient_name.".".$chr_ucsc_name.".$start.$end.bed";
			open(BED,">".$bed_file);
			print BED join("\n",@tmp);
			close(BED);
			my $fileoutemp = $dir_out_gvcf."/".$patient_name.".".$chr_ucsc_name.".$start.$end.data";
			
				warn "$chr_ucsc_name $start $end";
				my $cmd = qq{$javac -jar $gatk -T DepthOfCoverage -R $reference -I $filein -L $bed_file --minBaseQuality 20 --start 1 --stop 5000 --nBins 200  --countType COUNT_FRAGMENTS -o $fileoutemp -l off};
				system($cmd);
				unlink $bed_file;
				my $result;
				$result->{chromosome} = $chr_name;
				$result->{file} = $fileoutemp;
				$result->{start} = $start;
				$result->{end} = $end;
				$pm->finish(0,$result);
		}
}

$pm->wait_all_children;

open (OUT,">$fileout");
print OUT  $first_line;
foreach my $chr (@{$project->getChromosomes}){
	foreach my $key (sort {$a<=>$b} keys %{$all_data->{$chr->name}}){
		print OUT @{$all_data->{$chr->name}->{$key}};
	}
}
close(OUT);
system ("rm $dir_out_gvcf/$patient_name*");
warn $fileout;

