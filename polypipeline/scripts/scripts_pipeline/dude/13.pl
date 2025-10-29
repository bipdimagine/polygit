#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use GenBoNoSqlLmdb;

use Carp;
use strict;
use Set::IntSpan::Fast::XS ;
use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Carp;
 use JSON::XS;
 use List::MoreUtils qw(natatime uniq);
 use Statistics::Zscore; 
  use Statistics::Descriptive;
  my $buffer = new GBuffer;
my $project_name= "NGS2017_1534";
my $fork;
GetOptions(
	'project=s' => \$project_name,
	'fork=s' => \$fork,
);
die("man !!!  no fork ") unless $fork;
my $project = $buffer->newProjectCache( -name 			=> $project_name );
my $no =  $project->noSqlCnvs("r");
warn "start";
my $plexn_count = $no->get("raw_data","plexn_count");

#my $cai_count = $no->get("raw_data","cai_count");;

my $can_count = $no->get("raw_data","can_count");
$project->getPatients;
my $plexi_count = $no->get("raw_data","plexi_count");
my $controls = $no->get("raw_data","controls");
warn "end";
my $pm2 = new Parallel::ForkManager($fork);
my $chromosomes = $project->getChromosomes();
 delete $project->{noSqlCnvs};
 $project->buffer->dbh_deconnect();
 foreach my $chr (@{$chromosomes}){
 		warn "start ".$chr->name;
 		#next if $chr->name ne "X";
 	   my $pid      = $pm2->start() and next;	
 	   my $cai_count;
 	# $project->buffer->dbh_reconnect();
 	#warn $chr->name;
 	my $nb =0;
 	my $no = $chr->get_lmdb_cnvs("w");
 	my $no_raw =  $project->noSqlCnvs("r");
 	my $t = time;
 	while (my $primer = $chr->next_primer) {
 		$nb ++;
 		warn $nb." ".$chr->ucsc_name." ".abs(time-$t)  if $nb %2000 == 0;
 		 my $runs = $primer->getRuns();
		
		foreach my $run (@$runs) {
			
			my $run_id = $run->id;
			my $debug;

		#next unless $debug;
		my $patients = $primer->getPatientsByRun($run);
			
			if ($primer->id =~ /chr17_65943807/){
				warn " ".$run->id." ".scalar(@$patients);
		}
			foreach my $patient ( @{$primer->getPatientsByRun($run)}) {
			#	warn $patient->name()." ".$run_id;
				$primer->{cnv}->{$patient->id} = -1;
				$primer->{level}->{$patient->id} = -1;
				$primer->{zscore_data}->{$run->id} = [];
			}
			next if ($chr->name eq "Y");
			my $m1 = $can_count->{$primer->id}->{$run_id} / scalar(@$patients);
			my $limit_m = 5 ;
			$limit_m = 10 if $project->isDiagnostic;
			warn  $can_count->{$primer->id}->{$run_id}  if $debug;
			warn $m1 ." ==>  ".$limit_m if $debug;
			next if $m1 < $limit_m; ##
			my $den = ($can_count->{$primer->id}->{$run_id}/$plexn_count->{$run_id});
			my @data = ();
			my $hdata;
			my $cai_count->{$primer->id} = $no_raw->get("raw_data",$primer->id."_cai") unless exists $cai_count->{$primer->id} ;
			foreach my $c (@{$controls->{$run_id}}){
				die() unless exists $cai_count->{$primer->id}->{$c->{id}};
			 	my $score = ($cai_count->{$primer->id}->{$c->{id}}/$plexi_count->{$c->{id}}) /$den;
			 	$score = int($score*100)/100; 
			 	warn $score." ".$c->{id} if $debug;
				push(@data,$score*1.0)   if $score > 0.7 && $score < 1.3;
				$hdata->{$c->{id}} = $score   if $score > 0.7 && $score < 1.3;
			}
		
			next if scalar(@data) < 1;##
#				die() if $debug;
			$primer->{is_primer_ok} = 1;
			my $stats= new Statistics::Descriptive::Full;
			$stats->add_data(@data);
			my $sd = $stats->standard_deviation();
			$primer->{zscore_data}->{$run->id} = \@data;
			$primer->{sd}->{$run->id} = $sd;
			foreach my $patient ( @{$primer->getPatientsByRun($run)}) {
				my @data2 =sort {$a <=> $b} map{$hdata->{$_}} grep{$_ ne $patient->id} keys %$hdata;
				
			
				$primer->{sd}->{$patient->id} = $sd;
				my $score = 0;
				my $res = ($cai_count->{$primer->id}->{$patient->id}/$plexi_count->{$patient->id}) /$den;
		 		$score = int($res*100)/100; 
		 		if (scalar(@data2) > 5){
					pop(@data2);
					shift(@data2);
				}
		 		my $z = Statistics::Zscore->new;
				my $zscore1 = $z->standardize([$score,@data2]);
				my $zscore =$zscore1->[0];
				$primer->{cnv}->{$patient->id} = $score;
				
				$primer->{zscore}->{$patient->id} = $zscore;
				if ($debug){
					warn "level ".$patient->name;
					warn "level :: ".$primer->compute_level($patient,$debug);
					#die();
				}
				$primer->{level}->{$patient->id} = $primer->compute_level($patient);
				
			}
		}
		delete $primer->{project};
		delete $primer->{buffer};
 		 $no->put($primer->id,$primer);
 	}
 	$no->close();
 	 $pm2->finish;
	}
	warn "wait";
	$pm2->wait_all_children;
