#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";

#use Set::IntSpan;
use GenBoNoSqlLmdb;

use Carp;
use strict;
use Set::IntSpan::Fast::XS ;
use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Carp;
 use JSON::XS;
 use List::MoreUtils qw(natatime uniq samples  );
  use List::Util qw(sum );
 use Tabix;
 
#use GenBoNoSql;
#use GenBoNoSqlDejaVu;
#use Array::Diff;
#use UnQLite;
 # use Tie::LevelDB; 
 # use Devel::Size qw(size total_size);
  #use Devel::Size::Report qw/report_size/;
  #use Statistics::Descriptive;
  

 $| =1;               
#ENSG00000234585
my $buffer = new GBuffer;
my $project_name= "";
my $fork;
GetOptions(
	'project=s' => \$project_name,
	'fork=s' => \$fork,
);
die("he no fork ") unless $fork;
my $jobs;
my $project = $buffer->newProjectCache( -name 			=> $project_name );
my $runs = $project->getRuns();

my $no =  $project->noSqlCnvs("w");
my $control = $no->get("raw_data","controls");

my $hprojects;
my $no;
my $list_controls; 
my $patient_runs;
foreach my $run (keys %{$control}) {
	foreach my $hp (@{$control->{$run}}) {
		my $name = $hp->{patient};
		$list_controls->{$hp->{id}}++;
		$patient_runs->{$hp->{id}} = $run;
		my $dir = $hp->{file_depth};
		$hp->{no} = GenBoBinaryFile->new(name=>$name.".depth.lmdb",dir=>$hp->{dir_depth},mode=>"r");
		my $pr =$hp->{$run}->{project};
		push(@{$hprojects->{$run}->{$pr}},$hp);			
	}
}



foreach my $patient (@{$project->getPatients}){
			my $r = $patient->getRun();
			my $cov_file = $patient->fileNoSqlDepth;
			my $p;
			die() unless -e $cov_file;
			next if exists $list_controls->{$patient->id};
			$patient_runs->{$patient->id} = $r->id;
			$p->{file_depth} = $cov_file;
			$p->{dir_depth} = $project->getCoverageDir()."/lmdb_depth";
			$p->{bam} = $patient->getBamFile();
			$p->{bd_sex} = $patient->sex;
			$p->{sex} = $patient->compute_sex;
			$p->{id} = $patient->id;
			$p->{patient} = $patient->name;
			
			push(@{$control->{$r->id}},$p);
			
	
}
 my $pm2 = new Parallel::ForkManager($fork);
  my $cai_count;
  my $plexi_count;
  my $can_count;
  my $plexn_count;
  my $error;
  
 $pm2->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) = @_;
			$error ++  if $exit_code ne 0;
			warn "end process";
			my $count;
			my $cai = $data->{cai};
			my $no =  $project->noSqlCnvs("w");
			 #cai count
			 my $tcai_count = $data->{cai_count};
			 foreach my $primer_id (keys %$tcai_count){
			 	my $sum =0;
			 	$no->put("raw_data",$primer_id."_cai",$tcai_count->{$primer_id});
			 	my $covs;
			 	my $hrids;
			 	foreach my $patient_id (keys %{$tcai_count->{$primer_id}}){
			 		 	$plexi_count->{$patient_id} +=  $tcai_count->{$primer_id}->{$patient_id};
			 			$cai_count->{$primer_id}->{$patient_id} = $tcai_count->{$primer_id}->{$patient_id};
			 			
			 			if (exists $list_controls->{$patient_id}) {
			 				my $rid  = $patient_runs->{$patient_id};
			 				push(@{$covs->{$rid}}, $tcai_count->{$primer_id}->{$patient_id});
			 				
			 			
			 				$hrids->{$rid} ++;
			 				#$can_count->{$primer_id}->{$rid} +=  $tcai_count->{$primer_id}->{$patient_id};
			 				#$plexn_count->{$rid} += $tcai_count->{$primer_id}->{$patient_id};
			 				
			 			}
			 	}
			 foreach my $rid (keys %$hrids){	
			 	my @covs = @{$covs->{$rid}};
			 	@covs = sort{$a <=> $b} @covs;
				if (@covs >10){
					shift(@covs);
					pop(@covs);
				}
					my $sum = sum(@covs);
					$can_count->{$primer_id}->{$rid}  += $sum;
					$plexn_count->{$rid} +=  $sum;
			 	}
			 }
			 $no->close();
			 delete $project->{noSqlCnvs};
		}
	);
	
	my $hdepth;
	$project->getChromosomes;
	 $project->buffer->dbh_deconnect();
foreach my $chr (@{$project->getChromosomes}){
  my $pid      = $pm2->start() and next;	
   $project->buffer->dbh_reconnect();
	my $nb =0;
	warn "start ".$chr->name;
	my $res;
	while (my $primer = $chr->next_primer){
		warn $nb if $nb %50000 ==0;
		$nb++;
		
			my $runs = $primer->getRuns();
		
			foreach my $run (@$runs){
				my $run_id =  $run->id;
				foreach my $hp (@{$control->{$run_id}}) {
						my $name = $hp->{patient};
						my $dir = $hp->{dir_depth};
						$hdepth->{$name} =  GenBoBinaryFile->new(name=>$name.".depth.lmdb",dir=>$dir,mode=>"r") unless exists $hdepth->{$name};
						#warn $primer->start." ".$primer->end;
						my $array = $hdepth->{$name}->getDepth($chr->name,$primer->start,$primer->end);
						my $cov = sum(@$array)/scalar(@$array);
						if ($chr->name eq "X"){
							
							if ($hp->{sex} == 1 ){# && !($chr->isAutosomal($primer->start,$primer->end))) {
									$cov *= 2 unless $chr->isPseudoAutosomal($primer->start,$primer->end) ;
									
						 	} 
						
						}
						
						$res->{cai_count}->{$primer->id}->{$hp->{id}} = $cov;
						#warn $cov;	
				}
				
				
			}
		
	}	
	warn "end ".$chr->name;
  	$pm2->finish(0,$res);
	}
	$pm2->wait_all_children;
	$project->buffer->dbh_reconnect();
	my $no =  $project->noSqlCnvs("w");
	#$no->put("raw_data","cai_count",$cai_count);
	$no->put("raw_data","plexi_count",$plexi_count);
	$no->put("raw_data","plexn_count",$plexn_count);
	$no->put("raw_data","plexi_count",$plexi_count);
	$no->put("raw_data","can_count",$can_count);
	
	warn "verif";
	   $project->buffer->dbh_reconnect();
	my $can_count2 = $no->get("raw_data","can_count");
	foreach my $chr (@{$project->getChromosomes}){

		while (my $primer = $chr->next_primer){
			die($primer->id) unless exists $can_count2->{$primer->id};
		}
	}
	