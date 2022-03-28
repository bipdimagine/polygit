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
my $project_name= "NGS2017_1534";
my $fork;
GetOptions(
	'project=s' => \$project_name,
	'fork=s' => \$fork,
);
die("he no fork ") unless $fork;
my $jobs;
my $project = $buffer->newProjectCache( -name 			=> $project_name );
my $runs = $project->getRuns();

#choose control by run .
my $controls;
my $max_controls = 6;
foreach my $r (@$runs){
	my $hps =  $r->getAllPatientsInfos();
	warn $r->name;
	warn Dumper $hps if $r->id eq 1693;
	next if $r->id ne 1693;
	#only control
#	if (scalar(@$hps) < 8){
#		push(@{$controls->{$r->id}},@$hps);
#		next;
#	}
	my @hps2 =  grep{$_->{project} ne $project_name && $_->{status} eq 1} @$hps;

	if (@hps2>=$max_controls){
		@hps2 = samples $max_controls,@hps2;
		#push(@{$controls->{$r->id}},samples 6,@hps2);
	}
	
	if (@hps2<=$max_controls) {
		
		foreach my $h  (grep{$_->{project} ne $project_name && $_->{status} ne 1} @$hps){
			push(@hps2,$h);
			last if scalar(@hps2) >= $max_controls;
		}
			
		}

	if (@hps2<=$max_controls) {
		foreach my $h  (grep{$_->{project} eq $project_name && $_->{status} eq 1} @$hps){
			push(@hps2,$h);
			last if scalar(@hps2) >= $max_controls;
			
		}
	}
	#push(@{$controls->{$r->id}},@hps2);
	my %contr_projects;
	map {$contr_projects{$_->{project}} ++} @hps2;
	foreach my $pr (keys %contr_projects){
		my $buffer1 = new GBuffer;
		my $project1 = $buffer1->newProjectCache( -name 			=> $pr );
		foreach my $p (grep{$_->{project} eq $pr} @hps2){
			my $patient = $project1->getPatient($p->{patient});
			
			my $cov_file = $patient->fileNoSqlDepth;
			$p->{file_depth} = $cov_file;
			$p->{dir_depth} = $project1->getCoverageDir()."/lmdb_depth";
			warn 	$p->{dir_depth};
			#warn $patient->name." ".$pr  unless -e $cov_file;
			next unless -e $cov_file;
			$p->{file} = $cov_file;
			$p->{bam} = $patient->getBamFile();
			$p->{bd_sex} = $patient->sex;
			push(@{$controls->{$r->id}},$p);
		} 
		
	}
}


my $chr = $project->getChromosome("Y");

foreach my $r (keys %$controls){
	warn $r;
	foreach my $c (@{$controls->{$r}}){
		my $f = $c->{bam};
		
	#	warn Dumper $f;
	 	my $sambamba = $buffer->software("sambamba");
	 	my $nc = $project->getChromosome("Y")->fasta_name;
		my @vv =`$sambamba depth region $f -L $nc:2654896-2655740 2>/dev/null | cut -f 5`; #`tabix $f mean_chrY:1-15| cut -f 3`;
	 	#warn "$sambamba depth region $f -L $nc:2654896-2655740";
	 	my $v = $vv[-1];
		chomp($v);
		if ($c->{bd_sex} == 2) { 
			if ( $v < 3  ) {
				$c->{sex} = 2;
				}
				else {
					confess();
				}
		}
		else {
			if ( $v > 3  ) {
				$c->{sex} = 1;
			}
			else {
				warn $c->{bd_sex}."  $v";
				confess();
			}
		}
		
	}
	
}

warn Dumper $controls;
die();

foreach my $r (@$runs){
	warn $r->name();
	warn $r->id;
	die() if scalar(@{$controls->{$r->id}}) == 0;
}


my $primers;
foreach my $chr (@{$project->getChromosomes()}){
	#next if $chr->name eq "1";
	foreach my $p (sort{$a->start <=> $b->start} @{$chr->getPrimers}){
		push(@$primers,$p);
		$p->getPatients;
		$p->getRuns();
	}
}
my $nb2 = scalar(@$primers);
warn $nb2;
my $z=0;
my @name = map{$_->id} @$primers;
my $nb = int(scalar(@name)/($fork)) -1;
 my $it = natatime $nb, @name;
 my $pm2 = new Parallel::ForkManager($fork);
 #my $no =  $project->noSqlCoverage("c");
 my $global_res;
 my $cai_count;
 my $can_count;
 my $plexi_count;
 my $plexn_count;
 
foreach my $patient ( @{$project->getPatients}){
	my $run_id = $patient->getRun->id;
}
foreach my $primer (@$primers){
	$primer->getPatients();
}
 my $intspan_pseudo = Set::IntSpan::Fast::XS ->new("60001-2699520","154931044-155260560");
 my $error ;
 $pm2->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) = @_;
			$error ++  if $exit_code ne 0;
			warn "end process";
			my $count;
			 my $no =  $project->noSqlCnvs();
			 $project->buffer->dbh_reconnect();
			 #cai count
			 my $tcai_count = $data->{cai_count};
			 foreach my $primer_id (keys %$tcai_count){
			 	foreach my $patient_id (keys %{$tcai_count->{$primer_id}}){
			 			$cai_count->{$primer_id}->{$patient_id} = $tcai_count->{$primer_id}->{$patient_id};
			 	}
			 }
			 	 #can count
			  my $tcan_count =    $data->{can_count};
			  foreach my $primer_id (keys %$tcan_count){
			 	foreach my $run_id (keys %{$tcan_count->{$primer_id}}){
			 		
			 		$can_count->{$primer_id}->{$run_id} = $tcan_count->{$primer_id}->{$run_id}
			 	}
			  }
			  
			  #plexi_count 
			   my $tplexi_count =   $data->{plexi_count};
			   foreach my $patient_id (keys %$tplexi_count){
			   		$plexi_count->{$patient_id} +=  $tplexi_count->{$patient_id} ;
			    }
			  
			  #plexn_count
			  my $tplexn_count = $data->{plexn_count};
			   	foreach my $run_id (keys %{$tplexn_count}){
			   		$plexn_count->{$run_id} += $tplexn_count->{$run_id};
			   	}
			   	#samples
			
			  
			  $project->buffer->dbh_deconnect();
			   delete  $project->{noSqlCnvs};
		}
	);
 
 my $limit =0;
 $project->buffer->dbh_deconnect();
 delete  $project->{noSqlCnvs};
 my $process =0;
while (my @vals = $it->())
  { 
  	last if $error;
  	$process ++;
  my $pid      = $pm2->start() and next;	
   $project->buffer->dbh_reconnect();
   my $tabix = $buffer->software('tabix');
	my $h;
	map{$h->{$_} ++} @vals;
	my $primers = $project->myflushobjects($h,"primers");
	
	my $res;
	my $htabix;
	foreach my $primer (@$primers){
			my $debug;
			$debug =1 if $primer->id eq "primerchr4_55599161";
#		warn "$process ==>" .$limit."/".scalar (@$primers) if $limit%100 == 0;
		$limit ++;
		my @bams = map{$_->getBamFile} @{$primer->getPatients};
		foreach my $patient (@{$primer->getPatients}){
				
				 my $cov1 = return_cov_mean($patient,$primer,$patient->sex);
				# warn $cov1." ".$hash->{mean};
				$res->{mean}->{primers}->{$primer->id}->{$patient->id} = $cov1;
				$res->{mean}->{patients}->{$patient->id}->{$primer->id} = $cov1;
				#$res->{mean}->{patients}->{$patient->id}->{$primer->id} = $hash->{mean};
				$res->{plexi_count}->{$patient->id} += $cov1;
				#$res->{plexn_count}->{$patient->getRun->id} += $hash->{mean};
				
				$res->{cai_count}->{$primer->id}->{$patient->id} = $cov1;
				
		}
		my $st_pos = $primer->getChromosome->ucsc_name().":".$primer->start."-".$primer->end;
		my $len = abs($primer->start-$primer->end) +1;
		foreach my $run (@{$primer->getRuns}){
			my $ctrl = $controls->{$run->id};
			my @covs;
			foreach my $pc (@$ctrl){
			
				my $id = $pc->{id};
				my $cov =0;
				if (exists $res->{mean}->{patients}->{$id}->{$primer->id} ){
					$cov = $res->{mean}->{patients}->{$id}->{$primer->id};
				}
				else {
				my $file = $pc->{file};
				 $cov = return_cov_mean_depth($pc->{patient},$pc->{dir_depth},$primer,$pc->{sex});
				 
				 $res->{plexi_count}->{$id} += $cov;
				}
				$res->{mean}->{primers}->{$primer->id}->{$id}  = $cov;
				$res->{cai_count}->{$primer->id}->{$id} = $cov;
				
				warn $cov if $debug;
				
				push(@covs,$cov);
				
				#my $cov = `tabix $file $st_pos | awk '{sum += \$3} END {print sum}'`;
				#chomp($cov);
					
			#	$res->{plexn_count}->{$run->id} += $cov;
				#$res->{can_count}->{$primer->id}->{$run->id} += $cov;
			}
			@covs = sort{$a <=> $b} @covs;
			if (@covs >10){
				shift(@covs);
				pop(@covs);
			}
			
			if ($debug){
				  my $stat = Statistics::Descriptive::Full->new();
				   $stat->add_data(@covs);
				  my ($q, $m, $r, $rms) =  $stat->least_squares_fit();
				  warn "q:$q m:$m r:$r rms:$rms" . " ". join(";",@covs);
			}
			warn Dumper @covs  if $debug;  
			my $sum = sum(@covs);
			$res->{plexn_count}->{$run->id} += $sum;
				$res->{can_count}->{$primer->id}->{$run->id} += $sum;
			
		}
		
	}
	
	
	
  	$pm2->finish(0,$res);
	}
	$pm2->wait_all_children;
#	exit(0);
 $project->buffer->dbh_reconnect();
 die() if $error;
#warn Dumper $plexn_count;

#delete  $project->{noSqlCnvs};

my $no =  $project->noSqlCnvs("c");
$no->put("raw_data","cai_count",$cai_count);
$no->put("raw_data","can_count",$can_count);
$no->put("raw_data","plexn_count",$plexn_count);
$no->put("raw_data","plexi_count",$plexi_count);
$no->put("raw_data","controls",$controls);

exit(0);

my $htabix2;
sub return_cov_mean{
	my ($patient,$primer,$sex) = @_;
	#return return_cov_mean_tabix($patient->id,$patient->getCoverageFile(),$primer,$sex); 
	if ($patient->isNoSqlDepth){
		#my $sum =$patient->depth($primer->getChromosome()->name,$primer->start,$primer->end);
		
		#return $sum/(abs($primer->getChromosome()->start - $primer->getChromosome()->end) +1);

		my $array = $patient->depth($primer->getChromosome()->name,$primer->getChromosome()->start,$primer->getChromosome()->end);
		my $s = sum(@$array);
		my $cov = $s/scalar(@$array);
		if ($sex == 1  && !($primer->getChromosome()->isAutosomal($primer->start,$primer->end))) {
			$cov *= 2 ;
	 } 
	 return $cov;
	}
	
	return return_cov_mean_tabix($patient->id,$patient->getCoverageFile(),$primer,$sex); 
}
my $hdepth;
sub return_cov_mean_depth{
	my ($name,$dir,$primer,$sex) = @_;
	
	$hdepth->{$name} =  GenBoBinaryFile->new(name=>$name.".depth.lmdb",dir=>$dir,mode=>"r") unless exists $hdepth->{$name};
	my $array = $hdepth->{$name}->getDepth($primer->getChromosome()->name,$primer->start,$primer->end);
		my $s = sum(@$array);
		my $cov = $s/scalar(@$array);
		if ($sex == 1  && !($primer->getChromosome()->isAutosomal($primer->start,$primer->end))) {
			$cov *= 2 ;
	 } 
	 return $cov;
}

sub return_cov_mean_tabix{
	my ($id,$file,$primer,$sex) = @_;
	#warn $id ." ".$file;
	my $len = abs($primer->start-$primer->end) +1;
	$htabix2->{$id} = new Tabix(-data =>$file) unless exists $htabix2->{$id};
	my $query = $htabix2->{$id}->query($primer->getChromosome->ucsc_name,$primer->start,$primer->end);
		warn "$file "." ".$primer->getChromosome->ucsc_name.":".$primer->start."-".$primer->end  unless $query->get(); 
		my $cov =0;	
	while (my $line = $htabix2->{$id}->read($query)) {
		my ($thisChr, $thisPos, $thisCov) = split("\t", $line);
			$cov += $thisCov;
		}
	$cov = $cov/$len;				
	if ($sex == 1  && !($primer->getChromosome()->isAutosomal($primer->start,$primer->end))) {
		
		#		warn "coucou" ." ".$primer->getChromosome()->name;
			$cov *= 2 ;
	 } 
	
	return $cov;
}


