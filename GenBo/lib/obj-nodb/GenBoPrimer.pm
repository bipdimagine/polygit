package GenBoPrimer;

use strict;
use Moo;

use Data::Dumper;
use Config::Std;
use POSIX qw(ceil);
use Storable qw(store retrieve freeze thaw);
use Statistics::Descriptive;
use Statistics::Zscore; 
use Set::Intersection;
use List::Util qw(shuffle sum min max);

extends "GenBoGenomic";


has multiplex => (
	is		=> 'ro',
	required=> 1,
);
has intspan_pcr => (
	is		=> 'ro',
	default => sub {
		 return Set::IntSpan::Fast->new() ;
	}
);

has primersCacheFile =>(
is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
				 my $output   = $self->project->getCoverageDir . "/primers.".$self->id.".kct";
	 		 	return $output;
	
	}
);

has primersCache =>(
is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
				 my $output   = $self->primersCacheFile;
	 		 my $db1 = new KyotoCabinet::DB;
	 		 if ($self->project->no_cache){
	 		 	unlink $output if -e $output;
				if (!$db1->open($output, $db1->ONOLOCK | $db1->OCREATE |$db1->OWRITER| $db1->OREADER )){
					   
						printf STDERR ("open error: %s   => %s\n", $output,$db1->error);
						confess();
				}
					$db1->clear();	
	 		 }
	 		 else {
	 		 		if (!$db1->open($output, $db1->ONOLOCK | $db1->OREADER )){
						printf STDERR ("open error: %s\n", $db1->error);
							confess();
				}
	 		 }
	 		 return $db1;
	
	}
);

sub setPatients {
	my $self = shift;
	my %names;
	foreach my $c (@{$self->getCaptures}){
		my $patients = $c->getPatients();
		foreach my $p (@$patients){
			$names{$p->id} = undef;
		}
	}
	return \%names;
}

sub setRuns{
	my $self = shift;
	my %names;
	foreach my $p (@{$self->getPatients}) {
			
			$names{$p->getRun->id} = undef;
	}
	return \%names;
}

sub getPatientsByRun {
	my ($self,$run) = @_;
	my @patients = grep{$run->id eq $_->getRun->id} @{$self->getPatients};
	return \@patients;
}

sub setTranscripts {
	my $self = shift;
	my $debug;
#	$debug =1 if $self->id eq 'primerchr4_109108824';
	
#	foreach my $aids (@{$self->_transcriptsId()}){
#			my ($name,$chr) = split("_",$aids);
#			my $tr = $self->getProject->newTranscript($aids);	
#			if ($debug){
#				my %toto;
#				map {$toto{$_} ++ } @{$self->project->bundle_transcripts()} ;
#					
#			}
#		
#			#next unless $tr->ccds_name();
#			my $intspan = Set::IntSpan::Fast::XS->new($self->start."-".$self->end);
#				for my $ex (@{$tr->getExons}){
#			 		my $intersect = $self->getGenomicSpan->intersection($ex->getGenomicSpan);
#			 		#next if $intersect->is_empty();
#			 		
#					$hreturnTranscriptsId->{$tr->id} = undef;
#			 }
#			 
#	
#	}
	
my $hreturnTranscriptsId = $self->vannot;	
my @tt;
if ($self->project->bundle_transcripts()){
 @tt = get_intersection([keys %$hreturnTranscriptsId],$self->project->bundle_transcripts());
}
else {
	@tt = keys %$hreturnTranscriptsId;
}
if (@tt){
my $real_transcripts ;
foreach my $t (@tt){
	$real_transcripts->{$t} = undef;
}
	return $real_transcripts;
}

$self->{name} = "ctrl-".$self->{name};
return $hreturnTranscriptsId;

}
#sub setGenes {
#	my $self = shift;
#	my $hreturnGenesId ;#= $self->_transcriptsId();
#	foreach my $tr (@{$self->getTranscripts}){
#
#					$hreturnGenesId->{$tr->getGene->id} = undef;
#	
#	}
#
#	
#	return $hreturnGenesId;
#}

sub setExons {
	my $self = shift;
#	confess();
	my $hreturnExonsId ;#= $self->_transcriptsId();
	foreach my $tr (@{$self->getTranscripts}){
			my $intspan = Set::IntSpan::Fast::XS->new($self->start."-".$self->end);
				for my $ex (@{$tr->getExons}){
			 		my $intersect = $intspan->intersection($ex->getGenomicSpan);
			 		next if $intersect->is_empty();
	
					$hreturnExonsId->{$ex->id} = undef;
			 }
	
	}

	
	return $hreturnExonsId;
}



sub mean {
	my ($self,$p) =@_;
	#$self->getStatisticsCoverage($p);
	return $p->mean($self);

}




sub count {
	my ($self,$p,$compute) =@_;
	confess();
	return $self->{count_reads}->{$p->id} if exists  $self->{count_reads}->{$p->id};
	unless ($compute){
		$self->{count_reads}->{$p->id}  = $self->getData($p,"_count");
		if ($self->getChromosome->name eq "X"){
			confess();
			$self->{count_reads}->{$p->id} *=2 if $p->isMale();
		}
		return $self->{count_reads}->{$p->id}  if defined $self->{count_reads}->{$p->id} ;
	}
	
#	confess($p->name() ." ".$self->id." ".$p->transcriptsCoverageFile);
#	if  ($self->{count_reads}){
#	return $self->{count_reads}->{$p->id} if exists $self->{count_reads}->{$p->id};
#	#die() unless ($self->project->no_cache);
#	}
#	else {
#		$self->{count_reads}({});
#	}
		$self->{count_reads}->{$p->id} = 0;
	my $sam = $p->bio_db_sam();
	my @alignments = $sam->get_features_by_location(-seq_id => $self->getChromosome->ucsc_name,
                                                 -start  => $self->start,
                                                 -end    =>$self->end);
         			my @gal;
       				my $nb_start =0;
        				my $nb =scalar(@alignments);
        				
#         			foreach my $a (@alignments) {
#         					$nb++;
#         					my $start  = $a->start;
#   			if ($a->strand eq -1){
#             next if $a->end > $self->end+20;
#             next if $a->end < $self->start-20;
#   			}
#   			else {
#   				next if $start > $self->end+20;
#             	next if $start < $self->start-20;
#   			}
#             
#             $nb_start ++;
#         	
#         } 
                           
		$self->{count_reads}->{$p->id} =  $nb;
		$self->{count_start_end}->{$p->id} =  $nb_start;
		return 	$self->{count_reads}->{$p->id};

}






sub count_all {
	my ($self,$patient) = @_;
	
	my $run_id = $patient->getRun()->id;
	return $self->{count_all}->{$run_id} if exists $self->{count_all}->{$run_id};
	 $self->{count_all}->{$run_id} = 0;
	my $patients = $self->getProject->getPatients();
		my @data;
		my $nb = 0;
		foreach my $p (@$patients){
			next if $p->getRun()->id ne $run_id;
			#warn $p->name();
			$nb ++;
		#	last if $nb > 20;
			push(@data,$p->count($self->id,$self->getCaptures->[0]));
			#$self->{count_all}->{$run_id}  += $p->count($self->id,$self->getCaptures->[0]);
			
		}
		my $sum = sum @data;
#		if (scalar(@data) >10){
#		$sum -= min @data;
#		$sum -= max @data;
#		}
		$self->{count_all}->{$run_id} = $sum;
		
		return $self->{count_all}->{$run_id} ;
}


sub gene_id {
	my ($self,$patient) = @_;
	my $genes = $self->getGenes();
	my $gene_id = "none";
	$gene_id = $genes->[0]->id if scalar(@$genes);
	return $gene_id;
}

has cnv => ( 
	is		=> 'ro',
	required=>1,
	
);

sub compute_cnv_score{
	my ($self,$patient) = @_;
	my $score ;
#	warn "1";
	my $capture;
	foreach my $c (@{$self->getCaptures}){
		$capture= $c if $c->id eq $patient->getCapture->id();
	}
	return 0 unless ($capture);
	my $cai_count = $patient->count($self->id, $capture);
#	warn "2";
	my $genes = $self->getGenes();
	#warn $genes->[0]->name();
	my $gene_id = "none";
	$gene_id = $genes->[0]->id if scalar(@$genes);
		
	my $plexi_count = $patient->count_multiplex($self->multiplex,$gene_id);
	
#	warn "3";
#	warn $plexi_count;
	$plexi_count +=0;
	my $can_count = $self->count_all($patient);
	#warn "4";
	 $self->cnv->{$patient->id} = 0 unless $capture;
	#$can_count -= $cai_count;
	my $plexn_count = $capture->count_multiplex($self->multiplex,$gene_id,$patient);
#	warn "5";
	#$plexn_count -= $plexi_count;
		$score = 0;
	if ($plexi_count > 0 && $can_count >0 ){
		
		my $res = ($cai_count/$plexi_count) /($can_count/$plexn_count);
		$score = int($res*100)/100;
	}
	 $self->cnv->{$patient->id} = $score;
		return $score;
}

sub cnv_score_log2 {
	my ($self,$patient,$compute) = @_;
	return sprintf("%.2f", $patient->getProject->buffer->log2($self->cnv_score($patient,$compute)));
}

sub cnv_score {
	my ($self,$patient,$compute) = @_;
	return $self->cnv->{$patient->id} if exists $self->cnv->{$patient->id};
	#confess();
		#$self->cnv->{$patient->id}  = $self->getData($patient,"_cnv");
		#confess() if defined $self->cnv->{$patient->id} ;
		#die()	unless ($self->project->no_cache);
	
	 	$self->cnv->{$patient->id} = $self->compute_cnv_score($patient);
		return $self->cnv->{$patient->id};
}




has rawStatistics => (
	is		=> 'ro',
	lazy=>1,
	default => sub {
		my $self = shift;
		die();
		my $patients = $self->getProject->getPatients();
			my$stats= new Statistics::Descriptive::Full;
		my $score =0;
		foreach my $p (@$patients){
		
				 $stats->add_data($self->count($p));
		}
		
	
		
		return $stats;
	} 
);



sub getData {
	my ($self,$patient,$type) = @_;
	return undef;
	my $db1 = $patient->transcriptsCoverage();
	  if ($db1) {
			my $count  = $db1->get($self->id."$type");
			if (defined $count){
				#$self->{zscore}->{$patient->id} = $count;  
				return $count;
			}
	}
	return undef;
}

sub zscore {
	my ($self,$patient,$compute) = @_;
	confess() unless $patient;
	return $self->{zscore}->{$patient->id} if exists $self->{zscore}->{$patient->id};
	unless ($compute){
		$self->{zscore}->{$patient->id}  = $self->getData($patient,"_zscore");
		return $self->{zscore}->{$patient->id} if defined $self->{zscore}->{$patient->id};
	}
	
	my $patients = $self->getProject->getPatients();
	my $genes = $self->getGenes();
	my $gene_id = "none";
	$gene_id = $genes->[0]->id if scalar(@$genes);

	my $primers = $self->getCaptures->[0]->getListPrimers($self->multiplex,$gene_id);#$self->getProject->getCapture->getPrimersByMultiplex($self->multiplex());
	my $stats= new Statistics::Descriptive::Full;
	my $score =0;
	my @data;
	my @index;
	for(my $i =0 ;$i<@$patients;$i++ ){
		my $p = $patients->[$i];
		push(@index,$p->id);
		push(@data,$self->cnv_score($p));
	}
	



	$self->{zscore}->{$patient->id} = 1;
	my $z = Statistics::Zscore->new;
	my $zscore = $z->standardize( \@data );
	for(my $i =0 ;$i<@$patients;$i++ ){
		my $pid = $index[$i];
		$self->{zscore}->{$pid} = int($zscore->[$i]*100)/100;	
		}
	
	return $self->{zscore}->{$patient->id};
		
}



has statistics => (
	is		=> 'ro',
	lazy=>1,
	default => sub {
		my $self = shift;
		my $patients = $self->getPatients();
			my$stats= new Statistics::Descriptive::Full;
		my $score =0;
		foreach my $p (@$patients){
				 $stats->add_data($self->cnv_score($p));
		}
		

		
		return $stats;
	} 
);


has is_primer_ok => (
	is		=> 'ro',
	lazy=>1,
	default => sub {
		my $self = shift;
		my $patients = $self->getPatients();
		my $score =0;
		foreach my $p (@$patients){
				if ($self->cnv_score($p) > 1.5 || $self->cnv_score($p) <= 0.5){
				$score ++;
				}
		}
		return undef if $score * 4 > scalar(@$patients);
		return -1;
	} 
);

sub getStatisticsCoverage {
	my ($self,$p) =@_;
	return $self->{stats}->{$p->id} if exists $self->{stats}->{$p->id};
	my $db1 = $p->transcriptsCoverage();
	 if ( defined $db1){

			my $t = $db1->get($self->id);
			 $self->{count}->{$p->id} = $db1->get($self->id."_count");
			
			 
			if ($t){
			 my $v = thaw $t;
			
			 #die();
			  $self->{stats}->{$p->id} = new Statistics::Descriptive::Full;
		
  			   $self->{stats}->{$p->id}->add_data(@$v);
  			   
  			   return $self->{stats}->{$p->id};
  			  
			}
			
			die();
	 }	
	 else {
	 	die("restart coverage and cache");
	 }
}

has id_coverage =>(
is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->id;
	}
);

sub test_coverage {
		my ($self,$patient) = @_;
	 	my $hash ;
	 	$hash->{mean} = $patient->mean($self);
		 $hash->{min} = $patient->minimum($self);
		 return $hash;
}
sub cached_statistic_coverage{
	my ($self,$patient) = @_;
		my $no =  $self->project->noSqlCoverage();
		my $hash;
		eval{
		 $hash = $no->get($patient->name,$self->id_coverage);
		};
		$hash = {} unless $hash;
		
		
		return $hash if exists $hash->{mean} && exists $hash->{min} ;
		$hash->{mean} = $patient->mean($self);
		 $hash->{min} = $patient->minimum($self);
		 $no->set($patient->name,$self->id_coverage,$hash);
		 
		    warn $@ if $@;

		  return $hash;
		 
}

has id_cnv =>(
is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->id."_coverage";
	}
);

sub cached_cnv{
	my ($self,$patient,$debug) = @_;
		my $no =  $self->project->noSqlCoverage();
		my $hash;
		unless ($debug){
		eval{
		 	$hash = $no->get($patient->name."_cnv",$self->id);
		};
		}
		$hash = {} unless $hash;
		return $hash if exists $hash->{level_cnv} && exists  $hash->{cnv};
		$hash->{cnv} = $self->compute_cnv_score($patient);
		
		 $hash->{level_cnv} = $self->level($patient);
		 eval{
		 	$no->set($patient->name."_cnv",$self->id,$hash);
		 };
		
		  return $hash;
		 
}
sub sd {
	my ( $self , $patient) =  @_;
	return $self->{sd}->{$patient->id} if exists $self->{sd}->{$patient->id};
	$self->{sd}->{$patient->id} = int($self->statistics()->standard_deviation()*100)/100;
	return $self->{sd}->{$patient->id};
	
}

sub compute_level {
		my ( $self , $patient,$debug) =  @_;
		delete $self->{level}->{$patient->id};
		my $score = $self->cnv_score($patient);
		my $zscore = $self->zscore($patient);

#	if (!$patient->is_multiplex_ok($self->multiplex)){
#		return -2;
#	}
	if (!$self->is_primer_ok) {
		return -1;			
		}
		my $sdp= $self->sd($patient);#int($self->statistics()->standard_deviation()*100)/100;
		
		#my $zscore = 0;
		if ($debug){
			warn "score and sd ".$score." =+> ".$sdp;
			warn $score + $sdp/2 ;
			
		}
			if ($self->id =~ /chrX_41393905/){
			warn "\t".$patient->name." ".$score ." ".$sdp;
			#warn ($score + $sdp/2);
			#die();
		}
		if ($score - $sdp/2 >= 1.4  ){
					return 2;
		}	
		elsif ( ($score + $sdp/2   <= 0.6)  ){
			return 1;
		}
		elsif ($zscore>= 2.2 && ($score- $sdp/2) > 1.3){
			return 2;
		}
				#elsif ($score + $sdp  <= 0.6){
		elsif($zscore<= -2.2 && ($score + $sdp/2)   <= 0.7){
					return 1;
		}
	
		
		
		return 0;
}


sub level {
	my ( $self , $patient,$compute) =  @_;
	#warn $self->id;
		return $self->{level}->{$patient->id} if exists $self->{level}->{$patient->id};
	
#	unless ($compute){
#		$self->{level}->{$patient->id}  = $self->getData($patient,"_level");
#		warn "coucou"." ".$self->{level}->{$patient->id};
#		return $self->{level}->{$patient->id}  if defined $self->{level}->{$patient->id};
#		
#	}
#	
#	
#	delete $self->{level}->{$patient->id};
	$self->{level}->{$patient->id}  = $self->compute_level($patient);
	return $self->{level}->{$patient->id};
	
}

sub minimum {
	my ($self,$patient) = @_;
	return $patient->minimum($self);
	
}

sub return_raw_coverage_obj{
	my ($self,$p) = @_;
	return $self->getGenes()->[0]->get_coverage($p);
}
1;