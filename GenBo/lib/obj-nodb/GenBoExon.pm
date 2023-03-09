package GenBoExon;

use strict;
use Moose;
use MooseX::Method::Signatures;
use Data::Dumper;
use Config::Std;
use Storable qw(store retrieve freeze thaw);
use POSIX qw(ceil);
use List::Util qw(min sum);
extends "GenBoGenomic";

has isExon => (
	is		=> 'ro',
	default	=> 1,
);
has transcripts_object => (
	is		=> 'rw',
	required	=> 1,
);


has coverage_id =>(
	is		=> 'ro',
	lazy => 1,
	default => sub {
			my $self = shift;
			return join("-",$self->getChromosome->name(),$self->getTranscript()->name(),$self->start,$self->end);
	},
);

has utr => (
	is		=> 'ro',
	required =>1,
);
has intspan_no_utr => (
	is		=> 'ro',
	lazy  =>1,
	default => sub {
			my $self = shift;
			
			return  $self->getGenomicSpan()->diff($self->utr); 
	},
);


sub init_utr {
	my $self = shift;
	 my $string = $self->intspan_no_utr->as_string;
	
	 my($st,$en) = split("-",$string);
	 die() if $string =~/,/;
	 $self->{start_utr} = $st;
	 
	  $self->{end_utr} = $en;
	  $self->{end_utr} = $st unless $en;
	  $en++;
}

has start_utr => (
	is		=> 'ro',
	lazy =>1,
	default => sub {
			my $self = shift;
			return 0  if $self->intspan_no_utr->is_empty;
			$self->init_utr();
			return  $self->{start_utr};
	},
	
);
has is_start_utr => (
	is		=> 'ro',
	lazy =>1,
	default => sub {
			my $self = shift;
			return $self->start_utr ne $self->start;
	},
	
);

has end_utr => (
	is		=> 'ro',
	lazy =>1,
	default => sub {
			my $self = shift;
			return  0 if $self->intspan_no_utr->is_empty;
			$self->init_utr();
			return $self->{end_utr};
			#my $string = $self->intspan_no_utr->as_string();
			my @t = $self->intspan_no_utr->as_array();
			warn $t[-1];
			return  $t[-1];
			#return  $self->getGenomicSpan()->diff($self->utr); 
	},
);

has is_end_utr => (
	is		=> 'ro',
	lazy =>1,
	default => sub {
			my $self = shift;
			
				return $self->end_utr  ne $self->end;
		
			#return  $self->getGenomicSpan()->diff($self->utr); 
	},
	
);

has is_noncoding =>(
	is		=> 'ro',

	lazy =>1,
	default => sub {
		my $self = shift;
		my $intspan_noutr =  $self->getGenomicSpan()->diff($self->utr); 
		return 1 if $self->intspan_no_utr->is_empty;
		return ;
	}
);

method return_start_end ( Int :$padding){
	my $sstart = $self->start;
	my $send = $self->end;
	$sstart -= $padding; #unless $self->$self->start ;
	$send += $padding; #unless $self->$self->start ;
	$sstart =1 if $sstart <=0;
	#$send = $self->getChromosome->length -1 if $send > $self->getChromosome->length ;
	return {start=>$sstart,end=>$send};
}

method return_start_end_no_utr ( Int :$padding){
	
	return $self->{start_end_no_utr}->{$padding} if exists   $self->{start_no_utr}->{$padding};
	return undef  if $self->intspan_no_utr->is_empty;
		return (0,$self->intspan_no_utr) if $self->intspan_no_utr->is_empty;
	my $sstart = $self->start_utr;
	my $send = $self->end_utr;
	$sstart -= $padding unless $self->is_start_utr ;
	$send += $padding unless $self->is_end_utr ;
	$sstart =1 if $sstart <=0;
	#$send = $self->getChromosome->length -1 if $send > $self->getChromosome->length ;
	$self->{start_end_no_utr}->{$padding}->{start} = $sstart;
	$self->{start_end_no_utr}->{$padding}->{end} = $send;
	return $self->{start_end_no_utr}->{$padding};
}

method cached_statistic_coverage_coding (GenBoPatient :$patient, Int :$padding, Int :$limit, Int :$utr ){
		my $no =  $self->project->noSqlCoverage();

		my $t ={};
		eval {
			
		 $t  = $no->get($patient->name,$self->id."_".$padding."_".$utr);
		};
		return $t if defined $t;
		
		my  ($mean,$intspan,$min)  = $self->statistic_coverage_coding(patient=>$patient,padding=>$padding,limit=>$limit,utr=>$utr);
		my $h =  {mean=>$mean,intspan=>$intspan,min=>$min};
				
		eval{
		#	warn "save";
			$no->set($patient->name,$self->id."_".$padding."_".$utr,$h);
		};
	#	die();
		return $h;
		
}
method computed_statistic_coverage_coding (GenBoPatient :$patient, Int :$padding, Int :$limit, Int :$utr ){

		my  ($mean,$intspan,$min)  = $self->statistic_coverage_coding(patient=>$patient,padding=>$padding,limit=>$limit,utr=>$utr);
		my $h =  {mean=>$mean,intspan=>$intspan,min=>$min};
		
		return $h;
		
}

method statistic_coverage_coding (GenBoPatient :$patient, Int :$padding, Int :$limit, Int :$utr ){
	return $self->mean_intspan_coverage($patient,$padding,$limit) if $utr > 0;
	return $self->mean_intspan_coverage_coding($patient,$padding,$limit);
	
}
sub mean_intspan_coverage_coding{
 my ($self,$patient, $padding,$limit ) =@_;
	return @{$self->{coding_stats}->{$patient->id}->{$padding}} if exists $self->{coding_stats}->{$patient->id}->{$padding};
	return (0,$self->intspan_no_utr,-1) if $self->intspan_no_utr->is_empty;

			my $intspan = Set::IntSpan::Fast::XS->new();
			
	
	my $pos = $self->return_start_end_no_utr(padding=>$padding);

	my $sstart = $pos->{start};
	my $send = $pos->{end} ;
	
	return @{$self->{coding_stats}->{$patient->id}->{$sstart.$send}} if exists $self->{coding_stats}->{$patient->id}->{$sstart.$send};
	
	# $self->getTranscript->getGene->get_coverage($patient) unless exists  $self->getTranscript->{coverage_obj}->{$patient->id};
	my $res2 = $patient->depth($self->getChromosome->name,$sstart,$send);
	my $mean = sum(@$res2)/scalar(@$res2);
	my $min = min(@$res2);
	my $array = [$mean,$intspan,$min];

	#my $res2 =  $self->coverage_object($patient)->coverage($sstart,$send);#$self->getTranscript()->getGene->get_coverage($patient)->coverage($sstart,$send);

	$self->{coding_stats}->{$patient->id}->{$padding} = $array;
	$self->{coding_stats}->{$patient->id}->{$sstart.$send}  = $array;
	return @{$self->{coding_stats}->{$patient->id}->{$padding}};
}
sub  mean_intspan_coverage {
	
 my ($self,$patient,$padding,$limit ) =@_;
	
	return @{$self->{stats}->{$patient->id}->{$padding}} if exists $self->{stats}->{$patient->id}->{$padding};
	my $intspan = Set::IntSpan::Fast::XS->new();
	

	#my $pos = $self->return_start_end_no_utr(padding=>$padding);
	my $sstart = $self->{start};
	my $send = $self->{end};
	return @{$self->{coding_stats}->{$patient->id}->{$sstart.$send}} if exists $self->{coding_stats}->{$patient->id}->{$sstart.$send};
	 #$self->getTranscript->getGene->get_coverage($patient) unless exists  $self->getTranscript->{coverage_obj}->{$patient->id};
	
	my $res2 = $patient->depth($self->getChromosome->name,$sstart,$send);
	my $mean = sum(@$res2)/scalar(@$res2);
	my $min = min(@$res2);
	#$self->coverage_object($patient)->coverage($sstart,$send);
	my $array = [$mean,$intspan,$min];
	#my $res2 = $self->getTranscript->getGene->get_coverage($patient)->coverage($sstart,$send);
	$self->{stats}->{$patient->id}->{$padding} = $array;
	$self->{coding_stats}->{$patient->id}->{$sstart.$send}  = $array;
	
	return @{$self->{stats}->{$patient->id}->{$padding} };#($res2->{mean},$intspan,$res2->{min});
	
}

sub coverage_object {
	my ($self,$patient) = @_;
	return $self->{gc1}->{$patient->id} if exists $self->{gc1}->{$patient->id};
	my $titi = $self->getTranscript()->getGene->get_coverage($patient)->coverage($self->start-100,$self->end+100);
	
	my $toto = GenBoCoverageSamtools->new(chromosome=>$self->getChromosome, patient=>$patient, start=>$self->start-100, end=>$self->end+100);
	$toto->{array} = $titi->{array};
	$self->{gc1}->{$patient->id} = $toto; 
	return $self->{gc1}->{$patient->id};
	
}
sub compute_mean_min_intstpan{
	my ($self,$limit,$start,$end,$patient) = @_;
	my $res = $self->getTranscript()->get_coverage($patient)->get_range($start,$end);
	die();
	my $sum =0;
	my $intspan = Set::IntSpan::Fast::XS->new();
	my $intspan2 = Set::IntSpan::Fast::XS->new($start."-".$end);
	
	my $min = 9999999;
	my $nb;
	foreach my $a (@$res){
		$min = $a->[2] if $min > $a->[2];
		$intspan2->remove_range($a->[0],$a->[1]) if $a->[2] >= $limit;
		for (my$i=$a->[0];$i<=$a->[1];$i++){
			$nb ++;
			$sum += $a->[2];
		}
		
	}
	my $mean = int($sum / abs($end-$start+1));
	return ($mean,$intspan2,$min);
	
}


sub return_raw_coverage_obj{
	my ($self,$p) = @_;
	return $self->getTranscript()->getgene()->get_coverage($p);
}


sub getChromosome {
        my ($self) =@_;
        return $self->getTranscript->getChromosome();
}



1;