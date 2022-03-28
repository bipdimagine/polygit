package GenBoStructuralVariant;

use strict;
use Moose;
use MooseX::Method::Signatures;
use Parallel::ForkManager;
use Data::Dumper;
use Config::Std;
use GenBoCapture;
use GenBoGenomic;
use Position;
use List::Util;
use JSON::XS;
use Compress::Snappy;
 use List::Util qw( max min sum);
use Storable qw/thaw freeze/;
#use bytes;
extends "GenBoGenomic";



has isVariant => (
	is		=> 'ro',
	default	=> 1,
);
has isStructuralVariant => (
	is		=> 'ro',
	default	=> 1,
);

has info => (
	is =>'rw',
	#required => 1,
);


has intspan => (
	is		=> 'ro',
	reader	=> 'getGenomicSpan',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $intSpan = Set::IntSpan::Fast::XS->new();
    	$intSpan->add_range($self->start(), $self->end());
    	return $intSpan;
	},
);

sub vcf_line {
		my ($self, $patient) = @_;
		
		return $self->{vcf}->{$patient->id} ;
}

sub start {
	my ($self,$patient,$method) = @_;
	return $self->{start} unless $patient and $method;
	unless ($method){
		return undef unless exists $self->{vcf}->{$patient->id};
		return min( map{$_->{start}} values %{$self->{vcf}->{$patient->id}});
	}
	return $self->{vcf}->{$patient->id}->{$method}->{start};
	
}

sub end {
	my ($self,$patient,$method) = @_;
	return $self->{end} unless $patient and $method;
	unless ($method){
		return undef unless exists $self->{vcf}->{$patient->id};
		return max( map{$_->{end}} values %{$self->{vcf}->{$patient->id}});
	}
	return $self->{vcf}->{$patient->id}->{$method}->{end};
	
}
sub add_event {
	my ($self, $patient,$method,$hash) = @_;
	die() unless $hash;
	die() unless $method;
	
	$self->{vcf}->{$patient->id}->{$method} = $hash;
	$self->{start} = min ($self->{start},$hash->{start});
	$self->{end} = max ($self->{end},$hash->{end});
	delete $self->{intspan};
  $self->{length} = abs($self->{start}-$self->{end}); 
  $self->{vcf}->{$patient->id} = $hash;
 

}


sub isSameEvent {
	my ($self,$chr,$hash) = @_;
	return if $chr ->name ne  $self->getChromosome->name();
	return if  $hash->{end} < $self->start;
	return if  $hash->{start} > $self->end;
	my $start = $hash->{start};
	my $end = $hash->{end};
	confess() unless $end;
	my $intspan = Set::IntSpan::Fast::XS->new("$start-$end");
	my $l1 = abs($start-$end) +1;
	my $inter1 =  $intspan->intersection($self->getGenomicSpan);
	return if $inter1->is_empty;
	

	my $li = $self->buffer->Intspan_length($inter1);
	
	return  if ($li < 0.75 * $l1) and ($li < 0.75 * $self->length) ; 

	return 1;
	
}



1;
