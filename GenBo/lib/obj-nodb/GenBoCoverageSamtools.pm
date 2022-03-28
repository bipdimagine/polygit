package GenBoCoverageSamtools;
use strict;
use Moose;

use Data::Dumper;
use Config::Std;
use Storable qw(store retrieve freeze thaw fd_retrieve dclone);
use Clone 'clone';
use List::Util qw(first max maxstr min minstr reduce shuffle sum0);
extends "GenBoCoverage";

has chromosome => (
	is		=> 'rw',
	required =>1,
);

has patient => (
	is		=> 'rw',
	required =>1,
);

has raw => (
	is		=> 'rw',
	default=>0,
);

has array => (
	is		=> 'ro',
	lazy=>1,
	default => sub {
		my $self = shift;
		confess() unless  $self->patient;
		my $samtools = $self->patient->buffer->samtools;
		my $len = abs($self->start - $self->end) +1;
		my @v = ((0) x $len);
		my $region = $self->chromosome_name.":".$self->start."-".$self->end;
		my $bam = $self->patient->getBamFile();
		#warn "$samtools depth -Q 1 $bam -r $region -Q 1 | cut -f 2-3";
	#	my @t = `$samtools depth -Q 1 $bam -r $region -Q 1 | cut -f 2-3`;
		my $opt = "-Q 1";
	$opt ="" if $self->raw eq 1;
		my @t = `$samtools depth -d 50000 $opt  $bam -r $region $opt | cut -f 2-3` ;
	
	#	warn "$samtools depth $bam -r $region";
		chomp(@t);
		foreach my $tt (@t ){
			my ($p,$c) = split(" ",$tt);
			my $pv = $p - $self->start;
			next if $pv <0;
			$v[$pv] = $c;
			last if $pv > $self->end;
		
		}
			#my $span = Set::IntSpan::Fast::XS->new();
		#$self->span->add_range($self->start,$self->end);
		
		#return $$span;
		return \@v;
	#	my $gc  = GenBoCoverage->new(start=>$start,end=>$end,array=>\@v);
	},

);


sub add_value {
	my ($self,$pos_array,$pos_genomic) = @_;
	#confess("coucou  genomic:$pos_genomic array:$pos_array length:". $self->len." ");
	die();
	my $start = $pos_genomic;
	my $end = $pos_genomic;
	my $before;
	 if ($pos_genomic < $self->start()){
	 	$before =1;
	 	$start=$pos_genomic-10;
	 }
	 
	 
#	$start=$pos_genomic-10  if $pos_genomic < $self->start();
	$end=$pos_genomic+10  if $pos_genomic > $self->end();

	
	my $samtools = $self->patient->buffer->samtools;
	$self->array->[$pos_array] = 0;
	my $region = $self->chromosome_name.":".($pos_genomic-1)."-".($pos_genomic+1);
		my $bam = $self->patient->getBamFile();
		my @t = `$samtools depth -d 50000 -Q 1 $bam -r $region | cut -f 2-3`;
	#	warn "$samtools depth $bam -r $region";
		chomp(@t);
		foreach my $tt (@t ){
			my ($p,$c) = split(" ",$tt);
			
			next if $p ne $pos_genomic;
			my $pv = $p - $self->start;
			$self->array->[$pos_array] = $c;
			last;
		}
		$self->{start}=$pos_genomic  if $pos_genomic < $self->start();
		$self->{end}=$pos_genomic  if $pos_genomic > $self->end();
		$self->{len} = scalar(@{$self->{array}});
}

1;
