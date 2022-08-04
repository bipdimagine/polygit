package GenBoCoverageTabix;
use strict;
use Moose;
extends "GenBoCoverage";
use Data::Dumper;
use Config::Std;
use Storable qw(store retrieve freeze thaw fd_retrieve dclone);
use Clone 'clone';
use List::Util qw(first max maxstr min minstr reduce shuffle sum0);

has chromosome => (
	is		=> 'rw',
	required =>1,
);

has patient => (
	is		=> 'rw',
	required =>1,
);


has array => (
	is		=> 'ro',
		lazy=>1,
	default => sub {
		my $self = shift;
	my $tabix = $self->patient->tabix_coverage;
	my $len = abs($self->start - $self->end) +1;
	my @v = ((0) x $len);
	 my $res = $tabix->query_full($self->chromosome_name,$self->start-1,$self->end);
	  return \@v unless  $res;
		 my @data;
		
		 while(my $line = $res->next){
    		
				my($a,$p,$c) = split(" ",$line);
				confess() if $a ne $self->chromosome_name;
				my $pv = $p - $self->start;
				 $v[$pv] = $c;
			}
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

		my $tabix = $self->patient->tabix_coverage;
	my $samtools = $self->patient->buffer->samtools;
	$self->array->[$pos_array] = 0;
   my $res = $tabix->query_full($self->chromosome_name,$pos_genomic,$pos_genomic);

	#	warn "$samtools depth $bam -r $region";
		while(my $line = $res->next){
    		
				my($a,$p,$c) = split(" ",$line);
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