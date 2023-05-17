package GenBoNoSqlIntervalTree;
use strict;
use warnings;
use DBD::SQLite;
use Moo;

use Set::IntervalTree;
use Data::Dumper;
use Set::IntSpan::Fast::XS;
extends "GenBoNoSql";
has extension =>(
	is		=> 'rw',
default => sub {
		return "tree";
		
}
);

sub get {
	my ($self,$key1,$key2,$debug) = @_;
	my $data = $self->get_lite($key1,$key2);
	 my $tree = Set::IntervalTree->new;
	 foreach my $v (@{$data}) {
	 	#warn  $v->[0]  if $v->[1] == -5000;
	 	next if $v->[1] == -5000;
	 	$tree->insert(@$v);
	 }
	 
	return $tree;
}
sub get_intspan {
	my ($self,$key1,$key2,$debug) = @_;
	my $data = $self->get_data($key1,$key2);
	my $span3 = Set::IntSpan::Fast::XS->new();
	 foreach my $v (@{$data}) {
	 		next if $v->[1] == -5000;
	 	#warn Dumper $v if $v->[0] eq "ENSG00000203478_17";
	 	$span3->add_range($v->[1],$v->[2]-1);
	 	
	 }
	return $span3;
}

sub get_data  {
	my ($self,$key1,$key2,$debug) = @_;
	my $data = $self->get_lite($key1,$key2);
	
#	warn Dumper $data; 
	return $data;
}
1;
