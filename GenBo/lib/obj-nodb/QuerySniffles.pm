package QuerySniffles;
use strict;
use Moo;
use Data::Dumper;
use Carp;
use Clone 'clone';
use IPC::Open2;
use List::MoreUtils qw( uniq );
use Compress::Snappy;
use Storable qw/thaw freeze/;
use Align::NW;
use Bio::DB::HTS::Tabix;
use Bio::DB::HTS::VCF;
extends "QueryPbsv";


sub add_DP_AD {
	my ( $self, $patient_id, $hash, $x ) = @_;
	warn $x->{gt}->{DV};
	$hash->{annex}->{$patient_id}->{nb_all_mut} = $x->{gt}->{DV};
	$hash->{annex}->{$patient_id}->{nb_all_ref} = $x->{gt}->{DR};;
	$hash->{annex}->{$patient_id}->{dp} = $x->{gt}->{DV}+$x->{gt}->{DR};
}

sub add_SR_PR {
	my ( $self, $patient_id, $hash, $x ) = @_;
	my $a =  $x->{gt}->{DR};
	my $b =  $x->{gt}->{DR};
	$hash->{annex}->{$patient_id}->{pr} = "-1,-1";
	$hash->{annex}->{$patient_id}->{sr} = "-1,-1";
	$hash->{annex}->{$patient_id}->{pr} = 
	$hash->{annex}->{$patient_id}->{sr} ="$a,$b";
}
1;