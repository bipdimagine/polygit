package GenBoIndel;

use strict;
use Moose;
use MooseX::Method::Signatures;
use Data::Dumper;
use Config::Std;
use GenBoCapture;
use Position;
extends "GenBoVariation";

has type_object => (
	is		=> 'ro',
	default	=> "indels_object",
);

##
##
##
##
##has id => (
##	is		=> 'ro',
##	isa		=> 'Str',
##	required=> 1,
##	reader	=> 'id',
##);
##
##has project => (
##	is		=> 'ro',
##	isa		=> 'GenBoProject',
##	required=> 1,
##	weak_ref=> 1,
##	reader	=> 'getProject',
##);
##
##has chromosome => (
##	is		=> 'ro',
##	isa		=> 'Str',
##	lazy	=> 1,
##	reader	=> 'chromosome',
##	default => sub {
##		my $self = shift;
##		my $id = $self->id();
##		my @lFields = split("_", $id);
##		return $lFields[0];
##	},
##); 
##
##has position => (
##	is		=> 'ro',
##	isa		=> 'Position',
##	lazy	=> 1,
##	default => sub {
##		my $self = shift;
##		my $id = $self->id();
##		my @lFields = split("_", $id);
##		my $start = int($lFields[1]);
##		my $end = int($lFields[2]);
##		my $chrom = $self->chromosome();
##		my $posId = $chrom . '_' . $start . '_' . $end;
##		my $project = $self->getProject();
##		my $posObject = $project->flushObject('positions', $posId);
##		my $i = $start;
##		while ($i <= $end) {
##			push(@{$$project{'positions'}->{$chrom}->{$i}->{'indels'}}, $self);
##			$i++;
##		}
##		return $posObject;
##	},
##);
##
##has ref_allele => (
##	is		=> 'ro',
##	isa		=> 'Str',
##	lazy	=> 1,
##	reader	=> 'ref_allele',
##	default => sub {
##		my $self = shift;
##		my $id = $self->id();
##		my @lFields = split("_", $id);
##		return $lFields[2];
##	},
##);
##
##has var_allele => (
##	is		=> 'ro',
##	isa		=> 'Str',
##	lazy	=> 1,
##	reader	=> 'var_allele',
##	default => sub {
##		my $self = shift;
##		my $id = $self->id();
##		my @lFields = split("_", $id);
##		return $lFields[3];
##	},
##);
##
##has patients => (
##	is		=> 'rw',
##	isa		=> 'ArrayRef[Str]',
##	reader	=> 'patients',
##	lazy	=> 1,
##	default => sub { [] },
##);
#
1;