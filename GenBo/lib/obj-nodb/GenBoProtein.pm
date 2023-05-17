package GenBoProtein;

use strict;
use Vcf;
use Moo;

use Data::Dumper;
use Config::Std;
extends "GenBoGenomic";


has isProtein => (
	is		=> 'ro',
	default	=> 1,
);

has sequence => (
	is		=> 'ro',
);
	
has external_name => (
	is		=> 'ro',
);

has genomic_span => (
	is		=> 'ro',
	reader	=> 'getGenomicSpan',
	required	=> 1,
	
);


 has gene => (
	is		=> 'ro',
	required=> 1,
);
 has transcript => (
	is		=> 'ro',
	required =>1
	
);
 has gene_kyoto_id => (
	is		=> 'ro',
	required=> 1,
);



sub getSequence {
	my ($self,$start,$end) = @_;
	my $length =1;
	$length = ($end-$start) +1 if $end;
	die() if $length < 0;
	return '' if ($self->length() < ($start-1 + $length));
	return substr $self->sequence,$start-1, $length;
}
			
has length => (
	is		=> 'ro',
	reader	=> 'length',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $length = length($self->sequence());
		return $length;
	},
);
1;