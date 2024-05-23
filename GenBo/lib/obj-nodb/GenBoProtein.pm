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

#has sequence => (
#	is		=> 'ro',
#);
sub sequence {
	my ($self) = @_;
	return $self->{sequence_new} if exists $self->{sequence_new};
	foreach my $vg  (("HG38","HG19")) {
	my $index = $self->project->fastaProteinIndex($vg);
	if($index){
			my ($seq1, $length) = $index->get_sequence($self->getTranscript->name);
			if ($seq1){
				$self->{sequence_new} = $seq1;
				return $self->{sequence_new};
			} 
		}	
	 	
	}
	$self->{sequence_new} = $self->{sequence};
	return $self->{sequence_new};
}	
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
	return '' unless $self->length();
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