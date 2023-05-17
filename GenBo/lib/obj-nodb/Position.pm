=head1 NAME

Position : Specific GenBo for position 

=head1 SYNOPSIS

my $position = $buffer->new Position({-start => $start, -end => $end, -strand => $strand});

=head1 DESCRIPTION

Position provides a set of functions to get basics informations on positions and calculate positions between GenBo
The Position object is create between 2 GenBos (one passed to the method and the other one which the method is call)

=head1 METHODS

=cut

package Position;
use strict;
use Moo;
use Data::Dumper;
use Config::Std;


has start => (
	is		=> 'rw',
	required=> 1,
);

has end => (
	is		=> 'rw',
	required=> 1,
);

has strand => (
	is		=> 'rw',
	required=> 1,
	#default =>1,
);


has intspan => (
	is		=> 'ro',
	reader	=> 'getGenomicSpan',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $intSpan = Set::IntSpan::Fast->new();
    	$intSpan->add_range($self->start(), $self->end());
    	return $intSpan;
	},
);


=head2 toString
	Title   : toString
 	Usage   : $position->toString();
 	Function: Format the Position object in readable format
 	Returns : A string correponding to the Position
 	Args	: None
	Note    : [start-end](strand)
=cut

sub toString {
	my $self = shift;
	return "[".$self->start."-".$self->end."]" ."(".$self->length.") ". $self->strand();
}

=head2 length
	Title   : length
 	Usage   : $position->length();
 	Function: Get the length of the Position object
 	Returns : The length of the object (int)
 	Args	: None
	Note    : 
=cut

sub length {
	my $self = shift;
	return ($self->end - $self->start)+1;
}

1;
