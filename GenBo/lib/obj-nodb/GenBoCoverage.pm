package GenBoCoverage;

use strict;
use Moose;
use Data::Dumper;
use Config::Std;
use Storable qw(store retrieve freeze thaw fd_retrieve dclone);
use Clone 'clone';
use List::Util qw(first max maxstr min minstr reduce shuffle sum0 );

 
has start => (
	is		=> 'ro',
	required =>1,
);


has end => (
	is		=> 'ro',
	required =>1,
);


has array => (
	is		=> 'rw',
	
	#required =>1,

);



has chromosome_name =>(
is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->chromosome->fasta_name;
		return $self->chromosome->ucsc_name if $self->patient->is_ucsc_nomenclature();
		return  $self->chromosome->name();
	},
);

has len 	=> (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $len = abs($self->end-$self->start) +1;
		
		confess($len." ".scalar(@{$self->array})) if $len != scalar(@{$self->array});
		 return $len;
	},
);

has mean =>(
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		 return 0 if $self->sum ==0;
		 return  $self->sum/$self->len;
	},
);

has minimum =>(
is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return(min(@{$self->array}));
	},
);


has sum =>(
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		 return sum0 @{$self->array()};
	},
);


sub coverage {
	my ($self,$start,$end) = @_;
	#return $self->{gc}->{$start.$end} if exists $self->{gc}->{$start.$end};
	return [] if $start == -1;
	my @t;
	my $nb = 0 ;
	my $sum;
	my $min = 999999999999;
	my $array = $self->array;
	if ($start>=$self->start && $end <= $self->end){
		my $a = $self->translate($start);
		my $b = $self->translate($end);
		confess() if $a < 0;
		my @aa = (@$array[$a .. $b]);
		my $hash;
		 $hash->{array} = \@aa;
		 $hash->{sum} = sum0 (@aa);
		 $hash->{nb} = scalar(@aa);
		  $hash->{min} = min(@aa);
		  $hash->{mean} =$hash->{sum}/ $hash->{nb};
		  return $hash;
	}
	else {
		warn ($start-$self->start);
		confess($start."-".$self->start." end : ".$end."-> ".$self->end);
	}
	for (my $i=$start;$i<$end;$i++){
		my $x = $self->depth($i);
		$sum+= $x;
		$min = $x if $x< $min;
		$nb ++;
		push(@t,$x);
		#push(@t,$self->array->[$i]);
	}
	my $mean =0;
	 $mean = $sum/$nb if $nb>0;
	 my $hash;
	 $hash->{array} = \@t;
	
	  $hash->{min} = min (@t);
	  $hash->{nb} = $nb;
	  $hash->{sum} = $sum;
	   $hash->{mean} = $hash->{sum} / $nb;
	#   $self->{gc}->{$start.$end} = $hash;
	return  $hash;
}

sub sub_array {
	my ($self,$start,$end) = @_;
	return [] if $start == -1;
	my $array = $self->array;
	if ($start>=$self->start && $end <= $self->end){
			my $a = $self->translate($start);
		my $b = $self->translate($end);
		confess() if $a < 0;
		my @aa = (@$array[$a .. $b]);
		return \@aa;
	}
	confess($start."-".$end." ".$self->start."-".$self->end);
	
}

sub translate {
	my ($self,$pos) =@_;
	 my $vstart = ($pos - $self->start) ;
	 return $vstart;
}

sub depth {
	my ($self,$pos) =@_;
	 my $vstart = ($pos - $self->start) ;
	 confess("vstart: ".$vstart." pos :".$pos." ".$self->start." ".$self->end()." ".$self->len) if ($vstart <0);
	$vstart = 0 if $vstart < 0;
	confess() unless  defined $self->array->[$vstart];
	#return $self->array->[$vstart];
	#warn   $self->array->[$vstart];
	$self->add_value($vstart,$pos) unless  defined $self->array->[$vstart];
	confess("vstart: ".$vstart." pos :".$pos." ".$self->start." ".$self->end()." ".$self->len) if ($vstart > $self->len);
	#confess($pos." ".$self->start." ".$self->len) unless  defined $self->array->[$vstart];
	return $self->array->[$vstart]+0;
}

sub coverage_intspan {
	my ($self,$intspan) = @_;
	my @list = $intspan->as_array();
	my $vstart =  $self->start;
	confess() unless  defined $self->start;
	my @t;
	my $nb = 0 ;
	my $sum;
	my $min = 999999999999;
	my $array = $self->array;
	foreach my $pos (@list){
		my $x = $self->depth($pos);
		$sum+= $x;
		$min = $x if $x< $min;
		$nb ++;
		push(@t,$x);
	}
	 my $mean =0;
	 $mean = $sum/$nb if $nb>0;
	 my $hash;
	 $hash->{array} = \@t;
	 $hash->{mean} = $mean;
	  $hash->{min} = $min;
	  $hash->{nb} = $nb;
	  $hash->{sum} = $sum;
	return $hash;
	
}


1;