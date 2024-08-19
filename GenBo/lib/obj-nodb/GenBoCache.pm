package GenBoCache;
use strict;
use Storable qw(retrieve);
use Moo;
use Bit::Vector;
use Bit::Vector::Overload;
use Set::IntSpan::Fast;
use Data::Dumper;


has project => (
	is		=> 'ro',
	#reader	=> 'getProject',
);

has chromosome => (
	is		=> 'ro',
);

has variants => (
	is		=> 'rw',
);

has statistics => (
	is		=> 'rw',
	lazy    => 1,
	default => sub { return {};  }
);

has hash_convert => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		
	}
);

my $hash_convert;
$hash_convert->{'a'} = '0100';
$hash_convert->{'b'} = '01110';
$hash_convert->{'c'} = '0101';
$hash_convert->{'d'} = '01100';
$hash_convert->{'e'} = '01111';
$hash_convert->{'f'} = '11010';
$hash_convert->{'g'} = '01101';
$hash_convert->{'h'} = '11000';
$hash_convert->{'i'} = '11011';
$hash_convert->{'j'} = '11110';
$hash_convert->{'k'} = '11001';
$hash_convert->{'l'} = '11100';
$hash_convert->{'m'} = '11111';
$hash_convert->{'n'} = '10010';
$hash_convert->{'o'} = '11101';
$hash_convert->{'p'} = '10000';
$hash_convert->{'q'} = '10011';
$hash_convert->{'r'} = '10110';
$hash_convert->{'s'} = '10001';
$hash_convert->{'t'} = '10100';
$hash_convert->{'u'} = '10111';
$hash_convert->{'v'} = '0010';
$hash_convert->{'w'} = '10101';
$hash_convert->{'x'} = '0000';
$hash_convert->{'y'} = '0001';
$hash_convert->{'z'} = '0011';
$hash_convert->{'A'} = '0';
$hash_convert->{'B'} = '1';
$hash_convert->{'0100'} = 'a';
$hash_convert->{'01110'} = 'b';
$hash_convert->{'0101'} = 'c';
$hash_convert->{'01100'} = 'd';
$hash_convert->{'01111'} = 'e';
$hash_convert->{'11010'} = 'f';
$hash_convert->{'01101'} = 'g';
$hash_convert->{'11000'} = 'h';
$hash_convert->{'11011'} = 'i';
$hash_convert->{'11110'} = 'j';
$hash_convert->{'11001'} = 'k';
$hash_convert->{'11100'} = 'l';
$hash_convert->{'11111'} = 'm';
$hash_convert->{'10010'} = 'n';
$hash_convert->{'11101'} = 'o';
$hash_convert->{'10000'} = 'p';
$hash_convert->{'10011'} = 'q';
$hash_convert->{'10110'} = 'r';
$hash_convert->{'10001'} = 's';
$hash_convert->{'10100'} = 't';
$hash_convert->{'10111'} = 'u';
$hash_convert->{'0010'} = 'v';
$hash_convert->{'10101'} = 'w';
$hash_convert->{'0000'} = 'x';
$hash_convert->{'0001'} = 'y';
$hash_convert->{'0011'} = 'z';

##### METHODS #####



sub getNewVector {
	my ($self, $size) = @_;
	
	if ($size) { return Bit::Vector->new($size); }
	elsif (ref($self) eq 'GenBoGeneCache') { return Bit::Vector->new($self->len_subpart()); }
	return Bit::Vector->new($self->size_vector());
}

sub countVariants {
	my $self = shift;
	return $self->countThisVariants($self->getVariantsVector());
}

sub countThisVariants {
my ($self, $bitvector) = @_;
    return $bitvector->Norm();
}


sub count {
		my ($self, $bitvector) = @_;
		return $self->countThisVariants($bitvector);
}
sub setIdsBitOff {
	my ($self, $vector, $ids) = @_;
	foreach my $id (sort @$ids) {
		$vector->Bit_Off($id);
	}
}

sub setIdsBitOn {
	my ($self, $vector, $ids) = @_;
	foreach my $id (sort @$ids) {
		$vector->Bit_On($id);
	}
}

sub getIdsBitOn {
	my ($self, $vector) = @_;
	#warn $vector;
	#confess($self->name()) unless defined $vector;
	my @list = $vector->Index_List_Read();
	return \@list;
}

sub getListVarVectorIds {
	my ($self, $var) = @_;
	return $self->getIdsBitOn($var);
}

sub getListVarIds {
	my ($self, $var) = @_;
	my $hVarIds;
	my @lVarIds;
	return \@lVarIds unless (ref($var) eq 'Bit::Vector');
	return \@lVarIds if ($var->is_empty());
	foreach my $id (@{$self->getIdsBitOn( $var )}) {
		my $var_id = $self->chromosome->getVarId($id);
		$hVarIds->{$var_id} = undef;
	}
	@lVarIds = keys %$hVarIds;
	return \@lVarIds;
}

sub transformBitVectorToList {
	my ($self, $bitvector) = @_;
	my @lIds;
	unless (ref($bitvector)) { confess(); }
	my $set = Set::IntSpan::Fast::XS->new($bitvector->to_Enum());
	return [$set->as_array];

}

# set new variants vector
sub setVariantsVector {
	my ($self, $var) = @_;
	$self->variants($var);
}

# get variants vector
sub getVariantsVector {
	my $self = shift;
	return $self->variants();
}


# get variants vector in this category
sub getCategoryVariantsVector {
	my ($self, $cat) = @_;
	if ($self->getVariantsVector->is_empty()) {
		return $self->getNewVector();
	}
	if (exists $self->global_categories->{$cat}) {
		$self->global_categories->{$cat}->Intersection( $self->global_categories->{$cat}, $self->getVariantsVector() );
		return $self->global_categories->{$cat};
	}
	elsif ($cat eq 'ho' or $cat eq 'he') {
		my $vector = $self->getNewVector();
		foreach my $patient (@{$self->getPatients()}) {
			$vector += $self->patients_categories->{$patient->name().'_'.$cat};
		}
		return $vector;
	}
	elsif (exists $self->project->hash_ensembl_annotations->{$cat}) {
		unless (exists $self->categories->{$cat}) { return $self->getNewVector(); }
		$self->categories->{$cat}->Intersection( $self->categories->{$cat}, $self->getVariantsVector() );
		return $self->categories->{$cat};
	}
	return $self->getNewVector();
}

sub is_empty {
	my $self = shift;
	return $self->getVariantsVector($self->chromosome())->is_empty();
}

sub isTrue {
	my ($self,$vector,$index) = @_;
	return $vector->bit_test($index);
}

sub convert_intspan_to_vector {
	my ($self, $intspan, $chr) = @_;
	confess() unless ($chr);
	my @lNewCoord = $intspan->as_array();
	if (@lNewCoord) {
		return Bit::Vector->new_Enum($chr->size_vector(), join(',', @lNewCoord));
	}
	return $chr->getNewVector();
}

sub compress_vector_bin {
	my ($self, $value) = @_;
	my $text1 = $self->compress_bin_to_letters($value);
	my $text2 = $self->compress_letters($text1);
	return $text2;
}

sub compress_bin_to_letters {
	my ($self, $value) = @_;
	my $new_text = '';
	my $local_text = '';
	my @lCar = split('', $value);
    for (my $i = 0; $i < scalar(@lCar); $i++) {
    	$local_text .= $lCar[$i];
    	if (exists $hash_convert->{"$local_text"}) {
    		$new_text .= $hash_convert->{"$local_text"};
    		$local_text = '';
    	}
    }
    #cas reste des 0 ou 1
	my @lCar2 = split('', $local_text);
    for (my $i = 0; $i < scalar(@lCar2); $i++) {
    	if ($lCar2[$i] eq '0')    { $new_text .= 'A'; }
    	elsif ($lCar2[$i] eq '1') { $new_text .= 'B'; }
    }
	return $new_text;
}

sub compress_letters {
	my ($self, $value) = @_;
	my $new_text;
	my @lCar = split('', $value);
	my $previous_car = $lCar[0];
	my $nb = 1;
    for (my $i = 1; $i < scalar(@lCar); $i++) {
    	my $this_car = $lCar[$i];
    	if ($this_car eq $previous_car) { $nb++; }
    	else {
    		$new_text .= $previous_car;
    		if ($nb > 1) { $new_text .= $nb; }
    		$previous_car = $this_car;
    		$nb = 1;
    	}
    }
	$new_text .= $previous_car;
	if ($nb > 1) { $new_text .= $nb; }
	return $new_text;
}

sub decompress_vector_bin {
	my ($self, $value) = @_;
	my $text1 = $self->decompress_letters($value);
	my $text2 = $self->decompress_bin_to_letters($text1);
	return $text2;
}

sub decompress_letters {
	my ($self, $value) = @_;
	my $new_text;
	my @lCar = split('', $value);
	my $previous_car = $lCar[0];
	my $previous_nb = '';
	my $nb = 1;
    for (my $i = 1; $i < scalar(@lCar); $i++) {
    	if ($lCar[$i] =~ /^\d+?$/) {
		    $previous_nb .= $lCar[$i];
		}
    	else {
    		if ($previous_nb eq '') {
    			$new_text .= $previous_car;
    		}
    		else {
    			for (my $j = 0; $j < int($previous_nb); $j++) {
    				$new_text .= $previous_car;
    			}
    		}
    		$previous_car = $lCar[$i];
    		$previous_nb = '';
    	}
    	
    }
	if ($previous_nb eq '') {
		$new_text .= $previous_car;
	}
	else {
		for (my $j = 0; $j < int($previous_nb); $j++) {
			$new_text .= $previous_car;
		}
	}
	return $new_text;
}

sub decompress_bin_to_letters {
	my ($self, $value) = @_;
	my $new_text;
	my @lCar = split('', $value);
    for (my $i = 0; $i < scalar(@lCar); $i++) {
		$new_text .= $hash_convert->{$lCar[$i]};
    }
    return $new_text;
}

1;