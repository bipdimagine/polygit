package GenBoNoSqlRocksVector;
use Moo; 
use strict;
use warnings;
use Data::Dumper;
use JSON::XS;
use Carp;
use Storable qw/thaw freeze/;
extends "GenBoNoSqlRocks";

has has_config =>(
	is		=> 'ro',
default => sub {
		return 1;
}
);

has chromosome =>(
	is		=> 'ro',
	required =>1,
);
has size =>(
	is		=> 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		unless ($self->mode eq "c"){ 
			
			$self->load_config() ;
			return $self->{size};
		}
		return 0;
	
	}
);


sub encode {
		my ($self,$code) = @_;
		#return $code->to_Hex;
		return freeze($code);
}

sub decode {
	my ($self,$code) = @_;
	return undef unless $code;
	#return Bit::Vector->new_Hex($self->size,$code);
	my $z =  thaw($code);
	return $z;
}

sub load_config {
	my ($self) = @_;
		return if exists $self->{load};
	$self->{load} = 1;
	open(my $fh ,$self->json_file) or confess("can t open file ".$self->json_file);
	my $json_str = do { local $/; <$fh> };
	close($fh);
	my $h = decode_json($json_str);
	if($self->mode ne "c"){
		$self->{size} = delete $h->{size};
		$self->{date} = delete $h->{time};
		$self->{separator_category} = delete $h->{separator_category};
		$self->{vector_type_patient} = delete $h->{vector_type_patient};
		$self->{vector_type_transmission} = delete $h->{vector_type_transmission};
		$self->{vector_type_chromosome} = delete $h->{vector_type_chromosome};
		
		
	}
	return 1;
}
sub write_config {
	my ($self) = @_;
	open(my $fh,">".$self->json_file) or die("can t open file".$self->json_file);
	my $h;
	$h->{size} = $self->size;
	$h->{date} = time;
	$h->{separator_category} = $self->{separator_category};
	$h->{vector_type_patient} = $self->vector_type_patient;
	$h->{vector_type_transmission} = $self->vector_type_transmission;
	$h->{vector_type_chromosome} = $self->vector_type_chromosome;
	print $fh encode_json($h);
	close ($fh);
}

has separator_category =>(
	is		=> 'ro',
	default => sub {
		return "+";
	}
);



has separator_key =>(
	is		=> 'ro',
	default => sub {
		return "*";
	}
);


############################
#VECTOR transmission
#############################
has vector_type_transmission =>(
	is		=> 'rw',
	lazzy=>1,
	default => sub {
		my ($self) = @_;
		if ($self->mode eq "c"){
			return {};
		}
		else {
			$self->load_config;
			return $self->{vector_type_transmission};
		}
	}
);
sub id_vector_transmission {
	my ($self,$patient,$type) = @_;
	return $self->chromosome.$self->separator_category.$patient->id.$self->separator_category.$patient->getFamily->name.$self->separator_category.$type;
}
sub put_batch_vector_transmission {
	my ($self,$patient,$type,$v) =@ _;
	my $id = $self->id_vector_transmission($patient,$type);
	$self->vector_type_transmission->{$type} = 1;
	$self->put_batch($id,$v);
}
sub get_vector_transmission {
	my ($self,$patient,$type) =@ _;
	confess($type.":".join(";",keys %{$self->vector_type_transmission}) ) unless exists $self->vector_type_transmission->{$type};
	my $h = $self->cache_memory_patient($patient);
	my $id = $self->id_vector_transmission($patient,$type);
	die() unless $h;
	return $self->decode($h->{$id});
}

############################
#VECTOR CHROMOSOME
#############################
has vector_type_chromosome =>(
	is		=> 'rw',
	lazzy=>1,
	default => sub {
		my ($self) = @_;
		if ($self->mode eq "c"){
			return {};
		}
		else {
			$self->load_config;
			return $self->{vector_type_chromosome};
		}
	}
);
sub id_global_chromosome {
	my($self,$type) = @_;
	return $self->chromosome.$self->separator_category.$type;
}
sub put_batch_vector_chromosome {
	my ($self,$type,$v) =@ _;
	$self->vector_type_chromosome->{$type} = 1;
	$self->put_batch($self->id_global_chromosome($type),$v);
}
sub prepare_vector {
	my ($self,$types) = @_;
	my $ids = [];
	foreach my $t (@$types) {
		push(@$ids,$self->id_global_chromosome($t));
	}
	 $self->prepare($ids);
}


sub get_vector_chromosome {
	my ($self,$type) =@_;
	if ($type eq "all"){
		return Bit::Vector->new($self->size);
	}
	die($type."\n".Dumper $self->vector_type_chromosome) unless exists $self->vector_type_chromosome->{$type};
	my $v;
	my $id = $self->id_global_chromosome($type);
	if (exists $self->{buffer}->{$id} ){
		$v =  $self->{buffer}->{$id};
	}
	else {
		$v = $self->rocks->get($id);
	}
	confess($type." ".$self->chromosome) unless $v;
	return $self->decode($v);
}


sub get {
	my ($self,$key) = @_;
	confess();
}

sub get_vector_gene {
	my ($self,$key) = @_;
	return $self->SUPER::get($key);
}
sub count {
		my ($self, $bitvector) = @_;
		#if ($bitvector->is_empty());
		return ($bitvector->to_Bin() =~ tr/1//);
}
############################
#VECTOR PATIENT
#############################

sub id_global_patient {
	my($self,$patient,$type) = @_;
	return $self->chromosome.$self->separator_category.$patient->id.$self->separator_category.$type;
}

sub put_vector_patient_batch {
	my ($self,$patient,$type,$v) =@ _;
	$self->vector_type_patient->{$type} = 1;
	$self->put_batch($self->id_global_patient($patient,$type),$v);
}

has vector_type_patient =>(
	is		=> 'rw',
	lazzy=>1,
	default => sub {
			my ($self) = @_;
		if ($self->mode eq "c"){
			return {};
		}
		else {
			$self->load_config;
			return $self->{vector_type_patient};
		}
	}
);
sub cache_memory_patient{
	my ($self,$patient) = @_;
	confess() if $patient eq "55661";
	return $self->{hash}->{$patient->id} if exists  $self->{hash}->{$patient->id};
	my $search =$self->id_global_patient($patient,"");
	chop($search);
	my $x1 = $patient->id;
	my $iter = $self->rocks->new_iterator->seek($search);
	while (my ($key, $value) = $iter->each) {
		last if $key !~ /$x1/;
			$self->{hash}->{$patient->id}->{$key} = $value;
	}
	return $self->{hash}->{$patient->id};
}

sub get_vector_patient {
	my ($self,$patient,$search) = @_;
	die("problem ".$search) unless exists $self->vector_type_patient->{$search};
	my $h = $self->cache_memory_patient($patient);
	warn Dumper keys %$h  unless exists $h->{$self->id_global_patient($patient,$search)};
	confess($search." ".$self->id_global_patient($patient,$search)) unless exists $h->{$self->id_global_patient($patient,$search)};
	return $self->decode($h->{$self->id_global_patient($patient,$search)});
	
}


1;