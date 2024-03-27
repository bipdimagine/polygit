package GenBoNoSqlRocksVariation;
use Moo; 
use strict;
use warnings;
use Data::Dumper;
use JSON::XS;
use Set::IntSpan::Fast::XS;
extends "GenBoNoSqlRocks";
 
has has_config =>(
	is		=> 'ro',
default => sub {
		return 1;
}
);
 
has index =>(
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my ($self) = @_;
		return 0 if $self->mode eq "c";
		$self->load_config();
		$self->{index} = 0 unless $self->{index};
		return $self->{index};
}
); 

sub ranges{
	my ($self,$N) = @_;
	
	my $M = $self->index - 1;
	if ($M == -1){
		return [];
	} 
	if ($M < 1000){
		return([[0,$M]]);
	}
	my( $step, $s, @ranges ) = ( int $M / $N, 0 );
	push( @ranges, [ $s, $s+$step -1] ), $s+=$step for 1 ..$N;
	$ranges[ -1][1] = $M;
	return \@ranges;
	
}

sub size {
	my ($self) = @_;
	return $self->index;
}


sub raw_put_batch_variation {
	my ($self,$key,$value,$debug) = @_;
	my $index = $self->index;
	$self->put_batch_raw($key,$index);
	$self->put_batch_raw($index,$value);
	$self->index($index+1);
	return $index;	
}

sub put_batch_variation {
	my ($self,$key,$value,$debug) = @_;
	my $index = $self->index;
	$value->{vector_id} = $index;
	$self->put_batch_raw($key,$index);
	$self->put_batch_raw($value->id,$index);
	$self->put_batch($index,$value);
	$self->index($index+1);
	return $index;	
}

sub put_variation {
	my ($self,$key,$value,$debug) = @_;
	my $index = $self->index;
	#$value->{vector_id} = $index;
	$self->rocks->put($key,$index);
	
	$self->rocks->put($index,$self->encode($value));
	$self->index($index+1);
	return 1;	
}


sub get_variation  {
	my ($self,$key) = @_;
	warn "\t".$key;
	$key =  $self->get_raw($key) if ($key =~ /_/);
	warn "\t".$key;;
	return $self->get_index($key);
}
sub get_varid {
	my ($self,$key) = @_;
	return $key;
}

sub get_index {
	my ($self,$key) = @_;
	my $o = $self->get($key);
	return undef unless $o;
	$o->{vector_id} = $key;
	$o->{index_lmdb} = $key;
	return $o;
}
sub write_config {
	my ($self) = @_;
	open(my $fh,">".$self->json_file) or die("can t open file".$self->json_file);
	my $h;
	$h->{date} = time;
	$h->{version} = $self->version;
	$h->{index} = $self->index;
	print $fh encode_json($h);
	close ($fh);
}

sub load_config {
	my ($self) = @_;
	open(my $fh ,$self->json_file) or die("can t open file ".$self->json_file);
	my $json_str = do { local $/; <$fh> };
	close($fh);
	my $h = decode_json($json_str);
	if($self->mode ne "c"){
		$self->{index} = delete $h->{index};
	}
	return delete $h->{chunks};
}
1;
