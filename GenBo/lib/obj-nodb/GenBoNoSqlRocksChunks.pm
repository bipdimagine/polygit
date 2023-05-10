package GenBoNoSqlRocksChunks;
use Moose; 
use MooseX::Method::Signatures;
use strict;
use warnings;
use Data::Dumper;
extends "GenBoNoSqlRocks";

has start =>(
	is		=> 'ro',
	required=>1,

);
has genomic =>(
	is		=> 'ro',
	default => undef,

);

has pack =>(
	is		=> 'ro',
	required=>1,
); 
has compress =>(
	is		=> 'ro',
);
has integer =>(
	is		=> 'ro',
	default => undef,
);
has index =>(
	is		=> 'ro',
	default => undef,
);


sub put_raw_batch_pack {
	my ($self,$key,$value) = @_;
	confess();
	$self->add_index($key);
	$self->batch->put($key,$value);
}
has nb =>(
	is		=> 'rw',
	default => 1,
);
#has index_hash =>(
#	is		=> 'ro',
#	default => {},
#); 

sub write_batch {
	my ($self) = @_;
	$self->rocks->write($self->batch);
	delete $self->{batch};
}
sub put_batch {
	my ($self,$key,$value,$debug) = @_;

	$self->add_index($key);
	$key = $self->index_position($key);
	if ($self->nb%50_000 == 0){
		$self->write_batch();
		$self->nb(1);
		
		#$self->rocks->compact_range();
	}
	$self->nb($self->nb+1);
	if ($self->compress){
		
		return $self->batch->put($key,$self->encode($value));
	}
	elsif ($self->pack){
		return $self->batch->put($key,pack($self->pack,@$value));
	}
	else {
		return $self->batch->put($key,$self->encode($value));
	}
}

sub index_position {
	my ($self,$key) = @_;
	
	return $key unless $self->index;
	return sprintf("%010d", $key) if $self->index eq "genomic";
	if ($self->index eq "vcf") {
		my @t = split("_",$key);
		$t[0] = sprintf("%010d",$t[0]);
		my $l1  =length($t[1]);
		my $l2  =length($t[2]);
		if($l1==$l2){
			return $t[0]."_".$t[2]; 
		}
		elsif ($l1==1 && $l2 >1 ){
			return $t[0]."+".substr($t[2], 1);
		}
		elsif($l1==1 && $l2 >1 ){
			return $t[0]."-".($l1-1);
		}
		else {
			confess();
		}
		}
	confess();
}
override 'put' => sub {
	my ($self,$key,$value) = @_;
	$key = $self->index_position($key) ;
	if ($self->pack){
		confess();
		return unpack($self->pack,$self->rocks->get($key));
	}
	else {
		
		return  $self->SUPER::put($key,$value);
	}
};
override 'get' => sub {
	my ($self,$key) = @_;
	$key = $self->index_position($key) ;
	if ($self->pack){
		return unpack($self->pack,$self->rocks->get($key));
	}
	else {
		
		return  $self->SUPER::get($key);
	}
};

has intspan_keys => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $intspan = $self->SUPER::get("&intspan");
		unless ($intspan){
			 return Set::IntSpan::Fast::XS->new( );
		}
	},
);

sub add_index {
	my ($self,$key) = @_;
	if ($self->index eq "integer"){
		$self->intspan_keys->add($key);
		
	}
	else{
		my (@tab) = split("_",$key);
		$self->intspan_keys->add($tab[0]);
	} 
}

sub put_raw_integer {
	my ($self,$key,$value) = @_;
	confess();
	confess() unless $self->rocks;
	$self->add_index($key);
	$self->rocks->put($key,$value);
	#$self->_put_index($key) if ($self->is_index);
}

sub put_integer {
	my ($self,$key,$value) = @_;
	confess();
	confess() unless $self->rocks;
	
	$self->intspan_keys($key);
	$self->put($key,$value);
}
sub close {
	my ($self) =@_;
	if ($self->mode ne "r"){
	#if (exists $self->{intspan_keys}){
	$self->SUPER::put("&intspan",$self->{intspan_keys});
	#}
	
	if (exists $self->{batch}){
	#	warn "write ".$self->path_rocks;
		$self->rocks->write($self->batch);
	#	warn "end ".$self->path_rocks;
	}
	
		warn "\t\t compact ".$self->path_rocks;
		$self->rocks->compact_range();
		warn "\t\t end compact ".$self->path_rocks;
	}
	#$self->DESTROY();
	#$self->rocks->close();
	delete $self->{rocks};
	$self->{rocks} = undef;
	$self = undef;
}

1;