package GenBoNoSqlRocks;
use Moose;
use MooseX::Method::Signatures;
use strict;
use warnings;
use Data::Dumper;
use Carp qw(cluck longmess shortmess);
use Storable qw/thaw freeze/;
use JSON::XS;
use Sereal qw(sereal_encode_with_object sereal_decode_with_object);
use RocksDB;

has dir => (
	is		=> 'ro',
	required=> 1,
);

has name =>(
	is		=> 'ro',
	required=> 1,
);

has vmtouch => (
	is		=> 'rw',
default => sub {
		return undef;
		
}
);

has mode => (
	is		=> 'ro',
	required=> 1,
);

has  env => (
	is		=> 'rw',
);

has is_index =>(
	is		=> 'rw',
default => sub {
		return undef;
		
}
);

has is_compress =>(
	is		=> 'rw',
default => sub {
		return undef;
		
}
);
has path_rocks =>(
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $st = $self->dir."/".$self->name.".rocksdb";
		return $st;
		
}
);
has rocks =>(
	is		=> 'rw',
	lazy    => 1,
default => sub {
		my $self = shift;
	#	warn 
	#unless (-e $self->path_rocks){
	#	mkdir $self->path_rocks;
	#	system ("chmod a+rwx ".$self->path_rocks);
	#}
	if ($self->mode eq "r"){
#		confess();
		my $rocks = RocksDB->new($self->path_rocks,{IncreaseParallelism => 1,read_only=>1});
		#$rocks->IncreaseParallelism();
		return  $rocks;
	}
	elsif ($self->mode eq "c"){
		 RocksDB->destroy_db($self->path_rocks) if -e $self->path_rocks;
		 my $db =   RocksDB->new($self->path_rocks, {  log_file_time_to_roll=>1 ,IncreaseParallelism => 1,keep_log_file_num=>1,create_if_missing => 1,compression=>"zlib"});
		 system("chmod a+rwx ".$self->path_rocks);
		 return $db;
	}
	else {
		 return  RocksDB->new($self->path_rocks, {  log_file_time_to_roll=>1 ,IncreaseParallelism => 1,keep_log_file_num=>1,create_if_missing => 1,compression=>"zlib"});
	}
	
}
);

has sereal_encoder => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		#return Sereal::Encoder->new();
		return Sereal::Encoder->new({compress=>Sereal::SRL_ZSTD,compress_threshold=>0});
		return 0;
	},
);

sub encode {
	my ($self,$value)  =@_;
	return sereal_encode_with_object($self->sereal_encoder, $value);
}


has sereal_decoder => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return Sereal::Decoder->new({compress=>Sereal::SRL_ZSTD,compress_threshold=>0});
		return 0;
	},
);
sub decode {
	my ($self,$value)  = @_;
	return unless $value;
	return sereal_decode_with_object($self->sereal_decoder, $value);
}

sub get{
	my ($self,$key) = @_;
	my $v = $self->rocks->get($key);
	return unless $v;
	return $self->decode($v);
}
 
sub put{
	my ($self,$key,$value) = @_;
	
	confess() unless $self->rocks;
	$self->rocks->put($key,$self->encode($value));
	#$self->_put_index($key) if ($self->is_index);
} 
sub close {
	my ($self) =@_;
	if ($self->mode ne "r"){
		$self->rocks->compact_range();
	}
	#$self->DESTROY();
	$self->{rocks} = undef;
	$self = undef;
}
sub DESTROY {
	my ($self) =@_;
	#if ($self->mode ne "r"){
	#	$self->rocks->compact_range();
	#}
	system("rm -f ".$self->path_rocks()."/LOG*");
	system("rm -f ".$self->path_rocks()."/LOCK");
}
1;