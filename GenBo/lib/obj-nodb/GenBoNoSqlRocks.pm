package GenBoNoSqlRocks;
use Moo; 
#
use strict;
use warnings;
use Data::Dumper;
use Storable qw/thaw freeze/;
use JSON::XS;
use Sereal qw(sereal_encode_with_object sereal_decode_with_object);
use RocksDB;


has dir => (
	is		=> 'ro',
	required=> 1,
);
has nb =>(
	is		=> 'rw',
	default => 1,
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
has temporary =>(
	is		=> 'rw',
default => sub {
		return undef;
		
}
);

has pipeline =>(
	is		=> 'rw',
default => sub {
		return undef;
}
);
has tmp_dir =>(
	is		=> 'rw',
default => sub {
		return "/tmp/pipeline/";
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


sub exists_rocks  {
	my $self = shift;
	return -e ($self->path_rocks.'/CURRENT');
}


 sub rocks {
	my $self = shift;
	return 	$self->{rocks} if exists $self->{rocks};
	if($self->pipeline){
		my $zz =time + rand(5000);
		$self->path_rocks($self->tmp_dir."/".$self->name.".".$zz.".rocksdb");
	}
	if ($self->mode eq "r"){

		confess($self->path_rocks.'/CURRENT') unless ($self->exists_rocks());

		$self->{rocks} = RocksDB->new($self->path_rocks,{user_key_len=>10,read_only=>1});
		system("vmtouch -q -t ".$self->path_rocks);
		#$rocks->IncreaseParallelism();
		return  $self->{rocks};
	}
	elsif ($self->mode eq "m"){

		confess($self->path_rocks.'/CURRENT') unless ($self->exists_rocks());
		warn "MEMORY";
		#my $options = RocksDB::Options->new();
		my $bits_per_key = 10;
		#my $policy = RocksDB::BloomFilterPolicy->new($bits_per_key);
		$self->{rocks} = RocksDB->new($self->path_rocks,{user_key_len=>10,read_only=>1,allow_mmap_reads=>1,max_background_compactions=>0,IncreaseParallelism=>1});
		system("vmtouch -q -t ".$self->path_rocks);
		return  $self->{rocks};
	}
	elsif ($self->mode eq "c"){
		 RocksDB->destroy_db($self->path_rocks) if -e $self->path_rocks;
		  $self->{rocks} =   RocksDB->new($self->path_rocks, { log_file_time_to_roll=>1 ,IncreaseParallelism => 1,keep_log_file_num=>1,create_if_missing => 1,compression=>"zlib"});
		 system("chmod a+rwx ".$self->path_rocks);
		 return $self->{rocks};
	}
		elsif ($self->mode eq "t"){
		 RocksDB->destroy_db($self->path_rocks) if -e $self->path_rocks;
		  $self->{rocks} =   RocksDB->new($self->path_rocks, { user_key_len=>10,log_file_time_to_roll=>1 ,keep_log_file_num=>1,create_if_missing => 1,compression=>"zlib"});
		 system("chmod a+rwx ".$self->path_rocks);
		 return $self->{rocks};
	}
	else {
	
		 $self->{rocks} =   RocksDB->new($self->path_rocks, {  log_file_time_to_roll=>1 ,IncreaseParallelism => 1,keep_log_file_num=>1,create_if_missing => 1,compression=>"zlib"});
		  return $self->{rocks};
	}
}

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
 has batch => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		#return Sereal::Encoder->new();
		return RocksDB::WriteBatch->new;;
	},
);
sub write_batch {
	my ($self) = @_;
	warn "write";
	$self->rocks->write($self->batch);
	warn "end";
	delete $self->{batch};
}


sub put_batch {
	my ($self,$key,$value,$debug) = @_;
	if ($self->nb%50_000 == 0){
		$self->write_batch();
		$self->nb(1);
		#$self->rocks->compact_range();
	}
	$self->nb($self->nb+1);
	return $self->batch->put($key,$self->encode($value));

}
 sub put_raw{
	my ($self,$key,$value) = @_;
	confess() unless $self->rocks;
	$self->rocks->put($key,$value);
} 
 sub get_raw{
	my ($self,$key,$value) = @_;
	confess() unless $self->rocks;
	$self->rocks->get($key);
} 

sub put{
	my ($self,$key,$value) = @_;
	confess() unless $self->rocks;
	confess() unless $key;
	$self->rocks->put($key,$self->encode($value));
	#$self->_put_index($key) if ($self->is_index);
} 


sub close {
	my ($self) =@_;
	if ($self->mode ne "r"){
	if (exists $self->{intspan_keys} && $self->mode ne "r"){
		#$self->put("&intspan",$self->{intspan_keys});
	}
	if (exists $self->{batch}){
	#	warn "write ".$self->path_rocks;
		$self->rocks->write($self->batch);
	#	warn "end ".$self->path_rocks;
	}
	#confess() if $self->path_rocks =~ /_nosplit/;
#		warn "\t\t compact ".$self->path_rocks;
		$self->rocks->compact_range();
#		warn "\t\t end compact ".$self->path_rocks;
	}
	
	if ($self->pipeline){
		my $dir_prod = $self->dir."/".$self->name.".rocksdb";
		system("mkdir $dir_prod && chmod a+rwx $dir_prod") unless -e $dir_prod;
		system("find ".$dir_prod." -maxdepth 1 -type f -delete");
		system("rsync -rav --remove-source-files ".$self->path_rocks."/ $dir_prod/");
			system("rmdir ". $self->path_rocks);
		die() if $? ne 0;
		$self->pipeline(undef);
	}
	if ($self->temporary){
		system("find ".$self->path_rocks." -maxdepth 1 -type f -delete");
		system("rmdir ". $self->path_rocks);
	
		$self->temporary(undef);
	}
	#$self->DESTROY();
	#$self->rocks->close();
	delete $self->{rocks};
	$self->{rocks} = undef;
	#$self = undef;
}
sub DESTROY {
	my ($self) =@_;
	#if ($self->mode ne "r"){
	#	$self->rocks->compact_range();
	#}
	warn "DESTROY Rocks";
	system("rm -f ".$self->path_rocks()."/LOG*");
	system("rm -f ".$self->path_rocks()."/LOCK");
	if ($self->temporary && -e $self->path_rocks){
		system("find ".$self->path_rocks." -maxdepth 1 -type f -delete");
		system("rmdir ". $self->path_rocks);
	}
	if ($self->pipeline && $self->path_rocks =~ /tmp/) {
		warn "warning destruct object without closing delete ".$self->path_rocks();
		
		system("find ".$self->path_rocks." -maxdepth 1 -type f -delete");
		system("rmdir ". $self->path_rocks);
	}
	
}


1;