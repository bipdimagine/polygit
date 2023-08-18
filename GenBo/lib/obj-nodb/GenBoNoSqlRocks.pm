package GenBoNoSqlRocks;
use Moo; 
#
use strictures 2;
use warnings;
use Data::Dumper;
use Storable qw/thaw freeze/;
use JSON::XS;
use Sereal qw(sereal_encode_with_object sereal_decode_with_object write_sereal_file read_sereal_file);
use Carp;
use Compress::Zstd;
use RocksDB;

has version =>(
	is		=> 'ro',
	default => sub {
		return undef;
	}
);


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
sub  dictionary_array {
	my $self =@_;
	return $self->{dictionary_array} if exists $self->{dictionary_array};
	$self->{dictionary_array} =[];
	foreach my $v (sort{$self->dictionary->{$a} <=> $self->dictionary->{$b}} keys %{$self->dictonary}){
			push(@{$self->{dictionary_array}},$v);
		}
		return $self->{dictionary_array};
}	
has dictionary => (
is		=> 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $hash ={}; 
		if (-e $self->dictionary_file){
			return $self->sereal_decoder->decode_to_file($self->dictionary_file);
		}
		return $hash;
	}
);


sub id_dictionary {
	my ($self) = @_;
	$self->{id_dictionary} = -1 unless exists $self->{id_dictionary};
	$self->{id_dictionary} ++;
	return $self->{id_dictionary};
}

sub get_id_dictionary {
	my ($self,$text) = @_;
	delete $self->{dictionary_array};
	unless (exists $self->dictionary->{$text}) {
		$self->dictionary->{$text} = $self->id_dictionary;
	}
	return $self->dictionary->{$text} = $self->id_dictionary;
}
sub get_text_dictionary {
	my ($self,$id) = @_;
	return $self->dictionary_array->[$id];
}

has intspan_keys => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		if( -e $self->intspan_file){
			my $toto;
			$self->sereal_decoder->decode_from_file($self->intspan_file,$toto);
			return $toto;
		}
		return Set::IntSpan::Fast::XS->new( );
	},
);

has json_file =>(
is		=> 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		system("mkdir -p ".$self->path_rocks) unless -e $self->path_rocks;
		return $self->path_rocks."/configuration.db.json";
	}
);

has dictionary_file => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return  $self->path_rocks."/dictionary.sereal";
	},
);
has intspan_file => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return  $self->path_rocks."/intspan.sereal";
	},
);
sub exists_rocks  {
	my $self = shift;
	return -e ($self->path_rocks.'/CURRENT');
}


has pack =>(
	is		=> 'rw',
	default => sub {
		return undef;
	}
);

has description =>(
	is		=> 'rw',
	default => sub {
		return [];
	}
);

has hash_description =>(
	is		=> 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		 my $hash = { map {$self->description->[$_] => $_} 0 .. $#{$self->description} };
		return $hash;
	}
);
has factor =>(
	is		=> 'rw',
	default => sub {
		return [];
	}
);
sub write_config {
	my ($self) = @_;
	open(my $fh,">".$self->json_file) or die("can t open file");
	my $h;
	$h->{date} = time;
	$h->{version} = $self->version;
	$h->{pack} = $self->pack;
	$h->{description} = $self->description;
	$h->{factor} = $self->factor;
	print $fh encode_json($h);
	close ($fh);
}

sub load_config {
	my ($self) = @_;
	open(my $fh ,$self->json_file) or die("can t open file");
	my $json_str = do { local $/; <$fh> };
	close($fh);
	my $h = decode_json($json_str);
	if($self->mode ne "c"){
	$self->{pack} = delete $h->{pack};
	$self->{description} = delete $h->{description};
	$self->{factor} = delete $h->{factor};
	}
	return delete $h->{chunks};
}

 sub rocks {
	my $self = shift;
	return 	$self->{rocks} if exists $self->{rocks};
	if($self->pipeline){
		my $zz =time + rand(5000);
		$self->path_rocks($self->tmp_dir."/".$self->name.".".$zz.".rocksdb");
	}
	if ($self->mode eq "r"){
		#confess($self->json_file) if -e $self->json_file;
		$self->load_config() if -e $self->json_file;
		confess($self->path_rocks.'/CURRENT') unless ($self->exists_rocks());
		$self->{rocks} = RocksDB->new($self->path_rocks,{cache_index_and_filter_blocks=>1,verify_checksums=>0,user_key_len=>10,read_only=>1});
		system("vmtouch -q -t ".$self->path_rocks);
		#$rocks->IncreaseParallelism();
		return  $self->{rocks};
	}
	elsif ($self->mode eq "w"){
		$self->load_config() if -e $self->json_file;
		confess($self->path_rocks.'/CURRENT') unless ($self->exists_rocks());
		#my $options = RocksDB::Options->new();
		my $bits_per_key = 10;
		#my $policy = RocksDB::BloomFilterPolicy->new($bits_per_key);
		$self->{rocks} = RocksDB->new($self->path_rocks,{user_key_len=>10,allow_mmap_reads=>1,max_background_compactions=>0,IncreaseParallelism=>1});
		return  $self->{rocks};
	}
	elsif ($self->mode eq "m"){
		$self->write_config();
		confess($self->path_rocks.'/CURRENT') unless ($self->exists_rocks());
		#my $options = RocksDB::Options->new();
		my $bits_per_key = 10;
		#my $policy = RocksDB::BloomFilterPolicy->new($bits_per_key);
		$self->{rocks} = RocksDB->new($self->path_rocks,{user_key_len=>10,allow_mmap_reads=>1,max_background_compactions=>0,IncreaseParallelism=>1});
		return  $self->{rocks};
	}
	elsif ($self->mode eq "c"){
		 RocksDB->destroy_db($self->path_rocks) if -e $self->path_rocks;
		 $self->delete_files($self->path_rocks) if -e $self->path_rocks;
		 $self->write_config();
		  $self->{rocks} =   RocksDB->new($self->path_rocks, { log_file_time_to_roll=>1 ,IncreaseParallelism => 1,keep_log_file_num=>1,create_if_missing => 1,compression=>"zlib"});
		 system("chmod a+rwx ".$self->path_rocks);
		 return $self->{rocks};
	}
		elsif ($self->mode eq "t"){
		$self->write_config();
		 RocksDB->destroy_db($self->path_rocks) if -e $self->path_rocks;
		  $self->{rocks} =   RocksDB->new($self->path_rocks, { user_key_len=>10,log_file_time_to_roll=>1 ,keep_log_file_num=>1,create_if_missing => 1,compression=>"zlib"});
		 system("chmod a+rwx ".$self->path_rocks);
		 return $self->{rocks};
	}
	else {
		die();
		 $self->{rocks} =   RocksDB->new($self->path_rocks, {  log_file_time_to_roll=>1 ,IncreaseParallelism => 1,keep_log_file_num=>1,create_if_missing => 1,compression=>"zlib"});
		  return $self->{rocks};
	}
}
sub delete_files {
	my ($self,$dir) =@_;
	system("find ".$dir." -maxdepth 1 -type f -delete");
}

sub encode {
	my ($self,$value)  =@_;
	return sereal_encode_with_object($self->sereal_encoder, $value);
}

has sereal_encoder => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		#return Sereal::Encoder->new({compress=>Sereal::SRL_ZSTD,compress_threshold=>0});
		return Sereal::Encoder->new({compress=>Sereal::SRL_UNCOMPRESSED,compress_threshold=>0});
		return 0;
	},
);

has sereal_decoder => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return Sereal::Decoder->new({compress=>Sereal::SRL_UNCOMPRESSED,compress_threshold=>0});
		#return Sereal::Decoder->new({compress=>Sereal::SRL_ZSTD,compress_threshold=>0});
		return 0;
	},
);
sub decode {
	my ($self,$value)  = @_;
	return unless $value;
	return sereal_decode_with_object($self->sereal_decoder, $value);
}

sub get {
	my ($self,$key) = @_;
	confess() unless $key;
	my $v = $self->rocks->get($key);
	return unless $v;
	return $self->decode($v);
}
 has batch => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return RocksDB::WriteBatch->new;;
	},
);
sub write_batch {
	my ($self) = @_;
	$self->rocks->write($self->batch);
	delete $self->{batch};
}



sub put_batch_raw {
	my ($self,$key,$value,$debug) = @_;

	$self->batch->put($key,$value);
}
 sub put_batch_compress{
	my ($self,$key,$value) = @_;
	confess() unless $self->rocks;
	$self->batch->put($key,compress($value));
}
sub put_batch {
	my ($self,$key,$value,$debug) = @_;
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
		$self->sereal_encoder->encode_to_file($self->intspan_file,$self->intspan_keys);
	}
	if (keys %{$self->dictionary}){
		$self->sereal_encoder->encode_to_file($self->dictionary_file,$self->dictionary);
	}
	if (exists $self->{batch}){
		#warn "write ".$self->path_rocks;
		$self->rocks->write($self->batch);
		#warn "end ".$self->path_rocks;
	}
	#confess() if $self->path_rocks =~ /_nosplit/;
		#warn "\t\t compact ".$self->path_rocks;
		$self->rocks->compact_range();
		#warn "\t\t end compact ".$self->path_rocks;
	}
	
	if ($self->pipeline){
		my $dir_prod = $self->dir."/".$self->name.".rocksdb";
		system("mkdir $dir_prod && chmod a+rwx $dir_prod") unless -e $dir_prod;
		$self->delete_files($dir_prod);
		warn "rsync -rav --remove-source-files ".$self->path_rocks."/ $dir_prod/";
		system("rsync -rav --remove-source-files ".$self->path_rocks."/ $dir_prod/");
			system("rmdir ". $self->path_rocks);
		die() if $? ne 0;
		$self->pipeline(undef);
	}
	if ($self->temporary){
		$self->delete_files($self->path_rocks);
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
	system("rm -f ".$self->path_rocks()."/LOG*");
	system("rm -f ".$self->path_rocks()."/LOCK");
	if ($self->temporary && -e $self->path_rocks){
		$self->delete_files($self->path_rocks);
		system("rmdir ". $self->path_rocks);
	}
	if ($self->pipeline && $self->path_rocks =~ /tmp/) {
		warn "warning destruct object without closing delete ".$self->path_rocks();
		$self->delete_files($self->path_rocks);
		system("rmdir ". $self->path_rocks);
	}
	#warn "destroy rocks ".time;
}


1;