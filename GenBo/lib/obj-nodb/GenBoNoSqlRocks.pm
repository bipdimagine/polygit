package GenBoNoSqlRocks;
use Moo;
##
use strictures 2;
use warnings;
use Data::Dumper;
use Storable qw/thaw freeze/;
use JSON::XS;
use Sereal
  qw(sereal_encode_with_object sereal_decode_with_object write_sereal_file read_sereal_file);
use Carp;
use RocksDB;

has version => (
	is      => 'ro',
	default => sub {
		return undef;
	}
);

has dir => (
	is       => 'ro',
	required => 1,
);
has nb => (
	is      => 'rw',
	default => 1,
);
has name => (
	is       => 'ro',
	required => 1,
);

has vmtouch => (
	is      => 'rw',
	default => sub {
		return undef;

	}
);

has mode => (
	is       => 'ro',
	required => 1,
);

has env => ( is => 'rw', );

has compression => (
	is      => 'rw',
	default => sub {
		return "lz4hc";

	}
);
has temporary => (
	is      => 'rw',
	default => sub {
		return undef;

	}
);
has has_config => (
	is      => 'ro',
	default => sub {
		return undef;
	}
);

has pipeline => (
	is      => 'rw',
	default => sub {
		return undef;
	}
);
has tmp_dir => (
	is      => 'rw',
	default => sub {
		return "/tmp/pipeline/";
	}
);

has is_compress => (
	is      => 'rw',
	default => sub {
		return undef;

	}
);
has path_rocks => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $st   = $self->dir . "/" . $self->name . ".rocksdb";
		return $st;

	}
);


sub dictionary_array {
	my ($self) = @_;
	return $self->{dictionary_array} if exists $self->{dictionary_array};
	$self->{dictionary_array} = [];
	foreach my $v ( sort { $self->dictionary->{$a} <=> $self->dictionary->{$b} }
		keys %{ $self->dictionary } )
	{
		push( @{ $self->{dictionary_array} }, $v );
	}

	return $self->{dictionary_array};
}
has dictionary => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $hash = {};
		if ( -e $self->dictionary_file ) {
			return $self->sereal_decoder->decode_from_file(
				$self->dictionary_file );
		}
		return $hash;
	}
);

sub id_dictionary {
	my ($self) = @_;
	$self->{id_dictionary} = -1 unless exists $self->{id_dictionary};
	$self->{id_dictionary}++;
	return $self->{id_dictionary};
}

sub get_id_dictionary {
	my ( $self, $text ) = @_;
	return $self->dictionary->{$text} if exists $self->dictionary->{$text};
	delete $self->{dictionary_array}  if exists $self->{dictionary_array};

	$self->dictionary->{$text} = $self->id_dictionary;
	return $self->dictionary->{$text};

}

sub get_text_dictionary {
	my ( $self, $id ) = @_;
	confess($id) unless defined $id;
	return $self->dictionary_array->[$id];
}

has intspan_keys => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		if ( -e $self->intspan_file ) {
			my $toto;
			$self->sereal_decoder->decode_from_file( $self->intspan_file,
				$toto );
			return $toto;
		}
		return Set::IntSpan::Fast::XS->new();
	},
);

has json_file => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		system( "mkdir -p " . $self->path_rocks ) unless -e $self->path_rocks;
		my $file = $self->path_rocks . "/configuration.db.json";
		if ( $self->mode eq "c" && -e $file ) {
			unlink $file;
		}
		return $file;
	}
);

sub dictionary_file {
	my ($self) = @_;
	my $file = $self->path_rocks . "/dictionary.sereal";
	if ( $self->mode eq "c" && -e $file ) {
			unlink $file;
		}
	return $file;
}

has intspan_file => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $file = $self->path_rocks . "/intspan.sereal";
		if ( $self->mode eq "c" && -e $file ) {
			unlink $file;
		}
		return $file;
	},
);

sub exists_rocks {
	my $self = shift;
	return -e ( $self->path_rocks . '/CURRENT' );
}

has pack => (
	is      => 'rw',
	default => sub {
		return undef;
	}
);

has description => (
	is      => 'rw',
	default => sub {
		return [];
	}
);

has hash_description => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $hash = { map { $self->description->[$_] => $_ }
			  0 .. $#{ $self->description } };
		return $hash;
	}
);
has factor => (
	is      => 'rw',
	default => sub {
		return [];
	}
);


sub write_config {
	my ($self) = @_;
	open( my $fh, ">" . $self->json_file ) or die("can t open file");
	my $h;
	$h->{date}        = time;
	$h->{version}     = $self->version;
	$h->{pack}        = $self->pack;
	$h->{description} = $self->description;
	$h->{factor}      = $self->factor;
	print $fh encode_json($h);
	close($fh);
}

sub load_config {
	my ($self) = @_;
	return if exists $self->{load};
	$self->{load} = 1;
	my $file = $self->json_file;
	if ( $self->mode eq "c" ) {
		return;
	}
	open( my $fh, $self->json_file ) or die("can t open file");
	my $json_str = do { local $/; <$fh> };
	close($fh);
	my $h = decode_json($json_str);
	if ( $self->mode ne "c" ) {

		$self->{pack}        = delete $h->{pack};
		$self->{description} = delete $h->{description};
		$self->{factor}      = delete $h->{factor};
	}
	return delete $h->{chunks};
}

sub  create_mode {
	my ($self) = @_;
	return ($self->mode eq "c" or $self->mode eq "t");
}

sub rocks {
	my $self = shift;
	return $self->{rocks} if exists $self->{rocks};
	my $bits_per_key = 10;
	my $policy = RocksDB::BloomFilterPolicy->new($bits_per_key);
	if ( $self->pipeline ) {
		my $zz = time + rand(5000);
		$self->path_rocks(
			$self->tmp_dir . "/" . $self->name . "." . $zz . ".rocksdb" );
	}
	if ($self->create_mode){
		system("mkdir ".$self->dir." && chmod a+rwx ".$self->dir) unless -e $self->dir;
		
	}
	confess unless $self->mode;
	if ( $self->mode eq "r" ) {

		#confess($self->json_file) if -e $self->json_file;
		$self->load_config() if -e $self->json_file;
		confess( $self->path_rocks . '/CURRENT' )
		
		  unless ( $self->exists_rocks() );
		$self->{rocks} = RocksDB->new(
			$self->path_rocks,
			{
				cache_index_and_filter_blocks => 1,
				user_key_len               => 10,
				verify_checksums              => 0,
				read_only                     => 1,
				allow_mmap_reads           => 1,
				max_background_compactions => 3,
				IncreaseParallelism        => 1,
				filter_policy => $policy 
			}
		);

		#$rocks->IncreaseParallelism();
		return $self->{rocks};
	}
	elsif ( $self->mode eq "w" ) {
		
		$self->load_config() if -e $self->json_file;
		confess( $self->path_rocks . '/CURRENT' ) unless ( $self->exists_rocks() );

		#my $options = RocksDB::Options->new();
		my $bits_per_key = 10;

		#my $policy = RocksDB::BloomFilterPolicy->new($bits_per_key);
		$self->{rocks} = RocksDB->new(
			$self->path_rocks,
			{
				user_key_len               => 10,
				allow_mmap_reads           => 1,
				max_background_compactions => 0,
				IncreaseParallelism        => 1,
				filter_policy => $policy 
			}
		);
		return $self->{rocks};
	}
	elsif ( $self->mode eq "m" ) {
		$self->write_config();
		confess( $self->path_rocks . '/CURRENT' )
		  unless ( $self->exists_rocks() );

		#my $options = RocksDB::Options->new();
		my $bits_per_key = 10;

		#my $policy = RocksDB::BloomFilterPolicy->new($bits_per_key);
		$self->{rocks} = RocksDB->new(
			$self->path_rocks,
			{
				user_key_len               => 10,
				allow_mmap_reads           => 1,
				max_background_compactions => 0,
				IncreaseParallelism        => 1,
				filter_policy => $policy 
			}
		);
		return $self->{rocks};
	}
	elsif ( $self->mode eq "c" or $self->mode eq "t" ) {

		if ( -e $self->dictionary_file ) {
			unlink $self->dictionary_file;
			delete $self->{dictionary_file};
		}
		RocksDB->destroy_db( $self->path_rocks ) if -e $self->path_rocks;
		$self->delete_files( $self->path_rocks ) if -e $self->path_rocks;

		#$self->write_config();
		$self->{rocks} = RocksDB->new(
			$self->path_rocks,
			{
				log_file_time_to_roll => 1,
				IncreaseParallelism   => 1,
				keep_log_file_num     => 1,
				create_if_missing     => 1,
				filter_policy => $policy ,
				compression           => $self->compression
			}
		);
		system( "chmod a+rwx " . $self->path_rocks );
		return $self->{rocks};
	}
	elsif ( $self->mode eq "t" ) {
		$self->write_config();
		RocksDB->destroy_db( $self->path_rocks ) if -e $self->path_rocks;
		$self->{rocks} = RocksDB->new(
			$self->path_rocks,
			{
				user_key_len          => 10,
				log_file_time_to_roll => 1,
				keep_log_file_num     => 1,
				create_if_missing     => 1,
				filter_policy => $policy ,
				compression           => $self->compression
			}
		);
		system( "chmod a+rwx " . $self->path_rocks );
		return $self->{rocks};
	}
	else {
		die();
		$self->{rocks} = RocksDB->new(
			$self->path_rocks,
			{
				log_file_time_to_roll => 1,
				IncreaseParallelism   => 1,
				keep_log_file_num     => 1,
				create_if_missing     => 1,
				filter_policy => $policy ,
				compression           => "zlib"
			}
		);
		return $self->{rocks};
	}
}

sub delete_files {
	my ( $self, $dir ) = @_;
	system( "find " . $dir . " -maxdepth 1 -type f -delete" );
}

sub encode {
	my ( $self, $value ) = @_;
	return sereal_encode_with_object( $self->sereal_encoder, $value );
}

has sereal_encoder => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;

#return Sereal::Encoder->new({compress=>Sereal::SRL_ZSTD,compress_threshold=>0});
		return Sereal::Encoder->new(
			{ compress => Sereal::SRL_UNCOMPRESSED, compress_threshold => 0 } );
		return 0;
	},
);

has sereal_decoder => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return Sereal::Decoder->new(
			{ compress => Sereal::SRL_UNCOMPRESSED, compress_threshold => 0 } );

#return Sereal::Decoder->new({compress=>Sereal::SRL_ZSTD,compress_threshold=>0});
		return 0;
	},
);

sub decode {
	my ( $self, $value ) = @_;
	return unless $value;
	
	my $x;
		$x = sereal_decode_with_object( $self->sereal_decoder, $value );
	return $x;
}
sub get_raw {
	my ( $self, $key ) = @_;
	confess() unless $self->rocks;
	return  $self->{buffer}->{$key} if (exists $self->{buffer}->{$key});
	$self->rocks->get($key);
}
sub get {
	my ( $self, $key,$debug ) = @_;
	
	confess("-  -") unless defined $key;
	confess() unless $self->rocks;
	my $v = $self->get_raw($key);
	return unless $v;
	return $self->decode($v);
}

sub get_cached {
	my ( $self, $key,$debug ) = @_;
	confess() unless defined $key;
	confess() unless $self->rocks;
	my $v;
	if (exists $self->{buffer}->{$key}){
		$v =  $self->{buffer}->{$key}
	}
	else {
		confess();
	}
	return unless $v;
	return $self->decode($v);
}

has batch => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return RocksDB::WriteBatch->new;
	},
);

sub write_batch {
	my ($self) = @_;
	return unless $self->batch;
	$self->rocks->write( $self->batch );
	delete $self->{batch};
}

sub put_batch_raw {
	my ( $self, $key, $value, $debug ) = @_;

	$self->batch->put( $key, $value );
}

sub put_batch_compress {
	my ( $self, $key, $value ) = @_;
	confess() unless $self->rocks;
	$self->batch->put( $key, compress($value) );
}

sub put_batch {
	my ( $self, $key, $value, $debug ) = @_;
	return $self->batch->put( $key, $self->encode($value) );
}

sub put_raw {
	my ( $self, $key, $value ) = @_;
	confess() unless $self->rocks;
	$self->rocks->put( $key, $value );
}



sub put {
	my ( $self, $key, $value ) = @_;
	confess() unless $self->rocks;
	confess() unless $key;
	$self->rocks->put( $key, $self->encode($value) );

	#$self->_put_index($key) if ($self->is_index);
}

sub exists {
	my ( $self, $key ) = @_;
	return $self->rocks->exists($key);
}

sub close {
	my ($self) = @_;
	if ( $self->mode ne "r" ) {
		if ( $self->has_config() ) {
			$self->write_config();

		}
		if ( exists $self->{intspan_keys} && $self->mode ne "r" ) {
			$self->sereal_encoder->encode_to_file( $self->intspan_file,
				$self->intspan_keys );
		}
		if ( keys %{ $self->dictionary } ) {
			warn $self->dictionary_file;
			$self->sereal_encoder->encode_to_file( $self->dictionary_file,
				$self->dictionary );
		}
		if ( exists $self->{batch} ) {

			#warn "write ".$self->path_rocks;
			$self->rocks->write( $self->batch );

			#warn "end ".$self->path_rocks;
		}

		#confess() if $self->path_rocks =~ /_nosplit/;
		#warn "\t\t compact ".$self->path_rocks;

		$self->rocks->compact_range();

		#warn "\t\t end compact ".$self->path_rocks;
	}

	if ( $self->pipeline ) {
		my $dir_prod = $self->dir . "/" . $self->name . ".rocksdb";
		system("mkdir $dir_prod && chmod a+rwx $dir_prod") unless -e $dir_prod;
		$self->delete_files($dir_prod);
		warn "rsync -rav --remove-source-files "
		  . $self->path_rocks
		  . "/ $dir_prod/";
		system( "rsync -rav --remove-source-files "
			  . $self->path_rocks
			  . "/ $dir_prod/" );
		system( "rmdir " . $self->path_rocks );
		die() if $? ne 0;
		$self->pipeline(undef);
	}
	if ( $self->temporary ) {
		$self->delete_files( $self->path_rocks );
		system( "rmdir " . $self->path_rocks );

		$self->temporary(undef);
	}

	#$self->DESTROY();
	#$self->rocks->close();
	delete $self->{rocks};
	$self->{rocks} = undef;

	#$self = undef;
}


sub prepare {
	my ($self,$list) = @_;
	return  if scalar (@$list) == 0;
	$self->{buffer} = $self->rocks->get_multi(@$list);
	return 1;
}



sub return_rocks_id_from_gnomad_id {
	my ($self,$id) = @_;
	warn $id;
	my ($chr,$pos,$ref,$alt) = split("-",$id);
	warn "))))))))";
		warn $id." ".$ref." "."++++";
	return $self->return_rocks_id($pos,$ref,$alt);
}
sub return_genomic_rocks_id_from_gnomad_id {
	my ($self,$id) = @_;
	confess() unless $id;
	my ($chr,$pos,$ref,$alt) = split("-",$id);
	return $chr."!".$self->return_rocks_id($pos,$ref,$alt);
}

sub return_rocks_id_from_genbo_id {
	my ($self,$id) = @_;
	my ($chr,$pos,$ref,$alt) = split("_",$id);

	
	return $self->return_rocks_id($pos,$ref,$alt);
}

sub return_rocks_id {
	my ($self,$pos,$ref,$alt) = @_;
	 $pos  = $self->stringify_pos($pos);
#	 my $zid = $self->compress_vcf_position($ref,$alt);
#	 return $pos."!".$zid;
	my $l1 = length($ref);
	my $l2 = length($alt);
	return  ($pos."!".$alt) if ($l1 == 1 && $l2 ==1);
	my $seqid = $alt;
	if ($alt =~ /del/ ){
		return  ($pos."!".$alt);
	}
	elsif ($ref =~ /del/ ){
		return  ($pos."!".$ref);
	}
	elsif ($alt=~ /inv/ ){
		return  ($pos."!".$alt);
	}
	elsif ($alt=~ /ins/ ){
		return  ($pos."!".$alt);
	}
	elsif ($alt=~ /dup/ ){
		return  ($pos."!".$alt);
	}
	
	
	die($ref." ".$alt) if $alt =~ /INV/;
	#die(ref." ".$alt) if $alt =~ /DEL/;s
	die(ref." ".$alt) if $alt =~ /DUP/;
	
	
	if ($l1 ==1 && $l2 > 1){
		
		$seqid = "+".substr($alt, 1);
		return  ($pos."!".$seqid);
	}
	elsif ($l1 >1 && $l2 ==1){
		$ref = substr($ref, 1);
		$seqid = ($l1 -1);
		return  ($pos."!".$seqid);
	}
	 elsif ($l1 >1 && $l2 == $l1 && $l2>1 ){
	 	
		$ref = substr($ref, 1);
		$alt = substr($alt, 1);
		$seqid = "$ref*$alt";
		return  ($pos."!".$seqid);
	}
	
	else {
		return  ($pos."!".$ref."*".$alt);
		confess($l1." ".$l2);
	}
	
}
sub compress_vcf_position {
	my ($self,$ref,$alt) = @_;
	if ($alt =~/INV/){
		warn $alt;
		my @z = split("-",$alt);
		die("!"."@".$z[-1]);
		return "!"."@".$z[-1];
	}
	my $l1 = length($ref);
	my $l2 = length($alt);
	if ($l1 ==1 && $l2 ==1){
		return $alt;
	}
	if ($l1 ==1 && $l2 > 1){
		return "+".substr($alt, 1);
	}
	if ($l1 >1 && $l2 == 1){
		$ref = substr($ref, 1);
		return ($l1 -1);
	}
	if($l2 == 0 ){
		return ($ref."+");
	}
	confess($ref." ".$alt);
	
}

sub stringify_pos {
	my ($self,$pos) = @_;
	return ($pos,sprintf("%010d", $pos));
}

sub DESTROY {
	my ($self) = @_;
	#if ($self->mode ne "r"){
	#	$self->rocks->compact_range();
	#}
	system( "rm -f " . $self->path_rocks() . "/LOG*" );
	system( "rm -f " . $self->path_rocks() . "/LOCK" );
	if ( $self->temporary && -e $self->path_rocks ) {
		$self->delete_files( $self->path_rocks );
		system( "rmdir " . $self->path_rocks );
	}
	if ( $self->pipeline && $self->path_rocks =~ /tmp/ ) {
		warn "warning destruct object without closing delete "
		  . $self->path_rocks();
		$self->delete_files( $self->path_rocks );
		system( "rmdir " . $self->path_rocks );
	}

}

1;
