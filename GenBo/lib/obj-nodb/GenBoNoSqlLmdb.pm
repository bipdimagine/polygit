package GenBoNoSqlLmdb;



# ABSTRACT: An embedded, NoSQL SQLite database with SQL indexing
use Moo;

use strict;
use warnings;
use Data::Dumper;
use GenBoCursorLmdb;
use GenBoCursorLmdbNoIndex;
 use Set::IntervalTree;
#use Gzip::Faster;
use LMDB_File qw(:flags :cursor_op :error);
 use Digest::MD5 qw(md5 md5_hex md5_base64);
use Compress::Snappy;
#use Compress::Zstd;
use Carp;
use Digest::MD5 qw(md5 md5_hex md5_base64);
# use Compress::Zstd;
 # use Compress::LZ4;
use Compress::Zlib qw( zlib_version); 
use Storable qw/thaw freeze/;
use JSON::XS;
use Carp;
#use  Compress::Zstd;
use Compress::Zlib qw( zlib_version); 
has lmdb_extension =>(
	is		=> 'rw',
default => sub {
		return "lmdb";
		
}
);

has dir => (
	is		=> 'ro',
	required=> 1,
);
has vmtouch => (
	is		=> 'rw',
default => sub {
		return undef;
		
}
);
has name =>(
	is		=> 'ro',
	required=> 1,
);
has mode => (
	is		=> 'ro',
	required=> 1,
);
has  env => (
	is		=> 'rw',
);

has test => (
	is		=> 'rw',
	default => sub {
		return undef;
	}
);

has txn => (
	is		=> 'rw',
);
has is_compress =>(
	is		=> 'rw',
default => sub {
		return undef;
		
}
);
has is_index =>(
	is		=> 'rw',
default => sub {
		return undef;
		
}
);
sub exists_db {
	my ($self) = @_;
	return -e $self->filename();
}
has is_integer =>(
	is		=> 'rw',
default => sub {
		return undef;
		
}
);
has filename =>(
	is		=> 'rw',
	lazy=>1,
default => sub {
		my $self =shift;
		my $filename = $self->dir."/".$self->name;
		return $filename;
		
}
);


has filename_index =>(
	is		=> 'rw',
	lazy=>1,
default => sub {
		my $self =shift;
		my $filename = $self->dir."/".$self->name."_index";
		return $filename;
		
}
);

sub clean_files {
	my ($self) = @_;
	die() unless $self->mode eq "c";
	unlink $self->filename  if -e $self->filename;
	unlink $self->filename_index  if -e $self->filename_index;
}

sub create {
	my ($self) = @_;
	my $db_name = $self->name();
	
	$self->lmdb($db_name)->put("test","test");
	$self->lmdb($db_name)->del("test");
	system ("chmod a+rw ".$self->filename);
	if ($self->is_index){
		$self->lmdb($db_name."_index")->put(1,"test");
		$self->lmdb($db_name."_index")->del(1);
			system ("chmod a+rw ".$self->filename_index);
	}

}

sub lmdb {
	my ($self,$name) = @_;
	confess($name) unless $name;
	unless ($self->is_compress){
		$name.=".uc";
	}
	if (exists $self->{lmdb_file}->{$name}){
		return  $self->{lmdb_file}->{$name};
	}
	unless (-e $self->dir){
		mkdir $self->dir;
		system ("chmod a+rwx ".$self->dir);
	}
	

	my $filename = $self->dir."/".$name;
	
	#	my $filename = $self->filename;
		
		if ($self->mode eq "c"){
			#$self->clean_files();
  			unlink $filename  if -e $filename;
  		
 	 }
	
	
 	 if ($self->mode eq "r"){
  			confess("not find database $filename ".$self->dir."/".$self->name ) unless -e $filename;
 	 }
	 my $env;
	 	if ($self->mode eq "r"){
	 		my $filename =  $self->dir."/$name";
	 		#warn "read only ".$self->dir."/$name";
	 		
	 		if ($self->test  ){
#	 		my $newname =  md5_hex($self->dir."/$name");
#	 		my $filename = "/mnt/ramdisk/$newname";
#	 		system("rsync -auvz ".$self->dir."/$name $filename") unless -e "$filename" ;

	 		$env = LMDB::Env->new("$filename", {
      			mapsize => 100 * 1024 * 1024 * 1024, # Plenty space, don't worry
      			mode   => 0777,
      			maxreaders => 256,
      			#flags => LMDB_File::MDB_NOSUBDIR | MDB_RDONLY | MDB_NOLOCK | MDB_NOTLS | MDB_FIXEDMAP,
      			flags => LMDB_File::MDB_NOSUBDIR | MDB_RDONLY | MDB_NOLOCK | MDB_NOTLS | MDB_FIXEDMAP,
 	 			});
 			system("/software/bin/vmtouch -t $filename -q  ") if -e "/software/bin/vmtouch";
	 		}
	 		else {
	 	#	warn "rsync -rav ".$self->dir."/$name /tmp/tt/$newname";
	 		#warn "rsync -rav ".$self->dir."/$name /mnt/ramdisk/$newname";
	 		#$self->{dir} =  "/mnt/ramdisk/".$newdir; 
	 		#my $new_name = "/mnt/ramdisk/$name.".rand(10000);
	 		#
	 		
	 		$env = LMDB::Env->new($self->dir."/$name", {
      			mapsize => 100 * 1024 * 1024 * 1024, # Plenty space, don't worry
      			mode   => 0777,
      			maxreaders => 256,
      			#flags => LMDB_File::MDB_NOSUBDIR | MDB_RDONLY | MDB_NOLOCK | MDB_NOTLS | MDB_FIXEDMAP,
      			flags => LMDB_File::MDB_NOSUBDIR | MDB_RDONLY | MDB_NOLOCK | MDB_NOTLS | MDB_FIXEDMAP,
      # More options
 			 });
	 		}
 	 	#warn $self->dir."/$name";
 	 	my $sname = $self->dir."/$name";
 	  my $t = time;
 	  #system("/software/bin/vmtouch -t $sname   ");
 	   	system("/software/bin/vmtouch -t $sname -q  ") if $self->vmtouch;
 	   	#warn (abs(time -$t)." ".$sname) if $self->vmtouch;
 	   #	system("rsync -av $sname /tmp/toto.xxx");
 	   #
# 	   	 warn ($sname);
	 	}
		else {
			 
#			    warn "**** ".$self->dir."/$name";
	  			$env = LMDB::Env->new($self->dir."/$name", {
      			mapsize => 100 * 1024 * 1024 * 1024, # Plenty space, don't worry
      			mode   => 0777,
      			flags => LMDB_File::MDB_NOSUBDIR	,
      # More options
 	 });

		}
	
 	$self->{env}->{$name} = $env;
	my $txn = $env->BeginTxn();
	$txn->AutoCommit (1);
	$self->{txn}->{$name}   = $txn;

	 $self->{lockfile}->{$name} = $self->dir."/$name-lock";
	 if (-e $self->{lockfile}->{$name}){
	 	unlink $self->{lockfile}->{$name};
	 }
 	if ($name =~/index/){
 			$self->{lmdb_file}->{$name} = $txn->OpenDB( {    # Create a new database
      		flags => MDB_CREATE | MDB_INTEGERKEY 
  	});
  	
 	}
 	else {
 	if ($self->is_integer){

 	#$self->{txn}->{$name}->set_compare( sub { $a <=> $b } ) if $self->is_integer;
  	$self->{lmdb_file}->{$name} = $txn->OpenDB( {    # Create a new database
      #	dbname => "$name",
      flags => MDB_CREATE  | MDB_INTEGERKEY
  	});
  		#$$self->{lmdb_file}->{$name}->set_compare( sub { my($a,$b) =@_; $a <=> $b } ) if $self->is_integer;
 	}
 	else {	

 	$self->{lmdb_file}->{$name} = $txn->OpenDB( {    # Create a new database
      #	dbname => "$name",
      flags => MDB_CREATE
  	});
  	
 	}
 	}
 
  	return $self->{lmdb_file}->{$name} ;
	
}
sub _compress {
	my ($self,$data) = @_;
	  if ($self->is_compress == 3){
	  	require "Compress/Zstd.pm";
	#	return Compress::Zstd::compress($data);
	  	return  Compress::Zlib::compress($data);
	  }
	#return  Compress::Zlib::compress($data);
	require "Compress/Zstd.pm";
	return Compress::Zstd::compress($data);
#	return Compress::Snappy::compress($data);
}
sub _uncompress {
	my ($self,$data) = @_;
	 if ($self->is_compress == 3){
	  	return  Compress::Zlib::uncompress($data);
	  }
	my $val = unpack 'H8',$data;
	if ($val eq "28b52ffd" ){
		require "Compress/Zstd.pm";
		#my $obj =  thaw(Compress::Zlib::uncompress($code));
		return Compress::Zstd::uncompress($data);
	}
	return Compress::Snappy::uncompress($data);
}

sub decode {
	my ($self,$code) = @_;
	return undef unless $code;
	return $code unless $self->is_compress;
	my $obj = thaw ($self->_uncompress($code));
	return $obj->{data};
	
}


sub encode {
		my ($self,$code) = @_;
		return $code unless $self->is_compress;
		
	if  ($self->is_compress == 2){
		require "Compress/Zstd.pm";
		#return thaw(Compress::Zlib::compress($code));
		#my $string =  encode_json $code->{annex};
		#delete   $code->{annex};
		#$code->{annex_json} = $string;
		return  Compress::Zstd::compress(freeze ({data=>$code}));
		}
		
		#die();
#		return compress(freeze ({data=>$code}));
	#	return  Compress::Zstd::compress(freeze ({data=>$code}));
	#
		return $self->_compress(freeze ({data=>$code}));
}

sub lmdb_key {
	my ($self,$key) = @_;
	confess($key."-") unless defined $key;
	return $key if $self->is_integer;
	  if (length($key) > 500){
     	my (@s) = split ("_",$key);
     	#warn scalar(@s);
     	if (scalar(@s) eq 1){
     		my $m = md5_hex($key);
     		$key = $m;
     		
     	}
     	elsif (scalar(@s) eq 2){
     		my $m = md5_hex($s[1]);
     		$key = $s[0]."_".$m;
     	}
     	else {
     	my $a = shift @s;
     	my $b = shift @s;
     	my $m = md5_hex(join("_",@s));
     	$key = $a."_".$b."_".$m;
     	}
     }
     return $key;
}

sub _first {
	my ($self) = @_;
	my $db_name = $self->name();
	
	 $self->{$db_name}->{cursor} = $self->lmdb($db_name)->Cursor() unless exists  $self->{$db_name}->{cursor}; 
	my $data;
	$self->{$db_name}->{cursor}->get($self->{$db_name}->{last_id}, $data, MDB_LAST);
	$self->{$db_name}->{cursor}->get($self->{$db_name}->{first_id}, $data, MDB_FIRST);
	return ($self->{$db_name}->{first_id},$data);
}

#sub get_index {
#	my ($self,$key) = @_;
#	die() unless $self->is_index();
#	my $id = $self->index($db_name);
#	$self->lmdb($db_name."_index")->get($id,$key);
#	
#}


sub cursor {
	my ($self,$start,$end) = @_;
	my $db_name = $self->name();
	die() unless $self->is_index();
	$start = 0 unless ($start);
	$end = $self->nb_keys -1 unless $end;
	my $dbindex = $db_name."_index";
	my $cursor = GenBoCursorLmdb->new(lmdb_index=>$self->lmdb($dbindex),start=>$start,end=>$end);
	return $cursor;
	
}

sub _next {
	my ($self) = @_;
	my $db_name = $self->name();
	 my ($current_key,$current_value);
	 unless (exists  $self->{$db_name}->{first_id}){
	 	  ($current_key,$current_value) = $self->_first($db_name);
	 	   # $self->{$db_name}->{cursor}->get( $self->{$db_name}->{current_key}, $self->{$db_name}->{current_value},MDB_GET_CURRENT);
	 }
	 else {
	 	$current_key =  $self->{$db_name}->{next_key};
	 	$current_value =  $self->{$db_name}->{next_value};
	 	 $self->{$db_name}->{current_value} = $current_value;
	 }
	 if ($current_key eq "**"){
	 	return (undef,undef);
	 }
	elsif ($current_key eq $self->{$db_name}->{last_id}){
	 	 $self->{$db_name}->{next_key} ="**";
	 }
	 else {
	 		$self->{$db_name}->{cursor}->get($self->{$db_name}->{next_key},$self->{$db_name}->{next_value},MDB_NEXT);
	 }
	 
	 return ($current_key,$current_value);
	}
	

sub next_key {
	my ($self) = @_;
	my $db_name = $self->name();
	my @t = $self->_next($db_name);
	return $t[0];
}
sub next_value {
	my ($self) = @_;
	my $db_name = $self->name();
	my @t = $self->_next($db_name);
	return $self->decode($t[1]);
}
sub current_value {
	my ($self) = @_;
	my $db_name = $self->name();
	#my @t = $self->_next($db_name);
	return  $self->decode($self->{$db_name}->{current_value});
}
sub next_key_value {
	my ($self) = @_;
	my $db_name = $self->name();
	my @t = $self->_next($db_name);
	return ($t[0],$self->decode($t[1]));
}

sub get_keys {
	my ($self) = @_;
	my $db_name = $self->name();
	my $cursor = $self->lmdb($db_name)->Cursor();
	my $key;
	my $data;
	my @data;
	my $dh;
	my $key_end;

	$cursor->get($key_end, $data, MDB_LAST);
		my $h = $cursor->get($key, $data, MDB_FIRST);
		 push(@data,$key);
		$dh->{$key} ++;
	
	while( $key ne $key_end ){
		 $cursor->get($key, $data, MDB_NEXT);
		 push(@data,$key);
		$dh->{$key} ++;
	
	}
		 return [keys %$dh];
}




sub nb_keys {
	my ($self) = @_;
	my $db_name = $self->name();
	return 0 unless -e $self->filename();
	my $cursor = $self->lmdb($db_name)->stat();
	return $cursor->{entries};
}
sub stat {
	my ($self) = @_;
	my $db_name = $self->name();
	return  $self->lmdb($db_name)->stat();
}

sub index {
	my ($self) = @_;
	my $db_name = $self->name();
	unless (exists $self->{$db_name}->{index}){
		#$self->{$db_name}->{index} = $self->nb_keys($db_name) -1;
		$self->{$db_name}->{index} =0;
	}
	else {
		$self->{$db_name}->{index} ++;
	}
	return 	$self->{$db_name}->{index};
}


sub get_current_index {
	my ($self) = @_;
	my $db_name = $self->name();
	
	return 	$self->{$db_name}->{index};
}
sub _put_index{
	my ($self,$key) = @_;
	my $db_name = $self->name();
	my $id = $self->index($db_name);
	$self->lmdb($db_name."_index")->put($id,$key);
	return $id;
}
sub get_index{
	my ($self,$index) = @_;
	my $key = $self->get_key_index($index);
	return $self->get($key);
	
}

my $buffer_index;
sub get_varid{
	my ($self,$index) = @_;
	return $self->{key}->{$index} if exists $self->{key}->{$index};
	my $key = $self->get_key_index($index);
	$self->{key}->{$index} = $key;
	return $self->{key}->{$index};
	
}

sub get_key_index{
	my ($self,$index) = @_;
	my $db_name = $self->name();
	my $id = $self->index($db_name);
	my $key = $self->lmdb($db_name."_index")->get($index);
	return $key;
	
}

sub fast_put{
	my ($self,$key,$value) = @_;
	die();
	return $self->put($key,$value);

}

sub get_lmdb_cursor{
	my ($self,$start,$end) = @_;
	my $db_name = $self->name();
	my $cursor = GenBoCursorLmdbNoIndex->new(lmdb=>$self->lmdb($db_name));
	return $cursor;
}


sub array_get_next {
	my ($self,$cursor) = @_;
	$self->{no_index_cursor} = $self->get_lmdb_cursor() unless exists $self->{no_index_cursor};
	my ($key,$value) = $self->{no_index_cursor}->next_key();
	if ($key){
		return([$key,$self->decode($value)]);
	}
	else {
		return undef;
	}
}
sub get_next {
	my ($self,$cursor) = @_;
	$self->{no_index_cursor} = $self->get_lmdb_cursor() unless exists $self->{no_index_cursor};
	my ($key,$value) = $self->{no_index_cursor}->next_key();
	if ($key){
		return({$key=>$self->decode($value)});
	}
	else {
		return undef;
	}
}
sub get_first_index {
	my ($self) = @_;
	my $key;
	my $value;
	my $cursor =  $self->get_lmdb_cursor();
	 $cursor->get($key, $value, MDB_FIRST);
	 warn $key;
	 my $t = time;
	 my $nb;
	 my $previous = $key;
	 while ($key){
	 	
	 	$nb ++;
	 	if ($nb % 1000000 == 0 ){
	 		warn $nb." ".abs(time - $t);
	 		$t = time;
	 		
	 	}
	 $cursor->get($key, $value, MDB_NEXT);
	 }
	 die();
}

sub put{
	my ($self,$key,$value) = @_;
	my $db_name = $self->name();
	my $index = -1;
	
	 if ($self->is_index){
	 	unless (exists $value->{index_lmdb} ){
		 	$index = $self->_put_index($key);
			$value->{index_lmdb} = $index;
	 	}
	 	
	 }
	 
	$self->lmdb($db_name)->put($self->lmdb_key($key),$self->encode($value));
	return  $index;
	#$self->_put_index($key) if ($self->is_index);
} 
sub put_text{
	my ($self,$key,$value) = @_;
	my $db_name = $self->name();
	my $index = -1;
	
	# if ($self->is_index){
	 #	unless (exists $value->{index_lmdb} ){
	#	 	$index = $self->_put_index($key);
	#		$value->{index_lmdb} = $index;
	 #	}
	 	
	 #}
	 
	$self->lmdb($db_name)->put($self->lmdb_key($key),$self->encode($self->_compress($value)));
	return  $index;
	#$self->_put_index($key) if ($self->is_index);
} 
sub put_with_index{
	my ($self,$key,$value) = @_;
		my $db_name = $self->name();
	 if ($self->is_index){
	 	confess() unless exists $value->{index_lmdb};
		my $id = $value->{index_lmdb};
		$self->lmdb($db_name."_index")->put($id,$key);
	 }
	 else {
	 	confess();
	 }
	my $index = 	$value->{index_lmdb};
	$self->lmdb($db_name)->put($self->lmdb_key($key),$self->encode($value));
	return  $index;
	#$self->_put_index($key) if ($self->is_index);
} 


sub update{
	my ($self,$value) = @_;
	my $key =  $value->{id};
	die() unless defined $key;
	confess() unless exists $value->{index_lmdb};
	#warn $value->{index_lmdb};
	#confess();
	my $db_name = $self->name();	
	$self->lmdb($db_name)->put($self->lmdb_key($key),$self->encode($value));
	return  $value->{index_lmdb};
	#$self->_put_index($key) if ($self->is_index);
} 

sub get{
	my ($self,$key,$debug) = @_;
#	warn $self->lmdb_key($key) if $debug;
	my $db_name = $self->name();
	return  $self->decode($self->lmdb($db_name)->get($self->lmdb_key($key)));
}
sub del {
	my ($self,$key) = @_;
	my $db_name = $self->name();
	return  $self->decode($self->lmdb($db_name)->del($self->lmdb_key($key)));
}

sub get_cursor {
	my ($self,$key,$cursor) = @_;
	my $db_name = $self->name();
	my $value;
	#$cursor->get($key, $value, MDB_FIRST);
	$cursor->get($key, $value, MDB_SET_KEY);
	 $value = $self->get($key);
	#warn $key." ".Dumper $value;
	#warn $key." ".Dumper $self->decode($value);
	
	while ($value){   
		my $k = $key;
		 $cursor->get($key, $value, MDB_NEXT_DUP);
		
		#$cursor->get($key, $value, MDB_NEXT_DUP);
		warn $key." ".Dumper $self->decode($value);
		#die() if $k ne $key; 
	}    			
}

sub get_raw{
	my ($self,$key) = @_;
	my $db_name = $self->name();
	return $self->lmdb($db_name)->get($self->lmdb_key($key));
#	return  $self->decode($self->lmdb($db_name)->get($self->lmdb_key($key)));
}

sub put_raw{
	my ($self,$key,$value) = @_;
	my $db_name = $self->name();
	my $index = -1;
	if ($self->is_index){
	 	unless (exists $value->{index_lmdb} ){
		 	$index = $self->_put_index($key);
			$value->{index_lmdb} = $index;
	 	}
	 	
	 }
	
	$self->lmdb($db_name)->put($self->lmdb_key($key),$value);
	return  $index;
	#$self->_put_index($key) if ($self->is_index);
}

sub get_interval_tree {
	my ($self,$key) = @_;
	my $data = $self->get($key);
	my $tree = Set::IntervalTree->new;
	 foreach my $v (@{$data}) {
	 	$tree->insert(@$v);
	 }
	return $tree;
}

sub get_with_sequence{
	my ($self,$key,$key2,$debug) = @_;
	my $db_name = $self->name();
	my $hash = $self->decode($self->lmdb($db_name)->get($self->lmdb_key($key)));
	warn Dumper $hash if $debug;
	if ($hash){
		return $hash->{$key2};
	}
	return  undef;
}

sub is_lmdb_exists {
	my ($self, $name) = @_;
	my $filename = $self->dir."/".$name;
	return 1 if -e $filename;
	return;
}

sub exists {
	my ($self,$key) = @_;
	my $db_name = $self->name();
	return $self->lmdb($db_name)->get($self->lmdb_key($key));
}
sub size {
	my ($self) = @_;
	my $db_name = $self->name();
	
	return ($self->nb_keys($db_name));
}

sub ranges{
	my ($self,$N) = @_;
	my $db_name = $self->name();
	die() unless $db_name;
	
	my $M = $self->nb_keys($db_name) -1;
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

sub close {
	my ($self) = @_;
	foreach my $txn (values %{$self->{txn}}){
		next unless $txn;
		
		#delete $self->{txn}->{$txn};
		#$txn->abort();
		$txn->commit();
		
	}
		foreach my $txn (keys %{$self->{txn}}){
	#	warn $txn;
	#	$txn->commit();
		delete $self->{txn}->{$txn};
	}
	foreach my $f (values %{$self->{lockfile}}){
		unlink $f if -e $f;
	}
	foreach my $lmdb (keys %{$self->{lmdb_file}}){
		#$lmdb->close;
		delete $self->{lmdb_file}->{$lmdb};
	}
}

sub DESTROY {
	my ($self) = @_;
	foreach my $txn (keys %{$self->{txn}}){
	#	warn $txn;
	#	$txn->commit();
		delete $self->{txn}->{$txn};
	}
	foreach my $f (values %{$self->{lockfile}}){
		unlink $f if -e $f;
	}
	foreach my $lmdb (keys %{$self->{lmdb_file}}){
		#$lmdb->close;
		delete $self->{lmdb_file}->{$lmdb};
		
	}
	

   # print "In DESTROY\n";
}

sub get_values {
	my ($self) = @_;
	my $db_name = $self->name();
	return [] if $self->size == 0;
	my $cursor = $self->lmdb($db_name)->Cursor();
	my $key;
	my $data;
	my @data;
	my $dh;
	my $key_end;
	my $ref;
 eval {
		$cursor->get($key_end, $data, MDB_LAST);
	 };
	if ($@) {
		confess($self->dir."-".$self->name);
	}
	# or confess($self->dir."-".$self->name);
	my $h;
	

	 eval {
			$h = $cursor->get($key, $data, MDB_FIRST);# or confess($self->dir."-".$self->name);
	
	 };
	if ($@) {
		confess($self->dir."-".$self->name);
	}
	unless ($key){
		confess($self->dir."-".$self->name);
	}
		 push(@$ref,$self->decode($data));
		 #push(@data,$key);
		$dh->{$key} ++;
	while( $key ne $key_end ){
		
		 $cursor->get($key, $data, MDB_NEXT);
		 push(@$ref,$self->decode($data));
		$dh->{$key} ++;
	
	}
		 return $ref;
}


1;
