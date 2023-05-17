package GenBoNoSqlLmdbInteger;



# ABSTRACT: An embedded, NoSQL SQLite database with SQL indexing
use Moo;

use strict;
use warnings;
use Data::Dumper;
#use Gzip::Faster;

use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use LMDB_File qw(:flags :cursor_op :error);
 use Digest::MD5 qw(md5 md5_hex md5_base64);
 use Compress::Snappy;
use Storable qw/thaw freeze/;
has lmdb_extension =>(
	is		=> 'rw',
default => sub {
		return "lmdb";
		
}
);

has is_compress =>(
	is		=> 'rw',
default => sub {
		return undef;
		
}
);
has dir => (
	is		=> 'ro',
	required=> 1,
);

has mode => (
	is		=> 'ro',
	required=> 1,
);


has txn => (
	is		=> 'rw',
);


has env =>(
is		=> 'ro',
	lazy =>1,
default => sub {
		my $self = shift;
	
	
		unless (-d $self->dir."/".$self->lmdb_extension){
			
		mkdir $self->dir."/".$self->lmdb_extension unless -d $self->dir."/".$self->lmdb_extension;
		system("chmod a+rwx ".$self->dir."/".$self->lmdb_extension);
		}
  
 	 

	  my $env = LMDB::Env->new($self->dir."/".$self->lmdb_extension, {
      mapsize => 100 * 1024 * 1024 * 1024, # Plenty space, don't worry
      maxdbs => 30, # Some databases
      mode   => 0777,
  	  });
   unless ($self->mode eq "w"){
  	#system("chmod ")
	unlink  $self->dir."/".$self->lmdb_extension."/lock.mdb"  if -e  $self->dir."/".$self->lmdb_extension."/lock.mdb";
  }
  $env->sync (1);
  return $env;
}
);

#has lmdb =>(
#	is		=> 'ro',
#	lazy =>1,
#default => sub {
#		my $self = shift;
#		my $key1 = "lmdb_dejavu";
#	
#	
#		unless (-d $self->dir."/".$self->lmdb_extension){
#			
#		mkdir $self->dir."/".$self->lmdb_extension unless -d $self->dir."/".$self->lmdb_extension;
#		system("chmod a+rwx ".$self->dir."/".$self->lmdb_extension);
#		}
#	warn $self->dir();
#	  my $env = LMDB::Env->new($self->dir."/toto", {
#      mapsize => 100 * 1024 * 1024 * 1024, # Plenty space, don't worry
#      maxdbs => 30, # Some databases
#      mode   => 0777,
#      flags => LMDB_File::MDB_NOSUBDIR	,
#      # More options
#  });
#  $env->sync (1);
# my $txn = $env->BeginTxn(); # Open a new transactionif ()
#
#  unless ($self->mode eq "w"){
#  	#system("chmod ")
#  	$txn->AutoCommit ( [ 1 ] );
#	unlink  $self->dir."/".$self->lmdb_extension."/lock.mdb"  if -e  $self->dir."/".$self->lmdb_extension."/lock.mdb";
#  }
# my $db;
# if ($self->mode eq "w"){
# 	
# 	$db = $txn->OpenDB( {    # Create a new database
#      dbname => "dejavu",
#      flags => MDB_CREATE
#  });
#  $txn->AutoCommit (1);
#  $txn->AutoCommit ([1]);
# }
# else {
# 	$db = $txn->OpenDB( {    # Create a new database
#      dbname => "dejavu",
#      flags => MDB_RDONLY
#  });
##  warn "end";
# }
# warn $txn;
#  $self->txn($txn);
#  
# $self->{db2} =   $txn->OpenDB( {    # Create a new database
#      dbname => "test",
#      flags => MDB_CREATE
#  });
# #
## warn "end";
#  return $db;
#		
#}
#);

sub lmdb {
	my ($self,$name) = @_;
	die() unless $name;
	if (exists $self->{lmdb_file}->{$name}){
		return  $self->{lmdb_file}->{$name};
	}
	
	my $filename = $self->dir."/$name";
		if ($self->mode eq "c"){
			warn $filename;
  			unlink $filename  if -e $filename;
 	 }
	 my $env;
	 	if ($self->mode eq "r"){
	 		$env = LMDB::Env->new($self->dir."/$name", {
      			mapsize => 100 * 1024 * 1024 * 1024, # Plenty space, don't worry
      			mode   => 0777,
      			maxreaders => 256,
      			flags => LMDB_File::MDB_NOSUBDIR | MDB_RDONLY | MDB_NOLOCK,
      # More options
 	 });
	 	}
		else {
		
	  			$env = LMDB::Env->new($self->dir."/$name", {
      			mapsize => 100 * 1024 * 1024 * 1024, # Plenty space, don't worry
      			mode   => 0777,
      			flags => LMDB_File::MDB_NOSUBDIR | MDB_INTEGERKEY	,
      # More options
 	 });
			}
 	$self->{env}->{$name} = $env;
	my $txn = $env->BeginTxn();
	
	#$txn->AutoCommit (1);
	$self->{txn}->{$name}   = $txn;
	$self->{txn}->{$name}->set_compare( sub { $a <=> $b } ) if $self->is_integer;
	 $self->{lockfile}->{$name} = $self->dir."/$name-lock";
	 if (-e $self->{lockfile}->{$name}){
	 	unlink $self->{lockfile}->{$name};
	 }
 	
 	$self->{lmdb_file}->{$name} = $txn->OpenDB( {    # Create a new database
      #	dbname => "$name",
      flags => MDB_CREATE | MDB_INTEGERKEY
  	});
  	return $self->{lmdb_file}->{$name} ;
	
}


sub decode {
	my ($self,$code) = @_;
	return undef unless $code;
	return $code unless $self->is_compress;
	if  ($self->is_compress == 2){
		my $obj = thaw (gunzip($code));
		return $obj->{data};
	}
	
	my $obj = thaw (decompress($code));
	return $obj->{data};
}
sub encode {
		my ($self,$code) = @_;
		return $code unless $self->is_compress;
	
		return compress(freeze ({data=>$code}));
}

sub lmdb_key {
	my ($self,$key) = @_;
	  if (length($key) > 500){
     	my (@s) = split ("_",$key);
     	my $m = md5_hex($s[2]."_".$s[3]);
     	$key = $s[0]."_".$s[1]."_".$m;
     }
     return $key;
}
sub search {
	my ($self,$db_name,$key) = @_;
	my $cursor = $self->lmdb($db_name)->Cursor();
	my $data;
	 $cursor->get($key, $data, MDB_SET_RANGE);
	 my $key1;
	 my $data1;
	  $cursor->get($key1, $data1, MDB_GET_CURRENT);
	  
}

sub _first {
	my ($self,$db_name) = @_;
	
	 $self->{$db_name}->{cursor} = $self->lmdb($db_name)->Cursor() unless exists  $self->{$db_name}->{cursor}; 
	my $data;
	$self->{$db_name}->{cursor}->get($self->{$db_name}->{last_id}, $data, MDB_LAST);
	$self->{$db_name}->{cursor}->get($self->{$db_name}->{first_id}, $data, MDB_FIRST);
	return ($self->{$db_name}->{first_id},$data);
}

sub jump {
	my ($self,$db_name,$pos) = @_;
	
	 $self->{$db_name}->{cursor} = $self->lmdb($db_name)->Cursor() unless exists  $self->{$db_name}->{cursor}; 
	my $data;
	$self->{$db_name}->{cursor}->get($self->{$db_name}->{last_id}, $data, MDB_LAST);
	$self->{$db_name}->{cursor}->get($pos, $data, MDB_SET);
	$self->{$db_name}->{first_id} = $pos;
	$self->{$db_name}->{cursor}->get($self->{$db_name}->{next_key},$self->{$db_name}->{next_value},MDB_NEXT);
	return ($self->{$db_name}->{first_id},$data);
}



sub _next {
	my ($self,$db_name) = @_;
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
	my ($self,$db_name) = @_;
	my @t = $self->_next($db_name);
	return $t[0];
}
sub next_value {
	my ($self,$db_name) = @_;
	my @t = $self->_next($db_name);
	return $self->decode($t[1]);
}
sub current_value {
	my ($self,$db_name) = @_;
	#my @t = $self->_next($db_name);
	return  $self->decode($self->{$db_name}->{current_value});
}
sub next_key_value {
	my ($self,$db_name) = @_;
	my @t = $self->_next($db_name);
	return ($t[0],$self->decode($t[1]));
}

sub keys {
	my ($self,$db_name) = @_;
	my $cursor = $self->lmdb($db_name)->Cursor();
	my $key;
	my $data;
	my @data;
	my $dh;
	my $key_end;

	$cursor->get($key_end, $data, MDB_LAST);
		my $h = $cursor->get($key, $data, MDB_FIRST);
		 push(@data,$key);
		 warn $key;
	$dh->{$key} ++;
	
	while( $key ne $key_end ){
		 $cursor->get($key, $data, MDB_NEXT);
		 push(@data,$key);
	#	 warn $key;
		$dh->{$key} ++;
	
	}
		 return [keys %$dh];
}


sub nb_keys {
	my ($self,$db_name) = @_;
	my $cursor = $self->lmdb($db_name)->stat();
	return $cursor->{entries};
}

sub put{
	my ($self,$db_name,$key,$value) = @_;
	#$self->lmdb();
#	warn "again";
#	 $self->{db2}->put($self->lmdb_key($key),$self->encode($value));
	$self->lmdb($db_name)->put($key,$self->encode($value)) ;
	
	
    # $self->txn->put($self->lmdb->dbi,$self->lmdb_key($key),$self->encode($value));
   #  $self->txn->commit();# if $self->lmdb->Alive;
} 
sub get{
	my ($self,$db_name,$key) = @_;
	return  $self->decode($self->lmdb($db_name)->get($key));
}
sub exists {
	my ($self,$db_name,$key) = @_;
	return $self->lmdb($db_name)->get($self->lmdb_key($key));
}

sub DESTROY {
	my ($self) = @_;
	#warn "destroy";
	foreach my $txn (values %{$self->{txn}}){
		$txn->commit();
	#	warn "commit";
	}
	foreach my $f (values %{$self->{lockfile}}){
		unlink $f if -e $f;
	}
	foreach my $lmdb (values %{$self->{lmdb_file}}){
		$lmdb->close;
	}
	
	#warn "coucou";

   # print "In DESTROY\n";
}

1;
