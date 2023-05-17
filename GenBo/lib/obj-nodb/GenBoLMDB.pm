package GenBoNoSqlDejaVu;



# ABSTRACT: An embedded, NoSQL SQLite database with SQL indexing
use Moo;

use strict;
use warnings;
use Data::Dumper;
extends "GenBoNoSql";
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Tie::LevelDB; 
use LMDB_File qw(:flags :cursor_op);
 use Digest::MD5 qw(md5 md5_hex md5_base64);
 
 has 
 has lmdb_extension =>(
	is		=> 'rw',
default => sub {
		return "lmdb_dejavu";
		
}
);

 has mode => (
	is		=> 'ro',
	required=> 1,
);

 has lmdb =>(
	is		=> 'ro',
	lazy =>1,
default => sub {
		my $self = shift;
		my $key1 = "lmdb_dejavu";
		$self->{is_lmdb} = 1;
	
		unless (-d $self->dir."/".$self->lmdb_extension){
			
		mkdir $self->dir."/".$self->lmdb_extension unless -d $self->dir."/".$self->lmdb_extension;
		system("chmod a+rwx ".$self->dir."/".$self->lmdb_extension);
		}

	  my $env = LMDB::Env->new($self->dir."/".$self->lmdb_extension, {
      mapsize => 100 * 1024 * 1024 * 1024, # Plenty space, don't worry
      maxdbs => 20, # Some databases
      mode   => 0777,
      #flags => LMDB_File::MDB_RDONLY,
      # More options
  });
 my $txn = $env->BeginTxn(); # Open a new transactionif ()
  unless ($self->mode eq "w"){
  	#system("chmod ")
	unlink  $self->dir."/".$self->lmdb_extension."/lock.mdb"  if -e  $self->dir."/".$self->lmdb_extension."/lock.mdb";
  }
 my $db;
 if ($self->mode eq "w"){
 	$db = $txn->OpenDB( {    # Create a new database
      dbname => "dejavu",
      flags => MDB_CREATE
  });
 }
 else {
 	$db = $txn->OpenDB( {    # Create a new database
      dbname => "dejavu",
      flags => MDB_RDONLY
  });
#  warn "end";
 }
 #
# warn "end";
  return $db;
		
}
);
