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
 use Compress::Snappy;
 use Carp;
  use List::Util qw[min max];
has extension =>(
	is		=> 'rw',
default => sub {
		return "dejavu.lite";
		
}
);

has lmdb_extension =>(
	is		=> 'rw',
default => sub {
		return "lmdb_dejavu";
		
}
);
sub create_table {
	my ($self,$key1) = @_;
	return $self->{table}->{$key1}   if (exists $self->{table}->{$key1});
	 my $table_name =$self->name_table_data;
	 $self->{table}->{$key1} = $table_name;
	 return  $self->{table}->{$key1} unless  $self->write();
	$self->dbh($key1)->do("DROP TABLE IF EXISTS $table_name")  or die $DBI::errstr  if $self->mode eq 'c'  ;

	#$self->dbh($chr)->do("CREATE TABLE if not exists $table_name (_key TEXT PRIMARY KEY, _value BLOB) ;")  or die $DBI::errstr;;
	$self->dbh($key1)->do("CREATE TABLE if not exists $table_name (_key VARCHAR(250),ho INTEGER, projects INTEGER, _value BLOB,start INTEGER,end INTEGER, variation_type VARCHAR(2)) ")  or die $DBI::errstr;;
	#$self->dbh($key1)->do("CREATE UNIQUE INDEX if not exists _key_idx  on $table_name (_key); ")  or die $DBI::errstr;;
	my $T2 = 	$table_name."POSITION";
	$self->dbh($key1)->do("CREATE VIRTUAL TABLE  if not exists $T2  USING rtree( id ,start,end)")  or die $DBI::errstr;
#	$self->dbh($key1)->do("CREATE VIRTUAL TABLE  if not exists $table_name USING fts3(_key VARCHAR(250) NOT NULL ,_value blob ) ;")  or die $DBI::errstr;
	
	return 	$self->{table}->{$key1} ;
}
has lmdb =>(
	is		=> 'ro',
	lazy =>1,
default => sub {
		my $self = shift;
		my $key1 = "lmdb_dejavu";
		$self->{is_lmdb} = 1;
	
		unless (-d $self->dir."/".$self->lmdb_extension){
			warn  $self->dir."/".$self->lmdb_extension;
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
  if ($self->mode eq "c"){
  	unlink  $self->dir."/".$self->lmdb_extension."/data.mdb"  if -e  $self->dir."/".$self->lmdb_extension."/data.mdb";
  }
 
  unless ($self->mode eq "w"){
  	#system("chmod ")
	unlink  $self->dir."/".$self->lmdb_extension."/lock.mdb"  if -e  $self->dir."/".$self->lmdb_extension."/lock.mdb";
  }
 my $db;
 
  if ($self->mode eq "r") {
  	warn "read";
 	$db = $txn->OpenDB( {    # Create a new database
      dbname => "dejavu",
      flags => MDB_RDONLY
  });
#  warn "end";
 }
  else {
  	warn "create";
 	$db = $txn->OpenDB( {    # Create a new database
      dbname => "dejavu",
      flags => MDB_CREATE ,
  });
 }


 #
# warn "end";
  return $db;
		
}
);


sub lmdb_key {
	my ($self,$key) = @_;
	  if (length($key) > 500){
     	my (@s) = split ("_",$key);
     	my $m = md5_hex($s[2]."_".$s[3]);
     	$key = $s[0]."_".$s[1]."_".$m;
     }
     return $key;
}
sub put_lmdb{
	my ($self,$key,$value) = @_;
	
     $self->lmdb->put($self->lmdb_key($key),$value);
} 
sub get_lmdb{
	my ($self,$key) = @_;
	return  $self->lmdb->get($self->lmdb_key($key))
}


sub leveldb {
	my ($self,$key1) = @_;
	return $self->{leveldb}->{$key1}   if (exists $self->{leveldb}->{$key1});
	$self->{leveldb}->{$key1}= new Tie::LevelDB::DB($self->dir."/$key1");
	return 	$self->{leveldb}->{$key1} ;
}


sub delete_bulk {
	my ($self,$key1,$key2) = @_;
	my $id = $self->change_id($key1);
	my $table_name = $self->create_table($id);

	 $self->dbh($id)->do(qq{delete  from $table_name where _key match "$key2"});
	 

	 return 1;
}

sub sth_insert_dejavu_cached {
	my ($self,$chr) = @_;
	return $self->{sth_insert_dejavu_cached}->{$chr} if exists $self->{sth_insert_dejavu_cached}->{$chr};
	 my $table = $self->create_table($chr);
	$self->{sth_insert_dejavu_cached}->{$chr} = $self->dbh($chr)->prepare("insert  into $table (_key,ho,projects,_value) values(?,?,?,?)");
	return $self->{sth_insert_dejavu_cached}->{$chr};
}


sub sth_insert_position_cached {
	my ($self,$chr) = @_;
	return $self->{sth_insert_position_cached}->{$chr} if exists $self->{sth_insert_position_cached}->{$chr};
	 my $table = $self->create_table($chr)."POSITION";
	$self->{sth_insert_position_cached}->{$chr} = $self->dbh($chr)->prepare("insert  into $table (id,start,end) values(?,?,?)");
	return $self->{sth_insert_position_cached}->{$chr};
}

sub sth_select_position_cached {
	my ($self,$chr) = @_;
	return $self->{sth_select_position_cached}->{$chr} if exists $self->{sth_select_position_cached}->{$chr};
	 my $table_name = $self->create_table($chr);
	 my $T2 = $table_name."POSITION";
	 
	$self->{sth_insert_position_cached}->{$chr} = $self->dbh($chr)->prepare("select __DATA__._key from $T2,$table_name  where start>=? and end <= ? and $T2.id=$table_name.rowid;") ;
	return $self->{sth_select_position_cached}->{$chr};
}

sub _put_dejavu {
	my ($self,$chr,$hash,$debug) = @_;
	
	$self->dbh($chr)->begin_work;
	  foreach my $key (keys %$hash){
	  	
	  			$self->sth_insert_dejavu_cached($chr)->execute($key,$hash->{$key}->{ho},$hash->{$key}->{all},$self->encode($hash->{$key}->{data}) )  or die $DBI::errstr." $key";;
	  			my $id = $self->dbh($chr)->sqlite_last_insert_rowid();
	  			$self->sth_insert_position_cached($chr)->execute($id,$hash->{$key}->{start},$hash->{$key}->{end}) ;
	  
	 }
	 
	 $self->dbh($chr)->commit;
 return 1;
}

sub get_position {
	my ($self,$chr,$start,$end) = @_;
	 my $table_name = $self->create_table($chr);
	 my $T2 = $table_name."POSITION";
	# warn "select __DATA__._key from $T2,$table_name  where $T2.start>=$start and $T2.end<=$end and $T2.id=$table_name.rowid;";
	my $aray_ref = $self->dbh($chr)->selectcol_arrayref("select __DATA__._key from $T2,$table_name  where $T2.start<=$end and $T2.end>=$start and $T2.id=$table_name.rowid;",{ Columns=>[1] });
	#my $aray_ref = $self->dbh($chr)->selectcol_arrayref("select __DATA__._key from $T2,$table_name  where $T2.start>=$start and $T2.end<=$end and $T2.id=$table_name.rowid;",{ Columns=>[1] });
	
	
	return $aray_ref;
	
	warn Dumper $aray_ref;
	die();
	
}

sub prepare_cnv {
	my ($self,$chr) = @_;

	return $self->{prepare_cnv}->{$chr} if exists $self->{prepare_cnv}->{$chr};
		 my $table_name = $self->create_table($chr);

	$self->{prepare_cnv}->{$chr} = $self->dbh($chr)->prepare(qq{select $table_name.variation_type, $table_name.projects,$table_name.start,$table_name.end,$table_name._value from $table_name  where end>=? and start<=?  and variation_type=?   });
	return $self->{prepare_cnv}->{$chr};
}

sub get_cnv {
	my ($self,$chr,$start,$end,$type,$dejavu,$seuil) = @_;
	
	 my $table_name = $self->create_table($chr);
	 my $T2 = $table_name."POSITION";
	 
	 #my $max = int( ($end-$start)*0.75);
	 #my $lmax = ($end-$start) * 1.25;
	  #my $lmin = ($end-$start) * 0.75;

	  #$lmax = 1_000_000;
	 #my $sth = $self->dbh($chr)->prepare(qq{select $table_name.variation_type, $table_name.projects,$table_name.start,$table_name.end,$table_name._value from $table_name  where start<=$end and end>=$start   and variation_type="$type";});
	 #$self->sth_insert_dejavu_cached($chr)->execute($key,$hash->{$key}->{ho},$hash->{$key}->{all},$self->encode($hash->{$key}->{data}) )  or die $DBI::errstr." $key";;
	 
	# my $sth = $self->dbh($chr)->prepare(qq{select $table_name.variation_type, $table_name.projects,$table_name.start,$table_name.end,$table_name._value from $table_name, $T2  where $T2.start<=$end and $T2.end>=$start and $T2.id=$table_name.rowid  and $table_name.variation_type="$type";});
	#	my $sth = $self->dbh($chr)->prepare(qq{select ($table_name.start - $table_name.end),$table_name.variation_type,$table_name.projects,$table_name._value from $table_name, $T2  where $T2.start<=$end and $T2.end>=$start and $T2.id=$table_name.rowid   and $table_name.variation_type="$type";});
    #	warn qq{select $table_name.variation_type,$table_name.projects,$table_name._value from $table_name, $T2  where $T2.start<=$end and $T2.end>=$start and $T2.id=$table_name.rowid   and $table_name.variation_type="$type";};
	# $self->prepare_cnv->execute($start,$end,$type);
    # $self->{prepare_cnv}->{$chr} = $self->dbh($chr)->prepare(qq{select $table_name.variation_type, $table_name.projects,$table_name.start,$table_name.end,$table_name._value from $table_name  where end>=? and start<=?  and variation_type=? and  (min(end,$end)-max(start,$start))>=$lmin and  (max(end,$end) - min($start,start)<=$lmax) limit 20;});
	 
	  $self->prepare_cnv($chr)->execute($start,$end,$type);

	my $x;
	my$nb;
	while (my @row = $self->prepare_cnv($chr)->fetchrow_array)     # retrieve one row
	{ 
	 	my $start1 = $row[2];
	 	my $end1 = $row[3];
	 
	 	my $identity = $self->getIdentityBetweenCNV($start,$end,$start1,$end1);
	 	next if $identity <  $seuil;

	 	$nb++;
		my $z = $self->decode($row[-1]);
		foreach my $zz (keys %$z){
			foreach my $zzz (keys %{$z->{$zz}}) {
				foreach my $zzzz (keys %{$z->{$zz}->{$zzz}}) {
					$x->{$zz}->{$zzz}->{$zzzz} = $start1."_".$end1."_".$identity;
				}
			}
			unless ($dejavu eq "all") {last if scalar(keys %$x) > $dejavu+1};
		}
		unless ($dejavu eq "all") {last if scalar(keys %$x) > $dejavu+1};
	 }
	 return $x;
}

sub getIdentityBetweenCNV {
	my ( $self, $start1, $end1, $start2, $end2) = @_;
	
	#retourne le recouvrement en % de la longueur du plus long des deux evenements 
	
	my $overlap = min( $end1, $end2 ) - max( $start1, $start2 );
	confess if abs( $start1 - $end1 ) ==0;
	my $overlap1 = $overlap / abs( $start1 - $end1 );
	my $overlap2 = $overlap / abs( $start2 - $end2 );

	return min($overlap1*100,$overlap2*100);
}

sub areSameCNV {
	my ( $self, $start1, $end1, $start2, $end2, $seuil) = @_;
	confess("\n\nERROR: seuil not found\n\n") unless $seuil;
	return 1 if ($self->getIdentityBetweenCNV($start1, $end1, $start2, $end2) >= $seuil);
	return;
}


sub put_dejavu {
	my ($self,$key1,$key2,$data,$start,$end) = @_;
	my $h ;
	$h->{$key2}->{data} = $data;
	$h->{$key2}->{start} = $start;
	$h->{$key2}->{end} = $end;
	$self->_put($key1,$h);
	
}

sub put_dejavu_bulk {
	my ($self,$key1,$hash,$debug) = @_;
	$self->_put_dejavu($key1,$hash,$debug);
	
}

sub get_bulk_dejavu {
	my ($self,$key1,$keys,$max) = @_;
	 my $table_name = $self->create_table($key1);
	 my $aray_ref;
	if ($keys){
	my $ids = join(",",map{$_=qq{"$_"} unless $_ =~/\"/}@$keys);
	 $aray_ref = $self->dbh($key1)->selectall_arrayref("select _key,projects from $table_name where _key in($ids) ");
	}
	else {
			 $aray_ref = $self->dbh($key1)->selectall_arrayref("select _key,projects,all from $table_name ;");
	}
	
	my $h= [];
	foreach my $a (@$aray_ref){
		next if $a->[1] <= $max;
		push(@$h,$a->[0]);
		#$h->{$a->[0]} ++;
		#$h->{$a->[0]}->{all} =  $a->[2];
	}
return $h;
}
sub get_bulk_ho {
	my ($self,$key1,$keys,$max) = @_;
	 my $table_name = $self->create_table($key1);
	 my $aray_ref;
	if ($keys){
	my $ids = join(",",map{$_=qq{"$_"} unless $_ =~/\"/}@$keys);
	 $aray_ref = $self->dbh($key1)->selectall_arrayref("select _key,ho from $table_name where _key in($ids) ");
	}
	else {
			 $aray_ref = $self->dbh($key1)->selectall_arrayref("select _key,ho,all from $table_name ;");
	}
	
	my $h;
	foreach my $a (@$aray_ref){
		next if $a->[1] <= $max;
		push(@$h,$a->[0]);
		#$h->{$a->[0]} ++;
		#$h->{$a->[0]}->{all} =  $a->[2];
	}
return $h;
}
sub get_hohe_nbprojects {
	my ($self,$key1,$id) = @_;
	my $table_name = $self->create_table($key1);
 	my $aray_ref = $self->dbh($key1)->selectall_arrayref("select ho,projects from $table_name where _key ='$id' ");
	return ($aray_ref->[0]);
}
sub DESTROY {
	my ($self) = @_;
	unlink  $self->dir."/".$self->lmdb_extension."/lock.mdb"  if -e $self->dir."/".$self->lmdb_extension."/lock.mdb" ;
	#$self->lmdb->close();

	$self->close();
   # print "In DESTROY\n";
}

1;
