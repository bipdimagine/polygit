package GenBoNoSqlPosition;



# ABSTRACT: An embedded, NoSQL SQLite database with SQL indexing
use Moose;
use MooseX::Method::Signatures;
use strict;
use warnings;
use Data::Dumper;
extends "GenBoNoSql";
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;


has extension =>(
	is		=> 'rw',
default => sub {
		return "lite.pos";
		
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
	$self->dbh($key1)->do("CREATE  TABLE  if not exists __VALUE__ ( _key TEXT ,_id INTEGER primary key)")  or die $DBI::errstr;
	$self->dbh($key1)->do("CREATE VIRTUAL TABLE  if not exists $table_name USING rtree( _key ,_start,_end)")  or die $DBI::errstr;
#	$self->dbh($key1)->do("CREATE VIRTUAL TABLE  if not exists $table_name USING fts3(_key VARCHAR(250) NOT NULL ,_value blob ) ;")  or die $DBI::errstr;
	
	return 	$self->{table}->{$key1} ;
}

sub delete_bulk {
	my ($self,$key1,$key2) = @_;
	my $id = $self->change_id($key1);
	my $table_name = $self->create_table($id);

	 $self->dbh($id)->do(qq{delete  from $table_name where _key match "$key2"});
	 

	 return 1;
}
sub index {
	my ($self,$key1) = @_;
	$self->{index}->{$key1} = 0 unless exists $self->{index}->{$key1};
	$self->{index}->{$key1} ++;
	return $self->{index}->{$key1};
}
sub sth_insert_cached {
	my ($self,$chr) = @_;
	return $self->{sth_insert}->{$chr} if exists $self->{sth_insert}->{$chr};
	 my $table = $self->create_table($chr);
	$self->{sth_insert}->{$chr} = $self->dbh($chr)->prepare("insert  into $table (_key,_start,_end) values(?,?,?)");
	return $self->{sth_insert}->{$chr};
}
sub sth_insert_cached_value {
	my ($self,$chr) = @_;
	return $self->{sth_insert_cached_value}->{$chr} if exists $self->{sth_insert_cached_value}->{$chr};
	 my $table = $self->create_table($chr);
	$self->{sth_insert_cached_value}->{$chr} = $self->dbh($chr)->prepare("insert  into __VALUE__ (_key,_id) values(?,?)");
	return $self->{sth_insert_cached_value}->{$chr};
}

sub _put {
	my ($self,$chr,$hash) = @_;
	$self->dbh($chr)->begin_work;
	  foreach my $key (keys %$hash){
	  			my $index = $self->index($chr);
	  			$self->sth_insert_cached_value($chr)->execute($key,$index) ;
	  			$self->sth_insert_cached($chr)->execute($index,$hash->{$key}->{start},$hash->{$key}->{end}) ;
	  
	 }
	 
	 $self->dbh($chr)->commit;
 return 1;
}


sub put {
	my ($self,$key1,$key2,$start,$end) = @_;
	my $h ;
	$h->{$key2}->{start} = $start;
	$h->{$key2}->{end} = $end;
	$self->_put($key1,$h);
	
}

sub put_bulk {
	my ($self,$key1,$hash) = @_;
	$self->_put($key1,$hash);
	
}

1;
