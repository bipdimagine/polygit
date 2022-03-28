package GenBoNoSqlAnnotation;



# ABSTRACT: An embedded, NoSQL SQLite database with SQL indexing
use Moose;
use MooseX::Method::Signatures;
use strict;
use warnings;
use Data::Dumper;
use Storable qw/thaw freeze/;
extends "GenBoNoSqlText";

has extension =>(
	is		=> 'rw',
default => sub {
		return "annot.search";
		
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
	$self->dbh($key1)->do("CREATE VIRTUAL TABLE  if not exists $table_name USING fts4(_key ,_value blob ,notindexed='_value', tokenize=unicode61 'tokenchars=-_') ;")  or die $DBI::errstr;
#	$self->dbh($key1)->do("CREATE VIRTUAL TABLE  if not exists $table_name USING fts3(_key VARCHAR(250) NOT NULL ,_value blob ) ;")  or die $DBI::errstr;
	
	return 	$self->{table}->{$key1} ;
}

sub getIds {
	
	my ($self,$id) = @_;
	#$id =~s/xxx/_/g;
	return $id;
}



sub restore_id {
	my ($self,$id) = @_;
	#$id =~s/xxx/_/g;
	return $id;
}

sub change_id {
	my ($self,$id) = @_;
#	$id =~s/_/xxx/g;
	return $id;
}
1;