package GenBoNoSqlDejaVuSV;
# ABSTRACT: An embedded, NoSQL SQLite database with SQL indexing
use Moose;
use MooseX::Method::Signatures;
use strict;
use warnings;
use Data::Dumper;
extends "GenBoNoSqlDejaVu";

has toto =>(
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
	$self->dbh($key1)->do("CREATE TABLE if not exists $table_name (_key VARCHAR(250),ho INTEGER, projects INTEGER, _value BLOB,start INTEGER,end INTEGER, variation_type VARCHAR(2), length INTEGER, caller_infos VARCHAR(100), patients INTEGER)")  or die $DBI::errstr;;
	warn "coucou";
	#$self->dbh($key1)->do("CREATE UNIQUE INDEX if not exists _key_idx  on $table_name (_key); ")  or die $DBI::errstr;;
#	my $T2 = 	$table_name."POSITION";
#	$self->dbh($key1)->do("CREATE VIRTUAL TABLE  if not exists $T2  USING rtree( id ,start,end)")  or die $DBI::errstr;
#	$self->dbh($key1)->do("CREATE VIRTUAL TABLE  if not exists $table_name USING fts3(_key VARCHAR(250) NOT NULL ,_value blob ) ;")  or die $DBI::errstr;
	
	return 	$self->{table}->{$key1} ;
}





sub prepare_cnv {
	my ($self,$chr) = @_;

	return $self->{prepare_cnv}->{$chr} if exists $self->{prepare_cnv}->{$chr};
		 my $table_name = $self->create_table($chr);
	warn $self->dir." ".$chr;
	$self->{prepare_cnv}->{$chr} = $self->dbh($chr)->prepare(qq{select $table_name.variation_type, $table_name.projects,$table_name.start,$table_name.end,$table_name._value from $table_name  where end>=? and start<=?  and variation_type=? and length<? and length>? and projects<?  });
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
	 my $length = abs($start-$end) +1;
	 my $sp = int($length*abs(100-$seuil));
	  $self->prepare_cnv($chr)->execute($start,$end,$type,$length+$sp,$length-$sp,$dejavu+1);
	
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

1;