package GenBoNoSqlDejaVuSV;
# ABSTRACT: An embedded, NoSQL SQLite database with SQL indexing
use Moo;

use strict;
use warnings;
use Data::Dumper;
use Hash::Merge qw/ merge /;
extends "GenBoNoSqlDejaVu";

has toto =>(
	is		=> 'rw',
default => sub {
		return "lmdb_dejavu";
		
}
);
sub create_table {
	my ($self) = @_;
	my $key1 = "SV";
	return $self->{table}->{$key1}   if (exists $self->{table}->{$key1});
	 my $table_name =$self->name_table_data;
	 $self->{table}->{$key1} = $table_name;
	 return  $self->{table}->{$key1} unless  $self->write();
	$self->dbh($key1)->do("DROP TABLE IF EXISTS $table_name")  or die $DBI::errstr  if $self->mode eq 'c'  ;

	#$self->dbh($chr)->do("CREATE TABLE if not exists $table_name (_key TEXT PRIMARY KEY, _value BLOB) ;")  or die $DBI::errstr;;
	$self->dbh($key1)->do("CREATE TABLE if not exists $table_name (_key VARCHAR(250),ho INTEGER, _value BLOB,chromosome2 VARCHAR(2),chromosome1 VARCHAR(2) , start1 INTEGER, start2 INTEGER, variation_type VARCHAR(2), caller VARCHAR(20),nb_project INTEGER, nb_patient INTEGER)")  or die $DBI::errstr;;

	return 	$self->{table}->{$key1} ;
}

sub sth_insert {
	my ($self) = @_;
	return $self->{dbh_insert} if exists  $self->{dbh_insert};
	 $self->{dbh_insert}= $self->dbh("SV")->prepare(
		'insert into  __DATA__(_key,_value,chromosome1,start1,chromosome2,start2,caller,nb_project,nb_patient)  values(?,?,?,?,?,?,?,?,?) ;') or die $DBI::errstr;
	return $self->{dbh_insert};
}

sub insert_sv {
		my ($self,$id,$value,$caller) = @_;
		my ( $c1, $start1, $c2, $start2 ) = split( /_/,$id);
		my $new_hash ={};
		my $npa = 0;
		my $np =0;
		foreach my $p (keys %$value){
			push(@{$new_hash->{projects}},$p);
			$np ++;
			foreach my $pa (keys %{$value->{$p}}){
					$new_hash->{patients}->{$p."!".$pa} = $id;
					$npa ++;
			} 
			
		}
		$self->sth_insert->execute($id,$self->encode($new_hash),$self->change_chr($c1),$start1,$self->change_chr($c2),$start2,$caller,$np,$npa);
		return;
}


sub change_chr {
	my ($self,$chr) = @_;
	return "X" if $chr == 23;
	return "Y" if $chr == 24;
	return "MT" if $chr == 25;
	return $chr;
}

sub create_index {
	my ($self) = @_;
	my $chr = "SV";
	$self->dbh($chr)->do(qq{CREATE UNIQUE INDEX if not exists _key_idx  on __DATA__ (_key);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _start_idx  on __DATA__ (start1);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _end_idx  on __DATA__ (start2);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _type_idx  on __DATA__ (chromosome1,chromosome2);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _type_idx1  on __DATA__ (chromosome1);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _type_idx2  on __DATA__ (chromosome2);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _type_idx2  on __DATA__ (nb_patient);});
}


sub prepare_sv{
	my ($self,$chr) = @_;

	return $self->{prepare_cnv}->{$chr} if exists $self->{prepare_cnv}->{$chr};
	my $table_name = $self->create_table($chr);

	$self->{prepare_cnv}->{$chr} = $self->dbh($chr)->prepare(qq{select $table_name.chromosome1, $table_name.start1,$table_name.chromosome2,$table_name.start2,nb_patient,_value from $table_name
	  where chromosome1=? and chromosome2=? and start1 between ? and ?  and start2 between ? and ? order by nb_patient});
	
	return $self->{prepare_cnv}->{$chr};
}



sub get_sv {
		my ($self,$chr1,$start1,$chr2,$start2,$dejavu) = @_;
		
	 my $table_name = $self->create_table("SV");
	 
	  $self->prepare_sv("SV")->execute($chr1,$chr2,$start1-50,$start1+50,$start2-50,$start2+50);

	my $x = {};
	my $nb = 0;
	$x->{infos} = {};
	while (my @row = $self->prepare_cnv("SV")->fetchrow_array){ 
		 if ($row[4]  > $dejavu) {
		 	$x->{dv_patients}  = $row[4];
		 	last;
		 } 
		$x->{infos} = merge $x->{infos}, $self->decode($row[-1]);
		
		$x->{dv_patients} += keys %{$x->{infos}->{patients}};
		
		 if ($x->{dv_patients}  > $dejavu){
		 	last;
		 } 
		
	 }
	 return $x;
}

1;