package GenBoNoSqlDejaVuSV_agregate;
# ABSTRACT: An embedded, NoSQL SQLite database with SQL indexing
use Moo;

use strict;
use warnings;
use Data::Dumper;
use Hash::Merge qw/ merge /;
use List::Util qw(min max);
extends "GenBoNoSqlDejaVu";


sub create_table {
	my ($self) = @_;
	my $key1 = "SVa";
	return $self->{table}->{$key1}   if (exists $self->{table}->{$key1});
	 my $table_name =$self->name_table_data;
	 $self->{table}->{$key1} = $table_name;
	 return  $self->{table}->{$key1} unless  $self->write();
	$self->dbh($key1)->do("DROP TABLE IF EXISTS $table_name")  or die $DBI::errstr  if $self->mode eq 'c'  ;

	#$self->dbh($chr)->do("CREATE TABLE if not exists $table_name (_key TEXT PRIMARY KEY, _value BLOB) ;")  or die $DBI::errstr;;
	$self->dbh($key1)->do("CREATE TABLE if not exists $table_name ( _value BLOB,chromosome2 VARCHAR(3),chromosome1 VARCHAR(3) , start1_min INTEGER,  start1_max INTEGER,start2_min INTEGER,start2_max INTEGER, patients INTEGER, projects INTEGER)")  or die $DBI::errstr;;

	return 	$self->{table}->{$key1} ;
}

sub sth_insert {
	my ($self) = @_;
	return $self->{dbh_insert} if exists  $self->{dbh_insert};
	 $self->{dbh_insert}= $self->dbh("SVa")->prepare(
		'insert into  __DATA__(_value,chromosome1,chromosome2,start1_min,start1_max,start2_min,start2_max,projects,patients)  values(?,?,?,?,?,?,?,?,?) ;') or die $DBI::errstr;
	return $self->{dbh_insert};
}

sub insert_sv {
		my ($self,$obj) = @_;
		$self->sth_insert->execute($self->encode($obj->{value}),$obj->{chromosome1},$obj->{chromosome2},$obj->{start1_min},$obj->{start1_max},$obj->{start2_min},$obj->{start2_max},$obj->{nb_patients},$obj->{nb_projects});
		return;
}


sub create_index {
	my ($self) = @_;
	my $chr = "SVa";
	$self->dbh($chr)->do(qq{CREATE UNIQUE INDEX if not exists _key_idx  on __DATA__ (chromosome1,chromosome2,start1_min,start1_max,start2_min,start2_max);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _start_idx1  on __DATA__ (start1_min);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _start_idx2  on __DATA__ (start2_min);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _start_idx3  on __DATA__ (start1_max);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _start_idx4  on __DATA__ (start2_max);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _type_idx  on __DATA__ (chromosome1,chromosome2);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _type_idx1  on __DATA__ (chromosome1);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _type_idx2  on __DATA__ (chromosome2);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _type_stat  on __DATA__ (nb_patients);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _type_stat  on __DATA__ (nb_projects);});
}


sub prepare_sv{
	my ($self,$chr) = @_;
	$chr = "SVa";
	return $self->{prepare_sv}->{$chr} if exists $self->{prepare_sv}->{$chr};
	my $table_name = $self->create_table($chr);
	
	warn q{select $table_name.chromosome1, $table_name.start1_min, $table_name.start1_max,$table_name.chromosome2,$table_name.start2_min, $table_name.start2_max,projects,patients,_value from $table_name
	  where chromosome1=? and chromosome2=? and ($table_name.start1_min between ? or ?  and start1_max between ? and ?) and (start2_min between ? and ?  or start2_max between ? and ?)};
		$self->{prepare_sv}->{$chr} = $self->dbh($chr)->prepare(qq{select $table_name.chromosome1, $table_name.start1_min, $table_name.start1_max,$table_name.chromosome2,$table_name.start2_min, $table_name.start2_max,projects,patients,_value from $table_name
	  where  chromosome1=? and chromosome2=? and start1_min < ? and start1_max > ? and start2_min < ? and start2_max > ? });
	
	return $self->{prepare_sv}->{$chr};
}



sub get_sv {
	my ($self,$chr1,$start1,$chr2,$start2,$dejavu,$seuil) = @_;
		
	 my $table_name = $self->create_table("SVa");
	  $self->prepare_sv("SVa")->execute($chr1,$chr2,$start1,$start1,$start2,$start2);
	my $x;
	my$nb;
	my $h;
	
	 
	$h->{infos} = {};
	my $zz = 0;
	my @start1;
	my @start2;
	while (my @row = $self->prepare_sv("SVa")->fetchrow_array){ 
		next;
	$h->{chromosome1} = $row[0];
	push(@start1,$row[1]);
	push(@start1,$row[2]);
	push(@start2,$row[4]);
	push(@start2,$row[5]);
	
	$h->{chromosome2} = $row[3];

	$h->{patients} += $row[7];
	$h->{projects} += $row[6];
	$h->{infos} = merge $h->{infos}, $self->decode($row[8]);
	$h->{zz} ++;
	
	 }
	 $h->{start1_min} = min(@start1);
	 $h->{start1_max} = max(@start1);
	 $h->{start2_min} = min(@start2);
	 $h->{start2_max} = max(@start2);
	 return $h;
}

1;