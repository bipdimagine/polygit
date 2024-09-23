package GenBoNoSqlDejaVuSV;
# ABSTRACT: An embedded, NoSQL SQLite database with SQL indexing
use Moo;

use strict;
use warnings;
use Data::Dumper;
use Hash::Merge qw/ merge /;
use Carp;
extends "GenBoNoSqlDejaVu";


sub create_table {
	my ($self) = @_;
	my $key1 = "SV";
	return $self->{table}->{$key1}   if (exists $self->{table}->{$key1});
	 my $table_name =$self->name_table_data;
	 $self->{table}->{$key1} = $table_name;
	 return  $self->{table}->{$key1} unless  $self->write();
	$self->dbh($key1)->do("DROP TABLE IF EXISTS $table_name")  or die $DBI::errstr  if $self->mode eq 'c'  ;

	#$self->dbh($chr)->do("CREATE TABLE if not exists $table_name (_key TEXT PRIMARY KEY, _value BLOB) ;")  or die $DBI::errstr;;
	$self->dbh($key1)->do("CREATE TABLE if not exists $table_name (_key VARCHAR(250), _value BLOB,chromosome2 VARCHAR(2),chromosome1 VARCHAR(2) , start1 INTEGER, start2 INTEGER, variation_type VARCHAR(3), nb_project INTEGER, nb_patient INTEGER)")  or die $DBI::errstr;;
	$self->create_index();
	return 	$self->{table}->{$key1} ;
}

sub sth_insert {
	my ($self) = @_;
	return $self->{dbh_insert} if exists  $self->{dbh_insert};
	$self->create_table();
	 $self->{dbh_insert}= $self->dbh("SV")->prepare(
		'insert into  __DATA__(_key,_value,chromosome1,start1,chromosome2,start2,nb_project,nb_patient)  values(?,?,?,?,?,?,?,?) ;') or die $DBI::errstr;
	return $self->{dbh_insert};
}

sub insert_sv {
		my ($self,$id,$value,$caller) = @_;
		my ( $c1, $start1, $c2, $start2 ) = split( /_/,$id);
		my $new_hash ={};
		my $npa = 0;
		my $np =0;
		foreach my $p (keys %$value){
			push(@{$new_hash->{$p}},keys %{$value->{$p}});
			$np ++;
			foreach my $pa (keys %{$value->{$p}}){
					$npa ++;
			} 
		}
		$self->sth_insert->execute($id,$self->encode($new_hash),$c1,$start1,$c2,$start2,$np,$npa);
		return;
}
sub insert_local_sv {
		my ($self,$sv,$npa) = @_;
		my $new_hash ={};
		$self->sth_insert->execute($sv->{ID},$self->encode($sv),$sv->{CHROM1},$sv->{POS1},$sv->{CHROM2},$sv->{POS2},1,$npa);
		return;
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

	return $self->{prepare_sv}->{$chr} if exists $self->{prepare_sv}->{$chr};
	my $table_name = $self->create_table($chr);

	$self->{prepare_sv}->{$chr} = $self->dbh($chr)->prepare(qq{select start1 as start1, start2 as start2 , nb_project as nb_project,nb_patient as nb_patient,_value as hash from $table_name
	  where chromosome1=? and chromosome2=? and start1 between ? and ?  and start2 between ? and ? order by nb_patient});
	
	return $self->{prepare_sv}->{$chr};
}

sub get_sv_dejavu {
	my ($self,$sv,$limit,$project_name) = @_;
	die() unless $project_name;
	confess() unless $limit;
	return $self->get_sv($sv->{CHROM1},$sv->{POS1},$sv->{CHROM2},$sv->{POS2},$limit,$project_name);
}

sub get_sv {
	my ($self,$chr1,$start1,$chr2,$start2,$limit,$project_name) = @_;
		
	 my $table_name = $self->create_table("SV");
	  $self->prepare_sv("SV")->execute($chr1,$chr2,$start1-$limit,$start1+$limit,$start2-$limit,$start2+$limit);
	my $x = {};
	my $nb = 0;
	$x->{infos} = {};
	$x->{dv_patients} = 0;
	$x->{dv_projects} = 0;
	while (my $row = $self->prepare_sv("SV")->fetchrow_hashref){ 
		
		my $hashdv = $self->decode($row->{hash});
		 
		$hashdv->{identity1}  = abs($row->{start1} - $start1);
		
		$hashdv->{identity2} = abs($row->{start2} - $start2);
		my $nbp = 0;
		my $nbproject = 0;
		
		if (exists $hashdv->{$project_name}){
			$nbp = scalar(@{$hashdv->{$project_name}});
			delete  $hashdv->{project}->{$project_name};
			$nbproject =1;
		}
		delete  $hashdv->{project};
		
	#	warn Dumper $hashdv if $hashdv;
	#	die() if $hashdv;
		$x->{infos} = merge $x->{infos},$hashdv;
		$x->{dv_patients} = $row->{nb_patient} - $nbp;
		$x->{dv_projects} = $row->{nb_project} - $nbproject;
		 
		
	 }
	 return $x;
}








  sub prepare_select_global {
	my ($self) = @_;
		
	return $self->{prepare_select_global} if exists $self->{prepare_select_global};
	 my $table_name = $self->create_table();
	
	$self->{prepare_select_global} = $self->dbh("SV")->prepare(qq{select 
		$table_name._value as hash  																				
		from $table_name where nb_patient < ? });
	return $self->{prepare_select_global};
}

sub get_sv_by_position {
	my ($self,$dv,$patient) = @_;
	 my $table_name = $self->create_table("SV");
	  $self->prepare_select_global()->execute($dv);
		my $nb = 0;
		my @dj;
	while (my @row = $self->prepare_select_global()->fetchrow_array)     # retrieve one row
	{
		 my $hashdv = $self->decode($row[0]);
		 next unless exists $hashdv->{PATIENTS}->{$patient->id};
		 push(@dj,$hashdv);
	}
	return \@dj;
}
sub get_all_sv {
	my ($self,$dv,$patient) = @_;
	 my $table_name = $self->create_table("SV");
	  $self->prepare_select_global()->execute($dv);
		my $nb = 0;
		my @dj;
	while (my @row = $self->prepare_select_global()->fetchrow_array)     # retrieve one row
	{
		 my $hashdv = $self->decode($row[0]);
		 next unless exists $hashdv->{PATIENTS}->{$patient->id};
		 push(@dj,$hashdv);
	}
	return \@dj;
}

1;