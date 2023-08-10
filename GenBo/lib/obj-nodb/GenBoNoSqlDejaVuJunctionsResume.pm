package GenBoNoSqlDejaVuJunctionsResume;
# ABSTRACT: An embedded, NoSQL SQLite database with SQL indexing
use Moo;

use strict;
use warnings;
use List::Util qw[min max];
use Data::Dumper;
use Carp;
extends "GenBoNoSqlDejaVuJunctions";


has extension =>(
	is		=> 'rw',
		default => sub {
		return "dejavu_junctions.resume.lite";
	}
);



sub create_table {
	my ($self,$key1) = @_;
	return $self->{table}->{$key1}   if (exists $self->{table}->{$key1});
	my $table_name =$self->name_table_data;
	$self->{table}->{$key1} = $table_name;
	return  $self->{table}->{$key1} unless  $self->write();
	$self->dbh($key1)->do("DROP TABLE IF EXISTS $table_name")  or die $DBI::errstr  if $self->mode eq 'c'  ;
	$self->dbh($key1)->do("CREATE TABLE if not exists $table_name (_key VARCHAR(250),patients_ratios10 INTEGER,patients_ratios20 INTEGER, patients INTEGER, start INTEGER, end INTEGER, details BLOB, details_ratios10 BLOB, details_ratios20 BLOB)")  or die $DBI::errstr;;
	$self->dbh($key1)->do("CREATE UNIQUE INDEX if not exists _key_idx  on $table_name (_key); ")  or die $DBI::errstr;;
	return 	$self->{table}->{$key1} ;
}

sub prepare_junctions {
	my ($self,$chr) = @_;
	return $self->{prepare_junctions}->{$chr} if exists $self->{prepare_junctions}->{$chr};
	my $table_name = $self->create_table($chr);
	$self->{prepare_junctions}->{$chr} = $self->dbh($chr)->prepare(qq{select $table_name.start,$table_name.end,$table_name.patients,$table_name.details,$table_name.patients_ratios10,$table_name.details_ratios10,$table_name.patients_ratios20,$table_name.details_ratios20,$table_name._key  from $table_name  WHERE start=? OR end=? });
	return $self->{prepare_junctions}->{$chr};
}

sub getIdentityBetweenCNV {
	my ( $self, $start1, $end1, $start2, $end2) = @_;
	return 0 unless $start2;
	return 0 unless $end2;
	#retourne le recouvrement en % de la longueur du plus long des deux evenements 
	my $overlap = min( $end1, $end2 ) - max( $start1, $start2 );
	confess() if (abs( $start1 - $end1 ) ==0);
	my $overlap1 = $overlap / abs( $start1 - $end1 );
	my $overlap2 = $overlap / abs( $start2 - $end2 );
	return min($overlap1*100,$overlap2*100);
}

sub get_position {
	my ($self,$chr,$start,$end) = @_;
	my $table_name = $self->create_table($chr);
	my $T2 = $table_name."POSITION";
	$self->prepare_junctions($chr)->execute($start,$end);
	my $x;
	while (my @row = $self->prepare_junctions($chr)->fetchrow_array) {
	 	my $start1 = $row[0];
	 	my $end1 = $row[1];
	 	my $nbpat = $row[2];
	 	my $details = $row[3];
	 	my $nbpat_r10 = $row[4];
	 	my $details_r10 = $row[5];
	 	my $nbpat_r20 = $row[6];
	 	my $details_r20 = $row[7];
	 	$start =~ s/ //g;
	 	$end =~ s/ //g;
	 	$start1 =~ s/ //g;
	 	$end1 =~ s/ //g;
		next if not $start == $start1;
		next if not $end == $end1;
	 	$x->{$start1.'-'.$end1}->{nbpat} = $nbpat;
	 	$x->{$start1.'-'.$end1}->{details} = $details;
	 	$x->{$start1.'-'.$end1}->{nbpat_r10} = $nbpat_r10;
		$x->{$start1.'-'.$end1}->{details_r10} = $details_r10;
	 	$x->{$start1.'-'.$end1}->{nbpat_r20} = $nbpat_r20;
	 	$x->{$start1.'-'.$end1}->{details_r20} = $details_r20;
	}
	return $x;
}

sub get_junctions_generic {
	my ($self,$type,$chr,$start,$end,$seuil) = @_;
	my $table_name = $self->create_table($chr);
	my $T2 = $table_name."POSITION";
	$self->prepare_junctions($chr)->execute($start,$end);
	my $x;
	while (my @row = $self->prepare_junctions($chr)->fetchrow_array) {
	 	my $start1 = $row[0];
	 	my $end1 = $row[1];
	 	my $nbpat = $row[2];
	 	my $details = $row[3];
	 	my $nbpat_r10 = $row[4];
	 	my $details_r10 = $row[5];
	 	my $nbpat_r20 = $row[6];
	 	my $details_r20 = $row[7];
	 	
		next unless $start1;
		next unless $end1;
	 	$start =~ s/ //g;
	 	$end =~ s/ //g;
	 	$start1 =~ s/ //g;
	 	$end1 =~ s/ //g;
	 	next if ($start1 == $end1);
	 	
	 	next unless $start1;
	 	next unless $end1;
	 	
	 	my $identity;
	 	if ($start == $start1 and $end == $end1) { $identity = 100; }
	 	if ($start == ($start1 - 1) and $end == ($end1 + 1)) { $identity = 100; }
	 	else { $identity = $self->getIdentityBetweenCNV($start,$end,$start1,$end1); }
	 	next if $identity <  $seuil;
	 	
		if ($type eq 'ratio10') { 
		 	$x->{$start1.'-'.$end1}->{nbpat} = $nbpat_r10;
		 	$x->{$start1.'-'.$end1}->{details} = $details_r10;
		}
		elsif ($type eq 'ratio20') { 
		 	$x->{$start1.'-'.$end1}->{nbpat} = $nbpat_r20;
		 	$x->{$start1.'-'.$end1}->{details} = $details_r20;
		}
		else {
		 	$x->{$start1.'-'.$end1}->{nbpat} = $nbpat;
		 	$x->{$start1.'-'.$end1}->{details} = $details;
		}
		
	}
	return $x;
	
}

sub get_junctions {
	my ($self,$chr,$start,$end,$seuil) = @_;
	return $self->get_junctions_generic('all',$chr,$start,$end,$seuil);
}

sub get_junctions_ratio10 {
	my ($self,$chr,$start,$end,$seuil) = @_;
	return $self->get_junctions_generic('ratio10',$chr,$start,$end,$seuil);
}

sub get_junctions_ratio20 {
	my ($self,$chr,$start,$end,$seuil) = @_;
	return $self->get_junctions_generic('ratio20',$chr,$start,$end,$seuil);
}

sub get_nb_junctions_generic {
	my ($self,$hres,$patient_name) = @_;
	my $hp;
	foreach my $id (keys %{$hres}){
		next if $hres->{$id}->{nbpat} == 0;
		my @lPat = split(';', $hres->{$id}->{details});
		foreach my $p (@lPat) {
			if ($patient_name) {
				next if ($patient_name eq $p);
				if ($p =~ /[0-9]{2}_[0-9]{4}_/) {
					my @tmp = split('_', $p);
					shift(@tmp); shift(@tmp);
					$p = join('_', @tmp);
					next if ($patient_name eq $p);
				}
			}
			$hp->{$p} = undef;
		}
	}
	return scalar(keys %{$hp});
}

sub get_nb_junctions {
	my ($self,$chr,$start,$end,$seuil,$patient_name) = @_;
	my $h = $self->get_junctions($chr,$start,$end,$seuil);
	return $self->get_nb_junctions_generic($h,$patient_name);
}

sub get_nb_junctions_ratio10 {
	my ($self,$chr,$start,$end,$seuil,$patient_name) = @_;
	my $h = $self->get_junctions_ratio10($chr,$start,$end,$seuil);
	return $self->get_nb_junctions_generic($h,$patient_name);
}

sub get_nb_junctions_ratio20 {
	my ($self,$chr,$start,$end,$seuil,$patient_name) = @_;
	my $h = $self->get_junctions_ratio20($chr,$start,$end,$seuil);
	return $self->get_nb_junctions_generic($h,$patient_name);
}

1;