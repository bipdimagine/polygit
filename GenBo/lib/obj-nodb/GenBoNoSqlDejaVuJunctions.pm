package GenBoNoSqlDejaVuJunctions;
# ABSTRACT: An embedded, NoSQL SQLite database with SQL indexing
use Moo;

use strict;
use warnings;
use Data::Dumper;
extends "GenBoNoSqlDejaVu";


has extension =>(
	is		=> 'rw',
		default => sub {
		return "dejavu_junctions.lite";
	}
);

sub create_table {
	my ($self,$key1) = @_;
	return $self->{table}->{$key1}   if (exists $self->{table}->{$key1});
	my $table_name =$self->name_table_data;
	$self->{table}->{$key1} = $table_name;
	return  $self->{table}->{$key1} unless  $self->write();
	$self->dbh($key1)->do("DROP TABLE IF EXISTS $table_name")  or die $DBI::errstr  if $self->mode eq 'c'  ;
	$self->dbh($key1)->do("CREATE TABLE if not exists $table_name (_key VARCHAR(250),ratios INTEGER, patients INTEGER, projects INTEGER, _value BLOB,start INTEGER,end INTEGER, variation_type VARCHAR(2), length INTEGER, gene_name VARCHAR(100))")  or die $DBI::errstr;;
	$self->dbh($key1)->do("CREATE UNIQUE INDEX if not exists _key_idx  on $table_name (_key); ")  or die $DBI::errstr;;
	return 	$self->{table}->{$key1} ;
}

sub prepare_junctions_resume {
	my ($self,$chr) = @_;
	return $self->{prepare_junctions_resume}->{$chr} if exists $self->{prepare_junctions_resume}->{$chr};
	my $table_name = $self->create_table($chr);
	$self->{prepare_junctions_resume}->{$chr} = $self->dbh($chr)->prepare(qq{select $table_name.variation_type,$table_name.start,$table_name.end,$table_name._key,$table_name.patients,$table_name.projects,$table_name.ratios,$table_name.gene_name  from $table_name  WHERE start=? OR end=?  and variation_type=?   });
	return $self->{prepare_junctions_resume}->{$chr};
}

sub prepare_junctions {
	my ($self,$chr) = @_;
	return $self->{prepare_junctions}->{$chr} if exists $self->{prepare_junctions}->{$chr};
	my $table_name = $self->create_table($chr);
	$self->{prepare_junctions}->{$chr} = $self->dbh($chr)->prepare(qq{select $table_name.variation_type,$table_name.start,$table_name.end,$table_name._value,$table_name.gene_name from $table_name  WHERE start=? OR end=?  and variation_type=?   });
	return $self->{prepare_junctions}->{$chr};
}

sub get_junction_resume {
	my ($self,$chr,$start,$end,$type,$dejavu,$seuil) = @_;
	my $table_name = $self->create_table($chr);
	my $T2 = $table_name."POSITION";
	$self->prepare_junctions_resume($chr)->execute($start,$end,$type);
	my (@lRes);
	while (my @row = $self->prepare_junctions_resume($chr)->fetchrow_array) { 
	 	my $start1 = $row[1];
	 	my $end1 = $row[2];
	 	
	 	my $identity = $self->getIdentityBetweenCNV($start,$end,$start1,$end1);
	 	next if $identity <  $seuil;
	 	
#	 	my $identity = 100;
#	 	if ($start != $start1 or $end != $end1) {
#	 		my $len1 = ($end-$start);
#	 		my $len2 = ($end1-$start1);
#	 		if ($len1 > $len2) { $identity = ($len2/$len1)*100;  }
#	 		if ($len2 < $len1) { $identity = ($len1/$len2)*100;  }
#	 	} 
#	 	next if $identity <  $seuil;
	 	
	 	#warn join(' ', @row).$identity;
	 	
	 	my $x;
		$x->{type} = $row[0];
		$x->{id} = $row[-5];
		$x->{patients} = $row[-4];
		$x->{projects} = $row[-3];
		$x->{ratios} = $row[-2];
		$x->{ensg} = $row[-1];
		my @lTmp = split('_', $x->{id});
		if ($lTmp[1] == $start and $lTmp[2] == $end) {
			$x->{same_as} = '100%';
		}
		else {		
			$x->{same_as} = int($identity).'%';
		}
		push(@lRes, $x);
		unless ($dejavu eq "all") {last if scalar(keys %$x) > $dejavu+1};
	}
	return \@lRes;
}

sub get_position {
	my ($self,$chr,$start,$end,$type,$dejavu) = @_;
	my $table_name = $self->create_table($chr);
	my $T2 = $table_name."POSITION";
	$self->prepare_junctions($chr)->execute($start,$end,$type);
	my $x;
	while (my @row = $self->prepare_junctions($chr)->fetchrow_array) {
	 	my $type = $row[0];
	 	my $start1 = $row[1];
	 	my $end1 = $row[2];
	 	my $h_details = $self->decode($row[3]);
	 	my $gene_id = $row[4];
	 	$start =~ s/ //g;
	 	$end =~ s/ //g;
	 	$start1 =~ s/ //g;
	 	$end1 =~ s/ //g;
		next if not $start == $start1;
		next if not $end == $end1;
	 	$x->{$start1.'-'.$end1.'-'.$type}->{type} = $type;
	 	$x->{$start1.'-'.$end1.'-'.$type}->{start} = $start1;
	 	$x->{$start1.'-'.$end1.'-'.$type}->{end} = $end1;
	 	$x->{$start1.'-'.$end1.'-'.$type}->{details} = $h_details;
	 	$x->{$start1.'-'.$end1.'-'.$type}->{gene_id} = $gene_id;
	}
	return $x;
}

sub get_junction {
	my ($self,$chr,$start,$end,$type,$dejavu,$seuil) = @_;
	my $table_name = $self->create_table($chr);
	my $T2 = $table_name."POSITION";
	$self->prepare_junctions($chr)->execute($start,$end,$type);
	my $x;
	while (my @row = $self->prepare_junctions($chr)->fetchrow_array) {
	 	my $start1 = $row[1];
	 	my $end1 = $row[2];
	 	my $ensg = $row[-1];
	 	
	 	my $identity = $self->getIdentityBetweenCNV($start,$end,$start1,$end1);
	 	next if $identity <  $seuil;
	 	
#	 	my $identity = 100;
#	 	if ($start != $start1 or $end != $end1) {
#	 		my $len1 = ($end-$start);
#	 		my $len2 = ($end1-$start1);
#	 		if ($len1 > $len2) { $identity = ($len2/$len1)*100;  }
#	 		if ($len2 < $len1) { $identity = ($len1/$len2)*100;  }
#	 	}  
#	 	next if $identity <  $seuil;
	 	
		my $z = $self->decode($row[-2]);
		foreach my $zz (keys %$z){
			foreach my $zzz (keys %{$z->{$zz}}) {
				my $identite = int($identity);
				$identite = 100 if ($start1 == $start and $end1 == $end);
				foreach my $zzzz (keys %{$z->{$zz}->{$zzz}}) {
					$x->{$identite}->{$zz}->{$zzz}->{$zzzz} = $z->{$zz}->{$zzz}->{$zzzz};
				}
				$x->{$identite}->{$zz}->{$zzz}->{start} = $start1;
				$x->{$identite}->{$zz}->{$zzz}->{end} = $end1;
				$x->{$identite}->{$zz}->{$zzz}->{identity} = $identite;
			}
			unless ($dejavu eq "all") {last if scalar(keys %$x) > $dejavu+1};
		}
		unless ($dejavu eq "all") {last if scalar(keys %$x) > $dejavu+1};
	 }
#	 die;
	 return $x;
}


1;