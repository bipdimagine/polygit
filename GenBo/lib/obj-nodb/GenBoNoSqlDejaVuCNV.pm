package GenBoNoSqlDejaVuCNV;
# ABSTRACT: An embedded, NoSQL SQLite database with SQL indexing
use Moo;

use strict;
use warnings;
use Data::Dumper;
extends "GenBoNoSqlDejaVu";



sub create_table {
	my ($self,$key1) = @_;
	return $self->{table}->{$key1}   if (exists $self->{table}->{$key1});
	 my $table_name =$self->name_table_data;
	 $self->{table}->{$key1} = $table_name;
	 return  $self->{table}->{$key1} unless  $self->write();
	$self->dbh($key1)->do("DROP TABLE IF EXISTS $table_name")  or die $DBI::errstr  if $self->mode eq 'c'  ;
	$self->dbh($key1)->do("CREATE TABLE if not exists $table_name (_key VARCHAR(250),ho INTEGER, projects INTEGER, _value BLOB,start INTEGER,end INTEGER, variation_type VARCHAR(3), length INTEGER) ")  or die $DBI::errstr;;
	return 	$self->{table}->{$key1} ;
}

sub sth_insert {
	my ($self,$chr) = @_;
	return $self->{dbh_insert}->{$chr} if exists  $self->{dbh_insert}->{$chr};
	
	 $self->{dbh_insert}->{$chr}= $self->dbh($chr)->prepare('insert into  __DATA__(_key,_value,start,end,length,variation_type)  values(?,?,?,?,?,?) ;') or die $DBI::errstr;
	return $self->{dbh_insert}->{$chr};
}

sub insert_cnv {
		my ($self,$id,$value) = @_;
		my ( $type, $c1, $start, $end ) = split( /_/,$id);
		$self->sth_insert($c1)->execute($id,$self->encode($value),$start,$end,abs($start-$end)+1,$type);
		return;
}


sub create_index {
	my ($self,$chr) = @_;
	$self->dbh($chr)->do(qq{CREATE UNIQUE INDEX if not exists _key_idx  on __DATA__ (_key);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _start_idx  on __DATA__ (start);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _end_idx  on __DATA__ (end);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _type_idxle  on __DATA__ (length);});
	$self->dbh($chr)->do(qq{CREATE  INDEX if not exists _type_idx1  on __DATA__ (variation_type);});
	
}


###################
# READ DATA
###################

sub prepare_cnv {
	my ($self,$chr) = @_;
	confess() unless $chr;
	
	return $self->{prepare_cnv}->{$chr} if exists $self->{prepare_cnv}->{$chr};
	 my $table_name = $self->create_table($chr);

	$self->{prepare_cnv}->{$chr} = $self->dbh($chr)->prepare(qq{select $table_name.variation_type, $table_name.projects,$table_name.start,$table_name.end,$table_name._value from $table_name  where end >= ? and start<=?  and variation_type=?  and length between  ? and  ?});
	warn "$chr";
	return $self->{prepare_cnv}->{$chr};
}




sub get_cnv {
	my ($self,$chr,$start,$end,$type,$dejavu,$seuil) = @_;
	
	 my $l = (100 - $seuil)/100;
	 my $len = abs($start - $end);
	 	$l = abs($l*$len)-1;
	  my $table_name = $self->create_table($chr);
	  my $min = $len - $l;
	  my $max = $len + $l;
	  $self->prepare_cnv($chr)->execute($start,$end,$type,($len - $l),($len + $l));
		my $x;
		my$nb = 0;
	while (my @row = $self->prepare_cnv($chr)->fetchrow_array)     # retrieve one row
	{ 
		
	 	my $start1 = $row[2];
	 	my $end1 = $row[3];
	 $nb++;
	 	my $identity = $self->getIdentityBetweenCNV($start,$end,$start1,$end1);
	 	next if $identity <  $seuil;

	 	
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

 #
sub get_cnv_for_manue {
		my ($self,$chr,$start,$end,$type,$dejavu,$seuil,$TheProjectName) = @_;
	my $nbproject=0;
	my $nbpatient=0;
	my $nbDJV_Wisecondor=0;
	my $nbDJV_Canvas=0;
	my $nbDJV_Manta=0;
	
	 my $l = (100 - $seuil)/100;
	 my $len = abs($start - $end);
	 	$l = abs($l*$len)-1;
	  my $table_name = $self->create_table($chr);
	  my $min = $len - $l;
	  my $max = $len + $l;
	  my $scorecaller_evt=0;
	  $self->prepare_cnv($chr)->execute($start,$end,$type,($len - $l),($len + $l));
		my $nb = 0;
		my %dj;
	while (my @row = $self->prepare_cnv($chr)->fetchrow_array)     # retrieve one row
	{ 
	 	my $start1 = $row[2];
	 	my $end1 = $row[3];
	 	my $identity = $self->getIdentityBetweenCNV($start,$end,$start1,$end1);
	 	next if $identity <  $seuil;
	 	
		my $hashdv = $self->decode($row[-1]);
		
		foreach my $project_name (keys %{$hashdv})
	{
			next if ($project_name eq $TheProjectName);
			$nbproject++;
				unless ($dejavu eq "all") {last if $nbproject  > $dejavu+1};
			foreach my $pname (keys %{$hashdv->{$project_name}})
			{
					my $scorecaller=0;
					my $c;
					$nbpatient++;
					my $res ="";
					foreach my $caller (keys %{$hashdv->{$project_name}->{$pname}})
					{
						
							if ($caller eq "wisecondor") {$scorecaller +=4; $c="w";$nbDJV_Wisecondor++;}
							if ($caller eq "canvas") {$scorecaller +=2; $c="c";$nbDJV_Canvas++;}
							if ($caller eq "manta") {$scorecaller +=1; $c="m";$nbDJV_Manta++;}
							$res .= int($identity)."%_".$c." ";
					}
					$scorecaller_evt += $scorecaller;
					
					my $etiq = $project_name.":".$pname.":".$res;
					$dj{$etiq} ++;
				}
		}
			unless ($dejavu eq "all") {last if $nbproject  > $dejavu+1};
	}
	
	my $theliste;
	$theliste = join(",",sort keys %dj);
	$scorecaller_evt /= $nbpatient if $nbpatient;
	return ($nbproject,$nbpatient,$nbDJV_Wisecondor,$nbDJV_Canvas,$nbDJV_Manta,$scorecaller_evt,$theliste);
	
}


1;