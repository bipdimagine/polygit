package GenBoNoSqlDejaVuCNV;
# ABSTRACT: An embedded, NoSQL SQLite database with SQL indexing
use Moo;

use strict;
use warnings;
use Data::Dumper;
use Sereal qw(sereal_encode_with_object sereal_decode_with_object write_sereal_file read_sereal_file);
extends "GenBoNoSqlDejaVu";

sub encode {
	my ( $self, $value ) = @_;
	return sereal_encode_with_object( $self->sereal_encoder, $value );
}

has sereal_encoder => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;

	return Sereal::Encoder->new({compress=>Sereal::SRL_ZSTD});
	
	},
);

has sereal_decoder => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;

		return Sereal::Decoder->new({compress=>Sereal::SRL_ZSTD});
		return 0;
	},
);

sub decode {
	my ( $self, $value ) = @_;
	return unless $value;
	
	my $x;
	eval {
	$x = sereal_decode_with_object( $self->sereal_decoder, $value );
	
	};
if ($@){
	warn $value;
	confess();
}
return $x;
}



sub create_table {
	my ($self,$key1) = @_;
	return $self->{table}->{$key1}   if (exists $self->{table}->{$key1});
	 my $table_name =$self->name_table_data;
	 $self->{table}->{$key1} = $table_name;
	 return  $self->{table}->{$key1} unless  $self->write();
	$self->dbh($key1)->do("DROP TABLE IF EXISTS $table_name")  or die $DBI::errstr  if $self->mode eq 'c'  ;
	$self->dbh($key1)->do("CREATE TABLE if not exists $table_name (_key VARCHAR(250),ho INTEGER, projects INTEGER, _value BLOB, _value2 BLOB,dv_project INTEGER,dv_sample INTEGER,dv_sr INTEGER,dv_depth INTEGER,dv_cov INTEGER,start INTEGER,end INTEGER, variation_type VARCHAR(3), length INTEGER) ")  or die $DBI::errstr;;
	return 	$self->{table}->{$key1} ;
}

sub sth_insert {
	my ($self,$chr) = @_;
	return $self->{dbh_insert}->{$chr} if exists  $self->{dbh_insert}->{$chr};
	$self->create_table($chr); 
	 $self->{dbh_insert}->{$chr}= $self->dbh($chr)->prepare('insert into  __DATA__(_key,_value,_value2,start,end,length,variation_type,dv_project,dv_sample,dv_sr,dv_depth,dv_cov)  values(?,?,?,?,?,?,?,?,?,?,?,?) ;') or die $DBI::errstr;
	return $self->{dbh_insert}->{$chr};
}

sub insert_cnv {
		my ($self,$id,$value,$value2,$dv_project,$dv_sample,$dv_sr,$dv_depth,$dv_cov) = @_;
		my ( $type, $c1, $start, $end ) = split( /_/,$id);
		confess() unless $value2;
		$self->sth_insert($c1)->execute($id,$self->encode($value),$self->encode($value2),$start,$end,abs($start-$end)+1,$type,$dv_project,$dv_sample,$dv_sr,$dv_depth,$dv_cov);
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
	return $self->{prepare_cnv}->{$chr};
}




sub get_local_cnv {
	my ($self,$global_id,$seuil) = @_;
	my ($type,$chr,$start,$end) = split(/_/, $global_id);
		 my $l = (100 - $seuil)/100;
	 my $len = abs($start - $end);
	 	$l = abs($l*$len)-1;
	  my $table_name = $self->create_table($chr);
	  my $min = $len - $l;
	  my $max = $len + $l;
	  $self->prepare_cnv($chr)->execute($start,$end,$type,($len - $l),($len + $l));
	  	my $x;
		my$nb = 0;
	my $res;
	while (my @row = $self->prepare_cnv($chr)->fetchrow_array)     # retrieve one row
	{ 
		
	 	my $start1 = $row[2];
	 	my $end1 = $row[3];
	 	$nb++;
	 	my $identity = $self->getIdentityBetweenCNV($start,$end,$start1,$end1);
	 	next if $identity <  $seuil;
	 	my $z = $self->decode($row[-1]);
	 	foreach my $p  (keys %$z){
	 		map{$res->{$p}->{$_} ++ } keys %{$z->{$p}};
	 	}
	}
	return $res;
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



########±±±±#####
# get CNV local
################


######################
# for position
################

 sub get_cnv_project {
	my ($self,$type,$chr,$start,$end,$seuil,$pid) = @_;
	 my $l = (100 - $seuil)/100;
	 my $len = abs($start - $end);
	 $l = abs($l*$len)-1;
	  my $table_name = $self->create_table($chr);
	  my $min = $len - $l;
	  my $max = $len + $l;
	    $self->prepare_cnv_2($chr)->execute($start,$end,$type,($len - $l),($len + $l));
		my $nb = 0;
		my @dj;
	while (my $row = $self->prepare_cnv_2($chr)->fetchrow_hashref)     # retrieve one row
	{ 
	 	my $start1 = $row->{start};
	 	my $end1 = $row->{end};
	 	my $identity = $self->getIdentityBetweenCNV($start,$end,$start1,$end1);
	 	next if $identity <  $seuil;
	 	
	 	my $hashdv = $self->decode($row->{hash2});
	 	if($pid){
	 	next unless exists $hashdv->{$pid};
	 	}
	 	my $l = abs($start-$end);
	 	push(@dj,$self->decode($row->{hash}));
	}
	return \@dj;
	
 }





 #
 
 sub prepare_cnv_2 {
	my ($self,$chr) = @_;
	confess() unless $chr;
	
	return $self->{prepare_cnv2}->{$chr} if exists $self->{prepare_cnv2}->{$chr};
	 my $table_name = $self->create_table($chr);

	$self->{prepare_cnv2}->{$chr} = $self->dbh($chr)->prepare(qq{select 
		$table_name.variation_type as type , $table_name.projects as projects ,$table_name.start as start ,
		$table_name.end as end  ,$table_name._value2 as hash2 ,$table_name._value as hash
		from $table_name  where end >= ? and start<=?  and variation_type=?  and length between  ? and  ?});
	return $self->{prepare_cnv2}->{$chr};
}

 sub get_all_cnv {
 	my ($self,$chr,$dv,$patient) = @_;
 	  $self->prepare_select_global($chr)->execute($dv);
		my $nb = 0;
		my @dj;
	while (my @row = $self->prepare_select_global($chr)->fetchrow_array)     # retrieve one row
	{
		 my $hashdv = $self->decode($row[-1]);
		 next unless exists $hashdv->{$patient->id};
		 my $hash = $self->decode($row[0]);
		 push(@dj,$hash);
	}
	return \@dj;
	
 } 
 
 
  sub prepare_select_global {
	my ($self,$chr) = @_;
	confess() unless $chr;
	
	return $self->{prepare_select_global}->{$chr} if exists $self->{prepare_select_global}->{$chr};
	 my $table_name = $self->create_table($chr);
	
	$self->{prepare_select_global}->{$chr} = $self->dbh($chr)->prepare(qq{select 
		$table_name._value as hash  ,$table_name._value2 as hash 																				
		from $table_name where dv_sample < ? });
	return $self->{prepare_select_global}->{$chr};
}
 
 sub get_cnv_for_manue {
	my ($self,$global_id,$seuil,$TheProjectName) = @_;
	my ($type,$chr,$start,$end) = split(/_/, $global_id);
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
	  $self->prepare_cnv_2($chr)->execute($start,$end,$type,($len - $l),($len + $l));
		my $nb = 0;
		my @dj;
	while (my @row = $self->prepare_cnv_2($chr)->fetchrow_array)     # retrieve one row
	{ 
	 	my $start1 = $row[2];
	 	my $end1 = $row[3];
	 	my $identity = $self->getIdentityBetweenCNV($start,$end,$start1,$end1);
	 	next if $identity <  $seuil;
	 	$identity = int($identity);
		my $hashdv = $self->decode($row[-1]);
		delete $hashdv->{$TheProjectName};
		my $string ="";
		foreach my $pr (keys %$hashdv){
			#warn $hashdv->{$pr}->{string} ;
			$nbpatient += $hashdv->{$pr}->{patients};
			$nbproject ++;
			$string .= $hashdv->{$pr}->{string};
			$nbDJV_Wisecondor += $hashdv->{$pr}->{coverage} if exists $hashdv->{$pr}->{coverage};
			$nbDJV_Canvas += $hashdv->{$pr}->{depth} if exists  $hashdv->{$pr}->{depth};
			$nbDJV_Manta += $hashdv->{$pr}->{sr} if exists  $hashdv->{$pr}->{sr};
			
		}
		$string =~ s/!X!/$identity/g;
		push(@dj,$string) if $string;
	}
	my $theliste;
	$theliste = join(",",@dj);
	$scorecaller_evt /= $nbpatient if $nbpatient;
	#warn join (";",($nbproject,$nbpatient,$nbDJV_Wisecondor,$nbDJV_Canvas,$nbDJV_Manta,$scorecaller_evt,$theliste));
	#warn $self->dir();
	#warn "%%%%%%%%%%%%%%%%%%";
  #	my $magic_String_manue  = $global_id."+%NBPROJECT%;%LIST%;".$nbproject.";".$nbpatient.";".$theliste.";".$nbDJV_Wisecondor.";".$nbDJV_Canvas.";".$nbDJV_Manta.";".$flag;
	my $hash;
	$hash->{nb_patients} = $nbpatient;
	$hash->{nb_projects} = $nbproject;
	$hash->{caller_sr} =$nbDJV_Manta ;
	$hash->{caller_depth} = $nbDJV_Canvas ;
	$hash->{caller_coverage} =$nbDJV_Wisecondor ;
	$hash->{string} = $theliste;
	return  $hash;
	
}
sub get_cnv_for_manue_old {
	my ($self,$global_id,$dejavu,$seuil,$TheProjectName) = @_;
	my ($type,$chr,$start,$end) = split(/_/, $global_id);
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
	warn join (";",($nbproject,$nbpatient,$nbDJV_Wisecondor,$nbDJV_Canvas,$nbDJV_Manta,$scorecaller_evt,$theliste));
	warn $self->dir();
	return ($nbproject,$nbpatient,$nbDJV_Wisecondor,$nbDJV_Canvas,$nbDJV_Manta,$scorecaller_evt,$theliste);
	
}


1;