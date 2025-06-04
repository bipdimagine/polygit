package GenBoDuckDejaVuSv;


use Moo;

use strict;
use warnings;
use Data::Dumper;
use List::Util qw[min max];
	
has project =>(
	is		=> 'rw',
	required=>1,
);

has version =>(
	is		=> 'ro',
	lazy=>1,
	default => sub {
		my ($self)= @_;
		return 38 if $self->project->current_genome_version =~/38/ ;
		return 19;
}	
);

has dir =>(
	is		=> 'ro',
	lazy=>1,
default => sub {
		my ($self)= @_;
		return $self->project->buffer->config_path("dejavu","sv");;
}	
);


has colstart =>(
	is		=> 'ro',
	lazy=>1,
	default => sub {
		my ($self)= @_;
		return "start".$self->version;
}	
);
has colend =>(
	is		=> 'ro',
	lazy=>1,
	default => sub {
		my ($self)= @_;
		return "end".$self->version;
}
);

has colchr1 =>(
	is		=> 'ro',
	lazy=>1,
	default => sub {
		my ($self)= @_;
		return "chr1_".$self->version;
}
);
has colchr2 =>(
	is		=> 'ro',
	lazy=>1,
	default => sub {
		my ($self)= @_;
		return "chr2_".$self->version;
}
);
has colpos1 =>(
	is		=> 'ro',
	lazy=>1,
	default => sub {
		my ($self)= @_;
		return "position1_".$self->version;
}
);
has colpos2 =>(
	is		=> 'ro',
	lazy=>1,
	default => sub {
		my ($self)= @_;
		return "position2_".$self->version;
}
);
has dbh =>  (
	is		=> 'rw',
	lazy=>1,
default => sub {
	my $dbh =   DBI->connect("dbi:ODBC:Driver=DuckDB;Database=:memory:", "", "", { RaiseError => 1 , AutoCommit => 1});
	return $dbh;
});


sub create_table_project {
	my ($self) =@_;
	return $self->{create_table_project} if exists  $self->{create_table_project};
	my $colpos1 = $self->colpos1;
	my $colpos2 = $self->colpos2;
	my $colchr1 = $self->colchr1;
	my $colchr2 = $self->colchr2;
	my $dir = $self->dir;
	my $file = $self->project->name.".parquet";
	my $parquet_file = $dir."$file";
	die($parquet_file) unless -e $dir."/$file";
	my $query = qq{CREATE TABLE sv_call_project AS
                           SELECT project,type,patient,$colchr1 as chrom1,$colchr2 as chrom2,$colpos1 as pos1 ,$colpos2 as pos2,caller,
                           sr1,sr2,pr1,pr2,gq,qual
                           FROM '$parquet_file'
                           WHERE $colpos1 > -1  order by type,$colchr1,$colchr2,$colpos1,$colpos2;
	};
	$self->dbh->do($query);
	$self->{create_table_project} = "sv_call_project";
	return $self->{create_table_project};
}





has sth_query_project =>(
	is		=> 'rw',
	lazy=>1,
default => sub {
	my ($self) =@_;
	my $table = $self->create_table_project();
	my $project_id = $self->project->id;
	 my $sql = qq{
	 SELECT *
		FROM $table 
		WHERE project=? and patient = ? ;
	 };
	 my $sth = $self->dbh->prepare($sql);
	 return $sth;
}
);

sub get_sv_project {
	my ($self,$patient ) = @_;
	
	$self->sth_query_project->execute($self->project->id,$patient->id);
	#$self->sth_query2->execute($type,$chr, $start, $end, $minl, $maxl);
	# Récupération des résultats dans un tableau de hash
		my @results;
	
		while (my $row = $self->sth_query_project->fetchrow_hashref) {
			push(@results,$row);
	 		
		}
		return \@results;
}

sub return_sv_hash {
	my ($self,$row) = @_;
	my $h;
		$h->{CHROM1} = $row->{CHROM1};
		
			$h->{infos}->{sr} = [$row->{sr1},$row->{sr2}];
			$h->{infos}->{pr} = [$row->{pr1},$row->{pr2}];
			$h->{PATIENT_ID} = $row->{patient};
			$h->{TYPE} = $row->{type};
			
			$h->{ID} = join("_",$row->{CHROM1},$row->{POS1},$row->{CHROM2},$row->{POS2});
			return $h;
}

sub create_table_total {
	my ($self) =@_;
	return $self->{create_table_total} if exists  $self->{create_table_total};
	my $colpos1 = $self->colpos1;
	my $colpos2 = $self->colpos2;
	my $colchr1 = $self->colchr1;
	my $colchr2 = $self->colchr2;
	my $dir = $self->dir;
	my $file = "*.parquet";
	my $project_id = $self->project->id;
	my $parquet_file = $dir."$file";
	warn "start";
	my $query = qq{CREATE TABLE sv_call_all AS
                             SELECT project,type,patient,$colchr1 as CHROM1,$colchr2 as CHROM2,$colpos1 as POS1 ,$colpos2 as POS2,caller,sr1,sr2,pr1,pr2
                           FROM '$parquet_file'
                           WHERE $colpos1 > -1 order by type,$colchr1,$colchr2,$colpos1,$colpos2;
	};
	warn "end";
	$self->dbh->do($query);
	$self->{create_table_total} = "sv_call_all";
	return $self->{create_table_total};
}

sub create_table_partial {
	my ($self,$chr1,$start1,$chr2,$start2,$limit) = @_;
	
	
	my $a = $start1-$limit;
	my $b = $start1+$limit;
	my $c = $start2-$limit;
	my $d = $start2+$limit;
	
	my $colpos1 = $self->colpos1;
	my $colpos2 = $self->colpos2;
	my $colchr1 = $self->colchr1;
	my $colchr2 = $self->colchr2;
	my $dir = $self->dir;
	my $file = "*.parquet";
	my $project_id = $self->project->id;
	my $parquet_file = $dir."$file";
	my $query = qq{
		PRAGMA threads=4;
		CREATE TABLE sv_call_all AS
                             SELECT project,type,patient,$colchr1 as CHROM1,$colchr2 as CHROM2,$colpos1 as POS1 ,$colpos2 as POS2,caller,sr1,sr2,pr1,pr2
                           FROM '$parquet_file'
                           where   CHROM1=$chr1 and CHROM2=$chr2 and POS1 > $a and   POS1 < $b  and POS2>$c and  POS2 < $d
                           order by type,$colchr1,$colchr2,$colpos1,$colpos2;
	};
	$self->dbh->do($query);
	$self->{create_table_total} = "sv_call_all";
	return $self->{create_table_total};
}



sub prepare_sv_details{
	my ($self) = @_;

	return $self->{prepare_sv} if exists $self->{prepare_sv};
	my $colpos1 = $self->colpos1;
	my $colpos2 = $self->colpos2;
	my $colchr1 = $self->colchr1;
	my $colchr2 = $self->colchr2;
	my $table = $self->create_table_total;
	$self->{prepare_sv} = $self->dbh->prepare(qq{SELECT patient,project,CHROM1 as chr1,POS1 as pos1,CHROM2 as chr2,POS2 as pos2,caller from $table
	  where   CHROM1=? and CHROM2=? and POS1 > ? and   POS1 < ?  and POS2>? and  POS2 < ?});
	
	return $self->{prepare_sv};
}

sub prepare_sv{
	my ($self) = @_;

	return $self->{prepare_sv} if exists $self->{prepare_sv};
	my $colpos1 = $self->colpos1;
	my $colpos2 = $self->colpos2;
	my $colchr1 = $self->colchr1;
	my $colchr2 = $self->colchr2;
	my $table = $self->create_table_total;
	$self->{prepare_sv} = $self->dbh->prepare(qq{SELECT patient,project from $table
	  where   CHROM1=? and CHROM2=? and POS1 > ? and   POS1 < ?  and POS2>? and  POS2 < ?});
	
	return $self->{prepare_sv};
}

has caller_type_flag =>  (
	is		=> 'rw',
	lazy=>1,
	default => sub {
		return  {
			"caller_sr" => 1,
			"caller_depth" => 2,
			"caller_coverage" => 4,
		};
	}
);

sub get_dejavu_details {
	my ($self,$chr1,$start1,$chr2,$start2,$limit) = @_;
	my $table = $self->create_table_partial($chr1,$start1,$chr2,$start2,$limit);
	$self->prepare_sv_details("SV")->execute($chr1,$chr2,$start1-$limit,$start1+$limit,$start2-$limit,$start2+$limit);
	my $x;
	while (my $row = $self->prepare_sv("SV")->fetchrow_hashref){ 
		my $chr1 = $row->{chr1};
		my $pos1 = $row->{pos1};
		my $chr2 = $row->{chr2};
		my $pos2 = $row->{pos2};
		
		$x->{$row->{project}}->{$row->{patient}}->{chr1} = $chr1;
		$x->{$row->{project}}->{$row->{patient}}->{pos1} = $pos1;
		$x->{$row->{project}}->{$row->{patient}}->{chr2} = $chr2;
		$x->{$row->{project}}->{$row->{patient}}->{pos2} = $pos2;
		
 		my $caller1 = $row->{'caller'};
		$x->{$row->{project}}->{$row->{patient}}->{caller_coverage} ++ if $caller1 & $self->caller_type_flag->{caller_coverage};
		$x->{$row->{project}}->{$row->{patient}}->{caller_depth} ++  if $caller1 & $self->caller_type_flag->{caller_depth};
		$x->{$row->{project}}->{$row->{patient}}->{caller_sr} ++  if $caller1 & $self->caller_type_flag->{caller_sr};
		
	}
	return $x;
}


sub get_dejavu {
	my ($self,$chr1,$start1,$chr2,$start2,$limit) = @_;
	  $self->prepare_sv("SV")->execute($chr1,$chr2,$start1-$limit,$start1+$limit,$start2-$limit,$start2+$limit);
	  my @results;
	  my $x;
	  	 my $hproject;
	  	$x->{dv_patients} =0;
	  while (my $row = $self->prepare_sv("SV")->fetchrow_hashref){ 
	  	my $project = $row->{project};
	  	my $patient = $row->{patient};
		$x->{projects}->{$project}->{$patient} ++;
		
	 	 $x->{nb_patients} ++;
	  }
	 $x->{this_project} = delete $x->{projects}->{$self->project->id};
	 
	  $x->{nb_projects} = scalar(keys %{$x->{projects}});
	  
	return $x;
}



1;
