package GenBoDuckDejaVuCNV;

use Moo;

use strict;
use Carp;
use warnings;
use Data::Dumper;
use JSON::XS;
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
		
		warn $self->project->buffer->config_path("dejavu","cnv");;
		
		return $self->project->buffer->config_path("dejavu","cnv");;
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

has colchr =>(
	is		=> 'ro',
	lazy=>1,
	default => sub {
		my ($self)= @_;
		return "chr".$self->version;
}
);

has rawcol =>  (
	is		=> 'rw',
	lazy=>1,
default => sub {
	my $col = ["project","patient","chr38","start38","end38","chr19","start19","end19","type","callers","caller_type_flag",'sr1','sr2','sr_qual','pr1','pr2',"depth_qual",'depth_CN','coverage_zscore','coverage_ratio','elementary','cn'];
	
});

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
has type_by_caller =>  (
	is		=> 'rw',
	lazy=>1,
default => sub {
return  {
	"wisecondor" => 4,
	"manta" => 1,
	"pbsv" => 1,
	"dragen-sv" =>1,
	"hificnv" =>2,
	"canvas" =>2,
	"dragen-cnv" =>2,
	"cnvnator" =>2,
}
}
);

has cnv_callers =>  (
	is		=> 'rw',
	lazy=>1,
default => sub {
return  {
    "wisecondor"          => 1 << 0,  # 2^0 = 1
    "canvas"        => 1 << 1,  # 2^1 = 2
    "manta"        => 1 << 2,  # 2^2 = 4
    "pbsv"        => 1 << 3,  
    "dragen-sv"        => 1 << 4,  
    "hificnv"        => 1 << 5,  
    "dragen-cnv"        => 1 << 6, 
     "cnvnator"        => 1 << 7, 
}
}
);



#has sth_query_project =>(
#	is		=> 'rw',
#	lazy=>1,
#	default => sub {
#		my ($self) =@_;
#		my $table = $self->create_table_project();
#		my $project_id = $self->project->id;
#		 my $sql = qq{
#			PRAGMA threads=4; SELECT *
#			FROM $table 
#			WHERE project=? and patient = ?;
#		 };
#		 
#		my $duckdb = $self->project->buffer->software('duckdb');
#		my $cmd = qq{set +H | $duckdb -json -c "$sql"};
#		my $json_duckdb = `$cmd`;
#		my $h = decode_json $json_duckdb;
#		return $h;
#	}
#);
#sub create_table_project {
#	my ($self) =@_;
#	my $uid = "cnv_".$self->project->id;
#	return $self->{$uid} if exists  $self->{$uid};
#	my $colchr = $self->colchr;
#	my $colstart = $self->colstart;
#	my $colend = $self->colend;
#	my $dir = $self->dir;
#	my $file = $self->project->name."*.parquet";
#	my $parquet_file = $dir."$file";
#	my $query = qq{
#		PRAGMA threads=4; CREATE TABLE $uid AS SELECT *, abs($colend - $colstart) AS length,
#                           FROM '$parquet_file'
#                           WHERE $colstart > -1  order by $colchr,$colstart,$colend;
#	};
#	my $duckdb = $self->project->buffer->software('duckdb');
#	my $cmd = qq{set +H | $duckdb -json -c "$query"};
#	my $json_duckdb = `$cmd`;
#	$self->{$uid} = $uid;
#	return $self->{$uid};
#}

#sub get_cnvs_by_project {
#	my ($self,$patient ) = @_;
#	my $colchr = $self->colchr;
#	my $colstart = $self->colstart;
#	my $colend = $self->colend;
#	my $patient_id = $patient->id;
#	my $project_id = $self->project->id;
#	my $dir = $self->dir;
#	my $proj_id = $self->project->id();
#	my $file = $self->project->name.'.'.$proj_id.".parquet";
#	my $parquet_file = $dir."$file";
#	my $sql = qq{
#		PRAGMA threads=4;
#		SELECT *
#		FROM '$parquet_file' 
#		WHERE project=$project_id and patient=$patient_id;
#	 };
#	 
#	my $duckdb = $self->project->buffer->software('duckdb');
#	my $cmd = qq{set +H | $duckdb -json -c "$sql"};
#	my $json_duckdb = `$cmd`;
#	my $list = decode_json $json_duckdb;
#	my @res;
#	foreach my $row (@$list) {
#		$row->{chromosome} = delete $row->{$colchr};
#		$row->{start} = delete $row->{$colstart};
#		$row->{end} = delete $row->{$colend};
#		$row->{uid} = join("_",$row->{chromosome},$row->{start},$row->{end},$row->{patient});
#		$row->{id} = join("_",$row->{chromosome},$row->{start},$row->{end});
#		push(@res,$row);
#	}
#	return $list;
#}

#sub create_table_total {
#	my ($self) =@_;
#	return $self->{create_table_total} if exists  $self->{create_table_total};
#	my $colchr = $self->colchr;
#	my $colstart = $self->colstart;
#	my $colend = $self->colend;
#	my $dir = $self->dir;
#	my $file = "*.parquet";
#	my $parquet_file = $dir."$file";
#	my $create_view_cmd = qq{ CREATE TABLE cnv_call_project AS SELECT *, abs($colend - $colstart) AS length,
#                           FROM '$parquet_file'
#                           WHERE $colstart > -1  order by $colchr,$colstart,$colend;
#	};
#	my $dbfile = $self->dbfile;
#	my $duckdb = $self->project->buffer->software('duckdb');
#	system("$duckdb $dbfile \"$create_view_cmd\"") == 0 or die "Erreur création vue: $?";
#	$self->{create_table_total} = "cnv_call_project";
#	return $self->{create_table_total};
#}



sub dejavu_details {
	my ($self,$type,$chr,$start,$end,$seuil ) = @_;
	#my $project_id = $self->project->id;
	my $project_id = '1';
	my $p = (100 - $seuil)/100;
	my $len = abs($end - $start) +1;
	my $vv =   int($len * $p);
	my $minl = $len - $vv;
	my $maxl = $len + $vv;
	my $minx = $start - $vv;
	my $maxx = $start + $vv;
	my $miny = $end - $vv;
	my $maxy = $end + $vv;
	#$self->create_table_partial($chr, $minx, $maxx, $miny, $maxy);
	my $colchr = $self->colchr;
	my $colstart = $self->colstart;
	my $colend = $self->colend;
	#$self->create_table_total;
	my $dir = $self->dir;
	my $file = "*.parquet";
	my $parquet_file = $dir."$file";
	 my $sql = qq{
		PRAGMA threads=20;
	 SELECT $colstart as start ,$colend as end ,project,patient,callers,caller_type_flag, abs($colend - $colstart) AS length,
		FROM '$parquet_file' 
		WHERE type = '$type' AND $colchr = '$chr'
		AND $colstart > -1 
		AND $colstart BETWEEN $minx AND $maxx
  		AND $colend   BETWEEN $miny AND $maxy
  		AND length  BETWEEN $minl AND $maxl
  		AND  project != $project_id
  		order by $colchr,$colstart,$colend;
	 };
	my $duckdb = $self->project->buffer->software('duckdb');
	my $cmd = qq{set +H | $duckdb -json -c "$sql"};
	my $json_duckdb = `$cmd`;
	my $list = decode_json $json_duckdb;
	
	# Récupération des résultats dans un tableau de hash
	my @results;
	my $nb;
	my $hash;
	my %pp;
	my @dj;
	
	foreach my $row (@$list) {
 		my $start1 = $row->{start};
 		my $end1 = $row->{end};	
 		my $identity = $self->getIdentityBetweenCNV($start,$end,$start1,$end1);
 		my $project1 = $row->{project};
 		my $patient1 = $row->{patient};
 		
		$hash->{$project1}->{$patient1}->{caller_sr} = 0 ;
		$hash->{$project1}->{$patient1}->{caller_depth} =0;
		$hash->{$project1}->{$patient1}->{caller_coverage} = 0;
		$hash->{$project1}->{$patient1}->{start} = $start1;
		$hash->{$project1}->{$patient1}->{end} = $end1;
 		
 		#ici detail projet et patient
 		my $caller1 = $row->{callers};
		$hash->{$project1}->{$patient1}->{caller_coverage} ++ if $caller1 & $self->caller_type_flag->{caller_coverage};
		$hash->{$project1}->{$patient1}->{caller_depth} ++  if $caller1 & $self->caller_type_flag->{caller_depth};
		$hash->{$project1}->{$patient1}->{caller_sr} ++  if $caller1 & $self->caller_type_flag->{caller_sr};
 		$identity = int($identity);
		my $string ="$identity";
 		$hash->{$project1}->{$patient1}->{identity} = int($string);
 		
	}
	return  $hash;
}


has dbfile =>  (
	is		=> 'rw',
	lazy=>1,
	default => sub {
		my $self = shift;
		my $dir = $self->project->rocks_pipeline_directory().'cnv.mydb.duckdb';
		warn $dir;
		return $dir
	}
);


#	is		=> 'rw',
#	lazy=>1,
#default => sub {
#	my ($self) =@_;
#	my $colchr = $self->colchr;
#	my $colstart = $self->colstart;
#	my $colend = $self->colend;
#	$self->create_table_total;
#	 my $sql = qq{
#	 SELECT $colstart as start ,$colend as end ,project,patient,callers,caller_type_flag,
#		FROM cnv_call_project 
#		WHERE type = ? AND $colchr = ?
#		AND $colstart BETWEEN ? AND ?
#  		AND $colend   BETWEEN ? AND ?
#  		AND length  BETWEEN ? AND ?
#  		AND  project != ?;
#	 };

sub prepare_cnv{
	my ($self, $type, $chr, $minx, $maxx, $miny, $maxy, $minl, $maxl,$project_id) = @_;
	my $table = $self->create_table_total;
	my $colchr = $self->colchr;
	my $colstart = $self->colstart;
	my $colend = $self->colend;
	$self->create_table_total;
	 my $sql = qq{
	 SELECT $colstart as start ,$colend as end ,project,patient,callers,caller_type_flag,
		FROM cnv_call_project 
		WHERE type='$type' AND $colchr='$chr'
		AND $colstart BETWEEN $minx AND $maxx
  		AND $colend   BETWEEN $miny AND $maxy
  		AND length  BETWEEN $minl AND $maxl
  		AND  project != $project_id;
	 };
	return $sql;
}

sub run_duckdb_query {
    my ($self, $sql) = @_;
	my $dbfile = $self->dbfile;
	my $duckdb = $self->project->buffer->software('duckdb');
    my $cmd = qq{$duckdb -json $dbfile "$sql"};
    open my $fh, "-|", $cmd or die "Impossible d'exécuter '$sql': $!";
    my $json;
    while (my $line = <$fh>) {
        chomp $line;
        $json .= $line;
    }
    if ($json) {
	    my $list = decode_json $json;
	    return $list;
    }
    return \[];
}

##seuil: 90%
#sub dejavu {
#	my ($self,$type,$chr,$start,$end,$seuil ) = @_;
#	my $project_id = $self->project->id;
#	my $p = (100 - $seuil)/100;
#	my $len = abs($end - $start) +1;
#	my $vv =   int($len * $p);
#	my $minl = $len - $vv;
#	my $maxl = $len + $vv;
#	my $minx = $start - $vv;
#	my $maxx = $start + $vv;
#	my $miny = $end - $vv;
#	my $maxy = $end + $vv;
#	my $query = $self->prepare_cnv($type,$chr, $minx, $maxx, $miny, $maxy, $minl, $maxl,$project_id);
#	my $list = $self->run_duckdb_query($query);
#	# Récupération des résultats dans un tableau de hash
#	my @results;
#	my $nb;
#	my $hash;
#	$hash->{nb_patients} = 0;
#	$hash->{nb_projects} = 0;
#	$hash->{caller_sr} = 0 ;
#	$hash->{caller_depth} =0;
#	$hash->{caller_coverage} = 0;
#	$hash->{string} = "";
#	my %pp;
#	my @dj;
#	foreach my $row (@$list) {
#		
# 		my $start1 = $row->{start};
# 		my $end1 = $row->{end};	
# 		my $identity = $self->getIdentityBetweenCNV($start,$end,$start1,$end1);
# 		my $project1 = $row->{project};
# 		my $patient1 = $row->{patient};
# 		my $caller1 = $row->{callers};
# 		$pp{$project1} ++;
#	
#		$hash->{nb_patients} ++;
#		$hash->{caller_coverage} ++ if $caller1 & $self->caller_type_flag->{caller_coverage};
#		$hash->{caller_depth} ++  if $caller1 & $self->caller_type_flag->{caller_depth};
#		$hash->{caller_sr} ++  if $caller1 & $self->caller_type_flag->{caller_sr};
# 	
# 		$identity = int($identity);
#		my $string ="$identity";
#		push(@dj,$string) if $string;
#	}
#	$hash->{nb_projects} = scalar(keys %pp);
#	$hash->{string}  = join(",",@dj);
#	return  $hash;
#}


has dbh =>  (
	is		=> 'rw',
	lazy=>1,
default => sub {
	my $dbh =   DBI->connect("dbi:ODBC:Driver=DuckDB;Database=:memory:", "", "", { RaiseError => 1 , AutoCommit => 1});
	return $dbh;
});

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
		WHERE project=? and patient = ?;
	 };
	 
	 my $sth = $self->dbh->prepare($sql);
	 return $sth;
}
);

sub create_table_project {
	my ($self) =@_;
	my $uid = "cnv_".$self->project->id;
	return $self->{$uid} if exists  $self->{$uid};
	my $colchr = $self->colchr;
	my $colstart = $self->colstart;
	my $colend = $self->colend;
	my $dir = $self->dir;
	my $file = $self->project->name."*.parquet";
	my $parquet_file = $dir."$file";
	my $query = qq{CREATE TABLE $uid AS
                           SELECT *, abs($colend - $colstart) AS length,
                           FROM '$parquet_file'
                           WHERE $colstart > -1  order by $colchr,$colstart,$colend;
	};
	$self->dbh->do($query);
	$self->{$uid} = $uid;
	return $self->{$uid};
}

sub create_table_total {
	my ($self) =@_;
	return $self->{create_table_total} if exists  $self->{create_table_total};
	my $colchr = $self->colchr;
	my $colstart = $self->colstart;
	my $colend = $self->colend;
	my $dir = $self->dir;
	my $file = "*.parquet";
	my $parquet_file = $dir."$file";
	my $query = qq{CREATE TABLE cnv_call_project AS
                           SELECT *, abs($colend - $colstart) AS length,
                           FROM '$parquet_file'
                           WHERE $colstart > -1  order by $colchr,$colstart,$colend;
	};
	$self->dbh->do($query);
	$self->{create_table_total} = "cnv_call_project";
	return $self->{create_table_total};
}

has sth_query_identity =>(
	is		=> 'rw',
	lazy=>1,
default => sub {
	my ($self) =@_;
	my $colchr = $self->colchr;
	my $colstart = $self->colstart;
	my $colend = $self->colend;
	$self->create_table_total;
	 my $sql = qq{
	 SELECT $colstart as start ,$colend as end ,project,patient,callers,caller_type_flag,
		FROM cnv_call_project 
		WHERE type = ? AND $colchr = ?
		AND $colstart BETWEEN ? AND ?
  		AND $colend   BETWEEN ? AND ?
  		AND length  BETWEEN ? AND ?
  		AND  project != ?;
	 };
	 
	 my $sth = $self->dbh->prepare($sql);
	 
	 return $sth;
}
);

#seuil: 90%
sub dejavu {
	my ($self,$type,$chr,$start,$end,$seuil ) = @_;
	my $project_id = $self->project->id;
	my $p = (100 - $seuil)/100;
	my $len = abs($end - $start) +1;
	my $vv =   int($len * $p);
	 my $minl = $len - $vv;
	 my $maxl = $len + $vv;
	 my $minx = $start - $vv;
	 my $maxx = $start + $vv;
	  my $miny = $end - $vv;
	 my $maxy = $end + $vv;
	$self->sth_query_identity->execute($type,$chr, $minx, $maxx, $miny, $maxy, $minl, $maxl,$project_id);
	# Récupération des résultats dans un tableau de hash
	my @results;
	my $nb;
	my $hash;
	$hash->{nb_patients} = 0;
	$hash->{nb_projects} = 0;
	$hash->{caller_sr} = 0 ;
	$hash->{caller_depth} =0;
	$hash->{caller_coverage} = 0;
	$hash->{string} = "";
		my %pp;
		my @dj;
		while (my $row = $self->sth_query_identity->fetchrow_hashref) {
			
	 		my $start1 = $row->{start};
	 		my $end1 = $row->{end};	
	 		my $identity = $self->getIdentityBetweenCNV($start,$end,$start1,$end1);
	 		my $project1 = $row->{project};
	 		my $patient1 = $row->{patient};
	 		my $caller1 = $row->{callers};
	 		$pp{$project1} ++;
		
			$hash->{nb_patients} ++;
			$hash->{caller_coverage} ++ if $caller1 & $self->caller_type_flag->{caller_coverage};
			$hash->{caller_depth} ++  if $caller1 & $self->caller_type_flag->{caller_depth};
			$hash->{caller_sr} ++  if $caller1 & $self->caller_type_flag->{caller_sr};
	 	
	 		$identity = int($identity);
			my $string ="$identity";
			push(@dj,$string) if $string;
		}
	$hash->{nb_projects} = scalar(keys %pp);
	$hash->{string}  = join(",",@dj);
	return  $hash;
}

sub get_cnvs_by_project {
	my ($self,$patient ) = @_;
	$self->sth_query_project->execute($self->project->id,$patient->id);
	# Récupération des résultats dans un tableau de hash
		my @results;
		my $colchr = $self->colchr;
		my $colstart = $self->colstart;
		my $colend = $self->colend;
		while (my $row = $self->sth_query_project->fetchrow_hashref) {
			$row->{chromosome} = delete $row->{$colchr};
			$row->{start} = delete $row->{$colstart};
			$row->{end} = delete $row->{$colend};
			$row->{uid} = join("_",$row->{chromosome},$row->{start},$row->{end},$row->{patient});
			$row->{id} = join("_",$row->{chromosome},$row->{start},$row->{end});
			push(@results,$row);
		}
		return \@results;
}

sub getIdentityBetweenCNV {
	my ( $self, $start1, $end1, $start2, $end2) = @_;
	
	#retourne le recouvrement en % de la longueur du plus long des deux evenements 
	
	my $overlap = min( $end1, $end2 ) - max( $start1, $start2 );
	die() if abs( $start1 - $end1 ) ==0;
	my $overlap1 = $overlap / abs( $start1 - $end1 );
	my $overlap2 = $overlap / abs( $start2 - $end2 );

	return min($overlap1*100,$overlap2*100);
}

1;
