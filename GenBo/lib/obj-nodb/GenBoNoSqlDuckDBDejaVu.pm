package GenBoNoSqlDuckDBDejaVu;
#use lib "$Bin/";
use strict;
use warnings;
use Moo;
use Data::Dumper;
use JSON::XS;
use POSIX;
#use GenBoNoSqlRocksChunks;
use Carp;
use Fcntl ':flock';
use File::NFSLock qw(uncache);
our $json_chr_length = qq{{"HG19":{"6":"171115067","11":"135006516","9":"141213431","15":"102531392","14":"107349540","1":"249250621","8":"146364022","17":"81195210","7":"159138663","13":"115169878","MT":"16569","22":"51304566","Y":"59373566","3":"198022430","21":"48129895","18":"78077248","5":"180915260","X":"155270560","20":"63025520","16":"90354753","2":"243199373","12":"133851895","19":"59128983","4":"191154276","10":"135534747"},"HG38":{"21":"46709983","Y":"57227415","18":"80373285","3":"198295559","22":"50818468","MT":"16569","20":"64444167","16":"90338345","X":"156040895","5":"181538259","2":"242193529","10":"133797422","4":"190214555","19":"58617616","12":"133275309","11":"135086622","6":"170805979","14":"107043718","15":"101991189","9":"138394717","17":"83257441","8":"145138636","1":"248956422","13":"114364328","7":"159345973"}}};
sub getChromsosomes {
	my ($version) = @_;
	my $h = decode_json($GenBoNoSqlDuckDBDejaVu::json_chr_length);
	return [keys %{$h->{$version}}];
}
has genome =>(
	is		=> 'ro',
	required=>1,

);

has chr_length =>(
	is 	=> 'ro',
	lazy => 1,
	default => sub {
		my $self = shift;
		my $h = decode_json($json_chr_length);
		
		die(Dumper $h ) unless exists $h->{$self->genome}->{$self->chromosome_name};
		my $v = $h->{HG19}->{$self->{chromosome_name}};
		$v = $h->{HG38}->{$self->{chromosome_name}} if $h->{HG38}->{$self->{chromosome_name}} > $v;
		return  $v;
	}

 );
 

has pipeline =>(
	is		=> 'rw',
	default => sub {
		return undef;
	}

);

has compress =>(
	is		=> 'ro',
	default => sub {
		return undef;
	}

);
has index =>(
	is		=> 'ro',
	default => sub {
		return undef;
	}

);
has pack =>(
	is		=> 'ro',
	default => sub {
		return undef;
	}
);
has description =>(
	is		=> 'ro',
	default => sub {
		return undef;
	}
);
has factor =>(
	is		=> 'ro',
	default => sub {
		return undef;
	}
);
has dir =>(
	is		=> 'rw',
	required=>1,
);



sub dir_db {
		my ($self,$chr,$id) = @_;
		$chr = $self->chromosome_name unless $chr;
		my $dd = $self->dir."/".$chr."/";
		$dd .= $id."/" if $id;
		unless (-e $dd){
			if ($self->mode eq "r"){
				confess($dd)."-".$self->mode;
			}
			system("mkdir -p $dd && chmod a+w $dd");
		}
		return $dd;
		
}

has chunk_size =>(
	is		=> 'rw',
	lazy=>1,
	default => sub {
		return 5_000_000;
	}
);
has chromosome =>(
	is		=> 'rw',
	required=>1,
);
has chromosome_name =>(
	is		=> 'ro',
	lazy=>1,
	 default => sub {
		my $self = shift;
		my $chr = $self->chromosome;
		$chr =~ s/chr//;
		$chr="MT" if $chr eq "M";
		return $chr;
	 }
);
has mode =>(
	is		=> 'rw',
	required=>1,
);


has chunks =>(
	is		=> 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		if($self->mode eq "c" or $self->mode eq "t" or ($self->mode eq "w" && !(-e $self->json_file))){
			$self->clean() if -e $self->json_file;
			my $length = $self->chr_length;
#			warn $length;
			die() unless $length;
			my $c = $self->divide_by_chunks($self->chromosome_name,$length,$self->chunk_size);
			$self->write_config($c);
			return $c;
		}
		
		else {
			my $length = $self->chr_length;
			my $c = $self->divide_by_chunks($self->chromosome_name,$length,$self->chunk_size);
			$self->write_config($c);
			return $c;
			die() unless -e $self->json_file;
			return $self->load_config();
		}
		
}
);
sub clean {
	my ($self) =@_;
	
	#die();
}


has tree =>(
	is		=> 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
 		my $tree =  Set::IntervalTree->new;
 		for my $r (@{$self->regions}){
 			$tree->insert($r, $r->{start}, $r->{end});
 		}
 		return $tree;
	}
);


sub get_genomic_regions {
	my ($version,$size) = @_;
	warn $size;
	confess() unless $size;
	confess() unless $version;
	my $h = decode_json($GenBoNoSqlDuckDBDejaVu::json_chr_length);
	my $r;
	foreach my $chr (keys %{$h->{$version}}) {
	my $length = $h->{HG19}->{$chr};
	$length  = $h->{HG38}->{$chr} if $h->{HG38}->{$chr} > $length;
	die() unless $length;
	my $chunks = GenBoNoSqlDuckDBDejaVu::divide_by_chunks(undef,$chr,$length,$size);
	foreach my $id (sort{$chunks->{$a}->{start} <=> $chunks->{$b}->{end}} keys %{$chunks}){
			push(@$r,{chromosome=>$chr,id=>$id,start=>$chunks->{$id}->{start},end=>$chunks->{$id}->{end},tabix=>$chr.":".$chunks->{$id}->{start}."-".$chunks->{$id}->{end}});
			
		}
	}
	return $r;
}


has regions =>(
is		=> 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $r;
		foreach my $id (sort{$self->chunks->{$a}->{start} <=> $self->chunks->{$b}->{end}} keys %{$self->chunks}){
			
			push(@$r,{chromosome=>$self->chromosome_name,id=>$id,start=>$self->chunks->{$id}->{start},end=>$self->chunks->{$id}->{end},tabix=>$self->chromosome.":".$self->chunks->{$id}->{start}."-".$self->chunks->{$id}->{end}});
			
		}
 		return $r;
	}
);

has json_file =>(
is		=> 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->dir_db."chunk.json";
	}
);
has vector_index_region =>(
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return [];
	}
);


sub save_vector_index_region {
	my($self,$v) = @_;
	open(my $fh ,$self->json_file) or die("can t open file");
	my $json_str = do { local $/; <$fh> };
	close($fh);
	my $h = decode_json($json_str);
	$h->{vector_index_region} = $v;
	open(my $fh2,">".$self->json_file) or die("can t open file");
	print $fh2 encode_json($h);
	close ($fh);
}
sub write_config {
	my ($self,$chunks) = @_;
	open(my $fh,">".$self->json_file) or die("can t open file ".$self->json_file);
	my $h;
	$h->{chunks} = $chunks;
	$h->{compress} = $self->compress;
	$h->{pack} = $self->pack;
	$h->{description} = $self->description;
	$h->{factor} = $self->factor;
	$h->{index} = $self->index;
	$h->{vector_index_region} = $self->vector_index_region;
	print $fh encode_json($h);
	close ($fh);
}
sub load_config {
	my ($self) = @_;
	open(my $fh ,$self->json_file) or die("can t open file ".$self->json_file);
	my $json_str = do { local $/; <$fh> };
	close($fh);
	my $h = decode_json($json_str);
	if($self->mode ne "c"){
	$self->{compress} = delete $h->{compress};
	$self->{pack} = delete $h->{pack};
	$self->{description} = delete $h->{description};
	$self->{factor} = delete $h->{factor};
	$self->{index} = delete $h->{index};
	$self->{vector_index_region} = delete $h->{vector_index_region};
	}
	return delete $h->{chunks};
}
sub divide_by_chunks{
	my ($self,$chr,$length,$size) = @_;
	confess() unless $size;
		my $region;
		
		my $from =1;
    	my $to = $length;
  
    	my $start;
    	my $end;
    	
    	while ($from < $to){
        $start = $from;
        $end = $from + $size;
        if ($end > $to){
            $end = $from + $size;
        }
        my $id_chunk = $chr.".".$from.".".$end;
        $region->{$id_chunk}->{start} = $from;
        $region->{$id_chunk}->{end} = $end;
        
        $from = $end;
       }
       return $region;
}
sub _current_db {
	my ($self) = @_;
	#return 
	return $self->{_current_db};
}

sub change_id {
	my ($self,$key) = @_;
	my $start = $key;
	if ($key =~ /_/){
		my @t = split("_",$key);
		$start = shift(@t);
		my $chaine = sprintf("%010d", $start);
		return($start,$chaine."_".join("_",@t));
	}
	return ($start,sprintf("%010d", $key));
}

sub current {
	my ($self,$key) = @_;
	return $self->{current_db};
}

sub convert_rocksId_varId {
	my ($self, $rocksId) = @_;
	
}



sub _dbh {
	my $dbh =   DBI->connect("dbi:ODBC:Driver=DuckDB;Database=:memory:", "", "", { RaiseError => 1 , AutoCommit => 1});
	#$dbh->do("PRAGMA memory_limit='5GB'");
	return $dbh;
}

sub write_parquet {
	my ($self) = @_;
	
}
sub close {
	my ($self) = @_;
	return unless exists $self->{_currentdb};
	if ($self->{MODIFY}){
	if ($self->mode eq "w"){
		my $parquet_file = $self->{_currentdb}->{parquet};
		
		my $sql = "COPY (SELECT * FROM variants ORDER BY key) TO '$parquet_file'  (FORMAT 'parquet', COMPRESSION 'ZSTD');";
		$self->{_currentdb}->{dbh}->do($sql);
		$self->{_currentdb}->{lock}->unlock;
		uncache ($parquet_file);
	}
	elsif ($self->mode eq "c"){
		my $parquet_file = $self->{_currentdb}->{parquet};
		warn $parquet_file;
		my $sql = "COPY (SELECT * FROM variants ORDER BY key) TO '$parquet_file'  (FORMAT 'parquet', COMPRESSION 'ZSTD');";
		$self->{_currentdb}->{dbh}->do($sql);
		confess() unless -e $parquet_file;
	}
	delete $self->{MODIFY};
	}
	$self->{_currentdb}->{dbh}->disconnect();	
	delete $self->{_currentdb};
}

sub _create_table {
	my ($self) = @_;
	my $sql = qq{
    	CREATE TEMP TABLE IF NOT EXISTS variants (
        project INTEGER,
        start INTEGER,
        key TEXT,
        he INTEGER,
        ho INTEGER,
        value TEXT
    )
	};
	$self->{_currentdb}->{dbh}->do($sql);
}
sub get_dbh {
	my ($self,$region) = @_;
	if (exists $self->{_currentdb}){
		return $self->{_currentdb}->{dbh} if $self->{_currentdb}->{uniq_id} eq $region->{chromosome}.".".$region->{id};
		$self->close();
	}


	my $id = $region->{id};
	my $set = Set::IntSpan::Fast::XS->new();
 	$set->add_range($region->{start}, ($region->{end}-1));
	$self->{_currentdb}->{intspan} = $set;
	$self->{_currentdb}->{dbh} =  $self->_dbh();
	$self->{_currentdb}->{uniq_id} = $region->{chromosome}.".".$region->{id};
	# Charger le fichier Parquet dans une table temporaire
	my $d = $self->dir_db($region->{chromosome});
	my $parquet_file = "$d/$id.parquet";
	$self->{_currentdb}->{parquet} = $parquet_file;
	$self->{_currentdb}->{parquet_lock} = $parquet_file.".lock";
	
	if ($self->mode eq "w"){
		$self->mode("c") unless -e $parquet_file;
	}
	if ($self->mode eq "r") {
			$self->{_currentdb}->{dbh}->do(qq{CREATE TEMP TABLE variants AS SELECT *  FROM read_parquet("$d/$id.parquet") ;}) if (-e $parquet_file);
			$self->_create_table() unless (-e $parquet_file);
	}
	elsif ($self->mode eq "c") {
		unlink $parquet_file if -e $parquet_file;
		confess();
		$self->_create_table();
	}
	elsif ($self->mode eq "w") {
		$self->{_currentdb}->{lock} = new File::NFSLock {
  			file      => $self->{_currentdb}->{parquet_lock},
  			lock_type => LOCK_EX,
  			stale_lock_timeout => 5 * 60, # 30 min
		};
		if ($self->{_currentdb}->{lock}) {
			$self->{_currentdb}->{dbh}->do(qq{CREATE TEMP TABLE variants AS SELECT *  FROM read_parquet("$d/$id.parquet") ;});
		}
		else {
			confess();
		}
	}
	else {
		confess();
	}
	return $self->{_currentdb}->{dbh};
}






sub get_projects {
	my ($self,$region) = @_;
	my $dbh = $self->get_dbh($region);
	
	my $sth = $dbh->prepare("SELECT DISTINCT project FROM variants");
	$sth->execute();

# Récupération des résultats dans un tableau
			my @projects;
		while (my ($project) = $sth->fetchrow_array) {
    	push @projects, $project;
	}
	return \@projects;
}

sub insert {
	my ($self,$region,$pid,$start,$key,$he,$ho,$encode) =@_;
	my $dbh = $self->get_dbh($region);
	my $sth = $dbh->prepare("INSERT INTO variants (project,start, key,he,ho,value) VALUES (?, ?, ?, ?, ? , ?)");
	$sth->execute($pid,$start,$key,$he,$ho,$encode);
	
}
sub count_rows {
    my ($self,$file) = @_;
    my $query = "SELECT COUNT(*) FROM '$file'";
    my $sth = $self->_dbh->prepare($query);
    $sth->execute();
    my ($count) = $sth->fetchrow_array();
    return $count;
}
has temp_dir =>(
	is 	=> 'ro',
	lazy => 1,
	default => sub {
		return "/data-beegfs/tmp/";
	}
 );
 
 sub filetemp {
 	my ($self,$region) =@_;
 	warn $self->temp_dir."/".$region->{chromosome}.".".$region->{id};
 	return $self->temp_dir."/".$region->{chromosome}.".".$region->{id};
                        
 }
 sub parquet_file {
 	my ($self,$region)= @_;
 	my $id = $region->{id};
	my $d = $self->dir_db($region->{chromosome});
	return "$d/$id.parquet";
 }
  sub parquet_file_project {
 	my ($self,$region,$project_name)= @_;
 	my $id = $region->{id};
	my $d = $self->dir_db($region->{chromosome},$region->{id});
	return "$d/".$project_name.".parquet";
 }
 sub lock {
 	my ($self,$file) =@_;
 	return new File::NFSLock {
  			file      => $file.".lock",
  			lock_type => LOCK_EX,
  			stale_lock_timeout => 5 * 60, # 30 min
		};
 }
 

sub modify_parquet{
	my ($self,$parquet_file,$tmp_parquet,$query) = @_;
	my $lock = $self->lock($parquet_file);
	if ($lock) {
	uncache ($parquet_file);		
	my $dbh = $self->_dbh ;
	# Requête pour fusionner le Parquet et le CSV
	$dbh->do($query);
	system("mv $tmp_parquet $parquet_file");
	#$lock->unlock;
	uncache ($parquet_file);
	$dbh->disconnect();
	#}
	return ;
}
}

sub delete_project {
	my ($self,$region,$pid) = @_;
my $parquet_file =  $self->parquet_file($region);
my $tmp_parquet = $self->filetemp($region);

my $query = qq{
	COPY (SELECT * FROM '$parquet_file' WHERE project != '$pid') 
  TO '$tmp_parquet' (FORMAT PARQUET, COMPRESSION 'ZSTD');
};
my $dbh = $self->_dbh ;
$dbh->do("PRAGMA memory_limit='3GB'");
$dbh->do($query) or die "Failed to rename temp file: $!";;
rename($tmp_parquet, $parquet_file) or die "Failed to rename temp file: $!";
$dbh->disconnect;

#$self->modify_parquet($parquet_file,$tmp_parquet,$query);
}

sub add_csv {
	my ($self,$region,$csv_file) = @_;
	confess($csv_file) unless -e $csv_file;
	my $parquet_file =  $self->parquet_file($region);
	my $tmp_parquet = $self->filetemp($region);
	my $query = "
    	COPY (
        	SELECT * FROM '$parquet_file'
        	UNION ALL
        SELECT * FROM '$csv_file'
    )
    TO '$parquet_file' (FORMAT PARQUET , COMPRESSION 'ZSTD',OVERWRITE TRUE);";
    my $dbh = $self->_dbh ;
	$dbh->do($query);
	# Fermer la connexion
	#$self->modify_parquet($parquet_file,$tmp_parquet,$query);
	unlink $csv_file;
	$dbh->disconnect;
}

sub add_csv2 {
	my ($self,$region,$csv_file) = @_;
	confess($csv_file) unless -e $csv_file;
	my $parquet_file =  $self->parquet_file($region);
	my $query = "
    	COPY (
        SELECT * FROM '$csv_file' order by key
    )
    TO '$parquet_file' (FORMAT PARQUET , COMPRESSION 'ZSTD',PARTITION_BY ('project'),OVERWRITE TRUE);";
    my $dbh = $self->_dbh ;
	$dbh->do($query);
	# Fermer la connexion
	#$self->modify_parquet($parquet_file,$tmp_parquet,$query);
	unlink $csv_file;
	$dbh->disconnect;

}

sub add_csv_project {
	my ($self,$region,$csv_file,$project_name) = @_;
	confess($csv_file) unless -e $csv_file;
	my $parquet_file =  $self->parquet_file_project($region,$project_name);
	my $query = "
    	COPY (
        SELECT * FROM '$csv_file' order by key
    )
    TO '$parquet_file' (FORMAT PARQUET , COMPRESSION 'ZSTD',OVERWRITE TRUE);";
    my $dbh = $self->_dbh ;
	$dbh->do($query);
	unlink $csv_file;
	$dbh->disconnect;

}
sub get_dbh_by_position {
	my ($self,$pos,$pn) = @_;
	if (exists $self->{_currentdb}){
		return $self->{_currentdb}->{dbh} if $self->{_currentdb}->{intspan}->contains($pos);
	}
	my $results = $self->tree->fetch($pos,$pos+1);
	confess(@$results) if scalar(@$results) > 1;
	my $region = $results->[0];
	return $self->get_dbh($region);
	
}
sub get_region {
	my ($self,$pos) = @_;
	my $results = $self->tree->fetch($pos,$pos+1);
	confess(@$results) if scalar(@$results) > 1;
	my $region = $results->[0];
	return $region;
}
sub get_sth_select_dejavu {
	my ($self,$pos,$pn) = @_;
	my $dbh = $self->get_dbh_by_position($pos,$pn);
	return $self->{_currentdb}->{sth} if exists $self->{_currentdb}->{sth};
	$self->{_currentdb}->{sth}= $self->{_currentdb}->{dbh}->prepare("SELECT  key, he as he, ho as ho  FROM variants where key = ? ");
	return $self->{_currentdb}->{sth};
	
	
}
sub dejavu {
	my ($self,$id,$pn) = @_;
	my ($pos,$a) =split("!",$id);
	$pos *= 1;
	my $sth =  $self->get_sth_select_dejavu($pos,$pn);
	
	#my $sth = $dbh->prepare("SELECT  key, he as he, ho as ho,value as t  FROM variants where key = '$id' ");
	#my $sth = $dbh->prepare("SELECT   key, he as he, ho as ho,value as t FROM variants  where start>$pos-2 and start<$pos+2 ");
	$sth->execute($id);
	my $st;
	my $n =0;
	
	while (my $row = $sth->fetchrow_hashref) {
#	#	warn $row->{key}." ".$row->{he}." ".$row->{ho};
	#	$st.=$row->{t};
		$n += $row->{he} + $row->{ho};
	}
}
sub dejavu2 {
	my ($self,$id,$pn) = @_;
	my ($pos,$a) =split("!",$id);
	$pos *= 1;
	my $dbh =  $self->get_db($pos,$pn);
	my $sth = $dbh->prepare("SELECT  key, he as he, ho as ho,value as t  FROM variants where key='$id' ");
	#my $sth = $dbh->prepare("SELECT  key,SUM(he) as he, sum(ho) as ho FROM variants  where key='$id' group by key");
	$sth->execute();
	my $st;
	my $n =0;
	
	while (my $row = $sth->fetchrow_hashref) {
		
		$st.=$row->{t};
		$n += $row->{he} + $row->{ho};
	}
}

sub dejavu1 {
	my ($self,$pos) = @_;
	my $dbh =  $self->get_db($pos);
	my $sth = $dbh->prepare("SELECT  key, he as he, ho as ho,value as t  FROM variants where start = $pos ");
	#my $sth = $dbh->prepare("SELECT  key,SUM(he) as he, sum(ho) as ho FROM variants  where key='$id' group by key");
	$sth->execute();
	my $st;
	my $n =0;
	my $id;
	while (my $row = $sth->fetchrow_hashref) {
		$id->{$row->{key}} ++;
		
		#$st.=$row->{t};
		#$n += $row->{he} + $row->{ho};
	}
	warn Dumper $id;
}

sub stringify_pos {
	my ($self,$pos) = @_;
	return ($pos,sprintf("%010d", $pos));
}



sub spliceAI_id {
	my ($self,$id) = @_;
	my ($pos,$a) =split("!",$id);
	$pos *= 1;
	$self->get_db($pos);
	return $self->get_db($pos)->spliceAI_id($id);
}
sub spliceAI {
	my ($self,$pos,$id) = @_;
	$self->get_db($pos);
	return $self->get_db($pos)->spliceAI($id);
}





sub DESTROY {
	#warn "END rocks genome";
}

1;