package GenBoNoSql;



# ABSTRACT: An embedded, NoSQL SQLite database with SQL indexing

use strict;
use warnings;
use DBD::SQLite;
use Moo;

#use Sereal::Decoder qw(decode_sereal sereal_decode_with_object scalar_looks_like_sereal);
#use  Sereal::Encoder::Constants qw(SRL_PROTOCOL_ENCODING_SNAPPY);
#use Sereal::Encoder  qw(SRL_SNAPPY);
use Data::Dumper;
use Compress::Snappy;
use Storable qw/thaw freeze/;
use Tie::LevelDB;
use Carp;

has dir => (
	is		=> 'ro',
	required=> 1,
);

has extension =>(
	is		=> 'rw',
default => sub {
		return "lite";
		
}
);

has cache_limit =>(
	is		=> 'rw',
		lazy	=> 1,
	default => sub {
		return 10000;
	}
);

has no_lock =>(
	is		=> 'ro',
	default => sub {
		return undef;
	}
);
has mode => (
	is		=> 'ro',
	required=> 1,
);

has name_table_search => (
	is		=> 'ro',
		lazy	=> 1,
	default => sub {
		return "__SEARCH__"
	}
);
has name_table_data => (
	is		=> 'ro',
	default => sub {
		return "__DATA__"
	}
);

has new_database =>(
is		=> 'rw',
lazy=>1,
default => sub {
		my $self = shift;
	
		return 0 unless $self->mode eq 'c';
		
		 if ($self->mode eq 'c'){
		 	return 1;
		 }
		 unless (-e $self->dir){
		 	return 1;
		 }
		 return undef;
	
	},
);

has create =>(
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return 1 if  $self->mode =~/c/;
		return;
	}
);

has write =>(
		is		=> 'rw',
		lazy	=> 1,
		default => sub {
		my $self = shift;
		return 1 if  $self->mode =~/w/ or $self->mode =~/c/;
		return;
		
}
);


has integer=>(
	is		=> 'rw',
		lazy	=> 1,
default => sub {
		my $self = shift;
	
		return 1 if  $self->mode =~/i/;
		return;
		
}
);

has searchable=>(
	is		=> 'rw',
		lazy	=> 1,
	default => sub {
		return;
		
	}
);
has buffer=>(
	is		=> 'rw',
		lazy	=> 1,
default => sub {
		return {}
}
);
sub db {
	my ($self,$chr) = @_;
	return $self->{leveldb}->{$chr} if exists $self->{leveldb}->{$chr};
	 tie  %{$self->{leveldb}->{$chr}}, 'Tie::LevelDB', $self->dir()."/$chr.leveldb";
	return $self->{leveldb}->{$chr};
}


sub clear {
	my ($self,$key1) = @_;
	my $f = $self->db_name($key1);
	if (-e $f ){
		unlink $f ;
	}
	return 1;
}

sub db_name {
	my ($self,$key1) = @_;
	confess() unless $key1;
	return $self->dir."/$key1.".$self->extension;
}

sub exists_db {
	my ($self,$key1) = @_;

	return  $self->{file}->{$key1} if exists $self->{file}->{$key1} ;

	unless (-e $self->dir."/$key1.".$self->extension){
		return;
	}
	if (-z $self->dir."/$key1.".$self->extension){
		return;
	}
		$self->{file}->{$key1} =1;
	return 1;
}
sub exists_table {
	my ($self,$key,$table_name) = @_;
	 my $aray_ref = $self->dbh($key)->selectcol_arrayref(qq{SELECT count(*) FROM sqlite_master WHERE type='table' AND name="$table_name";});
	return 1 if $aray_ref->[0] >0;
	return undef;
}
sub dbh {
	my ($self,$key1) = @_;

	return $self->{dbh}->{$key1} if exists $self->{dbh}->{$key1};
	unless (-e $self->dir()){
		mkdir $self->dir() unless -e $self->dir();
		my $dir = $self->dir();
		system("chmod a+rwx $dir 2>/dev/null" );
	}
#	warn $key1." ".$self->db_name($key1);
	
	my $db = $self->db_name($key1);
	if ($self->mode eq 'c'){
		unlink $db if -e $db;
		confess("problem with $db") if -e $db;
	}
	
	 unless (-e $db){
	 	$self->{create}  =1;
	 	
	 	system("touch $db;chmod a+rw $db 2>/dev/null;")
	 }
		my $dsn = "dbi:SQLite:dbname=".$db;
	
	$self->{dbh}->{$key1} = DBI->connect(          
    				$dsn, 
    				"",                          
    				"",                   
		) or confess ($DBI::errstr);
		if ($self->write){
		$self->{dbh}->{$key1}->do("PRAGMA cache_size = 40000") or confess ($DBI::errstr." ".$key1);
		$self->{dbh}->{$key1}->do("PRAGMA synchronous = OFF");
		#PRAGMA journal_mode=WAL;
		
		$self->{dbh}->{$key1}->do("PRAGMA journal_mode = OFF;");
		$self->{dbh}->{$key1}->do("PRAGMA locking_mode = EXCLUSIVE;") unless $self->no_lock;
		$self->{dbh}->{$key1}->do("PRAGMA temp_store = MEMORY;");
		$self->{dbh}->{$key1}->do("PRAGMA auto_vacuum = NONE;");
		$self->{dbh}->{$key1}->do("PRAGMA count_changes=OFF;");
		#$self->{dbh}->{$key1}->do("PRAGMA threads = 2");
		
		}
		else {
		$self->{dbh}->{$key1}->do("PRAGMA cache_size = 40000") or confess ($DBI::errstr." ".$key1);
		$self->{dbh}->{$key1}->do("PRAGMA synchronous = OFF");
		$self->{dbh}->{$key1}->do("PRAGMA journal_mode = DELETE;");
		$self->{dbh}->{$key1}->do("PRAGMA locking_mode = EXCLUSIVE;")  unless $self->no_lock;;
		$self->{dbh}->{$key1}->do("PRAGMA temp_store = MEMORY;");
		$self->{dbh}->{$key1}->do("PRAGMA auto_vacuum = NONE;");
		$self->{dbh}->{$key1}->do("PRAGMA count_changes=OFF;");
		$self->{dbh}->{$key1}->do("PRAGMA threads = 8");
			#$self->{dbh}->{$key1}->do("PRAGMA cache_size = 4000000") or confess ($DBI::errstr." ".$key1);
			#$self->{dbh}->{$key1}->do("PRAGMA temp_store = MEMORY;");
			
#			$self->{dbh}->{$key1}->do("PRAGMA query_only = 1; ");
		}
		return $self->{dbh}->{$key1};
		return $self->{dbh}->{$key1};
}

sub is_searchable {
	my ($self,$key) = @_;
	
	return 1 if $self->searchable;
	return $self->{$key}->{searchable} if exists $self->{$key}->{searchable};
	my $table_name = $self->name_table_search();
	
	 my $aray_ref = $self->dbh($key)->selectcol_arrayref(qq{SELECT count(*) FROM sqlite_master WHERE type='table' AND name="$table_name";});
	 
	$self->{$key}->{searchable} = 1 if $aray_ref->[0] >0;
	$self->{$key}->{searchable} = undef;
	return $self->{$key}->{searchable};	
} 


sub create_table {
	my ($self,$chr) = @_;
	return $self->{table}->{$chr}   if (exists $self->{table}->{$chr});
	#if ($self->mode eq 'c'){
	#	$self->clear($chr);
	#}
	 my $table_name =$self->name_table_data;
	 	
	 $self->{table}->{$chr} = $table_name;
	 	 my $aray_ref = $self->dbh($chr)->selectcol_arrayref(qq{SELECT count(*) FROM sqlite_master WHERE type='table' AND name="$table_name";});
	return $self->{table}->{$chr} if ($self->exists_table($chr,"__DATA__")) &&  $self->mode ne 'c'  ;;

	 if ($self->write){
	$self->dbh($chr)->do("DROP TABLE IF EXISTS $table_name")  or confess($DBI::errstr." ".." ".$table_name)  if $self->mode eq 'c'  ;
	
	if ($self->integer){
			$self->dbh($chr)->do("CREATE TABLE if not exists $table_name (_key INTEGER PRIMARY KEY, _value BLOB) ;")  or confess($DBI::errstr);;
	}
	else{
			$self->dbh($chr)->do("CREATE TABLE if not exists $table_name (_key VARCHAR(250), _value BLOB) ")  or confess($DBI::errstr." ".$table_name);;
			$self->dbh($chr)->do("CREATE UNIQUE INDEX if not exists _key_idx  on $table_name (_key); ")  or confess($DBI::errstr." ".$table_name);;
			
	}
	 	}
	 	
	#system("chmod  a+rw  ".$self->dir()."/$chr.lite 2>/dev/null") ;
	return 	$self->{table}->{$chr} ;
	
}

sub decode {
	my ($self,$code) = @_;
	return undef unless $code;
	my $obj = thaw (decompress($code));
	return $obj->{data};
}
sub encode {
		my ($self,$code) = @_;
		return compress(freeze ({data=>$code}));
}

sub get_ids {
	my ($self,$chr) = @_;
	my $table_name = $self->create_table($chr);
	my $toto;
	my $array = $self->dbh($chr)->selectall_hashref(
      "SELECT _key as key FROM $table_name ",1
  );
 
  return [ keys %$array];
}
#sub get_ids1 {
#	my ($self,$chr) = @_;
#	my $table_name = $self->create_table($chr);
#	my $toto;
#	my $array = $self->dbh($chr)->selectall_hashref(
  #    "SELECT _key as key FROM $table_name ",1
  #);
 
  #return [ keys %$array];
#}
sub get {
	my ($self,$key1,$key2,$debug) = @_;
	my $t = $self->get_lite($key1,$key2);
	return $t;
}
sub get_like{
	my ($self,$key1,$key2) = @_;
	return $self->get_like_lite($key1,$key2);
}

sub get_bulk{
	my ($self,$key1,$key2) = @_;
	return $self->get_bulk_lite($key1,$key2);
}
sub count_bulk{
	my ($self,$key1,$key2) = @_;
	return $self->count_bulk_lite($key1,$key2);
}
sub get_all{
	my ($self,$key1) = @_;
	return $self->get_all_lite($key1);
}
sub put {
	my ($self,$key1,$key2,$value) = @_;
	return $self->put_lite ($key1,$key2,$value);
}
sub set {
	my ($self,$key1,$key2,$value) = @_;
	return $self->put_lite ($key1,$key2,$value);
}
sub exists {
		my ($self,$key1,$key2) = @_;
		#confess() if exists $self->{exist}->{$key1};
		# $self->{exist}->{$key1} = 1;
		return $self->{exist}->{$key1.$key2} if exists  $self->{exist}->{$key1.$key2};
		#warn $key1;
		#warn  $self->{exist}->{$key1};
		#return undef unless ($self->exists_db($key1));
		if ( $self->exists_lite ($key1,$key2)){
			$self->{exist}->{$key1.$key2} = 1;
			return 1;
		}
		return undef;
}


sub exists_leveldb {
	my ($self,$key1,$key2) = @_;
	return defined $self->db($key1)->{$key2};
}

sub get_buffer {
	my ($self,$key1,$key2) = @_;
	#return $self->decode($self->buffer->{$key1} ->{$key2}) if exists $self->buffer->{$key1} ->{$key2};
	return $self->buffer->{$key1} ->{$key2} if exists $self->buffer->{$key1} ->{$key2};
	return $self->get_lite($key1,$key2);
}

sub put_buffer {
	my ($self,$key1,$key2,$value) = @_;
	my $hbuffer = $self->buffer->{$key1} ;
	$hbuffer->{$key2}= $value;#$self->encode($value);
	 my $table = $self->create_table($key1);
	 
	if (scalar (keys %{$self->buffer->{$key1}}) > $self->cache_limit ) {
		my $t =time;
		$self->save_buffer($key1);
		die() if scalar (keys %{$self->buffer->{$key1}}) >0;
		delete $self->buffer->{$key1};
	}
	
}

sub save_buffer{
	my($self,$key1) = @_; 
	my $hbuffer = $self->buffer->{$key1} ;
	$self->dbh($key1)->begin_work;
		foreach my $k (keys %{$hbuffer}){
			$self->sth_insert_cached($key1)->execute($k,$self->encode($hbuffer->{$k})) ;
			delete $self->buffer->{$key1}->{$k};
		}
		$self->dbh($key1)->commit;
		
		return;
}

sub purge_buffer {
	my ($self) = @_;
	foreach my $key1 (keys %{$self->buffer}){
		$self->save_buffer($key1);
	}
}

sub get_leveldb {
	my ($self,$key1,$key2) = @_;
	return $self->decode($self->db($key1)->{$key2});
	
}

sub put_leveldb {
	my ($self,$chr,$key,$value) = @_;
	my $z = $self->encode($value);
	
	$self->db($chr)->{$key} = $z;
	
}

sub get_lite {
	my ($self,$chr,$key) = @_;
	my $rs = $self->_get($chr,[$key]);
	return $rs->{$key};
}

sub get_bulk_lite {
	my ($self,$key1,$keys) = @_;
	 my $table_name = $self->create_table($key1);
	 my $aray_ref;
	my $h;
	if ($keys){
		my $z;
		foreach my $k (@$keys){
			$k =qq{"$k"} unless  $k=~/\"/;
			push(@$z,$k);
		}
	#confess() unless $z;
	return $h unless ($z);
	my $ids = join(",",@$z);

	 $aray_ref = $self->dbh($key1)->selectall_arrayref("select _key,_value from $table_name where _key in($ids)");
	}
	else {
			 $aray_ref = $self->dbh($key1)->selectall_arrayref("select _key,_value from $table_name ;");
	}
	foreach my $a (@$aray_ref){
		$h->{$a->[0]} = $self->decode($a->[1]);
	}

return $h;
}

sub count_bulk_lite {
	my ($self,$key1,$keys) = @_;
	 my $table_name = $self->create_table($key1);
	 my $aray_ref;
	if ($keys){
		my $z;
		foreach my $k (@$keys){
			#$k =qq{"$k"};# unless  $k=~/\"/;
			push(@$z,qq{"$k"});
		}
		confess() unless $z;
	my $ids = join(",",@$z);

	 $aray_ref = $self->dbh($key1)->selectcol_arrayref("select count(*) from $table_name where _key in($ids)");
	}
	else {
			 $aray_ref = $self->dbh($key1)->selectcol_arrayref("select count(*) from $table_name ;");
	}
	

	return  $aray_ref->[0];

}

sub delete_bulk {
	my ($self,$key1,$key2) = @_;
	confess();
	my $table_name = $self->create_table($key1);
	$self->dbh($key1)->do(qq{delete  from $table_name where _key =  "$key1"});
	
	 return 1;
}

sub get_like_lite {
	my ($self,$key1,$key2) = @_;
	 my $table_name = $self->create_table($key1);
	 my $aray_ref;
	 $aray_ref = $self->dbh($key1)->selectall_arrayref("select _key,_value from $table_name where _key like \'$key2\'   ;");
#	 warn "select _key,_value from $table_name where _key like \'$key2\'   ;";
	my $h;
	foreach my $a (@$aray_ref){
		$h->{$a->[0]} = $self->decode($a->[1]);
	}
return $h;
}
sub get_like2 {
	my ($self,$key1,$type,$type2) = @_;
	
	 my $table_name = $self->create_table($key1);
	 my $aray_ref;
	 my $type3 = $type."{";
	 $aray_ref = $self->dbh($key1)->selectall_arrayref("select _key,_value from $table_name where _key between \'$type\'   and \'$type3\' ;");
	 
#	warn "select _key,_value from $table_name where _key between \'$type\'   and \'$type3\' ;";
	my $h;
	foreach my $a (@$aray_ref){
		next unless  $a->[0] =~  /$type2/;
		$h->{$a->[0]} = $self->decode($a->[1]);
	}
return $h;
}


sub get_all_lite {
	my ($self,$key1) = @_;
	#my $sth = $self->dbh($chr)->prepare("select _key,_value from $table_name");
	 my $table_name = $self->create_table($key1);
	my $aray_ref = $self->dbh($key1)->selectall_arrayref("select _key,_value from $table_name");
	my $h;
	foreach my $a (@$aray_ref){
		$h->{$a->[0]} = $self->decode($a->[1]);
	}
return $h;
}

sub get_all_for_cached {
	my ($self,$key1) = @_;
	#my $sth = $self->dbh($chr)->prepare("select _key,_value from $table_name");
	 my $table_name = $self->create_table($key1);
	my $aray_ref = $self->dbh($key1)->selectall_arrayref("select _key,_value from $table_name");
	my $h;
	foreach my $a (@$aray_ref){
		$h->{$a->[0]}->{compress} = $a->[1];
	}
return $h;
}

sub put_bulk_lite {
	my ($self,$chr,$hash)= @_;
	$self->dbh($chr)->begin_work;
	foreach my $key (keys %$hash){
		next unless $key;
		$self->sth_insert_cached($chr)->bind_param( 1, "$key");
		$self->sth_insert_cached($chr)->bind_param( 2, $self->encode($hash->{$key}) );
		
	}
	$self->sth_insert_cached($chr)->execute();
	#$self->sth_insert_cached($chr)->execute($key,$self->encode($value));
		$self->dbh($chr)->commit;
}
sub put_lite {
	my ($self,$chr,$key,$value)= @_;
	$self->dbh($chr)->begin_work;
	$self->sth_insert_cached($chr)->bind_param( 1, "$key");
	$self->sth_insert_cached($chr)->bind_param( 2, $self->encode($value));
	$self->sth_insert_cached($chr)->execute();
	$self->dbh($chr)->commit;
	
	#$self->sth_insert_cached($chr)->execute($key,$self->encode($value));
	
	
}

sub exists_lite {
	my ($self,$key1,$key2) = @_;
	 $self->sth_select_cached($key1)->execute($key2);
	return defined $self->sth_select_cached($key1)->fetchrow_array();
}
sub _get {
	my ($self,$chr,$keys) = @_;
	
	#my $sth = $self->dbh($chr)->prepare_cached("select _value from $table_name where _key =?");
	my $res= {};
	foreach my $key (@$keys ){
		 $res->{$key} = undef;
		 $self->sth_select_cached($chr)->execute($key);
		
		 
		my ($s) = $self->sth_select_cached($chr)->fetchrow_array();
		
		if ($s){
			
			$res->{$key} =$self->decode($s);
		}
	}
	# $self->sth_select_cached($chr)->finish();
	
	return $res;
}

sub put_bulk {
	my ($self,$chr,$hash,$no_decode)= @_;
	 #my $table = $self->create_table($chr);
	#$self->dbh($chr)->do("DROP  INDEX  $table._key_idx ")  or die $DBI::errstr;;
	my $t = time; 
		#$self->dbh($chr)->begin_work;
		$self->_put($chr,$hash,$no_decode);
			#$self->dbh($chr)->commit;
		 $t = time; 
}

sub sth_select_cached {
	my ($self,$chr) = @_;
	return $self->{sth_select}->{$chr} if exists $self->{sth_select}->{$chr};
	 my $table_name = $self->create_table($chr);
	$self->{sth_select}->{$chr} = $self->dbh($chr)->prepare("select _value from $table_name where _key =?")  or confess();;
	return $self->{sth_select}->{$chr};
}
sub sth_insert_cached {
	my ($self,$chr) = @_;
	return $self->{sth_insert}->{$chr} if exists $self->{sth_insert}->{$chr};
	
	 my $table = $self->create_table($chr);
	$self->{sth_insert}->{$chr} = $self->dbh($chr)->prepare("insert or replace into $table (_key,_value) values(?,?)");
	return $self->{sth_insert}->{$chr};
}
sub _put {
	my ($self,$chr,$hash,$no_decode) = @_;
	$self->dbh($chr)->begin_work;
	  foreach my $key (keys %$hash){
	  	if ($no_decode){
	  			$self->sth_insert_cached($chr)->execute($key,$hash->{$key}) ;
	  	}
	  	else {
	  			$self->sth_insert_cached($chr)->execute($key,$self->encode($hash->{$key}) ) ;
	  	}
	  
	 }
	 
	 $self->dbh($chr)->commit;
 return 1;
}

sub close {
	my ($self) = @_;
	return if ($self->is_closed());
	if (exists $self->{is_lmdb}){
		$self->lmdb->Txn->commit;
	}
	$self->purge_buffer();
	foreach my $key (keys %{$self->{sth_select}}){
		delete $self->{sth_select}->{$key};
		
	}
	foreach my $key (keys %{$self->{sth_insert}}){
		delete $self->{sth_select}->{sth_insert};
		
	}
	foreach my $key (keys %{$self->{dbh}}){
		next unless $self->{dbh}->{$key};
		$self->{dbh}->{$key}->disconnect();
		delete $self->{dbh}->{$key};
	}
	$self->is_closed(1);
}

has is_closed =>(
	is		=> 'rw',
	default => 0,
);

sub DESTROY {
	my ($self) = @_;
	$self->close();
	$self->is_closed(1);
#	foreach my $key (keys %{$self->{create}}){
#		system("chmod  a+rw  ".$self->dir()."/$key.lite 2>/dev/null") ;
#	}
   # print "In DESTROY\n";
}
1;
