package GenBoNoSqlText;



# ABSTRACT: An embedded, NoSQL SQLite database with SQL indexing
use Moo;

use strict;
use warnings;
use Data::Dumper;
extends "GenBoNoSql";
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;


has extension =>(
	is		=> 'rw',
default => sub {
		return "search";
		
}
);

sub create_table {
	my ($self,$key1) = @_;
	return $self->{table}->{$key1}   if (exists $self->{table}->{$key1});
	 my $table_name =$self->name_table_data;
	 $self->{table}->{$key1} = $table_name;
	 return  $self->{table}->{$key1} unless  $self->write();
	   sub compressor {
        my $in = shift;
        my $out;
        #return compress($in);
        gzip \$in => \$out;
        return ($out);
    }
    sub uncompressor {
        my $in = shift;
        my $out;
        gunzip \$in => \$out;
       # return decompress($in);
        return ($out);
    }	
	 
	$self->dbh($key1)->do("DROP TABLE IF EXISTS $table_name")  or die $DBI::errstr  if $self->mode eq 'c'  ;
	#$self->dbh($chr)->do("CREATE TABLE if not exists $table_name (_key TEXT PRIMARY KEY, _value BLOB) ;")  or die $DBI::errstr;;
	$self->dbh($key1)->do("CREATE VIRTUAL TABLE  if not exists $table_name USING fts4(_key ,_value blob ,notindexed=_value, tokenize=unicode61 'tokenchars=-') ;")  or die $DBI::errstr;
#	$self->dbh($key1)->do("CREATE VIRTUAL TABLE  if not exists $table_name USING fts3(_key VARCHAR(250) NOT NULL ,_value blob ) ;")  or die $DBI::errstr;
	
	return 	$self->{table}->{$key1} ;
}

sub delete_bulk {
	my ($self,$key1,$key2) = @_;
	my $id = $self->change_id($key1);
	my $table_name = $self->create_table($id);
	confess() unless $key2;
	 $self->dbh($id)->do(qq{delete  from $table_name where _key match "$key2"});
	 

	 return 1;
}


sub get_bulk{
	my ($self,$key1,$keys) = @_;
	
	my $table_name = $self->create_table($key1);
	 my $aray_ref;
	if ($keys){
		my $ids = join(" OR ",map{$self->change_id($_)}@$keys);
		#warn $ids;
	 	$aray_ref = $self->dbh($key1)->selectall_arrayref("select _key,_value from $table_name where _key match \"$ids\"" );
	}
	else {
		confess();
			 $aray_ref = $self->dbh($key1)->selectall_arrayref("select _key,_value from $table_name ;");
	}
	
	my $h;
	foreach my $a (@$aray_ref){
		 my $id = $self->restore_id($a->[0]);
		 my (@find) = grep {$id =~/$_/} @$keys;
		 confess("") if scalar (@find) ne 1;
		   
		$h->{$find[0]} = $self->decode($a->[1]);
	}
	
	return $h;
	
}


sub get_like_lite {
	my ($self,$key1,$key2) = @_;
	$key2 = $self->change_id($key2);
	confess() unless $key2;
	 my $table_name = $self->create_table($key1);
	 my $aray_ref;
	 	#$key2 =~s/%/\*/;
	 	 $aray_ref = $self->dbh($key1)->selectall_arrayref("select _key,_value from $table_name where _key match \'$key2 \'    ;") or confess($key1." ".$key2." ".$self->dir()." - "."select _key,_value from $table_name where _key match \'$key2 \'    ;" );
		my $h;
		foreach my $a (@$aray_ref){
			my $id = $self->restore_id($a->[0]);
			$h->{$id} = $self->decode($a->[1]);
		}
return $h;
}

sub get_for_update {
	my ($self,$key1,$key2) = @_;
	$key2 = $self->change_id($key2);
	 my $table_name = $self->create_table($key1);
	 my $aray_ref;
	 	$key2 =~s/%/\*/;
	 	 $aray_ref = $self->dbh($key1)->selectall_arrayref("select rowid,_key,_value from $table_name where _key match \' $key2 \' ;") or confess();
	my $h;
	foreach my $a (@$aray_ref){
		my $id = $self->restore_id($a->[1]);
		$h->{data}->{$id} = $self->decode($a->[2]);
		$h->{rowid} = $a->[0];
	}
return $h;
}


sub get_text {
	my ($self,$key1,$key2) = @_;
	 #my $table_name = $self->create_table($key1);
	 	
	 	$key2 = $self->change_id($key2);
	 	#my $ary_ref = $dbh->selectcol_arrayref("select id, name from table", { Columns=>[1,2] });
	 	 $self->sth_select_text($key1)->execute($key2);
	 	 my @t;
		while ( my $ref = $self->sth_select_text($key1)->fetchrow_arrayref() ) {
   			push(@t,$self->restore_id($ref->[0]));
			}
	
	#	die();
		#warn Dumper $s;
	 	#my $ary_ref = $dbh->selectcol_arrayref("select id, name from table", { Columns=>[1,2] });
	 	#my $aray_ref = $self->dbh($key1)->selectcol_arrayref("select _key from __DATA__  where _key match \'$key2\'  ;", { Columns=>[1] }) or confess();
	 #	warn Dumper $aray_ref;
	 	#warn "------";
	 	#die();
	 	return \@t ;

}

sub sth_select_text {
	my ($self,$chr) = @_;
	return $self->{sth_text}->{$chr} if exists $self->{sth_text}->{$chr};
	 my $table_name = $self->create_table($chr);
	$self->{sth_text}->{$chr} = $self->dbh($chr)->prepare("select _key from $table_name where _key match ?  ;") ;
	return $self->{sth_text}->{$chr};
}

sub sth_select_count {
	my ($self,$chr) = @_;
	return $self->{sth_count}->{$chr} if exists $self->{sth_count}->{$chr};
	 my $table_name = $self->create_table($chr);
	$self->{sth_count}->{$chr} = $self->dbh($chr)->prepare("select count(*) from $table_name where _key  match ?") ;
	return $self->{sth_count}->{$chr};
}

sub get_count {
	my ($self,$key1,$key2) = @_;
	# my $table_name = $self->create_table($key1);
	 	$key2 = $self->change_id($key2);
	 	#my $ary_ref = $dbh->selectcol_arrayref("select id, name from table", { Columns=>[1,2] });
	 	 $self->sth_select_count($key1)->execute($key2);
		
		 
		my ($s) = $self->sth_select_count($key1)->fetchrow_array();
		
	 	#my $aray_ref = $self->dbh($key1)->selectall_arrayref("select count(_key) from $table_name where _key match \'$key2\'  ;") or confess();
	 	return  $s;

}

sub get_key_value {
	my ($self,$chr,$id) = @_;
	if ($id =~/dloop-1$/){
		$id = "genboid:".$id."_MT";
	}
	my $key = $self->change_id($id);
	#warn $id;
	my $res = $self->get_like_lite($chr,$key);
	return unless $res;
	my @values = values %$res;
	my @keys = keys  %$res;
	
	if   (scalar(@values) == 2  ) {
			return ($keys[0],$values[0]) if $values[0] eq $values[1]; 
		 
		#here for par X and Y chromosome same name so by default return transcript on X chromosome;
		my $k = join(";",@keys);
	
		#warn $k." ".$id;
		my $i1 = $id.'_X';
		my $i2 = $id.'_Y';
		if ($k=~/$i1/ && $k=~/$i2/){
			my $toto = \@values;
			 if (ref(\$values[0]) eq "SCALAR"){
			 	return ($keys[0],$values[0]) if $values[0] =~/_X/;
			 		return ($keys[1],$values[1]);
			 }
			return ($keys[0],$values[0]) if $values[0]->{chromosome} eq "X";
			return ($keys[1],$values[1]);
		}
		else {
			my @v = sort {$a cmp $b} @values;
			my @k = sort {$a cmp $b} @keys;
			if ($v[0] =~/_X/ && $v[1] =~/_Y/){
				return ($k[0],$v[0]);
			}
			
		}
		
		
	}
	confess("more than one answer for $key : ".Dumper($res)."\n".$res) if scalar(@values)>1;
	return if  scalar(@values) == 0;
	return ($keys[0],$values[0]);
} 

sub get {
	my ($self,$chr,$id) = @_;
 	my ($key,$value) = $self->get_key_value($chr,$id);
 	return $value;

	#return $rs->{$key};
}


sub get_key_values {
	my ($self,$chr,$id) = @_;
	if ($id =~/dloop-1$/){
		$id = "genboid:".$id."_MT";
	}
	my $key = $self->change_id($id);
	#warn $id;
	my $res = $self->get_like_lite($chr,$key);
}

sub get_all {
	my ($self,$chr,$id) = @_;
	if ($id =~/dloop-1$/){
		$id = "genboid:".$id."_MT";
	}
	my $key = $self->change_id($id);
	#warn $id;
	my $res = $self->get_like_lite($chr,$key);
	return unless $res;
	my @values = values %$res;
	
	return if  scalar(@values) == 0;
	return \@values;
	#return $rs->{$key};
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
	 $self->sth_select_cached($chr)->finish();
	return $res;
}


sub _put {
	my ($self,$chr,$hash) = @_;
	$self->dbh($chr)->begin_work;
	  foreach my $key (keys %$hash){
	  	my $id = $self->change_id($key);
	  	 if  (defined $hash->{$key}){
	  			$self->sth_insert_cached($chr)->execute($id,$self->encode($hash->{$key}) );
	  	 }
	  	 else {
	  	 		$self->sth_insert_cached($chr)->execute($id,undef) ;
	  	 }
	 }
	 
	 $self->dbh($chr)->commit;
 return 1;
}

sub save_buffer{
	my($self,$key1) = @_; 
	my $hbuffer = $self->buffer->{$key1} ;
		$self->dbh($key1)->begin_work;
		foreach my $k (keys %{$hbuffer}){
				my $id = $self->change_id($k);
			$self->sth_insert_cached($key1)->execute($id,$self->encode($hbuffer->{$k})) ;
			delete $self->buffer->{$key1}->{$k};
		}
			 $self->dbh($key1)->commit;
		return;
}

sub get_bulk_lite {
	my ($self,$key1,$keys) = @_;
	 my $table_name = $self->create_table($key1);
	 my $aray_ref;
	if ($keys){
	my $ids = join(",",map{$_=qq{"$_"} unless $_ =~/\"/}map{$self->change_id($_)} @$keys);
	 $aray_ref = $self->dbh($key1)->selectall_arrayref("select _key,_value from $table_name where _key in($ids)");
	}
	else {
			 $aray_ref = $self->dbh($key1)->selectall_arrayref("select _key,_value from $table_name ;");
	}
	my $h;
	foreach my $a (@$aray_ref){
		my $id = $self->restore_id($a->[0]);
		$h->{$id} = $self->decode($a->[1]);
	}
return $h;
}

sub sth_update_cached {
	my ($self,$chr) = @_;
	return $self->{sth_update}->{$chr} if exists $self->{sth_update}->{$chr};
	 my $table = $self->create_table($chr);
	$self->{sth_update}->{$chr} = $self->dbh($chr)->prepare("update $table set _value=? where rowid =? ");
	return $self->{sth_update}->{$chr};
}

sub update_hash {
	my ($self,$key1,$key2,$hash) = @_;
	$key2 = $self->change_id($key2);
	 my $table_name = $self->create_table($key1);
	 my $aray_ref;
	 	$key2 =~s/%/\*/;
	 	 $aray_ref = $self->dbh($key1)->selectall_arrayref("select rowid,_key,_value from $table_name where _key match \'$key2\'  ;") or confess();
		
	die(scalar(@$aray_ref)) if scalar(@$aray_ref) > 1;
	if (scalar(@$aray_ref) == 0 ){
		$self->put($key1,$key2,$hash);
		return;
	}
	
	#else {
	my	$h1 = $self->decode($aray_ref->[0]->[2]);
	#}
	foreach my $v (keys %$hash){
		$h1->{$v} = $hash->{$v};
	}
	$self->sth_update_cached($key1)->bind_param( 1,$self->encode($h1));
	$self->sth_update_cached($key1)->bind_param( 2, $aray_ref->[0]->[0]  );
	$self->sth_update_cached($key1)->execute();
	return $self->get($key1,$key2);
}
sub put_bulk_lite {
	my ($self,$chr,$hash)= @_;
	$self->dbh($chr)->begin_work;
	foreach my $key (keys %$hash){
		next unless $key;
		
		  	my $id = $self->change_id($key);
	
		$self->sth_insert_cached($chr)->bind_param( 1, "$id");
		$self->sth_insert_cached($chr)->bind_param( 2, $self->encode($hash->{$key}) );
		$self->sth_insert_cached($chr)->execute();
	}
	#$self->sth_insert_cached($chr)->execute($key,$self->encode($value));
	$self->dbh($chr)->commit;
}
sub put_lite {
	my ($self,$chr,$key,$value)= @_;
	my $id = $self->change_id($key);
	$self->sth_insert_cached($chr)->bind_param( 1, "$id");
	$self->sth_insert_cached($chr)->bind_param( 2, $self->encode($value));
	$self->sth_insert_cached($chr)->execute();
	#$self->sth_insert_cached($chr)->execute($key,$self->encode($value));
	
	
}


sub restore_id {
	my ($self,$id) = @_;
	$id =~s/xxx/_/g;
	return $id;
}

sub change_id {
	my ($self,$id) = @_;
	$id =~s/_/xxx/g;
	return $id;
}

1;
