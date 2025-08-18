package GenBoNoSqlLmdbCache;
use strict;
use warnings;
use DBD::SQLite;
use Moo;
use Try::Tiny;
use Data::Dumper;
use Set::IntSpan::Fast::XS;
use Module::Load;
use Compress::Snappy;
use JSON::XS;

#use Gzip::Faster;

extends "GenBoNoSqlLmdb";

#use Compress::Zlib qw( zlib_version); 
#use Storable qw/thaw freeze/;
has extension =>(
	is		=> 'rw',
default => sub {
		return "lmdb.cache";
		
}
);



sub put_cache {
	my ($self,$key1,$key2,$limit,$debug) = @_;
	require "Compress/Zstd.pm";
	$limit = 48 unless $limit;
	$limit = $limit *60 *60;
 	my $end = time + $limit;
 	$self->put($key1,{expiration=>$end,cache=>$key2,creation=>time});
 	return;
}
sub save_cache_hash {
	my ($self,$key1,$key2,$limit) = @_;
	my $end =-1;
	if ($limit){
		$limit = $limit *60 *60;
		$end = time + $limit;
	}
	my $text = encode_json ($key2);
	require "Compress/Zstd.pm";
	$self->put($key1,{expiration=>$end,cache=>Compress::Zstd::compress($text),snappy_html=>3});
}

sub put_cache_hash {
	my ($self,$key1,$key2,$limit,$debug) = @_;
#	require "Compress/Zstd.pm";
	my $v = $self->get_cache($key1);
	if($v){
		my $kk = $key1."::".rand(100)."::".time;
		$self->save_cache_hash($kk,$v);
	}
	$self->save_cache_hash($key1,$key2,$limit);
 	return;
}

sub save_cache_text {
	my ($self,$key1,$key2,$limit,$debug) = @_;
	my $end =-1;
	if ($limit){
		$limit = $limit *60 *60;
		$end = time + $limit;
	}
	require "Compress/Zstd.pm";	
	$self->put($key1,{expiration=>$end,cache=>Compress::Zstd::compress($key2),snappy_html=>4,date=>time});
	#$self->put($key1,{expiration=>$end,cache=>Compress::Zlib::compress($key2),snappy_html=>1,date=>time});
}
sub put_cache_text {
	my ($self,$key1,$key2,$limit,$debug) = @_;
#	require "Compress/Zstd.pm";
	my $v = $self->get_cache($key1);
	if($v){
			my $kk = $key1."::".rand(100)."::".time;
			$self->save_cache_text($kk,$v);
			
	}
	$self->save_cache_text($key1,$key2,$limit);
 	return;
}

sub get_cache {
	my ($self,$key,$debug) = @_;
 	my $h = $self->get($key);
 	warn $x;
 	return undef unless $h;
 	warn $h->{snappy_html};
 #	if ($h->{expiration}<time && $h->{expiration} > 0){
 #		$self->del($key);
 #		return "";
 #	}
 	if  ( $h->{snappy_html} == 1 ){
 		return Compress::Zlib::uncompress($h->{cache});
 	}
 	elsif ( $h->{snappy_html} == 2 ) {
 		my $x =  decode_json(Compress::Zlib::uncompress($h->{cache}) );
 		return $x;
 	}
 	 elsif ( $h->{snappy_html} == 3 ) {
 		my $x =  decode_json(Compress::Zstd::uncompress($h->{cache}) );
 		return $x;
 	}
 	elsif ( $h->{snappy_html} == 4 ) {
 		my $x;
			
 		 $x =  Compress::Zstd::uncompress($h->{cache});
 		 warn "ici";
 		 warn length($x);
 		return ($x);
 	}
 	else {
 		return $h->{cache};
 	}
 	die();
 	
 	return;
}
1;
