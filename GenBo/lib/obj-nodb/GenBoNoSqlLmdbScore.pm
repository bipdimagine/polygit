package GenBoNoSqlLmdbScore;
use strict;
use warnings;
use DBD::SQLite;
use Moo;

use Data::Dumper;
 use JSON::XS;
 use POSIX;
 use LMDB_File qw(:flags :cursor_op :error);
extends "GenBoNoSqlLmdb";


has score_description => (
	is		=> 'ro',
	lazy=> 1,
	default => sub {
		my $self = shift;
		$self->parse_description_file();
		return $self->{score_description};
	}
	
	
);
has pack_string => (
	is		=> 'rw',
	lazy=>1,
	default => sub {
		my $self = shift;
		$self->parse_description_file();
		return $self->{pack_string};
	}
	#required=> 1,
	
);
has factor => (
	is		=> 'rw',
	lazy=>1,
	default => sub {
		my $self = shift;
		$self->parse_description_file();
		return $self->{pack_string};
	}
);

has coding_type => (
	is		=> 'rw',
	lazy=>1,
	default => sub {
		my $self = shift;
		$self->parse_description_file();
		return "hash_alleles" unless exists $self->{type};
		return $self->{type};
	}
);

sub parse_description_file {
	my ($self) = @_;
	confess("no description file for".$self->dir) unless -e $self->description_file; 
	open(JSON, $self->description_file) or confess();
 	my $desc = decode_json <JSON>;
 	confess($self->description_file) unless $desc->{pack_string};
 	die() unless $desc->{score_description};
 	die() unless $desc->{factor};
 	$self->{pack_string} = $desc->{pack_string};
 	$self->{score_description} = $desc->{score_description};
 	$self->{factor} = $desc->{factor};
 	$self->{type} = $desc->{type};
 	
  	close JSON;
 	
}
has description_file => (
	is		=> 'rw',
	#required=> 1,
	lazy=>1,
	default => sub {
		my $self = shift;
		my $file = $self->dir."/description.json";
		return $file;
	}
	
);

has nb_score => (
	is		=> 'ro',
	lazy=>1,
	#required=> 1,
	default => sub {
		my $self = shift;
		confess() unless $self->score_description;
		return scalar(@{$self->score_description});
	}
);



sub put_score {
	my ($self,$key,$value,$debug) = @_;
	confess() if $key == 0;
	die() unless @{$self->score_description};
	my @t;
	for (my $i=0; $i<@{$self->score_description};$i++ ){
		my $k = $self->score_description->[$i];
		push(@t,int($value->{$k}*$self->factor->[$i]));
	}
	my $rec = pack($self->pack_string,@t);
	$self->put($key,$rec);
}

sub get_allele_score {
	my ($self,$key,$allele,$debug) = @_;
	my $rec = $self->get($key);
	
	return undef unless $rec;
	return undef unless exists  $rec->{$allele};
	return $self->unpack_score($rec->{$allele});
}

sub get_score {
	my ($self,$key,$debug) = @_;
	my $rec = $self->get($key);
	return undef unless $rec;
	return $self->unpack_score_hash_alleles($rec) if (ref($rec));
	return $self->unpack_score($rec);
}

sub unpack_score_hash_alleles{
	my ($self,$rec,$debug) = @_;
	my $final;

	foreach my $k (keys %$rec){
			$final->{$k} = $self->unpack_score($rec->{$k});
	}
	return $final;
}
	
sub pack_score{
	my ($self,$hscore) = @_;
	
	my @val;# = unpack($self->pack_string,$rec);
	
	foreach my $k (@{$self->score_description}){
		confess($k) unless exists $hscore->{$k};
		push(@val,$hscore->{$k});
	}
	
	return pack($self->pack_string,@val);
}

sub unpack_score{
	my ($self,$rec,$debug) = @_;

	my @val = unpack($self->pack_string,$rec);
	my $i=0;
	my $h;
	foreach my $k (@{$self->score_description}){
		$h->{$k} = $val[$i];
		$h->{$k} = $h->{$k}/$self->factor->[$i] if  $self->factor->[$i] > 1 ;
		$i++;
	}
	return $h;
}



has first => (
	is		=> 'ro',
	lazy=>1,
	#required=> 1,
	default => sub {
		my $self = shift;
		my $fid;
		my $data;
		$self->last();
		$self->cursor->get($fid, $data,MDB_FIRST);
		$self->current({id=>$fid,data=>$self->unpack_score($data)});
		   
		return $fid;
	}
);
     
has last => (
	is		=> 'ro',
	lazy=>1,
	#required=> 1,
	default => sub {
		my $self = shift;
		my $fid;
		my $data;
		$self->cursor->get($fid, $data, MDB_LAST);
		#test if position eq "-1" , j'ajoutais ca pour fichier fide avec une value de  "coucou"
		if ($fid eq -1){
			$self->cursor->get($fid, $data,MDB_PREV);
			
		}
		my $cid;
		$self->cursor->get($cid, $data,MDB_FIRST);
		return $fid;
	}
);

has cursor => (
	is		=> 'ro',
	lazy=>1,
	#required=> 1,
	default => sub {
		my $self = shift;
		my $db_name = $self->name();
		my $data;
		my $cursor = $self->lmdb($db_name)->Cursor();
		return $cursor;
	}
);
has current => (
	is		=> 'rw',
	
	
);

sub next {
	my ($self) = @_;
	return unless  $self->current();
	my $key =  $self->current()->{id};
	my $res = $self->current()->{data};
	
	return unless $key;
	 my $cid;
	my $data;
	
	if ($key eq $self->last){
		$self->current(undef,undef);	
	}
	else {
		$self->cursor->get($cid, $data, MDB_NEXT);
		
		#warn $data;
		#warn Dumper $self->transform_score($data);
		$self->current({id=>$cid,data=>$self->unpack_score($data)});
	}
	
	return($key,$res);	
}

1;

