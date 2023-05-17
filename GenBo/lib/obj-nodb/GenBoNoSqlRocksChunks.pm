package GenBoNoSqlRocksChunks;
use Moo; 
use strict;
use warnings;
use Data::Dumper;
extends "GenBoNoSqlRocks";

has start =>(
	is		=> 'ro',
	required=>1,

);
has genomic =>(
	is		=> 'ro',
	default => undef,

);

has compress =>(
	is		=> 'ro',
);
has integer =>(
	is		=> 'ro',
	default => undef,
);
has index =>(
	is		=> 'ro',
	default => undef,
);




#has index_hash =>(
#	is		=> 'ro',
#	default => {},
#); 


sub put_batch {
	my ($self,$key,$value,$debug) = @_;
	$key = $self->index_position($key);
	if ($self->nb%50_000 == 0){
		$self->write_batch();
		$self->nb(1);
		
		#$self->rocks->compact_range();
	}
	$self->nb($self->nb+1);
	if ($self->compress){
		
		return $self->batch->put($key,$self->encode($value));
	}
	else {
		return $self->batch->put($key,$self->encode($value));
	}
}

sub put_batch_raw {
	my ($self,$key,$value,$debug) = @_;
	#$key = $self->index_position($key);
	if ($self->nb%50_000 == 0){
		$self->write_batch();
		$self->nb(1);
		
		#$self->rocks->compact_range();
	}
	$self->nb($self->nb+1);
	$self->batch->put($key,$value);
}
sub index_position {
	my ($self,$key) = @_;
	
	return $key unless $self->index;
	return sprintf("%010d", $key) if $self->index eq "genomic";
	if ($self->index eq "vcf") {
		my @t = split("_",$key);
		$t[0] = sprintf("%010d",$t[0]);
		my $l1  =length($t[1]);
		my $l2  =length($t[2]);
		if($l1==$l2){
			return $t[0]."_".$t[2]; 
		}
		elsif ($l1==1 && $l2 >1 ){
			return $t[0]."+".substr($t[2], 1);
		}
		elsif($l1==1 && $l2 >1 ){
			return $t[0]."-".($l1-1);
		}
		else {
			confess();
		}
		}
	confess();
}
 sub put {
	my ($self,$key,$value) = @_;
	$key = $self->index_position($key) ;
		return  $self->SUPER::put($key,$value);
}

 sub get {
	my ($self,$key) = @_;
	$key = $self->index_position($key) ;
	return  $self->SUPER::get($key);
};

sub get_raw {
	my ($self,$key) = @_;
	return  $self->rocks->get($key);
}

has intspan_keys => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $intspan = $self->SUPER::get("&intspan");
		unless ($intspan){
			 return Set::IntSpan::Fast::XS->new( );
		}
	},
);





sub close {
	my ($self) =@_;
	if ($self->mode ne "r"){
	#if (exists $self->{intspan_keys}){
	#$self->SUPER::put("&intspan",$self->{intspan_keys});
	#}
	
	if (exists $self->{batch}){
	#	warn "write ".$self->path_rocks;
		$self->rocks->write($self->batch);
	#	warn "end ".$self->path_rocks;
	}
	
		$self->rocks->compact_range();
	}
	warn "close chunks ".$self->path_rocks;
	if ($self->pipeline){
		my $dir_prod = $self->dir."/".$self->name.".rocksdb";
		system("rsync -rav --remove-source-files ".$self->path_rocks." $dir_prod");
			system("rmdir ". $self->path_rocks);
		die() if $? ne 0;
		$self->pipeline(undef);
	}
	#$self->DESTROY();
	#$self->rocks->close();
	delete $self->{rocks};
	$self->{rocks} = undef;
}


1;