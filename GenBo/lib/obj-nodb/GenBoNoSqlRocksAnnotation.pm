package GenBoNoSqlRocksAnnotation;
use Moo; 
use strict;
use warnings;
use Data::Dumper;
use JSON::XS;
use Set::IntSpan::Fast::XS;
use Carp;
use Sereal qw(sereal_encode_with_object sereal_decode_with_object write_sereal_file read_sereal_file);
extends "GenBoNoSqlRocks";





has has_config =>(
	is		=> 'ro',
default => sub {
		return 1;
}
);

# sub rocks {
#	my $self = shift;
#	return 	$self->{rocks} if exists $self->{rocks};
#	warn "coucou ".$self->json_file;
#	if($self->mode ne "r" && !(-e $self->json_file)){
#			$self->write_config();
#			die(); 
#		}
#		
#		else {
#			die() unless -e $self->json_file;
#			$self->load_config();
#		}
#		die();
#		return  $self->SUPER::rocks();
#}

has sereal_encoder => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		#return Sereal::Encoder->new();
		return Sereal::Encoder->new({compress=>Sereal::SRL_UNCOMPRESSED});
		return 0;
	},
);

has sereal_decoder => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return Sereal::Decoder->new({compress=>Sereal::SRL_UNCOMPRESSED});
		return 0;
	},
);


sub buffer {
	my ($self,$array) = @_; 
	$self->{buffer} = $self->rocks->get_multi(@$array);
}

sub compress_gnomad_position {
	my ($self,$id) = @_;
	my ($chr,$pos,$ref,$alt) = split("-",$id);
	my $l1 = length($ref);
	my $l2 = length($alt);
	
	return  ($chr,$pos,$alt) if ($l1 == 1 && $l2 ==1);
	my $seqid = $alt;
	if ($l1 ==1 && $l2 > 1){
		$seqid = "+".substr($alt, 1);
	}
	elsif ($l1 >1 && $l2 ==1){
		$ref = substr($ref, 1);
		$seqid = ($l1 -1);
	}
	
	else {
		confess($l1." ".$l2);
	}
	return ($chr,$pos,$seqid);
}

sub put_cadd {
	my ($self,$pos,$allele,$value) = @_;
	my $a2 = $self->compress_gnomad_position($allele);
	$self->batch->put($pos."!".$a2,pack("C",$value));
}

sub put_cosmic {
	my ($self,$id,$value) = @_;
	my $rid = $self->return_rocks_id_from_gnomad_id($id);
	$self->batch->put($rid,join(":",@$value));
}

sub _get_no_sereal_pack {
	my ($self,$value) = @_;
#	$aa = "!".$aa;
	 my $res = {};
	 return undef unless $value;
	 my @t = unpack($self->pack,$value);
	 my $i =-1;
	 
		foreach my $k (@{$self->description}){
				$i++;
				$res->{$k} = $t[$i];
				unless (defined $res->{$k}){
				warn scalar(@t);
				warn $value;
					
				
				 warn Dumper(@t)." :: ".$value unless defined 
				warn $i." - $k -  ".$res->{$k}." - ".$self->factor->[$i];
				die();
				}
				$res->{$k} *= $self->factor->[$i]  if ($self->factor->[$i]);
				
		}
		return $res;	
}

sub value {
		my ($self,$id) = @_;
		 my $v =  $self->_get_no_sereal_pack($self->get_raw($id));
		 return undef unless $v;
		  return $v;
}

sub cadd_score_gnomad_id {
	my ($self,$id) = @_;
	
	 my $v =  $self->_get_no_sereal_pack($self->get_raw($self->return_rocks_id_from_gnomad_id($id)));
	 return undef unless $v;
	 return $v->{cadd_score};
}

sub dbscSNV_score {
	my ($self,$id) = @_;
	 my $v =  $self->_get_no_sereal_pack($self->get_raw($id));
	 
	 return undef unless $v;
	 warn Dumper $v;
	 return $v;
}
sub dbscSNV_score_gnomad_id {
	my ($self,$id) = @_;
	
	 my $v =  $self->_get_no_sereal_pack($self->get_raw($self->return_rocks_id_from_gnomad_id($id)));
	 return $v;
}
sub cadd_score {
	my ($self,$id) = @_;
	
	 my $v =  $self->_get_no_sereal_pack($self->get_raw($id));
	 return undef unless $v;
	 return $v->{cadd_score};
}

sub polyphen_score {
	my ($self,$id,$pos,$aa) = @_;
	my $hash = $self->prediction_score($id,$pos,$aa);
	return $hash->{polyphen};
}

sub alphamissense_score {
	my ($self,$enst,$pos,$aa) = @_;
	return "-" unless $aa;
	my $array = $self->get($enst);
	return "-" unless $array;
	my $pos_aa_array = $self->hash_description->{$aa};
	return $array->[$pos]->[$pos_aa_array];
}

sub clinvar {
	my ($self,$rocksid) = @_;
	#warn $self->return_genomic_rocks_id_from_gnomad_id($id);
	my $array = $self->get($rocksid);
	return unless $array;
	my $res;
	my $i=0;
	foreach my $k (@{$self->description}){
		my $v = shift(@$array);
		if ($k eq "clinvar_id" or $k eq "score"){
			$res->{$k} =  $v;
			
		}
		elsif ($k eq "genes" ){
			foreach my $g (@$v){
				$res->{$k}->{$self->get_text_dictionary($g) } ++;
			}
		}
		else {
		$res->{$k} = $self->get_text_dictionary($v);
		}
	}
	return $res;
}


sub prediction_score {
	my ($self,$id,$pos,$aa) = @_;
	$id = $id."!".$pos;
	#warn $self->get_raw($id);
	#next;
	my $vraw = $self->get_raw($id);
	#warn $id." ".$vraw;
	#die();
	return  {sift=>"-",polyphen=>"-"} unless $vraw;
	 my $key1 = "polyphen!".$aa;
	 my $key2 = "sift!".$aa;
	 my $pos1 = $self->hash_description->{$key1};
	 my $pos2 = $self->hash_description->{$key2};
	 my @t = unpack($self->pack,$vraw);
	 
	 return ({sift=>$t[$pos2]*$self->factor->[$pos2],polyphen=>$t[$pos1]*$self->factor->[$pos1]*0.1});

}
sub sift_score {
	my ($self,$id,$pos,$aa) = @_;
	my $hash = $self->prediction_score($id,$pos,$aa);
	return $hash->{sift};
}

sub gnomad_gnomad_id {
	my ($self,$id) = @_;
	 my $v =  $self->_get_no_sereal_pack($self->get_raw($self->return_rocks_id_from_gnomad_id($id)));
	 return $v;
	 return undef unless $v;
}

sub gnomad {
	my ($self,$id) = @_;
	 my $v =  $self->_get_no_sereal_pack($self->get_raw($id));
	 return undef unless $v;
	 return $v;
}

sub revel_decode {
	my ($self,$value) = @_;
	my $res = {};
	 foreach my $l (@$value){
	 	$res->{$l->[1]} = $l->[0];
	
	 }
	return $res;
}
sub revel {
	my ($self,$id) = @_;
	 my $value =  $self->get($id);
	 return undef unless $value;
	 return $self->revel_decode($value);
}

sub cosmic {
	my ($self,$id) = @_;
	 my $value =  $self->get_raw($id);
	  return undef unless $value;
	  return $value;
}

sub spliceAI {
	my ($self,$id) = @_;
	 my $value =  $self->get($id);
	 return undef unless $value;
#	 return $self->array_score($value);
	# return $value;
	 my $end;
	 my $res = {};
	 foreach my $l (@$value){
	 	my @t = unpack($self->pack,$l);
	 	next unless @t;
		
		my $gene = $t[-1];
		my $i=0;
		my $max = -10 ;
		my $max_cat;
			foreach my $k (@{$self->description}){
				$res->{$gene}->{$k} = shift @t;
		
				$res->{$gene}->{$k} *=$self->factor->[$i]  if ($self->factor->[$i]);
				if ($i <4 ){
					if ($res->{$gene}->{$k} > $max){
						$max = $res->{$gene}->{$k}; 
						$max_cat = $k;
					}
					
					
				}
				$i++;
			}
			
			$res->{$gene}->{max} = $max;
			$res->{$gene}->{max_cat} = $max_cat;
		}
	return $res;
}
sub spliceAI_gnomad_id {
	my ($self,$id) = @_;
	 my $value =  $self->get($self->return_rocks_id_from_gnomad_id($id));
	# return $value;
	 my $end;
	 my $res = {};
	 foreach my $l (@$value){
	 	my @t = unpack($self->pack,$l);
	 	next unless @t;
		
		my $gene = $t[-1];
		my $i=0;
			foreach my $k (@{$self->description}){
				$res->{$gene}->{$k} = shift @t;
				$res->{$gene}->{$k} *= $self->factor->[$i]  if ($self->factor->[$i]);
				$i++;
			}
		}
	return $res;
}
	 



1;