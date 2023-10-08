package GenBoNoSqlRocksPolyviewerVariant;
use FindBin qw($RealBin);
use Moo; 
use strict;
use Data::Dumper;
use JSON::XS;
use Digest::MD5 qw(md5_hex);
use Carp;
use File::Basename;
extends "GenBoNoSqlRocks";



has has_config =>(
	is		=> 'ro',
default => sub {
		return 1;
}
);
sub get_id_dictionary {
	my ($self,$value) =@ _;
	return $value;
}
sub get_text_dictionary {
	my ($self,$value) =@ _;
	return $value;
}

sub array_keys {
	my($self,$type) = @_;
	unless (exists $self->{array_keys} ){
		$self->load_config() if $self->mode ne "c";
	}
	
	confess() unless $type;
	return $self->{array_keys}->{$type} if exists $self->{array_keys}->{$type};
	$self->{array_keys}->{$type} = [];
	return $self->{array_keys}->{$type};
}


sub hash_keys {
	my ($self,$type) = @_;
	my $array = $self->array_keys($type);
	my $hash ={};
	for (my $i=0; $i<@$array;$i++){
		$hash->{$array->[$i]} = $i;
	}
	return $hash;
}


sub md5_array_keys {
	my($self,$type) = @_;
	return $self->{md5_array_keys}->{$type} if exists $self->{md5_array_keys}->{$type};
	 $self->{md5_array_keys}->{$type} = md5_hex(join(",",@{$self->array_keys($type)}));
	 return $self->{md5_array_keys}->{$type};
}

sub test_keys {
	my ($self,$pv,$type) = @_;
	my $hash = $self->hash_keys($type);
	my @new;
	foreach my $k (keys %$pv){
		push(@new,$k) unless exists $hash->{$k};
		#push(@{$self->array_keys($type)},$k) unless exists $hash->{$k};
	}
	foreach my $s (sort {$a cmp $b} @new){
		push(@{$self->array_keys($type)},$s);
	}
	return 1;
}

sub index_string {
	my ($self,$index,$chr) =@_;
	
	my $id =  sprintf("%010d", $index);
	if ($chr){
		 if( $chr eq "MT"){
		 	$chr = "MT>" ;
		 }
		 elsif( $chr eq "X"){
		 	$chr = "0X>"
		 }
		  elsif( $chr eq "Y"){
		 	$chr = "0Y>"
		 }
		 else {
		 		$chr = sprintf("%02d", $chr).">";
		 }
		$id = $chr.$id;
	}
	return $id;
}


sub PolyviewerVariant {
	my $dir = dirname(__FILE__);
	#require ($dir."/polyviewer/PolyviewerVariant.pm");
	return PolyviewerVariant->new();
}
has index_text_key => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return  {"name"=>1,"ccds"=>1,"appris"=>1,"prot"=>1,"nm"=>1,"id"=>1,"index_text_key"=>1,"impact_score_text"=>1};
	},
);

sub put_batch_gene {
	my ($self,$index,$hgenes) = @_;
	foreach my $gid (keys %$hgenes){
		my $tr = delete $hgenes->{$gid}->{tr};
		my $score = delete $hgenes->{$gid}->{score};
		my @all_score;
		$self->test_keys($score,"score");
		foreach my $k (@{$self->array_keys("score")}){
			my $value = delete $score->{$k};
			 push(@all_score,$value);
		}
	
	my @all;
	foreach my $l (@$tr){
		my @avalues;
		$self->test_keys($l,"transcripts");
		foreach my $k (@{$self->array_keys("transcripts")}){
		
			my $value = delete $l->{$k};
			if (exists $self->index_text_key->{$k}){
				$value = "-" unless defined $value ;
				$value = $self->get_id_dictionary($value);
			}
		 	push(@avalues,$value);
		}
		confess(Dumper $l) if keys %$l;
		push(@all,\@avalues);
	}
	 
	my $vgid =  $self->get_id_dictionary($gid);
	$self->put_batch($index."*".$vgid,[\@all_score,\@all]);
	#$self->put_batch($self->index_string($index,$chr)."*".$vgid,[\@all_score,\@all]);
	}
	return 1;
}

sub put_batch_patient {
	my ($self,$index,$hv) = @_;

	foreach my $patient_id (keys %{$hv}){
		my @all;
		my $array = $hv->{$patient_id};
		foreach my $l (@$array){
			my @avalues;
			$self->test_keys($l,"patients");
			my $debug;
			$debug =1 if $l->{pr};
			foreach my $k (@{$self->array_keys("patients")}){
				my $value = undef;
				if (exists  $l->{$k}){
					 $value = delete $l->{$k};
				}
					if (exists $self->index_text_key->{$k}){
						$value = $self->get_id_dictionary($value);
						}
		 			push(@avalues,$value);
		 			last unless keys %$l;
			}
			die(Dumper $l) if keys %$l;
			push(@all,\@avalues);
		}
		
		$self->put_batch($index."=".$patient_id,\@all);
		#$self->put_batch($self->index_string($index,$chr)."=".$patient_id,\@all);
	}
	
}
sub put_batch_PolyviewerVariant {
	my ($self,$pv) = @_;
	$self->test_keys($pv,"global");
	my $debug;
	$self->put_raw($pv->{global_vector_id},$pv->id);
	my @avalues;
	my $ids = delete $pv->{genes_id};
	$pv->{genes_id} = [];
	foreach my $gid (@{$ids}){
		my $id_dict = $self->get_id_dictionary($gid);
		push(@{$pv->{genes_id}},$id_dict);
	}
	
	foreach my $k (@{$self->array_keys("global")}){
		 push(@avalues, $pv->{$k});
	}
	$self->put_batch($pv->id,\@avalues);
}

sub encode_pvg {
	my ($self,$value) =@_;
	return encode_json $value;
}
sub decode_pvg {
	my ($self,$value) =@_;
	return decode_json $value;
}

sub put_batch_patient_gene {
	my ($self,$patient,$gid,$pv) = @_;
	$self->put_batch($patient->id."_gene_".$gid,$pv);
	return 1;
}


sub get_patient_gene {
	my ($self,$patient,$gid) = @_;
	my $array_global = $self->get($patient->id."_gene_".$gid);
	return $array_global;
	return 1;
}


sub get_id_patient_variant_genes{
	my ($self,$pid,$vid) =@_;
	return $pid."/".$vid;
}

sub put_batch_patient_variant_genes {
	my ($self,$pid,$vid,$pv) = @_;
	my $aig = my $ain = delete $pv->{array};
	#die(Dumper $pv) if scalar(keys %$pv);
	#$self->batch->put($vid."/".$pid,$self->encode_pvg($pv));
	$self->put_batch($self->get_id_patient_variant_genes($pid,$vid),$aig);
	return 1;
	
}

sub _get_patient_variant_genes {
	my ($self,$pid,$vid) = @_;
	my $array_global = $self->get($vid."/".$pid);
	return {array=>$array_global};
}

sub cache_variant {
	my ($self,$pid) = @_;
	return $self->{cached} if exists $self->{cached};
	my $search = $self->get_id_patient_variant_genes($pid,"");
		my $iter = $self->rocks->new_iterator->seek($search);
		my $vp = {};
		
		while (my ($key, $value) = $iter->each) {
    		last if $key !~ /$search/;
    		$self->{cached}->{$key} = $value;
		}
		return $self->{cached};
}

sub get_patient_variant_genes {
	my ($self,$pid,$vid) = @_;
	my $cache = $self->cache_variant($pid);
	my $id = $self->get_id_patient_variant_genes($pid,$vid);
	confess unless exists $self->cache_variant($pid)->{$id};
	return {array=>$self->decode($self->cache_variant($pid)->{$id})};
	my $array_global = $self->get($self->get_id_patient_variant_genes($pid,$id));
	return {array=>$array_global};
	
}

sub _get_gene {
	my ($self,$array_global,$vp,$debug) = @_;
	my $array_score = $array_global->[0];
	for (my $i=0;$i<scalar(@{$self->array_keys("score")});$i++){
		my $k = $self->array_keys("score")->[$i];
		my $value =  $a->[$i];
		$vp->{$k} = $value;
	}
	
	my $array = $array_global->[1];
	confess() unless $array;

	my $all;
	foreach my $a (@$array){
		my $h1;
	for (my $i=0;$i<scalar(@{$self->array_keys("transcripts")});$i++){
			my $k = $self->array_keys("transcripts")->[$i];
			my $value =  $a->[$i];
			if (exists $self->index_text_key->{$k}){
				$value = $self->get_text_dictionary($value);
			}
			$h1->{$k} = $value;
		}
		$h1->{enst} = $h1->{name};
		$h1->{codons_AA} ="XXX";
		push(@$all,$h1);
	}
	$vp->{transcripts} = $all;
}
sub getGene {
	my ($self,$index,$gene_id,$vp) = @_;
	my $r_gene_id = $self->get_id_dictionary($gene_id);
	my $array_global = $self->get($self->index_string($index)."*".$r_gene_id);
	$self->_get_gene($array_global,$vp);
	return 1;
}

sub testPolyviewerVariant{
	my ($self,$key) = @_;
	my $value = $self->rocks->get($key);
	
	my $array = $self->decode($value);
	warn Dumper $array;
	my $vp = {};
    $self->_getPolyviewerVariant($array,$vp);
   return bless $vp , 'PolyviewerVariant';
}

sub test {
	my ($self,$var_id,$patient_id,$gene_id) = @_;
	my $pid = $patient_id;
		my $debug;
		$debug =1 if $gene_id =~ /ENSG00000171759/;
		my $toto =  $self->get($var_id);
		my $index = $var_id;
		#my $index =  $self->get_raw($var_id);
		#my $sindex = $self->index_string($index);
		my $iter = $self->rocks->new_iterator->seek($index);
		$patient_id =$index."=".$patient_id;
		$gene_id =$index."*".$gene_id;
		my $vp = {};
		while (my ($key, $value) = $iter->each) {
    		last if $key !~ /$index/;
    		if($key eq $index){
    				my $array = $self->decode($value);
    				$self->_getPolyviewerVariant($array,$vp);
    				next;
    		}
    		
    		if($key  eq $patient_id) {
    			my $array = $self->decode($value);
    			$self->_get_patient($pid,$array,$vp);
    			next;
    		}
    		if($key  eq $gene_id) {
    			my $array = $self->decode($value);
    			$self->_get_gene($array,$vp,$debug);
    			next;
    		}
		}
		confess($var_id ) unless %$vp;  
		return bless $vp , 'PolyviewerVariant';
	
}

sub _get_patient {
	my ($self,$pid,$array,$vp) = @_;
	
	my $h;
	foreach my $a (@$array){
		my $h1;
	for (my $i=0;$i<scalar(@{$self->array_keys("patients")});$i++){
			my $k = $self->array_keys("patients")->[$i];
			next unless $k;
			my $value = $a->[$i];
				if (exists $self->index_text_key->{$k}){
						$value = $self->get_text_dictionary($value);
				}
			$h1->{$k} = $value;
		}
		$h->{$h1->{id}} = $h1;
	}
	$vp->{patients_calling} = $h;
	$vp->{text_caller} = $h->{$pid}->{array_text_calling};
	return 1;
	
} 
sub getPatient {
	my ($self,$index,$patient_id,$vp) = @_;
	my $array = $self->get($self->index_string($index)."=".$patient_id);
	confess($self->index_string($index)."=".$patient_id) unless $array;
	$self->_get_patient($patient_id,$array,$vp);
	
}


sub _getPolyviewerVariant {
	my ($self,$array,$vp) = @_;
	
	for (my $i=0;$i<scalar(@{$self->array_keys("global")});$i++){
			my $k = $self->array_keys("global")->[$i];
			$vp->{$k} = $array->[$i];
	}
	for (my $i=0;$i<@{$vp->{genes_id}};$i++){  
		my $value = $vp->{genes_id}->[$i];
		$vp->{genes_id}->[$i] = $self->get_text_dictionary($value);
		$vp->{genes}->{$vp->{genes_id}->[$i]} = 1;
	}
	
}



sub getPolyviewerVariant {
	my ($self,$key) = @_;
	
	#$key =  $self->get_raw($key) if ($key =~ /!/);
	my $vector_id =  $self->get_raw($key);
	
	my $array =  $self->get_raw($self->index_string($vector_id));
	my $vp = {};
	confess($vector_id." ".$key) unless $array;
	
	$self->_getPolyviewerVariant($array,$vp);
	$vp->{vector_id} = $vector_id;
	return bless $vp , 'PolyviewerVariant';
}

sub load_config {
	my ($self) = @_;
	open(my $fh ,$self->json_file) or confess("can t open file");
	my $json_str = do { local $/; <$fh> };
	close($fh);
	my $h = decode_json($json_str);
	
	
	if( $self->mode ne "c"){
		foreach my $k (keys %$h){
			$self->{$k} = delete $h->{$k};
		}	

	}
	return delete $h->{chunks};
}

sub write_batch {
	my ($self) = @_;
	return unless $self->batch;
	$self->rocks->write($self->batch);
	delete $self->{batch};
}
sub write_config {
	my ($self) = @_;
	open(my $fh,">".$self->json_file) or die("can t open file");
	my $h;
	$h->{date} = time;
	$h->{version} = $self->version;
	$h->{pack} = $self->pack;
	$h->{description} = $self->description;
	$h->{factor} = $self->factor;
	$h->{array_keys} = $self->{array_keys};
	$h->{index_text_key} = $self->index_text_key;
	print $fh encode_json($h);
	close ($fh);
}

1;
