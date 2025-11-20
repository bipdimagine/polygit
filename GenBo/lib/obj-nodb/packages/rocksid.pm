package rocksid;
use strict;




sub return_rocks_id_from_gnomad_id {
	my ($id) = @_;
	my ($chr,$pos,$ref,$alt) = split("-",$id);
	return return_rocks_id($pos,$ref,$alt);
}
sub return_genomic_rocks_id_from_gnomad_id {
	my ($id) = @_;
	confess() unless $id;
	my ($chr,$pos,$ref,$alt) = split("-",$id);
	return $chr."!".return_rocks_id($pos,$ref,$alt);
}

sub return_rocks_id_from_genbo_id {
	my ($id) = @_;
	my ($chr,$pos,$ref,$alt) = split("_",$id);
	
	return return_rocks_id($pos,$ref,$alt);
}
sub stringify_pos {
	my ($pos) = @_;
	return ($pos,sprintf("%010d", $pos));
}


sub return_rocks_id {
	my ($pos,$ref,$alt) = @_;
	my $l1 = length($ref);
	my $l2 = length($alt);
	return  (stringify_pos($pos)."!".$alt) if ($l1 == 1 && $l2 ==1);
	my $seqid = $alt;
	if ($alt =~ /del/ ){
		return  (stringify_pos($pos)."!".$alt);
	}
	elsif ($ref =~ /del/ ){
		return  (stringify_pos($pos)."!".$ref);
	}
	elsif ($alt=~ /inv/ ){
		return  (stringify_pos($pos)."!".$alt);
	}
	elsif ($alt=~ /ins/ ){
		return  (stringify_pos($pos)."!".$alt);
	}
	elsif ($alt=~ /dup/ ){
		return  (stringify_pos($pos)."!".$alt);
	}
	
	
	die($ref." ".$alt) if $alt =~ /INV/;
	#die(ref." ".$alt) if $alt =~ /DEL/;s
	die(ref." ".$alt) if $alt =~ /DUP/;
	
	
	if ($l1 ==1 && $l2 > 1){
		
		$seqid = "+".substr($alt, 1);
		return  (stringify_pos($pos)."!".$seqid);
	}
	elsif ($l1 >1 && $l2 ==1){
		$ref = substr($ref, 1);
		$seqid = ($l1 -1);
		return  (stringify_pos($pos)."!".$seqid);
	}
	 elsif ($l1 >1 && $l2 == $l1 && $l2>1 ){
	 	
		$ref = substr($ref, 1);
		$alt = substr($alt, 1);
		$seqid = "$ref*$alt";
		return  (stringify_pos($pos)."!".$seqid);
	}
	
	else {
		return  (stringify_pos($pos)."!".$ref."*".$alt);
		confess($l1." ".$l2);
	}
	
}
1;
