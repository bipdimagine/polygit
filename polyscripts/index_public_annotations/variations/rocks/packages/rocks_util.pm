package rocks_util;

sub compress_vcf_position {
	my ($ref,$alt) = @_;
	my $l1 = length($ref);
	my $l2 = length($alt);
	if ($l1 ==1 && $l2 ==1){
		return $alt;
	}
	elsif ($l1 ==1 && $l2 > 1){
		return "+".substr($alt, 1);
	}
	elsif ($l1 >1 && $l2 ==1){
		$ref = substr($ref, 1);
		return ($l1 -1);
	}
	confess();
	
}

sub rocks_id {
	my ($pos,$ref,$alt) = @_;
	my $l1 = length($ref);
	my $l2 = length($alt);
	if ($l1 ==1 && $l2 ==1){
		return stringify_pos($pos)."!".$alt;
	}
	elsif ($l1 ==1 && $l2 > 1){
		return stringify_pos($pos+1)."!"."+".substr($alt, 1);
	}
	elsif ($l1 >1 && $l2 ==1){
		$ref = substr($ref, 1);
		return stringify_pos($pos+1)."!".($l1 -1);
	}
	else {
		$ref = substr($ref, 1);
		$alt  = substr($ref, 1);
		return stringify_pos($pos+1)."!".($l1 -1)."=$ref"."x"."$alt";
	}
	
	die($ref." ".$alt);
}

sub stringify_pos {
	my ($self,$pos) = @_;
	return ($pos,sprintf("%010d", $pos));
}
1;