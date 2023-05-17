package GenBoBitVector;
use Moo;
use Bit::Vector;
use Bit::Vector::Overload;
use Set::IntervalTree;
use Set::IntSpan::Fast::XS ;
use Array::IntSpan;
use Data::Dumper;

has length =>(
	is		=> 'ro',
	required =>1,
);
has size =>(
	is		=> 'rw',
	default => sub {
		return 1_000_000;
	}
);

has nb_chunks =>(
	is		=> 'ro',
	lazy =>1,
	default => sub {
		my $self =shift;
			my $from =1;
    		my $to = $self->length;
    	my $start;
    	my $end;
    	my $part = 0;
       while ($from < $to){
        $start = $from;
        $end = $from + $self->size;
        if ($end > $to){
            $end = $to;
        }
		
		$part ++;
 
        $from = $end;
       }
	
	#warn scalar(@$regions);
	return $part;
		
	}
);
has chunks =>(
	is		=> 'ro',
	lazy =>1,
	default => sub {
		my $self =shift;
		my $tree_AC = Set::IntervalTree->new;
			my $from =1;
    		my $to = $self->length;
    		my $regions =[];
    	my $start;
    	my $end;
    	my $part = 0;
       while ($from < $to){
        $start = $from;
        $end = $from + $self->size;
        if ($end > $to){
            $end = $to;
        }
		#my $inter = $intspan->intersection($intspan_region);
		#warn $inter->as_string;
		#warn $inter->as_string;
		#warn "coucou" if $inter->is_empty();
		#$self->array_intspan->set_range($from,$end-1,$part);
		$tree_AC->insert($part,$from,$end+1);
		#push(@$regions,$from);
		#push(@$regions,{start=>$from,end=>$end-1,part=>$part});
		$part ++;
      #  push(@$res,{start=>$from,end=>$end,intspan=>$intspan_region,ext_gvcf=>$self->ucsc_name.".$from.$end.g.vcf",chromosome=>$self->ucsc_name}) unless $intspan_region->is_empty;
      
        #print chrom_name + ":" + str(region_start) + "-" + str(end)
        $from = $end;
       }
	
	#warn scalar(@$regions);
	return $regions;
		
		
	}
);

has array_intspan =>(
	is		=> 'ro',
	#lazy =>1,
	#weak_ref => 1,
	default => sub {
	my $self =shift;
	
	return Array::IntSpan->new();
	}
);
has tree =>(
is		=> 'ro',
	lazy =>1,
	#weak_ref => 1,
	default => sub {
	my $self =shift;
	my $foo = Array::IntSpan->new();
	
	my $tree_AC = Set::IntervalTree->new;

	foreach my $chunk (@{$self->chunks}){
		$foo->set_range($chunk->{start},$chunk->{end},$chunk->{part});
		$tree_AC->insert($chunk,$chunk->{start},$chunk->{end}+1);
	}
	return $foo;
	}
);

has intspan_chunk =>(
is		=> 'ro',
	lazy =>1,
	#weak_ref => 1,
	default => sub {
	return Set::IntSpan::Fast::XS->new;
	}
);

has hchunk =>(
is		=> 'ro',
	lazy =>1,
	#weak_ref => 1,
	default => sub {
	return {};
	}
);
sub chunk_vector {
	my($self,$nb) = @_;
	
	return $self->hchunk->{$nb}->{bitvector} if exists $self->hchunk->{$nb}->{bitvector};
	$self->hchunk->{$nb}->{bitvector} =  Bit::Vector->new($self->size); 
	$self->intspan_chunk->add($nb);
	return $self->hchunk->{$nb}->{bitvector};
}


sub add {
	my ($self,$pos) = @_;
	$pos ++;
	#my $res = $self->tree->fetch($pos,$pos+1);
	#die() if scalar(@$res) ne 1;
	my $posa = int($pos / $self->size);
	my $tpos = $pos - ($posa*$self->size);
	$self->chunk_vector($posa)->Bit_On($tpos);	
}

sub union {
	my ($self,$v) = @_;
	die() if scalar(@{$self->chunks}) ne scalar(@{$v->chunks});
	
	my $a1 = $self->intspan_chunk->intersection($v->intspan_chunk);
	my @a = $a1->as_array;
	warn scalar(@a) if (scalar(@a)>0);
	foreach my $i (@a){
		my $ca = $self->hchunk->{$i};#$i];
		my $cb = $v->hchunk->{$i};
		$ca->{bitvector} += $cb->{bitvector};
	}
	my $a2 = $v->intspan_chunk->diff($self->intspan_chunk);
	foreach my $i ($a2->as_array){
		my $cb = $v->hchunk->{$i};
		$self->hchunk->{$i}->{bitvector} = $cb->{bitvector}->Clone;
		$self->intspan_chunk->add($i);
	}
	
}

sub intersect {
	my ($self,$v) = @_;
	die() if scalar(@{$self->chunks}) ne scalar(@{$v->chunks});
	my $a1 = $self->intspan_chunk->intersection($v->intspan_chunk);
	foreach my $i ($a1->as_array){
		my $ca = $self->hchunk->{$i};
		my $cb = $v->hchunk->{$i};
		$ca->{bitvector} &= $cb->{bitvector};
		my $nb = 	 unpack("%32b*", $ca->Block_Read);
		if ($nb ==0 ){
			delete $ca->{bitvector};
			$self->intspan_chunk->remove($ca->{part});
		}
	}
	
	my $a2 = $self->intspan_chunk->diff($v->intspan_chunk);
	foreach my $i ($a2->as_array){
		my $ca = $self->hchunk->{$i};
		delete $ca->{bitvector};
		$self->intspan_chunk->remove($ca->{part});
	}
	
#	for (my $i=0;$i<@{$self->chunks};$i++){
#		my $ca = $self->chuncks->[$i];
#		my $cb = $v->chuncks->[$i];
#		if (exists  $ca->{bitvector} && exists $cb->{bitvector}){
#			$ca->{bitvector} &= $cb->{bitvector};
#		}
#		elsif (exists $ca->{bitvector}){
#				delete $ca->{bitvector};
#				$self->intspan_chunk->remove($ca->{part});
#		}
#		
#	}
}

sub vector {
	my ($self) = @_;
	my @cs;
	for (my $i=0;$i<$self->nb_chunks;$i++){
		if (exists $self->hchunk->{$i}->{bitvector}){
			push(@cs,$self->hchunk->{$i}->{bitvector});
		}
		else {
			push(@cs,Bit::Vector->new($self->size));	
		}
	}
	return   Bit::Vector->Concat_List(reverse @cs);
}
sub intspan {
	my ($self) = @_;
	my $string = $self->vector->to_Enum();
	return Set::IntSpan::Fast->new($string);
}

sub nb {
	my ($self) = @_;
	my $nb =0;
	foreach my $c (values %{$self->hchunk}){
		next unless exists ($c->{bitvector});
		$nb += 	 unpack("%32b*", $c->{bitvector}->Block_Read);
		
	}
	return $nb;
	
}
sub del {
	my ($self) = @_;
	confess();
	my $nb =0;
	foreach my $c (@{$self->chunks}){
		next unless exists ($c->{bitvector});
		delete $c->{bitvector};
		#$nb += 	 unpack("%32b*", $c->{bitvector}->Block_Read);
		
	}
	#return $nb;
	
}
sub to_Enum {
	my ($self) = @_;
	my $nb =0;
	my $vector =  $self->vector;
	return $self->vector->to_Enum();
	 
	
}
sub to_machin {
	my ($self) = @_;
	my $start = 0;
	my $vector = $self->vector();
	my @z;
while (($start < $self->length) &&
    (my ($min,$max) = $vector->Interval_Scan_inc($start)))
{
	push(@z,$min."-".$max);
    $start = $max + 1;
 	
    # do something with $min and $max
}
warn join(",",@z);
}


sub to_array{
	my ($self) = @_;
	return $self->intspan->as_array;
}
1;