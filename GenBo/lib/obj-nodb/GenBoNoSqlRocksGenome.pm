package GenBoNoSqlRocksGenome;
#use lib "$Bin/";
use strict;
use warnings;
use Moo;
use Data::Dumper;
use JSON::XS;
use POSIX;
use GenBoNoSqlRocksChunks;

my $json_chr_length = qq{{"HG19":{"6":"171115067","11":"135006516","9":"141213431","15":"102531392","14":"107349540","1":"249250621","8":"146364022","17":"81195210","7":"159138663","13":"115169878","MT":"16569","22":"51304566","Y":"59373566","3":"198022430","21":"48129895","18":"78077248","5":"180915260","X":"155270560","20":"63025520","16":"90354753","2":"243199373","12":"133851895","19":"59128983","4":"191154276","10":"135534747"},"HG38":{"21":"46709983","Y":"57227415","18":"80373285","3":"198295559","22":"50818468","MT":"16569","20":"64444167","16":"90338345","X":"156040895","5":"181538259","2":"242193529","10":"133797422","4":"190214555","19":"58617616","12":"133275309","11":"135086622","6":"170805979","14":"107043718","15":"101991189","9":"138394717","17":"83257441","8":"145138636","1":"248956422","13":"114364328","7":"159345973"}}};

has genome =>(
	is		=> 'ro',
	required=>1,

);

has chr_length =>(
	is 	=> 'ro',
	lazy => 1,
	default => sub {
		my $self = shift;
		my $h = decode_json($json_chr_length);
		die(Dumper $h) unless exists $h->{$self->genome}->{$self->chromosome_name};
		return  $h->{$self->genome}->{$self->{chromosome_name}};
	}

 );
has pipeline =>(
	is		=> 'rw',
	default => sub {
		return undef;
	}

);

has compress =>(
	is		=> 'ro',
	default => sub {
		return undef;
	}

);
has index =>(
	is		=> 'ro',
	default => sub {
		return undef;
	}

);
has pack =>(
	is		=> 'ro',
	default => sub {
		return undef;
	}
);
has description =>(
	is		=> 'ro',
	default => sub {
		return undef;
	}
);
has dir =>(
	is		=> 'rw',
	required=>1,
);

has dir_db =>(
	is		=> 'ro',
	lazy => 1,
	default => sub {
		my $self = shift;
		my $dd = $self->dir."/".$self->chromosome_name.".g.rocks/";
		unless (-e $dd){
			if ($self->mode eq "r"){
				confess($dd);
			}
			system("mkdir -p $dd && chmod a+w $dd");
		}
		return $dd;
		
	}
);

has chunk_size =>(
	is		=> 'rw',
	lazy=>1,
	default => sub {
		return 5_000_000;
	}
);
has chromosome =>(
	is		=> 'rw',
	required=>1,
);
has chromosome_name =>(
	is		=> 'ro',
	lazy=>1,
	 default => sub {
		my $self = shift;
		my $chr = $self->chromosome;
		$chr =~ s/chr//;
		$chr="MT" if $chr eq "M";
		return $chr;
	 }
);
has mode =>(
	is		=> 'rw',
	required=>1,
);


has chunks =>(
	is		=> 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		if($self->mode eq "c" or ($self->mode eq "w" && !(-e $self->json_file))){
			$self->clean() if -e $self->json_file;
			my $length = $self->chr_length;
			die() unless $length;
			my $c = $self->divide_by_chunks($length);
			$self->write_config($c);
			return $c;
		}
		
		else {
			die() unless -e $self->json_file;
			return $self->load_config();
		}
		
}
);
sub clean {
	my ($self) =@_;
	my $tchunks = $self->load_config();
	foreach my $id (keys %{$tchunks}){
		my $dir = $self->dir_db."$id";
		#warn $dir;
		RocksDB->destroy_db($dir) if -e $dir;
	}
	foreach my $id (keys %{$tchunks}){
		system("rmdir ".$self->dir_db."$id") if -e $self->dir_db."$id";
	}
	#die();
}

has tree =>(
is		=> 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		
 		my $tree =  Set::IntervalTree->new;
 		for my $r (@{$self->regions}){
 			  $tree->insert($r->{id}, $r->{start}, $r->{end});
 		}
 		return $tree;
	}
);

sub nosql_rocks {
	my ($self,$region) = @_;
	my $id = $region->{id};
	return $self->_chunks_rocks($id);
};
sub _pile {
	my ($self,$id) = @_;
	$self->{array_rocksdb} = [] unless exists $self->{array_rocksdb};
	if (scalar(@{$self->{array_rocksdb}} > 1)){
		my $sid = shift(@{$self->{array_rocksdb}});
		$self->chunks->{$sid}->{rocksdb}->close();
		$self->chunks->{$sid}->{rocksdb} = undef;
		delete $self->chunks->{$sid}->{rocksdb};
	}
	push(@{$self->{array_rocksdb}},$id);
}
sub _chunks_rocks {
	my ($self,$id) = @_;

	 unless (exists $self->chunks->{$id}->{rocksdb}){
	 	$self->_pile($id);
	 	$self->chunks->{$id}->{rocksdb} = GenBoNoSqlRocksChunks->new(dir=>$self->dir_db,mode=>$self->mode,name=>$id,start=>$self->chunks->{$id}->{start},index=>$self->index,compress=>$self->compress,pipeline=>$self->pipeline);
	 #	system("vmtouch -q -t ".$self->chunks->{$id}->{rocksdb}->path_rocks);
	 }
	 return $self->chunks->{$id}->{rocksdb};
};


#my $iter ; 
#my $iter = $no->rocks->new_iterator->seek_to_first;
#while (my ($key, $value) = $iter->each) {
#    #printf "%s: %s\n", $key, $value;
#    warn $key." ".join(";",$no->get1($key));
#}
has regions =>(
is		=> 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $r;
		foreach my $id (sort{$self->chunks->{$a}->{start} <=> $self->chunks->{$b}->{end}} keys %{$self->chunks}){
			push(@$r,{id=>$id,start=>$self->chunks->{$id}->{start},end=>$self->chunks->{$id}->{end},tabix=>$self->chromosome.":".$self->chunks->{$id}->{start}."-".$self->chunks->{$id}->{end}});
		}
 		return $r;
	}
);

has json_file =>(
is		=> 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->dir_db."chunk.json";
	}
);

sub write_config {
	my ($self,$chunks) = @_;
	open(my $fh,">".$self->json_file) or die("can t open file");
	my $h;
	$h->{chunks} = $chunks;
	$h->{compress} = $self->compress;
	$h->{pack} = $self->pack;
	$h->{description} = $self->description;
	$h->{index} = $self->index;
	print $fh encode_json($h);
	close ($fh);
}
sub load_config {
	my ($self) = @_;
	open(my $fh ,$self->json_file) or die("can t open file");
	my $json_str = do { local $/; <$fh> };
	close($fh);
	my $h = decode_json($json_str);
	if($self->mode ne "c"){
	$self->{compress} = delete $h->{compress};
	$self->{pack} = delete $h->{pack};
	$self->{description} = delete $h->{description};
	$self->{index} = delete $h->{index};
	}
	return delete $h->{chunks};
}
sub divide_by_chunks{
	my ($self,$length) = @_;
		my $region;
		my $size = $self->chunk_size;
		
		my $from =1;
    	my $to = $length;
  
    	my $start;
    	my $end;
    	
    	while ($from < $to){
        $start = $from;
        $end = $from + $size;
        if ($end > $to){
            $end = $to;
        }
        my $id_chunk = $self->chromosome_name.".".$from.".".$end;
        $region->{$id_chunk}->{start} = $from;
        $region->{$id_chunk}->{end} = $end;
        $from = $end;
       }
       return $region;
}
sub _current_db {
	my ($self) = @_;
	#return 
	return $self->{_current_db};
}

sub change_id {
	my ($self,$key) = @_;
	my $start = $key;
	if ($key =~ /_/){
		my @t = split("_",$key);
		$start = shift(@t);
		my $chaine = sprintf("%010d", $start);
		return($start,$chaine."_".join("_",@t));
	}
	return ($start,sprintf("%010d", $key));
		
}

sub put_batch {
	my ($self,$key,$value) = @_;
	my ($start,$key1) = $self->change_id($key);
	my $results = $self->tree->fetch($start,$start+1);
	confess(@$results) if scalar(@$results) > 1;
	 $self->_chunks_rocks( $results->[0])->put($key1,$value);
}

sub compress_gnomad_position {
	my ($self,$id) = @_;
	my ($chr,$pos,$ref,$alt) = split("-",$id);
	my $l1 = length($ref);
	my $l2 = length($alt);
	my $seqid = $alt;
	#if ($l1 ==1 && $l2 ==1){
	#	$nalt = $alt;
	#}
	if ($l1 ==1 && $l2 > 1){
		$seqid = "+".substr($alt, 1);
	}
	if ($l1 >1 ){
		$ref = substr($ref, 1);
		$seqid = ($l1 -1);
	}
	return ($chr,$pos,$seqid);
}
sub get {
	my ($self,$id) = @_;
	my ($chr,$key,$seq_id) = $self->compress_gnomad_position($id);
	my ($start,$key1) =  $self->change_id($key);
	my $results = $self->tree->fetch($start,$start+1);
	confess(@$results) if scalar(@$results) > 1;
	my $region_id = $results->[0];
	my $d = $self->_chunks_rocks($region_id)->get($key);
	return $d;
	return undef unless $d;
	return undef unless exists $d->{$seq_id};
	if ($self->pack){
		my @t = unpack($self->pack,$d->{$seq_id});
		#warn Dumper(@t);
		if ($self->description){
		my $res = {};
			foreach my $k (@{$self->description}){
				$res->{$k} = shift @t;
			}
			return $res;
		}
		return \@t;
	}
	return $d->{$seq_id};
}

sub get_raw {
	my ($self,$start,$id) = @_;
	
	my $results = $self->tree->fetch($start,$start+1);
	confess(@$results) if scalar(@$results) > 1;
	my $region_id = $results->[0];
	my $d = $self->_chunks_rocks($region_id)->get_raw($id);
	return $d;
	
}
sub get_position {
	my ($self,$key) = @_;
	my ($start,$key1) =  $self->change_id($key);
	my $results = $self->tree->fetch($start,$start+1);
	confess(@$results) if scalar(@$results) > 1;
	my $id = $results->[0];
	return $self->_chunks_rocks($id)->get($key);
}


sub close {
	my ($self) = @_;
	foreach my $id (keys %{$self->chunks}){
		$self->chunks->{$id}->{rocksdb}->close() if exists $self->chunks->{$id}->{rocksdb};
		delete $self->chunks->{$id}->{rocksdb} if exists $self->chunks->{$id}->{rocksdb};;
	}
	#if (exists $self->chunks->{$id}->{rocksdb}){
	#	$self->chunks->{$id}->{rocksdb}->close;
	#	delete $self->chunks->{$id}->{rocksdb};
	#}
}

sub DESTROY {
	warn "END rocks genome";
}

1;