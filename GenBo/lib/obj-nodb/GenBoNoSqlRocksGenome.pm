package GenBoNoSqlRocksGenome;
#use lib "$Bin/";
use strict;
use warnings;
use Moo;
use Data::Dumper;
use JSON::XS;
use POSIX;
#use GenBoNoSqlRocksChunks;
use Carp;
use GenBoNoSqlRocksAnnotation;

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
		
		die(Dumper $h ) unless exists $h->{$self->genome}->{$self->chromosome_name};
		my $v = $h->{HG19}->{$self->{chromosome_name}};
		$v = $h->{HG38}->{$self->{chromosome_name}} if $h->{HG38}->{$self->{chromosome_name}} > $v;
		return  $v;
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
has factor =>(
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
#		warn $self->mode;
		unless (-e $dd){
			if ($self->mode eq "r"){
				confess($dd)."-".$self->mode;
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
		if($self->mode eq "c" or $self->mode eq "t" or ($self->mode eq "w" && !(-e $self->json_file))){
			$self->clean() if -e $self->json_file;
			my $length = $self->chr_length;
#			warn $length;
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

has tree_db =>(
is		=> 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
 		my $tree =  Set::IntervalTree->new;
 		for my $r (@{$self->regions}){
 			#$self->chunks->{$id}->{rocksdb};
 			my $id = $r->{id};
 			my $db  = GenBoNoSqlRocksAnnotation->new(dir=>$self->dir_db,mode=>$self->mode,name=>$id,pipeline=>$self->pipeline,pack=>$self->pack,description=>$self->description,factor=>$self->factor);
 			 $tree->insert($db, $r->{start}, $r->{end});
 		}
 		return $tree;
	}
);

sub nosql_rocks {
	my ($self,$region) = @_;
	my $id = $region->{id};
	return $self->_chunks_rocks($id);
};

sub nosql_rocks_tmp {
	my ($self,$region) = @_;
	my $id = $region->{id};
	return  GenBoNoSqlRocksAnnotation->new(dir=>$self->dir_db,mode=>$self->mode,name=>$id,pack=>$self->pack,description=>$self->description,factor=>$self->factor,pipeline=>1);
	
};
sub nosql_rocks_read {
	my ($self,$region) = @_;
	my $id = $region->{id};
	return  GenBoNoSqlRocksAnnotation->new(dir=>$self->dir_db,mode=>"r",name=>$id,pack=>$self->pack,description=>$self->description,factor=>$self->factor,pipeline=>1);
	
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
	 	#$self->chunks->{$id}->{rocksdb} = GenBoNoSqlRocksChunks->new(dir=>$self->dir_db,mode=>$self->mode,name=>$id,start=>$self->chunks->{$id}->{start},index=>$self->index,compress=>$self->compress,pipeline=>$self->pipeline);
		$self->chunks->{$id}->{rocksdb} = GenBoNoSqlRocksAnnotation->new(dir=>$self->dir_db,mode=>$self->mode,name=>$id,pipeline=>$self->pipeline,pack=>$self->pack,description=>$self->description,factor=>$self->factor);
	 #	system("vmtouch -q -t ".$self->chunks->{$id}->{rocksdb}->path_rocks);
	 }
	 return $self->chunks->{$id}->{rocksdb};
};


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
has vector_index_region =>(
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return [];
	}
);

sub save_vector_index_region {
	my($self,$v) = @_;
	open(my $fh ,$self->json_file) or die("can t open file");
	my $json_str = do { local $/; <$fh> };
	close($fh);
	my $h = decode_json($json_str);
	$h->{vector_index_region} = $v;
	open(my $fh2,">".$self->json_file) or die("can t open file");
	print $fh2 encode_json($h);
	close ($fh);
}
sub write_config {
	my ($self,$chunks) = @_;
	open(my $fh,">".$self->json_file) or die("can t open file ".$self->json_file);
	my $h;
	$h->{chunks} = $chunks;
	$h->{compress} = $self->compress;
	$h->{pack} = $self->pack;
	$h->{description} = $self->description;
	$h->{factor} = $self->factor;
	$h->{index} = $self->index;
	$h->{vector_index_region} = $self->vector_index_region;
	print $fh encode_json($h);
	close ($fh);
}
sub load_config {
	my ($self) = @_;
	open(my $fh ,$self->json_file) or die("can t open file ".$self->json_file);
	my $json_str = do { local $/; <$fh> };
	close($fh);
	my $h = decode_json($json_str);
	if($self->mode ne "c"){
	$self->{compress} = delete $h->{compress};
	$self->{pack} = delete $h->{pack};
	$self->{description} = delete $h->{description};
	$self->{factor} = delete $h->{factor};
	$self->{index} = delete $h->{index};
	$self->{vector_index_region} = delete $h->{vector_index_region};
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
            $end = $from + $size;
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

sub current {
	my ($self,$key) = @_;
	return $self->{current_db};
}

sub convert_rocksId_varId {
	my ($self, $rocksId) = @_;
	
}
sub decode_dejavu {
	my ($self,$value) = @_;
	return undef unless $value;
	my @tab = split("!-x-!",$value);
	#warn $value;
	my $hash;
	foreach my $z (@tab){
		next unless $z;
	#	warn $z." ***";
		my ($p,$he,$ho,@a)=  unpack("w*",$z);
	#	warn "ok";
		$hash->{$p}->{he} = $he;
		$hash->{$p}->{ho} = $ho;
		
		$hash->{$p}->{patients} = \@a;
	} 
	return $hash;
}
sub dejavu {
	my ($self,$id) = @_;
	if ($self->current) {
		my $res =  $self->{current_db}->get_raw($id);
		return $self->decode_dejavu($res);
	}
	my ($pos,$a) =split("!",$id);
	$pos *= 1;
	$self->{current_db} = $self->get_db($pos);
#	warn $self->{current_db};
#	my $iter = $self->current->rocks->new_iterator->seek_to_first;
###	warn $id;
#my $x;
#	while (my ($key, $value) = $iter->each) {
#		$x ++;
#		next unless $key =~ /0082461765/;
#    	printf "%s  \n", $key;
#    	#warn Dumper($self->decode_dejavu($value));
#    	#warn $id;
#	}
#	die($x);
#	#$self->{current_db}->rocks->compact_range();
	my $h = $self->current->get_raw($id);
	return $self->decode_dejavu($h);
}


sub stringify_pos {
	my ($self,$pos) = @_;
	return ($pos,sprintf("%010d", $pos));
}

sub dejavu_hg19_id {
	my ($self, $rocks_id) = @_;
	my $var_id_hg19 = $self->current->get_raw('#'.$rocks_id);
	return $var_id_hg19;
}

sub get_dbs_interval {
	my ($self,$start,$end) = @_;
	my $results = $self->tree_db->fetch($start,$end);
	foreach my $db (@$results){
		my ($a,$b,$c)  = split(/\./,$db->name);
		$db->{start} = $b;
	}
	return $results;
}


sub dejavu_interval {
	my ($self, $start, $end) = @_;
	
	my $pos = $self->stringify_pos($start);
	my $dbs = $self->get_dbs_interval($start,$end);
	my $h_res;
	foreach my $db (sort{$a->{start} <=> $b->{start}} @$dbs){
		warn $db->{start};
		$self->{current_db} = $db;
		my $iter = $db->rocks->new_iterator->seek($pos);
		while (my ($key, $value) = $iter->each) {
			my ($this_pos, $id) = split('!', $key);
			next if int($this_pos) < $start;
			if ($this_pos > $end) {
				last;
			}
			$h_res->{$key} = $self->decode_dejavu($value);
		}
	}
	
	return $h_res;
}

sub spliceAI_id {
	my ($self,$id) = @_;
	my ($pos,$a) =split("!",$id);
	$pos *= 1;
	$self->get_db($pos);
	return $self->get_db($pos)->spliceAI_id($id);
}
sub spliceAI {
	my ($self,$pos,$id) = @_;
	$self->get_db($pos);
	return $self->get_db($pos)->spliceAI($id);
}

sub get_db {
	my ($self,$pos) = @_;
	my $results = $self->tree_db->fetch($pos,$pos+1);
	confess(@$results) if scalar(@$results) > 1;
	return $results->[0];
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
	#warn "END rocks genome";
}

1;