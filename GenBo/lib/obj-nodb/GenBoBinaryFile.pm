package GenBoBinaryFile;
use strict;
use warnings;
use Moose;
use MooseX::Method::Signatures;
use Data::Dumper;
use Carp;
use IO::File;
use List::Util qw(max sum min);
use Storable qw(store retrieve freeze thaw);
use Carp qw(cluck longmess shortmess);


has dir => (
	is		=> 'ro',
	required=> 1,
);


has name => (
	is		=> 'ro',
	required=> 1,
);

has filename =>(
is		=> 'rw',
lazy=>1,
default => sub {
		my $self = shift;
		confess() unless $self->dir();
		return $self->dir."/".$self->name;
		
	
	},
);

has treeindexname =>(
is		=> 'rw',
lazy=>1,
default => sub {
		my $self = shift;
		return $self->dir."/".$self->name.".tree.idx";
		
	
	},
);

	
has no =>(
is		=> 'ro',
lazy=>1,
default => sub {
		my $self = shift;
		if ($self->mode eq "r"){
			confess("++".$self->filename) unless -e $self->filename;
		}
		
		my $no2 = GenBoNoSqlLmdb->new(dir=>$self->dir,mode=>$self->mode,name=>$self->name,is_compress=>1);
	
		return $no2;
	},
);


sub tree {
	my ($self,$chr) = @_;
	return $self->{tree}->{$chr} if exists $self->{tree}->{$chr} ;
	$self->{tree}->{$chr} = Set::IntervalTree->new();

	foreach my $a (@{$self->array_index($chr)}){
		$self->{tree}->{$chr}->insert(@$a);
	}
	
	return $self->{tree}->{$chr};

}

sub array_index {
	my ($self,$chr) = @_;
	return $self->{tree_array}->{$chr} if exists $self->{tree_array}->{$chr};
	$self->{tree_array}->{$chr} = $self->no->get("index_".$chr);

	unless ($self->{tree_array}->{$chr}) {
		warn "\n\n###\nPB with file ".$self->dir().'/'.$self->name()."\nGenBoBinaryFile::array_index->get('index_'.$chr) == UNDEF...\n###\n\n";
	}
	
	return $self->{tree_array}->{$chr};
}




has mode => (
	is		=> 'ro',
	required=> 1,
	default => sub {
		return "r";
	}
);

has array_cached => (
	is		=> 'ro',
	required=> 1,
	default => sub {
		return [];
	}
);

sub get_cached {
	my ($self,$id) = @_;
	#return $self->no->get($id);
	return $self->{cached}->{$id} if exists $self->{cached}->{$id};
	$self->{cached}->{$id} = $self->no->get($id);
	#$self->{cached}->{$id} = $self->no->get($id);
	
	if (scalar(@{$self->array_cached}) < 10){
		push(@{$self->array_cached},$id);
	} 
	else {
		my $zid = shift(@{$self->array_cached});
		delete 	$self->{cached}->{$zid};
	}
	#$self->{cached}->{$id} = $self->no->get($id);
	return $self->{cached}->{$id};
}

sub getDepthForIntspan {
	my($self,$chr,$sintspan) = @_;
	my ($start,$end) = $self->start_end($sintspan);
	return $self->_getDepth($chr,$start,$end,$sintspan);
	
	
	
}

sub start_end {
	my ($self,$intspan) = @_;
	my (@cpl) = split(",",$intspan->as_string);
	my ($start,$st1) = split("-",$cpl[0]);
	 my ($e1,$end) = split("-",$cpl[-1]);
	 return ($start,$end) if $end;
	 return ($start,$e1);
	
}

sub _getDepth {
		my($self,$chr,$start,$end,$sintspan) = @_;
		my $debug;
		$debug =1 if $start == 14849758;
		$end+=1 if $start == $end;
		my $z = $self->tree($chr)->fetch($start,$end);
		
		warn Dumper $z if $debug;
#		die();
		my @final;
		foreach my $sid (sort{$a->{start} <=> $b->{start}} @$z){
			warn "-------" if $debug;
			warn Dumper $sid if $debug;
			my $intspan = $sid->{intspan}->intersection($sintspan);
			my $iter = $intspan->iterate_runs();
			warn $intspan->as_string if $debug;
			my $c = $self->get_cached($sid->{name});
			while (my ( $a, $b ) = $iter->()) {
				
				$b= $a unless $b;
				my $st = $a - $sid->{start};
				warn "\t $a $b  $st".$sid->{type} if $debug; 
				my $l = $b-$a +1;
				warn "\t $a $b  $st $l ".$sid->{type} if $debug;  
				if ($sid->{type} eq "0"){
					push(@final,(0)x$l);	
		 			next;
				}		
			
				 my $multi = $self->getMulti($sid->{type});
				
				my $forward="";
				if ($st > 0){
					$st *= $multi;
					$forward = "x".$st." ";
				}
			#
		 	#die() if $debug;
		 	warn "len :".length($c)."  " if $debug;
		 	
		 	if ($st+($multi * $l) > length($c)){
		 		#warn $chr." ".$start." ".$end." st: $st multi".$forward." ".$sid->{type}." $l  ==> "." length ".length($c)." ".$sid->{name};# if $start eq 14849758;
		 		my $ll = ($st /4)+$l;
 		 		push(@final,unpack($forward.$sid->{type}."$l",$c)) ;
 		 	#	warn length($c) ;
 		 	#	die();
		 		#push(@final,(0)x$st);	
		 	#	cluck $chr." ".$start." ".$end." ".$forward." ".$sid->{type}." $l  ==> "." ".length($c) ;
		 	#		die();
		 		#push(@final,unpack($forward.$sid->{type}."$l",$c)) ;
		 		#	#if $start eq 14849758;
	
		 		#die();
		 		next;
		 	}
		 	push(@final,unpack($forward.$sid->{type}."$l",$c)) ;
			}

	}
	return \@final;
}

sub _unpacksum {
	my ($self,$a,$b,$c,$l) = @_;
	
	#if ($b == "S"){
		my $nb=0;
		my $sum = 0;
		foreach my $x (unpack($a.$b."$l",$c)){
			if ($x == 65535 && $b eq "S") {
				$x = -1;
			}
			if ($x < 0){
				next;
			}
			$sum += $x;
			$nb ++;
		}
		return($sum,$nb);
	#} 
}

sub _getSum {
		my($self,$chr,$start,$end,$sintspan) = @_;
		confess() if $end <= $start;
		my $z = $self->tree($chr)->fetch($start,$end);
		my $sum =0;
		my $nb;
		my @final;
		my $toto;
		foreach my $sid (sort{$a->{start} <=> $b->{start}} @$z){
			my $intspan = $sid->{intspan}->intersection($sintspan);
			my $iter = $intspan->iterate_runs();
			my $c =  $self->get_cached($sid->{name});
			while (my ( $a, $b ) = $iter->()) {
				$b= $a unless $b;
				my $st = $a -$sid->{start};
				my $l = $b-$a +1;
				if ($sid->{type} eq "0"){
					push(@final,(0)x$l);	
		 			next;
				}		
			
			 my $multi = $self->getMulti($sid->{type});
				my $forward="";
				if ($st > 0){
					$st *= $multi;
					$forward = "x".$st." ";
				}
			my @t;
			
			eval {
			my ($s1,$n1) = $self->_unpacksum($forward,$sid->{type},$c,$l);
			$sum += $s1;
			$nb += $n1;
			# @t =  unpack($forward.$sid->{type}."$l",$c);
			 
			};
			if ($@) {
			#			die();
			if ($st>length($c)){
				$st = length($c) - $l;
				$forward = "x".$st." ";
			}
		
			my ($s1,$n1) = $self->_unpacksum($forward,$sid->{type},$c,$l);
			$sum += $s1;
			$nb += $n1;
			#@t = unpack($forward.$sid->{type}."$l",$c);
			} 
		#	@t = @t;
	#	my @x = @t;
#		foreach my $p (@t){
#			if ($p == 65535){
#				warn $sid->{type};
#				die();
#			}
#			if ($p>0 && $p < 65535){
#				push (@x,$p);
#			}
#			else {
#				push (@x,0);
#			}
#		}
			#my @x = map {$_ = 0 if $_ < 0} @t;
#			if (@x){
#				$sum += sum(@x);
#				$nb += scalar(@x);
#			}
		 	#push(@final,unpack($forward.$sid->{type}."$l",$c));
			}

	}
	#die();
	return ($sum,$nb);
}
sub getDepth {
	my($self,$chr,$start,$end) = @_;
	die() unless $chr;
	die('ici') unless $start;
	die() unless $end;
	my $sintspan = Set::IntSpan::Fast::XS->new("$start-$end" );
	return $self->_getDepth($chr,$start,$end,$sintspan);
#	die() unless exists $self->index->{$name};
}
sub getMean {
	my($self,$chr,$start,$end) = @_;
	die() unless $chr;
	die('ici') unless $start;
	die() unless $end;
	my $sintspan = Set::IntSpan::Fast::XS->new("$start-$end" );
	my ($sum,$nb) = $self->_getSum($chr,$start,$end,$sintspan);
	return $sum/(abs($start-$end)+1);
#	die() unless exists $self->index->{$name};
}
sub getSum {
	my($self,$chr,$start,$end) = @_;
	confess() unless $chr;
	die('ici') unless $start;
	die() unless $end;
	my $sintspan = Set::IntSpan::Fast::XS->new("$start-$end" );
	my ($sum,$nb) = $self->_getSum($chr,$start,$end,$sintspan);
	return $sum;
	warn $sum.' '.abs($start-$end);
	return $sum/abs($start-$end);
#	die() unless exists $self->index->{$name};
}
sub no_index {
	my ($self) = @_;
	$self->{index} = 0 unless exists $self->{index};
	$self->{index} ++;
	return $self->{index};
}
sub putValue {
	my($self,$key,$value) = @_;
	$self->no->put($key,$value);
}
sub getValue {
	my($self,$key,$value) = @_;
	return $self->no->get($key);
}

has bytes =>(
is		=> 'rw',
lazy=>1,
default => sub {
		my $self = shift;
		my $h;
		$h->{c} = 1;
		$h->{C} =1;
		$h->{s} = 2;
		$h->{S} = 2;
		$h->{l} = 4;
		$h->{L} = 4;
		$h->{q} = 8;
		$h->{Q} = 8;
		$h->{f} = 4;
		return $h;
	},
);

sub getMulti {
	my($self,$type) = @_;
	return $self->bytes->{$type} if exists $self->bytes->{$type};
	confess("error with $type");
}
sub putDepth {
	my($self,$chr,$start,$end,$data) = @_;
	die() unless $chr;
	die() unless $start;
	die() unless $end;
	my $vlength = abs($end -$start) +1;
	#warn $self->no_index;
	my $id = $chr."_".$self->no_index;
	my $sum  = sum(@$data) ;
	$sum = 0 unless $sum;
	my $max = max(@$data);
	my $min = min(@$data);
	if ($sum == 0){
		$self->no->put($id,"0");
		push(@{$self->{tree_array}->{$chr}},[{chr=>$chr,name=>$id,start=>$start,end=>$end,intspan=>Set::IntSpan::Fast::XS->new("$start-$end" ),type=>0},$start,$end+1]);
		return ;
	}
	
	
	my $type = "S";

	if ( abs($max - int($max)) > 0 ||  abs($min - int($min)) > 0 ){
	 	$type = "f";
	}
	elsif ($min < 0)  {
		if ($max < 127 && $min > -127){
			$type = "c";
		}
		elsif($max< 32764 && $min > -32764){
			$type = "s";
		}
		elsif($max< 2147483647 && $min > -2147483647){
			$type = "l";
		}
		else {
			$type = "q";
		}
		
	}
	else {
		
		if ($max < 255 ){
			$type = "C";
		}
		elsif($max < 65535){
			$type = "S";
		}
		elsif($max< 4294967295 ){
			$type = "l";
		}
		else {
			$type = "q";
		}
		
	}
	
	my $multi = $self->getMulti($type);
	#warn $type." ".$max." ".scalar(@$data); 
	
	
	#warn min(@$data);
	#die();
	my $x = pack($type."*",@$data);
	confess($vlength." ".length($x)) if length($x) != ($multi*scalar(@$data)) or $vlength != (length($x)/$multi);
	$self->no->put($id,$x);

	push(@{$self->{tree_array}->{$chr}},[{chr=>$chr,name=>$id,start=>$start,end=>$end,intspan=>Set::IntSpan::Fast::XS->new("$start-$end" ),type=>$type},$start,$end+1]);
}

sub save_index {
	my ($self) = @_;
	foreach my $chr (keys %{$self->{tree_array} }){
			$self->no->put("index_$chr",$self->{tree_array}->{$chr});
	}
}
sub length {
	my($self,$data) = @_;
	my $fh = $self->fh;
	$fh->seek(0,SEEK_END);
	return $fh->tell;
}

sub close  {
	my($self,$data) = @_;
	#warn "close";
	
	 $self->no->close() if  exists $self->{no};
	 delete $self->{no};
}
sub unlink {
	my($self) = @_;
	unlink $self->filename;
	unlink $self->treeindexname;
}
1;