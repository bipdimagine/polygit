package GenBoPcrMultiplex;

use Moose;

use Data::Dumper;
use Config::Std;
use FindBin qw($Bin);
use lib "$Bin/";
use lib "$Bin/../GenBoDB";
use Set::IntSpan::Fast::XS;
use Tabix;
use Storable qw(store retrieve freeze thaw);
use List::Util qw( shuffle sum min max);
extends "GenBoCapture";


has multiplex =>(
	is		=> 'ro',
	lazy	=> 1,
	reader	=> 'getMultiplex',
	default => sub {
		my $self = shift;
		my %t;
		foreach my $line (@{$self->primers_lines()}){
		my $hpos;
		chomp($line);
		my (@toto) = split(" ",$line);
		$t{$toto[-1]} ++; 
		}
		
		my @tt = keys %t;
		return \@tt;
	},

);

has isPcr =>(
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		return 1;
	},

);
has isCapture =>(
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		return undef;
	},

);

has files =>(
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $project = $self->getProject();
		my $version = $project->getVersion();
		$version = "HG19" if $version =~/HG19/;
		my  $dir = $self->project->capture_dir;
		my $return = {};
		my $files;
		my $file = $dir . "/" . $self->type . "/" . $self->file_name . ".gz";
		confess("Unable to find capture file :" .$file."\n") unless -e $file;
		$files->{gz} = $file;
		my $file2 = $dir . "/" . $self->type . "/" . $self->file_name;
		$files->{bed} = $file2;
		
		$files->{hotspot} = $file2;
		
		$files->{hotspot} =~ s/bed/hotspot.bed/;
		my $fp =  $dir . "/" . $self->type . "/" . $self->primers_filename;
		
		if (-s $fp){
			$files->{primers} =  $dir . "/" . $self->type . "/" .$self->primers_filename;
		}
		else {
			confess();
		}
		return $files;
	}
);




sub getPrimersByMultiplex {
	my ($self,$multi) = @_;
	return $self->{multiplexes}->{$multi}->{primers} if exists $self->{multiplexes}->{$multi}->{primers};
	$self->{multiplexes}->{$multi}->{primers}  = $self->parsePrimersForMultiplex($multi);
	return $self->{multiplexes}->{$multi}->{primers} ;
}

sub getPrimersIdByMultiplex {
	my ($self,$multi) = @_;
	return  $self->{primers_ids}->{$multi} if exists  $self->{primers_ids}->{$multi};
	my %rprimers;
	foreach my $line (@{$self->primers_lines()}){
		my (@sp) = split(" ",$line) ;
		my ($chr,$start,$end,$startr,$endr2,$m);# = split(" ",$line) ;
		$m = $sp[-1];
		$chr = $sp[0];
		$start = $sp[1];
		my $id;
		if (scalar(@sp) ==4){
			my $sid = $start-15;
			 $id = "primer".$chr."_".$sid;
			$startr = $sp[1];
			$end =  $sp[2];
			$self->{primer_size}->{$id} = abs($sp[2]-$sp[1])-1; 
		}
		else {
			 $id = "primer".$chr."_".$start;
			$self->{primer_size}->{$id} = abs($sp[3]-$sp[2])-1; 
		}
		next if $m ne $multi;
		
		$rprimers{$id} = $line; 
		#$self->{primer_size}->{$id} =1; 
		 
	}
	$self->{primers_ids}->{$multi} = \%rprimers;
	return $self->{primers_ids}->{$multi} ;
}

has primers_lines =>(
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
	
		my $multiplex_file = $self->multiplexFile();
			die("unable to find $multiplex_file, I'm dyiiiiiinnnng." ) unless -e $multiplex_file;
		my @lines = `cat $multiplex_file`;
		chomp(@lines);
		my $line = $lines[0];
		my @test = split(" ",$line);
		if (scalar(@test) < 4){
			confess("you defined Pcr amplification method but I don't have enough column in primer file : $line ;")
		}
		return \@lines;
	},
);

sub getListPrimers {
	my ($self,$multi) = @_;
	return $self->{multiplex}->{$multi}->{list} if exists  $self->{multiplex}->{$multi}->{list};
	$self->{multiplex}->{$multi}->{list} =  $self->getPrimersIdByMultiplex($multi);
	return  $self->{multiplex}->{$multi}->{list} 
}

sub set_count_multiplex{
	my ($self,$multi,$value) = @_;
	$self->{multiplex}->{$multi} = $value;
} 

sub count_multiplex {
	my ($self,$multi,$gene_id,$patient) = @_;
	my $run_id = $patient->getRun->id;
	return $self->{multiplex}->{$multi}->{$run_id}->{sum} if exists  $self->{multiplex}->{$multi}->{$run_id}->{sum};
	my $primers = $self->getListPrimers($multi);

	$self->{multiplex}->{$multi}->{sum} =0;
	my $nn =0;
	my @patients = grep {$_->getRun->id eq $run_id} @{$self->getPatients()};
	
	my $sum =0;
		foreach my $primer (keys %$primers){
	
		my @data;
		foreach my $p (@patients){
			push(@data,$p->count($primer));
		}
		 $sum += sum @data;

	}
	  $self->{multiplex}->{$multi}->{$run_id}->{sum} = $sum;
	return  $self->{multiplex}->{$multi}->{$run_id}->{sum} ;
}

sub parsePrimers{
	my ($self,$list) = @_;
	my $primers;
	foreach my $line (@$list){
		my $hpos;
		chomp($line);
		my (@toto) = split(" ",$line);
		my ($chrname,$startf,$endf,$startr,$endr,$plex) ;
		my $start_id;
		if (scalar(@toto) == 6){
		 ($chrname,$startf,$endf,$startr,$endr,$plex) = split(" ",$line);
		 $start_id = $startf;
		}
		elsif (scalar(@toto) == 4){ 
			 $chrname = $toto[0];
			 $startf =  $toto[1]-15;
			  $endf=  $toto[1];
			  $startr = $toto[2];
			  $endr = $toto[2]+15;
			  $plex = $toto[3];
			   $start_id = $startf;
			
		}
		elsif (scalar(@toto) == 3){ 
			 $chrname = $toto[0];
			 $startf =  $toto[1];
			  $endf=  $toto[1];
			  $startr = $toto[2];
			  $endr = $toto[2];
			  $plex = 1;
			   $start_id = $startf;
			
		}
		else{
			confess("you defined Pcr amplification method but I don't have enought column in primer file : $line ;")
		}
		my $chromosome = $self->getProject->getChromosome($chrname);
		$hpos->{chromosomes_object}->{$chromosome->id} =undef;
		$hpos->{gstart} = $startf;
		$hpos->{gend} = $endr;
		$hpos->{start} = $endf+1;
		$hpos->{end}   = $startr-1;
		$hpos->{start_forward}   = $startf;
		$hpos->{end_forward}   = $endf;
		$hpos->{start_reverse}   = $startr;
		$hpos->{end_reverse}   = $endr;
		$hpos->{intspan_pcr} =  Set::IntSpan::Fast::XS->new( );
		$hpos->{intspan_pcr}->add_range($startf,$endf,$startr,$endr);
		$hpos->{id}="primer".$chrname."_$start_id";
		$hpos->{name}= $plex."_".$chrname."_$startf";
		$hpos->{multiplex}= $plex;
		$hpos->{length}   = abs($endf-$startr)+1;
		$hpos->{intspan} =  Set::IntSpan::Fast::XS->new($hpos->{start}."-".$hpos->{end});
		$hpos->{cnv} ={};
		$self->{primer_size}->{$hpos->{id}} = abs($endf-$startr)-1; 
		#$self->{primer_size}->{$hpos->{id}} = 1; 
		push(@$primers,$hpos);
		
	
	}
	return $primers;
} 

sub parsePrimersForMultiplex {
	my ($self,$multiplex) = @_;
	my $data;
	foreach my $line (@{$self->primers_lines()}){
	
		my (@toto) = split(" ",$line);
		next if $toto[-1] ne $multiplex;
		push(@$data,$line);
	}
	
	my $primers = $self->parsePrimers($data);
	my $objs = $self->getProject()->flushObjects("primers",$primers);
	foreach my $o (@$objs){
	
		$o->{$self->type_object}->{$self->id} = undef;
	}
	return $objs;
}


sub parsePrimersForChromosome {
	my ($self,$chromosome) = @_;
		my $no = $self->project->noSqlCoverage();
	if ($no->exists_db("primers",$chromosome->name) && $no->exists("primers",$chromosome->name)){
		return $self->restore_save_primers($no,$chromosome);
	}
	my $n = $chromosome->name;
	my $u = $chromosome->ucsc_name;
	my @lines = grep{$_=~/^$n|^$u\s/} @{$self->primers_lines()};
	my $primers = $self->parsePrimers(\@lines);
	my $objs = $self->getProject()->flushObjects("primers",$primers);
	foreach my $o (@$objs){
		$o->{$self->type_object}->{$self->id} = undef;
	}
	return $objs;
}
1;