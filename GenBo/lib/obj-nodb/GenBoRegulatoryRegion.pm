package GenBoRegulatoryRegion;
use strict;
use Vcf;
use Moo;
use Data::Dumper;
use Config::Std;
use List::Util qw( max );
extends "GenBoGenomic";

has isRegulatoryRegion => (
	is		=> 'ro',
	default	=> 1,
);

has genomic_span => (
	is		=> 'rw',
	reader	=> 'getGenomicSpan',
	lazy=>1,
	default => sub {
		my $self = shift;
		my $set = Set::IntSpan::Fast->new($self->start."-".$self->end);
		return $set;
		
	},
	
);
has target_genes => (
	is		=> 'ro',
);

has type => (
	is		=> 'ro',
);
sub setGenes {
	my $self = shift;
	return {};
	return [$self];
}

sub setTranscripts {
	my $self = shift;
	return {};
	#return [$self];
}
sub setProteins {
	my $self = shift;
	my $returnObj = {};
	#my $prot = $self->protein();
	#$returnObj->{$prot} = undef if $prot;
	#return $returnObj;
}

sub setExons {
	my ($self) = @_;
	my $span_genomic = $self->genomic_span();
	my $exons;
	my $hpos;
	my $iter     = $span_genomic->iterate_runs();
	my $hreturn;
	my $pos= [];
	while ( my ( $from, $to ) = $iter->() ) {
		my $hpos;
		$hpos->{start} = $from;
		$hpos->{end} = $to;
	
		push(@$pos,$hpos);
	}

	my $num_exon = 1;
	my $start_cds =1;
	$num_exon = scalar(@$pos) if  $self->strand == -1;

	foreach my $hp (@$pos){
		my $from = $hp->{start};
		my $to = $hp->{end};
		
		my $hpos;
		my $ps = new Set::IntSpan::Fast::XS($from."-".$to);
		$hpos->{chromosomes_object}->{$self->getChromosome->id} =undef;
		$hpos->{transcripts_object}->{$self->id}= undef;
		$hpos->{gstart} = $from;
		$hpos->{gend} = $to;
		$hpos->{id}= $num_exon;
		$num_exon+= $self->strand;
		$hpos->{ext} ="ex";
		$hpos->{start} = $from;
		$hpos->{end}   = $to;
		$hpos->{name}= $hpos->{ext}.$hpos->{id};#.$hpos->{ext2};
		$hpos->{id}= $self->id.$hpos->{ext}.$hpos->{id};#.$hpos->{ext2};
		$hreturn->{$hpos->{id}} = undef;
		$hpos->{length}   = ($to-$from)+1;
		$hpos->{strand}   = $self->strand();
		$hpos->{intspan} = $ps;
		$hpos->{utr} = 0;
		$hpos->{name} .= "NC" ;#if $hpos->{utr}->equals($hpos->{intspan});
		
		my $len = $hpos->{start}-$hpos->{end}+1;
		push( @$exons, $hpos );
	}
	#if (scalar(@$exons) == 1){
		
	#	 $exons->[0]->{utr}->empty ;
		
	#}
	my @temp = sort {$a->{start} <=> $b->{start}} @$exons;

	my $objs = $self->getProject()->flushObjects("exons",$exons);
	return $hreturn;

}



1;