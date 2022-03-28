=head1 NAME

GBuffer : A buffer for the GenBo API

=head1 SYNOPSIS


=head1 DESCRIPTION

GBuffer provides a set of functions to create a buffer and get informations from the buffer 

=head1 METHODS

=cut
package BioTools;
use Carp;
use Moose;
use strict;

my $hcomplement = {
A=>'T',
C=>'G',
G=>'C',
T=>'A',
'-'=>'-',
'N'=>'N',
'B'=>'V',
'D'=>'H',
'H'=>'D',
'M'=>'K',
'N'=>'N',
'R'=>'Y',
'S'=>'S',
'U'=>'A',
'V'=>'B',
'W'=>'W',
'X'=>'X',
'Y'=>'R',
};

my  %IUB = ( A => [qw(A)],
	     C => [qw(C)],
	     G => [qw(G)],
	     T => [qw(T)],
	     U => [qw(U)],
	     M => [qw(A C)],
	     R => [qw(A G)],
	     W => [qw(A T)],
	     S => [qw(C G)],
	     Y => [qw(C T)],
	     K => [qw(G T)],
	     V => [qw(A C G)],
	     H => [qw(A C T)],
	     D => [qw(A G T)],
	     B => [qw(C G T)],
	     X => [qw(G A T C)],
	     N => [qw(G A T C)]
	     );
has iupac =>(
	is		=> 'ro',
	isa		=> 'HashRef',
	lazy	=> 1,	     
	default => sub { 
		return {
				A => qw(A),
				C => qw(C),
				G => qw(G),
				T => qw(T),
				U => qw(U),
				M => qw(AC),
				R => qw(AG),
				W => qw(AT),
				S => qw(CG),
				Y => qw(CT),
				K => qw(GT),
				V => qw(ACG),
				H => qw(ACT),
				D => qw(AGT),
				B => qw(CGT),
				X => qw(ACGT),
				N => qw(ACGT)	
		};
	}
	);
has rev_iupac =>(
	is		=> 'ro',
	isa		=> 'HashRef',
	lazy	=> 1,
	default => sub {
		return {
			A => "A",
			C => "C",
			G => "G",
			T => "T",
			U => "U",
			AC => "M",
			AG => "R",
			AT=>"W",
			CG=>"S",
			CT => "Y",
			GT=>"K",
			ACG=>"V",
			ACT => "H",
			AGT => "D",	
			CGT => "B",
			ACGT => "N",
	};
 
	},
);

has codons_table =>(
	is		=> 'ro',
	isa		=> 'HashRef',
	lazy	=> 1,
	default => sub {

	my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');
 	return \%g;
	},
);
has codons_table_mito =>(
	is		=> 'ro',
	isa		=> 'HashRef',
	lazy	=> 1,
	default => sub {

	my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F',,'TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'W','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'M','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'*','AGG'=>'*','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');
 
 	return \%g;
	},
);   
sub getAA {
	my ($self,$codon,$mt) =@_;
	unless ($mt){
	return "X" unless exists $self->codons_table->{$codon};
	return $self->codons_table->{$codon};
	}
	else {
	#die();
	return "X" unless exists $self->codons_table_mito->{$codon};
	return $self->codons_table_mito->{$codon};
	}
}

sub complement {
	my ($base) = @_;
	return $hcomplement->{uc($base)} if exists $hcomplement->{uc($base)};
	return "N";
	confess($base);
}

sub complement_sequence{
	my ($seq) = @_;
	my $seqr = scalar reverse($seq);
	my (@tt) = split ("",$seqr);
	my $res;
	foreach my $aa (@tt){
		my $z = complement($aa);
		unless ($z ){
			warn $seq;
			warn join('',@tt);
			warn $aa;
			confess();
		}
		$res .= $z;
	}
	return $res;
}
sub getIUPAC {
	my ($self,$chars) = @_;
	my %all;
	foreach my $c (@$chars){
		$all{$c}++;
	}
	my @bases = keys %all;
	
	return  $bases[0] if scalar(@bases) == 1;
	my $seq = uc(join("",sort{$a cmp $b} @bases));
	confess($seq) unless exists $self->rev_iupac()->{$seq};
	my $iupac = $self->rev_iupac()->{$seq};
	
}
sub translate {
	my ($self,$seq,$mt) = @_;
	unless (length($seq)%3 ==0) {
		for (my $i =0;$i<length($seq)%3;$i++){
			$seq.="N";
		}	
	}
	my @array = split("",$seq);
	my $protseq;
	for (my $i=0;$i<@array;$i+=3){
		my $codon = $array[$i].$array[$i+1].$array[$i+2];
		$protseq .= $self->getAA($codon,$mt);
	} 
	return $protseq;
}

1;