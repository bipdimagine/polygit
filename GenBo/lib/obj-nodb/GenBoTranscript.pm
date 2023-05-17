package GenBoTranscript;
use strict;
use Moo;
use GenBoCoverage;
use Data::Dumper;
use Config::Std;
use Storable qw(store retrieve freeze thaw);
#use Bio::Tools::CodonTable;
use POSIX qw(ceil floor);
use List::Util qw(min sum);
use List::MoreUtils qw(any bsearch_index  );
 use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);
 use Scalar::Util qw(looks_like_number);
extends "GenBoGenomic";

has isTranscript => (
	is		=> 'ro',
	default	=> 1,
);

has biotype => (
	is		=> 'ro',
	required=> 1,
);

has remap_status => (
	is		=> 'ro',
	required=> 1,
);

has is_partial_transcript => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return if ($self->getProject->gencode_version() eq '943');
		#return 1 if $self->remap_status() and $self->remap_status() eq 'partial';
		return 1 if $self->hash_partial_infos();
		return;
	}
);

has hash_partial_infos => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		#my $no = $self->getProject->rocksPartialTranscripts();
		my $no = $self->getProject->lmdbPartialTranscripts();
		return unless $no;
		my $h = $no->get($self->id);
		return $h if $h;
		return;
	}
);


has sequence_from_hg38 => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->hash_partial_infos->{seq38};
	}
);

has is_omim => (
        is              => 'ro',
        lazy    => 1,
        default => sub {
                my $self = shift;
                my ($t_id, $chr_id) = split('_', $self->id());
                return 1 if (exists $self->getProject->buffer->getOmimTranscriptsNames->{$t_id});
                return;
        }
);

has gene_kyoto_id => (
	is		=> 'ro',
	required=> 1,
);

has pLI => (
	is		=> 'ro',
	lazy	=> 1,
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		 return $self->project->lmdbPLI->get($self->name());
		 },
);
has isMain => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		#return 1;
		return 1 if  $self->project->lmdbMainTranscripts->exists($self->id);
	}
);


has gene_id => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		my @t = keys %{$self->genes_object()};
		confess() if scalar(@t) ne 1;
		return $t[0];
		
	},
	
);


has gene_external_name => (
	is		=> 'ro',
	required=> 1,
);

has ccds_name => (
	is		=> 'ro',
	required=> 1,
);
has gene => (
	is		=> 'ro',
	required=> 1,
);

has strand => (
	is		=> 'ro',
	required=> 1,
);

has genomic_span => (
	is		=> 'ro',
	required=> 1,
);

has essential_splice_site_span => (
	is		=> 'ro',
	reader =>'getSpanEssentialSpliceSite',
	default => sub {
		return Set::IntSpan::Fast::XS->new();
	}
);

has span_transcript => (
	is		=> 'ro',
	required=> 1,
);

has span_mature => (
	is		=> 'ro',
	default => sub {
		return Set::IntSpan::Fast::XS->new();
	}
);
has splice_site_span => (
	is		=> 'ro',
	reader => 'getSpanSpliceSite',
	default => sub {
		return Set::IntSpan::Fast::XS->new();
	}
);

has intspan => (
	is		=> 'ro',
	reader	=> 'getGenomicSpan',
	lazy	=> 1,
	default => sub {
		my $self = shift;
    	return $self->genomic_span();
	},
);
has intronic_span => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $span = Set::IntSpan::Fast::XS->new($self->start."-".$self->end);
		$span =  $span->diff($self->getGenomicSpan);
		$span =  $span->diff($self->getSpanSpliceSite);
		$span =  $span->diff($self->getSpanEssentialSpliceSite);
    	return $span;
	},
);

has tag => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		return {};
	},
);



has appris_type => (
        is              => 'ro',
        lazy    => 1,
        default => sub {
                my $self= shift;
                my ($app) = grep {$_ =~ /appris/} keys %{$self->{tag}};
                return "-" unless $app;
                my ($a,$t,$n) = split("_",$app);
                $n ="" unless $n;
                 if ($t =~ /princ/){
                        $t = "P"
                 }
                 else {
                        $t ="ALT"
                 }


                return $t.$n;
        },
);
has appris_level => (
        is              => 'ro',
        lazy    => 1,
        default => sub {
                my $self= shift;
                my ($app) = grep {$_ =~ /appris/} keys %{$self->{tag}};
                return 20 if  !($app) && $self->ccds_name;
                return "100" unless $app;
                my ($a,$t,$n) = split("_",$app);
                my $sc = 0;

                $n =0 unless $n;
                 if ($t =~ /princ/){
                         $sc= $n;
                 }
                 elsif ($n =~ /long/){
                 	$sc = 10;
                 }
                 else {
                  $sc = 10+$n;
                 }

                return $sc;
        },
);
has protein => (
	is		=> 'ro',
);

sub getOrfSequence {
	my $self = shift;
	my $seq = $self->sequence();
	if ($self->is_partial_transcript) {
		$seq = $self->sequence_from_hg38();
	}
	my $cs = substr($seq,$self->orf_start-1,abs($self->orf_start-$self->orf_end)+1);
	$self->{coding_sequence} = $cs;
	return $self->{coding_sequence} if ($self->{coding_sequence});
	
	#TODO: region pseudo autosomal
	if ($self->getChromosome->isPseudoAutosomal($self->start(), $self->end())) {
		my $new_id = $self->id();
		if ($self->getChromosome->id() eq 'Y') { $new_id =~ s/_Y/_X/; }
		elsif ($self->getChromosome->id() eq 'X') { $new_id =~ s/_X/_Y/; }
		return $self->getProject->lmdbGenBo->get($new_id)->{coding_sequence};
	}
}

has coding_sequence => (
	is		=> 'rw',
	#reader =>'getOrfSequence',
);

#has orf_start => (
#	is		=> 'ro',
#	
#);

sub orf_start {
	my $self = shift;
	return $self->{orf_start_new} if exists  $self->{orf_start_new};
	 $self->{orf_start_new} = scalar(@{$self->array_5_utr})+1;
	
	return $self->{orf_start_new};
}
sub orf_end {
	my $self = shift;
	return $self->{orf_end_new} if exists  $self->{orf_end_new};
	 $self->{orf_end_new} =  $self->length - scalar(@{$self->array_3_utr});
	return $self->{orf_end_new};
}
#has orf_end => (
#	is		=> 'ro',
#);

has sequence => (
	is		=> 'ro',
);
	
has external_name => (
	is		=> 'ro',
);

has external_protein_name => (
	is		=> 'ro',
);

has span_coding => (
	is		=> 'ro',
	reader =>'getSpanCoding',
	default => sub {
		return new Set::IntSpan::Fast::XS;
	}
);

has  genomic_orf_start => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my @array = $self->getSpanCoding()->as_array();
		return $array[0];
	}
);
has  genomic_orf_end => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my @array = $self->getSpanCoding()->as_array();
		return $array[-1];
	}
);
has spanCodonStartEnd => (
	is		=> 'ro',
	reader =>'getSpanCodonStartEnd',
);

has length => (
	is		=> 'ro',
	reader	=> 'length',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
	
		my $iter = $self->getGenomicSpan->iterate_runs();
		 my $len=0;
			while (my ( $from, $to ) = $iter->()) {
				$len += abs($to-$from) +1;
				#push(@t,[$from, $to]);
			}
			
			return $len;
	}
);

has cds_length => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
	
		my $iter = $self->getSpanCoding->iterate_runs();
		 my $len=0;
			while (my ( $from, $to ) = $iter->()) {
				$len += abs($to-$from) +1;
				#push(@t,[$from, $to]);
			}
			
			return $len;
	}
);		
		



sub setGenes {
	my $self = shift;
	my $returnObj = {};
	my $geneId = $self->gene;
	my $chr_name =  $self->getChromosome()->name;
	$returnObj->{$geneId."_".$chr_name} = undef;
	return $returnObj;
}

sub setProteins {
	my $self = shift;
	my $returnObj = {};
	
	my $prot = $self->protein();
	
	$returnObj->{$prot} = undef if $prot;
	return $returnObj;
}

sub isncRNA {
	my $self = shift;
	return $self->biotype() !~ /pseudo/ && $self->biotype() =~ /RNA/; 
	return undef;	
}

sub ispseudoGene{
	my $self = shift;

	return 1 unless  $self->protein();
	return undef;
}

sub codonsConsequenceForVariations {
	my ($self,$var,$startg,$endg) = @_;
	return $self->codonsConsequenceForDuplication($var,$startg,$endg) if ($var->isLargeDuplication());
	return $self->codonsConsequenceForLargeInsertion($var,$startg,$endg) if ($var->isLargeInsertion() or $var->isMei());
	return $self->codonsConsequenceForDeletion($var,$startg,$endg) if ($var->isDeletion() or $var->isLargeDeletion());
	return $self->codonsConsequenceForMnp($var) if $var->isMnp();
	return $self->codonsConsequenceForInsertion($var) if $var->isInsertion();
	return $self->codonsConsequenceForSubstitution($var); 
}



sub codonsConsequenceForDuplication {
		my ($self,$var) = @_;
	my $span =  Set::IntSpan::Fast::XS->new(($var->start()-5)."-".($var->end()+5));
	my @tt = $self->getGenomicSpan->intersection($span)->as_array();
	my $real_start = $tt[0];
	my $real_end = $tt[-1];
	my $pos_transcript = $self->translate_position($real_start);
	my $pos_orf = ($pos_transcript - $self->orf_start()) + 1;
	my $pos_orf_end = $pos_orf + 1; 
		my $codon1 = $self->getCodon($pos_orf);
		
	my $results = {
		transcript_position => $pos_transcript,
		orf_position => $pos_orf,
		orf_end => $pos_orf_end,
		seq_orf =>$var->sequence,
		prot_position => ceil($pos_orf/3),
		codon => $codon1,
		codon_mut => "",
		aa => "DUP",#$self->getProject->biotools->translate($codon1,$self->isMT),
		aa_mut => "dup",
	};
	return $results;
}

sub codonsConsequenceForLargeInsertion {
	my ($self,$var) = @_;
	my $span =  Set::IntSpan::Fast::XS->new(($var->start()+3)."-".($var->start()+3));
	my @tt = $self->getGenomicSpan->intersection($span)->as_array();
	my $real_start = $tt[0];
	my $real_end = $tt[-1];
	my $pos_transcript = $self->translate_position($real_start);
	my $pos_orf = ($pos_transcript - $self->orf_start()) + 1;
	my $pos_orf_end = $pos_orf + 1; 
	my $codon1 = $self->getCodon($pos_orf);
		
	my $results = {
		transcript_position => $pos_transcript,
		orf_position => $pos_orf,
		orf_end => $pos_orf_end,
		seq_orf =>$var->sequence,
		prot_position => ceil($pos_orf/3),
		codon => $codon1,
		codon_mut => "?",
		aa => "?",#$self->getProject->biotools->translate($codon1,$self->isMT),
		aa_mut => "?",
	};
	return $results;
}

sub codonsConsequenceForInsertion {
	my ($self,$var) = @_;
	my $posvar = $var->start();
	my $mt;
	
	
	my $pos_transcript = $self->translate_position($posvar);
	
	confess($var->start." ".$var->end) if $pos_transcript == -1;
	my $pos_orf = ($pos_transcript - $self->orf_start()) + 1;
	my $seq_orf = $var->sequence();
	my $codon1 = $self->getCodon($pos_orf);
	my $codonter ="";
	
	 $codonter = $self->getCodon($pos_orf+3) ;
	
	my $codon2 = $codon1;
	my $splice = ($pos_orf+2) % 3;
	#substr($codon2,$splice,1,$var->sequence());
	if ($self->{strand} == -1) {
		my $rseq_orf = reverse $seq_orf;
		$seq_orf ="";
		foreach my $a (split("",$rseq_orf)){
			warn $self->sequence() if $a eq "-";
			$seq_orf .= reverse_base($a);
		} 
	}
	substr($codon2,$splice,0) = $seq_orf;
	my $toto =$codon2.$codonter;
	my $nb_codons = ceil($var->length /3 );
	my $new_length = $nb_codons*3+3;

	$codon2 = substr ($toto,0,$new_length); 
	#my $codonTable = Bio::Tools::CodonTable->new();

	my $pos_orf_end = $pos_orf + 1; 
		my $results = {
		transcript_position => $pos_transcript,
		orf_position => $pos_orf,
		orf_end => $pos_orf_end,
		seq_orf =>$seq_orf,
		prot_position => ceil($pos_orf/3),
		codon => $codon1,
		codon_mut => $codon2,
		aa => $self->getProject->biotools->translate($codon1,$self->isMT),
		aa_mut => $self->getProject->biotools->translate($codon2,$self->isMT)
		
	};
	return $results;
}

sub codonsConsequenceForComplex {
	my ($self,$var) = @_;
	my $span =  Set::IntSpan::Fast::XS->new($var->start()."-".$var->end());
	my @tt = $self->getSpanCoding->intersection($span)->as_array();
	
	my $start = $tt[0];
	my $end = $tt[-1];
	if ($self->strand == -1 ){
		($start,$end) = ($end,$start);
	}


	my $var_translate_start =  $self->translate_position($start);
#	$var_translate_start = $self->orf_start+1;

	my $var_translate_end =  $self->translate_position($end);
	#$var_translate_end = $self->orf_start + 1 ;
	my $var_orf_start = $var_translate_start - $self->orf_start;
	my $var_orf_end = ($var_translate_start+$var->alt_len) - $self->orf_start;
	my $ref_orf_end =  ($var_translate_start+$var->ref_len) - $self->orf_start;
	my $seq1 = $self->getOrfSequence();
	my $len_delete = abs($var_orf_end - $var_orf_start);

	
	substr($seq1,$var_orf_start,$len_delete) = $var->var_allele;
	my $relative_codon_start = floor($var_orf_start /3) *3;
	my $relative_codon_end = floor($var_orf_end /3 +0.5) *3;
	my $ref_relative_codon_end = floor($ref_orf_end /3 +0.5) *3;
	#warn $seq1;
	my $ref_small_seq =  substr($self->getOrfSequence(),$relative_codon_start,abs($ref_relative_codon_end-$relative_codon_start ));
	my $small_seq = substr($seq1,$relative_codon_start,abs($relative_codon_end-$relative_codon_start ));
	my $results = {
		transcript_position => $var_translate_start,
		orf_position => $var_orf_start,
		orf_end => $var_orf_end,
		seq_orf =>$ref_small_seq,
		prot_position => ceil($var_orf_start/3),
		codon => $ref_small_seq,
		codon_mut => $small_seq,
		aa => $self->getProject->biotools->translate(uc($ref_small_seq),$self->isMT),
		aa_mut =>$self->getProject->biotools->translate($small_seq,$self->isMT),
	};
die();
}



sub codonsConsequenceForDeletion {
	my ($self,$var,$startg,$endg) = @_;
	my $posvar = $startg ;
	if ($self->strand == -1 ){
		$posvar = $endg;
	}
	
	my $ldeletion = abs($startg - $endg) +1;
	
	my $pos_transcript = $self->translate_position($posvar);
	return {} if $pos_transcript < 0;
	confess($var->name) if $pos_transcript < 0;
	my $pos_orf = ($pos_transcript - $self->orf_start()) + 1;
	$pos_orf = 1 if $pos_orf < 0;
	
	my $pos_orf_end = $pos_orf + $ldeletion -1; 
	
	my $seq_orf = substr($self->getOrfSequence,$pos_orf-1,$ldeletion);
	my $nb_codons =  abs (ceil($pos_orf_end/3) - ceil($pos_orf/3)) +1 ;# - ceil($ldeletion /3 );
	
	my $codon1;
	for (my $i = 0;$i< $nb_codons;$i++){
		$codon1 .= $self->getCodon($pos_orf+($i*3));
	}  
	my $length_total = $nb_codons*3 ;
	my $codon2 = "-";
	 	$codon2 = $codon1;
 	$codon1 =~ s/($seq_orf)/lc($1)/e;

 	$codon2 =~ s/($seq_orf)//;
	my $results = {
		transcript_position => $pos_transcript,
		orf_position => $pos_orf,
		orf_end => $pos_orf_end,
		seq_orf =>$seq_orf,
		prot_position => ceil($pos_orf/3),
		codon => $codon1,
		codon_mut => $codon2,
		aa => $self->getProject->biotools->translate(uc($codon1),$self->isMT),
		aa_mut =>$self->getProject->biotools->translate($codon2,$self->isMT),
	};
	return $results;
}

sub codonsConsequenceForSubstitution {
	my ($self,$var) = @_;
	my $posvar = $var->start;
	
	my $pos_transcript = $self->translate_position($posvar);
	confess() if $pos_transcript == -1;
	
	my $pos_orf = ($pos_transcript - $self->orf_start()) + 1;
	my $codon1 = $self->getCodon($pos_orf);
	#warn $codon1." ".$pos_orf." ".$self->name." ".length($self->getOrfSequence);
	my $seq_orf = substr($self->getOrfSequence,$pos_orf-1,1);
	
	my $seq= $var->getSequence();
	
	
	$seq = reverse_base($seq) if ($self->strand == -1);
	my $codon2 = $codon1;
	my $splice = ($pos_orf+2) % 3;
	
	substr($codon2,$splice,1,$seq);
	my $results = {
		transcript_position => $pos_transcript,
		orf_position => $pos_orf,
		seq_mut =>$seq,
		seq_orf =>$seq_orf,
		prot_position => ceil($pos_orf/3),
		pos_codon=>$splice,
		codon => $codon1,
		codon_mut => $codon2,
		aa =>$self->getProject->biotools->translate($codon1,$self->isMT),
		aa_mut => $self->getProject->biotools->translate($codon2,$self->isMT),
	};
	return $results;
	
}


sub codonsConsequenceForMnp {
	my ($self,$var) = @_;
	my @seq = split("",$var->sequence);
	my $codons;
     my @pos_trans;
     my @pos_orf;
	for (my $i=0;$i<@seq;$i++) {
		my $posvar = $var->start+$i;
		my $pos_transcript = $self->translate_position($posvar);
		next if $pos_transcript == -1;
		my $pos_orf = ($pos_transcript - $self->orf_start()) + 1;
		my $num_codon = int(($pos_orf-1) / 3);
		my $splice = ($pos_orf+2) % 3;
		push(@pos_orf,$pos_orf);
		push(@pos_trans,$pos_transcript);
		unless (exists $codons->{$num_codon}){
			$codons->{$num_codon}->{sequence_ref} = $self->getCodon($pos_orf);
			my $codon2 = $codons->{$num_codon}->{sequence_ref};
			my $splice = ($pos_orf+2) % 3;
			my $seq = $seq[$i];
			$seq = reverse_base($seq) if ($self->strand == -1);
			substr($codon2,$splice,1,$seq);
			$codons->{$num_codon}->{sequence_mut} =$codon2;
			$codons->{$num_codon}->{pos_orf} =$pos_orf;
			$codons->{$num_codon}->{seq_orf} = $seq; 
		}
		else {
			my $codon2 = $codons->{$num_codon}->{sequence_mut} ;
			my $splice = ($pos_orf+2) % 3;
			my $seq = $seq[$i];
			$seq = reverse_base($seq) if ($self->strand == -1);
			substr($codon2,$splice,1,$seq);
			$codons->{$num_codon}->{seq_orf} .= $seq; 
			$codons->{$num_codon}->{sequence_mut} =$codon2;
			$codons->{$num_codon}->{pos_orf} =$pos_orf if $pos_orf < $codons->{$num_codon}->{pos_orf};
		}
		
		
	}
	my $aa ="";
	my $aa_mut ="";
	my $seqcodon="";
	my $seqcodon_mut="";
	my $seq_orf = "";
	foreach my $cp (sort{$a <=> $b} keys %$codons) {
		my $codon = $codons->{$cp};
		$seqcodon .= $codon->{sequence_ref};
		$seqcodon_mut .= $codon->{sequence_mut};
		$aa .= $self->getProject->biotools->translate($codon->{sequence_ref},$self->isMT);
		$aa_mut .= $self->getProject->biotools->translate($codon->{sequence_mut},$self->isMT);
		$seq_orf .= $codon->{seq_orf};
	}
	my $pos_orf = min @pos_orf;
	
	my $moin = min @pos_trans;
	my $results = {
		transcript_position =>  $moin,
		orf_position => $pos_orf,
		seq_mut =>$var->sequence,
		seq_orf =>$seq_orf,
		prot_position => ceil($pos_orf/3),
		pos_codon=>"-",
		codon => $seqcodon,
		codon_mut => $seqcodon_mut,
		aa =>$aa,
		aa_mut => $aa_mut,
	};

	return $results;
	
}


sub getCodon{
	my ($self,$pos_tr) = @_;
	$pos_tr = 1 if $pos_tr ==0;
	confess() unless $pos_tr;
	return "" if $pos_tr == -1;
	return "" if  $pos_tr > length($self->getOrfSequence);
	my $splice = ($pos_tr+2) % 3;

	#confess() if (length($self->getOrfSequence) < $pos_tr- $splice - 1);
	#warn length($self->getOrfSequence)." ".($pos_tr- $splice - 1) if (length($self->getOrfSequence) < abs(($pos_tr- $splice - 1)+3));
	#warn length($self->getOrfSequence)." ".($pos_tr- $splice - 1) if ( ($pos_tr- $splice - 1) <0);
	my $codon1  = substr($self->getOrfSequence,$pos_tr- $splice - 1,3);
	#confess() if (length($self->getOrfSequence) < abs(($pos_tr- $splice - 1)+3));
	return $codon1;
	
	
	
}
sub search {
  my ($compare, $elem, $arr) = @_;
  my ($start, $end) = (0, $#{$arr});
  while( $start <= $end ) {
    my $middle = int(($start + $end) / 2);
    my $cmp_result = do {
      local($a,$b) = ($elem,$arr->[$middle]);
      $compare->()
    };
    if( $cmp_result == 0 ) {
      return $middle;
    }
    elsif( $cmp_result < 0 ) {
      $end = $middle - 1;
    }
    else {
      $start = $middle + 1;
    }
  }
  return -1;
}
 
sub reverse_base {
	my $base = shift;
	return "T" if $base eq "A" ;
	return "A" if $base eq "T" ;
	return "C" if $base eq "G" ;
	return "G" if $base eq "C" ;
	return "N" if $base eq "N" ;
	confess("base ".$base);
}

#here you have an array with all genomic_position use to trabnslate genomic position to transcript position
has array_pos => (
	is      => 'ro',
	lazy    => 1,
	default => sub { 
	my $self = shift;
	my @all;

	 @all = $self->getGenomicSpan->as_array;
	
	#if ($self->strand == - 1){
	#	@all = reverse ($self->getGenomicSpan->as_array);
	#}
	return \@all;
	
	}
);




#same than array_pos but on the coding part of the transcript 

has array_coding_pos => (
is      => 'ro',
	lazy    => 1,
	default => sub { 
	my $self = shift;
	my @all;
	 @all = $self->getSpanCoding->as_array;
	
	#if ($self->strand == - 1){
	#	@all = reverse ($self->getGenomicSpan->as_array);
	#}
	return \@all;
	
	}
);



has array_5_utr => (
is      => 'rw',
	lazy    => 1,
	default => sub { 
	my $self = shift;
		my $span_genomic = $self->genomic_span();
		my $span_utr = $span_genomic->diff($self->getSpanCoding);
		if ($self->strand == 1 ){
			my $end_coding  = $self->array_coding_pos->[-1];
			$span_utr->remove_range($end_coding,$end_coding+50000000);
			
		}
		else {
		my $start_coding = $self->array_coding_pos->[0];
		$span_utr->remove_range(0,$start_coding);
		}
		return [$span_utr->as_array];
	}
);
has array_3_utr => (
is      => 'rw',
	lazy    => 1,
	default => sub { 
		my $self = shift;
		my $span_genomic = $self->genomic_span();
		my $span_utr = $span_genomic->diff($self->getSpanCoding);
		if ($self->strand == 1 ){
			my $start_coding = $self->array_coding_pos->[0];
			$span_utr->remove_range(0,$start_coding);
		}
		else {
				my $end_coding  = $self->array_coding_pos->[-1];
			$span_utr->remove_range($end_coding,$end_coding+50000000);
		}
			return [$span_utr->as_array];
	}
);
#generic method a genomic_position and return the index in the array 

sub search_position {
	my ($self,$pos,$array) = @_;
	confess() unless $pos;
	my $xx = bsearch_index{$_ <=> $pos} @{$array};
	return -1 if $xx == -1;
	return -1 if $array->[$xx] ne $pos;
	return $xx if $xx == -1;
	if ($self->strand == -1){
		$xx = $#{$array} - $xx;
	}
	return $xx+1;
	
}
has span_genomic_intervaltree => (
is      => 'rw',
	lazy    => 1,
	default => sub { 
		my $self = shift;
		my @t ;
		my $iter = $self->getGenomicSpan->iterate_runs();
		 my $tree = Set::IntervalTree->new;
		 my $index = 0;
		 my $len=0;
			while (my ( $from, $to ) = $iter->()) {
				
				$tree->insert([$from,$len],$from,$to+1);
				$len += abs($to-$from) +1;
				#push(@t,[$from, $to]);
			}
			return $tree;
	}
);


sub get_correct_translate_position_hg38 {
	my ($self, $pos) = @_;
	return 0 if not $self->is_partial_transcript();
	return 0 unless $self->hash_partial_infos;
	$self->getSpanSpliceSite->empty();
	$self->getSpanSpliceSite->add_from_string( $self->hash_partial_infos->{splice_site_span}->as_string() );
	$self->getSpanEssentialSpliceSite->empty();
	$self->getSpanEssentialSpliceSite->add_from_string( $self->hash_partial_infos->{essential_splice_site_span}->as_string() );
	$self->{orf_start_new} = $self->hash_partial_infos->{HG38}->{cds}->{start};
	$self->{orf_end_new} = $self->hash_partial_infos->{HG38}->{cds}->{end};
	foreach my $nt (sort {$b <=> $a} keys %{$self->hash_partial_infos->{intspan}}) {
		if ($self->hash_partial_infos->{intspan}->{$nt}->contains($pos)) {
			return $nt;
		}
	}
	warn "\n";	
	warn "Transcript: ".$self->id();
	warn "Checking position : ".$pos;
	warn Dumper $self->hash_partial_infos();
	confess("\n\nProblem transcript ".$self->id()." in get_correct_translate_position_hg38 method. Die\n\n");
}

#translate genomic position to transcript position call to search_position on the corresponding array
sub translate_position{
	my ($self,$pos,$debug) = @_;
	confess()  unless  $self->getGenomicSpan;	
	
	return -1 unless ($self->getGenomicSpan->contains($pos));
	my $r1 = $self->span_genomic_intervaltree->fetch($pos,$pos+1);
	my $r = $r1->[0];
	my $rpos = abs($r->[0]-$pos)+1 + $r->[1];
 	$rpos = abs($rpos - $self->length) +1 if $self->strand == -1;
 	
# 	warn $rpos;
 	$rpos += $self->get_correct_translate_position_hg38($rpos);
 	
	return $rpos;
	
#			$self->project->start_timer();
#	my $t = $self->search_position($pos,$self->array_pos);
#	$self->project->add_timer();
#	die($t." ".$rpos) if $t ne $rpos;
#	return $t;
	#return $self->search_position($pos,$self->array_pos);
}
#translate genomic position to orf position call to search_position on the corresponding array

has span_coding_intervaltree => (
is      => 'rw',
	lazy    => 1,
	default => sub { 
		my $self = shift;
		my @t ;
		my $iter = $self->getSpanCoding->iterate_runs();
		 my $tree = Set::IntervalTree->new;
		 my $index = 0;
		 my $len=0;
			while (my ( $from, $to ) = $iter->()) {
				
				$tree->insert([$from,$len],$from,$to+1);
				$len += abs($to-$from) +1;
				#push(@t,[$from, $to]);
			}
			return $tree;
	}
);



sub translate_coding_position{
	my ($self,$pos,$debug) = @_;
	return -1 unless ($self->getSpanCoding->contains($pos));
	return -1 if $self->getGenomicSpan->is_empty;
	my $r1 = $self->span_coding_intervaltree->fetch($pos,$pos+1);
	my $r = $r1->[0];
	my $rpos = abs($r->[0]-$pos)+1 + $r->[1];
	
 	$rpos = abs($rpos - $self->cds_length) +1 if $self->strand == -1;
 	return $rpos;
 	
	my $t =  $self->search_position($pos,$self->array_coding_pos);
	die($rpos." ".$t) if $rpos ne $t; 
	return $t;
}

#HGVS to genomic position



sub translate_coding_position_on_genomic {
	my ($self,$pos,$debug) = @_;
	if ($pos =~/^[-*]/){
		return $self->translate_position_on_genomic($pos);
		
	}
	my ($pos1,$pos2);
	$pos2 =0;
	if ($pos =~/\+/){
	($pos1,$pos2) = split(/\+/,$pos) ;
	$pos = $pos1;
	
	}
	if ($pos =~/\-/){
	($pos1,$pos2) = split("-",$pos) ;
	$pos = $pos1;
	$pos2 *=-1;
	
	}
	unless ($pos ){
		$pos =1;
		
	}
	$pos--;
	
	if ($self->strand == -1){
		my @t  = reverse(@{$self->array_coding_pos});
		
		 $pos2*=-1;
		 die($self->name) if $pos > scalar(@t);
		return ($t[$pos]+$pos2)
	}
	
	return ($self->array_coding_pos->[$pos]+$pos2)
}

sub translate_position_on_genomic {
	my ($self,$pos,$debug) = @_;
	#$pos.="+10";
	my ($pos1) = $pos =~/^([\-\*]\d+)/;
	
	my $pos2 = $pos;
	my $pos3 = "\\".$pos1;
	$pos2 =~s/$pos3//;
	$pos2=0 unless $pos2;
	my $array;
	if ($pos1 =~/^\*/){
		$array = $self->array_3_utr();
		$pos1 =~ s/\*//;
	}
	else {
		$array = $self->array_5_utr();
	}
	
	

	#([\-\+]\d+)
	if ($self->strand == -1){
		$pos1 *= -1;
		$pos2 *= -1;
	}
	if ($pos1>0){
		$pos1 -=1;
	}
	if (abs($pos1)>scalar(@$array)){
		confess($pos);

	}
	return ($array->[$pos1]+$pos2);
	
	
	
}




sub translate_hgvs_vcf {
	my ($self,$hgvs,$debug) = @_;
	my $type ="";
	 $type = "snp"  if $hgvs=~/\>/;
	$type .= "del" if $hgvs=~/del/;
		$type .= "dup" if $hgvs=~/dup/;
	$type .= "ins" if $hgvs=~/ins/;
	
	if ($type eq "del"){
		my ($c,$pos,$t,$all) = $hgvs =~/(.)\.(.*)(del)(.*)/;
		my ($s,$e) = split("_",$pos);
		my $a = $self->translate_coding_position_on_genomic($s);
		my $b = $a +$self->strand*(length($all)-1);
		$all = BioTools::complement_sequence($all) if ($self->strand == -1);
		($a,$b) = sort ($a,$b);
		my $seq = $self->getChromosome()->sequence($a,$b);
		my $seq_ref = $self->getChromosome()->sequence($a-1,$a-1);
		die($seq." DIE $hgvs ==> $a ".$all." $a $b".$self->name) if $seq ne $all;
		my $v;
		$v->{type} = "deletion";
		$v->{genbo}->{start} = $a;
		$v->{genbo}->{ref} = $seq;
		$v->{genbo}->{alt} = "-";
		$v->{vcf}->{start} = $a-1;
		$v->{vcf}->{ref} = $seq_ref.$seq;
		$v->{vcf}->{alt} = $seq_ref;
		
		return ($v);
	}
	elsif ($type eq "dup" or $type eq "ins"){
		my ($c,$pos,$t,$all) = $hgvs =~/(.)\.(.*)($type)(.*)/;
		my ($s,$e) = split("_",$pos);
		my $a = $self->translate_coding_position_on_genomic($s);
		$all =BioTools::complement_sequence($all) if ($self->strand == -1);
		my $seq_ref=$self->getChromosome()->sequence($a-1,$a-1);
		my $v;
		$v->{type} = "insertion";
		$v->{genbo}->{start} = $a;
		$v->{genbo}->{ref} = $self->getChromosome()->sequence($a,$a);
		$v->{genbo}->{alt} = $all;
		$v->{vcf}->{start} = $a-1;
		$v->{vcf}->{ref} = $seq_ref;
		$v->{vcf}->{alt} = $seq_ref.$all;
		return($v);
		
	}
	elsif ($type eq "snp") {
	my ($c,$pos,$all) = $hgvs =~/(.)\.(.*)([ATCG]>[ATCG])/;
	my $a = $self->translate_coding_position_on_genomic($pos);
	my ($ref,$alt) = split(">",$all);
	$alt = BioTools::complement_sequence($alt) if ($self->strand == -1);
	$ref = BioTools::complement_sequence($ref) if ($self->strand == -1);
	my $seq_ref = $self->getChromosome()->sequence($a,$a);
	die("DIE :".$hgvs."  pos :$pos genomic:$a : $seq_ref $ref ".$self->name) if $seq_ref ne $ref;
	my $v;
	$v->{type} = "snp";
	$v->{genbo}->{start} = $a;
	$v->{genbo}->{ref} = $ref;
	$v->{genbo}->{alt} = $alt;
	$v->{vcf}->{start} = $a;
	$v->{vcf}->{ref} = $ref;
	$v->{vcf}->{alt} = $alt;
	return ($v);
	}
	else {
		warn "****==>".$hgvs;
		return;
	}
	
}
#test construct an annotations tree
 
 
 sub transform_instpan_to_interval_tree {
 	my ($self,$intspan,$tree,$type) = @_;
 	my $iter = $intspan->iterate_runs();
 	 while (my ( $from, $to ) = $iter->()) {
 	 	$tree->insert($type,$from,$to+1);
 	 }
 }

 has annotations_tree => (
 is      => 'ro',
	lazy    => 1,
	default => sub { 
	my $self = shift;
	 my $tree = Set::IntervalTree->new;
	 	foreach my $a (@{$self->annotations_array}){
	 		$tree->insert(@$a);
	 	}
	 return $tree;
	
	
	 
	}
	
 );





has exons_introns_simple_tree => (
is      => 'ro',
	lazy    => 1,
	default => sub { 
			my $self = shift;
	 my $tree = Set::IntervalTree->new;
	 #my $hh;
	foreach my $ex (@{$self->getExons}){
		my $h = {start=>$ex->start,end=>$ex->end,type=>"exon",name=>$ex->name};
		$tree->insert($ex,$ex->start,$ex->end+1) ;
	}
	foreach my $ex (@{$self->getIntrons}){
		my $h = {start=>$ex->start,end=>$ex->end,type=>"intron",name=>$ex->name};
		next if $ex->start == $ex->end;
		$tree->insert($ex,$ex->start,$ex->end+1) ;
	}

    return $tree;
	}
);

#ok construct an intervaltree on each exon intron position 
# the value for each position of the tree are 
#tip : I divide the intron in two equal part to directly have the name of the nearest exon 

# type =>"NC" non coding exon,  "exon",  "intron" , "5" 5 prime of the transcript "3" 3 prime of the transcript 
# nb => the exon number (for exon it's trivial) for intron the nearest exon
#   position => "+" or "-" the nearest exon is at left or right always for nomenclature t>C"-"50 or T>C"+"50 see constructnomenclaure on variation object
# pos_ref => useful for nomenclature 
#	-> exon first base depending of the strand
# -> intron first or last base of the nearest exon 
# use pos_ref to determine the distance of your genomic position and the nearest exon 
# from start of the exon
# to end of the exon

has exons_introns_tree => (
	is      => 'ro',
	lazy    => 1,
	default => sub { 
	my $self = shift;
	 my $tree = Set::IntervalTree->new;
	 #my $hh;
	 my $span_genomic = $self->genomic_span();
	 my $span_utr = $span_genomic->diff($self->getSpanCoding);
	my $span1 =  new Set::IntSpan::Fast::XS($self->start."-".$self->end);
	my $span_intronic = $span1->diff($span_genomic);
	
	my $iter = $span_genomic->iterate_runs();
	my $exon_nb =1;
	my @apos;
	#get all exon from to for the exon position 
    while (my ( $from, $to ) = $iter->()) {
    		my $type = "";
    		#add NC if it's NC coding exon 
    		$type = "NC" if $span_utr->contains_all_range($from,$to); 
    		
    		my $h = {from=>$from,to=>$to,type=>"exon",complement=>$type};
    		#pos_ref is from or to depending of the strand 
    		$h->{pos_ref} =$from if $self->strand == 1;
    		$h->{pos_ref} = $to if $self->strand == -1;
    		push(@apos,$h);
    }
 #die($self->name);
 # now I have a array with all exon in hashtable   
    
#now the introns 
# initialise the first exon it could be the 1 or the last of the array depending of the strand
 
     my $nb_exon = 1;
     if ($self->strand == -1){
     	  $nb_exon = scalar(@apos) ;
     }
     

     
     
     #my $start = $apos[0]->{from};
     $apos[0]->{nb} = $nb_exon;
     $apos[0]->{exon_name} = "ex".$apos[0]->{nb}."".$apos[0]->{complement};
  
   
   #ok let's start with all @apos
     	
    for  (my $i=0;$i<scalar(@apos);$i++) {
		#apos is an exon so directly insert in the tree the exon
		#	push(@$hh,["exon",$apos[$i],$apos[$i]->{from},$apos[$i]->{to}+1]);
			$apos[$i]->{rpos} = $apos[$i]->{to};
    		$tree->insert($apos[$i],$apos[$i]->{from},$apos[$i]->{to}+1);
    		#and for one exson I will have 2 pseudo-introns 
    		
    		#next exon is ++ or -- depending of the strand
    		$nb_exon +=  $self->strand;
    		
    		# semi intron before exon 
    		 if ($i <  scalar(@apos)-1)  {
    		 
    		#semi intron after exon 
    			my $h = dclone($apos[$i]);
    			$h->{type} = "intron";
    			$h->{position} = "+";
    			$h->{pos_ref} = $apos[$i]->{to}; #the neraest exon position is the end position of the previous exon
    			#$h->{pos_ref} = $apos[$i]->{to} if $self->strand == -1;
    			
    			#end position of this semi intron at the middle of the intron 
    			my $mean = int((abs(($apos[$i]->{to}+1)-($apos[$i+1]->{from}-1)))/2+1);
    			$h->{rpos} = $apos[$i]->{to}+$mean+1;
    			
   			#insert semi intron in the tree
   				#push(@$hh,["1",$apos[$i]->{to}+1,$apos[$i]->{to}+$mean+1,$mean,$h]);
    			$tree->insert($h,$apos[$i]->{to}+1,$apos[$i]->{to}+$mean+1);
    			
    			#second semi intron the nereast exon is now the next exon 
    			
    			$apos[$i+1]->{nb} = $nb_exon;
    			$apos[$i+1]->{exon_name} = "ex".$apos[$i+1]->{nb}."".$apos[$i+1]->{complement};;
    			my $h2 = dclone $apos[$i+1];
    			$h2->{pos_ref} = $apos[$i+1]->{from};
    			#$h->{pos_ref} = $apos[$i+1]->{from} if $self->strand == -1;
    			$h2->{type} = "intron";
    			$h2->{position} = "-";
    			#insert it but this time the start is the end position of the previoous semi intron +1 
    			#push(@$hh,["2",$apos[$i]->{to}+$mean+1,$apos[$i+1]->{from},$mean,$h2]) ;
    			next if ($apos[$i]->{to}+$mean+1 >= $apos[$i+1]->{from});
    			$h->{rpos} = $apos[$i+1]->{from};
    			$tree->insert($h2,$apos[$i]->{to}+$mean+1,$apos[$i+1]->{from}) ;
    		}
    		
    
    		
			
    	}
    	
    	#ok this one is special for the 5' and 3' position on the transcript limit is +/- 211 bases 
  		my $hstart = dclone $apos[0];
    		$hstart->{type} = "5";
    		$hstart->{position} = "-";
    		$hstart->{pos_ref} =  $apos[0]->{from};
    		my $tpos = $hstart->{from};
    		$tpos =   $hstart->{to};
    	# push(@$hh,["3",$hstart->{from}-250,$hstart->{from}]);
    		$tree->insert($hstart,$hstart->{from}-250,$hstart->{from});
    		my $hend = dclone $apos[-1];
    		$hend->{pos_ref} =  $hend->{to};
    		$hend->{type} = "3";
    		$hend->{position} = "+";
    		#push(@$hh,["4",$hend,$hend->{to}+1,$hend->{to}+250]);
    		$tree->insert($hend,$hend->{to}+1,$hend->{to}+250);
    	#	$self->{tree_test} = $hh;
	 return $tree;
		
		 }
);

# find on which exon or intron you are 
# I use the exons_introns_tree check that before 

sub find_exon_intron {
	my ($self,$start,$end) = @_;
	confess() unless  $self->getGenomicSpan;
	#warn $self->getGenomicSpan->as_string." $pos" unless $self->getGenomicSpan->contains($pos);
	#return  unless $self->getGenomicSpan->contains($pos);
	
	my $r1 = $self->exons_introns_tree()->fetch($start,$start+1);
	 if ( scalar(@$r1) == 0 && defined $end){
	 	$r1 = $self->exons_introns_tree()->fetch($end,$end+1);
	 }
	return if  scalar(@$r1) == 0;
		if (scalar(@$r1) ne 1){
			#warn $self->genomic_span->as_string();
			confess($start." ".$end);
		}
		return $r1->[0];
}


#return the exact number of the exon at this position if you are in an exon return the {nb} if not return -1 see exons_introns_tree 
sub findExonNumber {
	my ($self,$start,$end) = @_;
	my $r = $self->find_exon_intron($start,$end);
	return -1 unless defined $r;
	if ($r->{type} eq "exon"){
		return $r->{nb};
	}
	return -1;
}



# try to find the nearest exon form one positon 
#return distance from the nearest start exon , exon_name , and reference_position the nearest exon position 
sub computeNearestExon {
		my ($self,$start,$end) = @_;
		my $r = $self->find_exon_intron($start,$end);
		unless ($r){
			$r = $self->find_exon_intron($start,$end);
			return $r;
			warn $self->name." ".$self->getGene->external_name();
			confess($start." ".$end." ".$self->start." ".$self->end);
			
		}
		
	#	warn $start;
		
		die() unless $r->{pos_ref};
		my $n = $r->{exon_name};
		my $p = ($start - $r->{pos_ref})*$self->strand;
		return ($p,$n,$r->{pos_ref});
}

# tokeep compatiblity  
#return string  with distance."_".exon_number

sub findNearestExon {
	my ($self,$start,$end) = @_;
	my ($dist_min,$nearest) = $self->computeNearestExon($start,$end);
#	if ($r->{type} eq "exon")
	 return "-" unless $nearest;
	return $dist_min."_".$nearest;
	
}

sub setIntrons {
	my ($self) = @_;
	my $span_genomic = $self->genomic_span();
	my $span1 =  new Set::IntSpan::Fast::XS($self->start."-".$self->end);
	my $span_intronic = $span1->diff($span_genomic);
	#my $span_coding = $self->span_coding();
	
	my $exons;
	
	
	my $hpos;
	my $iter     = $span_intronic->iterate_runs();
	return {} if $span_intronic->is_empty; 
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
	$num_exon = scalar(@$pos)+1 if  $self->strand == -1;
	
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
		$hpos->{ext} ="intron";
		$hpos->{start} = $from;
		$hpos->{end}   = $to;
		$hpos->{name}= $hpos->{ext}."[".$hpos->{id}."-".($num_exon)."]";#.$hpos->{ext2};
		$hpos->{id}= $self->id.$hpos->{ext}.$hpos->{id};#.$hpos->{ext2};
		$hreturn->{$hpos->{id}} = undef;
		$hpos->{length}   = ($to-$from)+1;
		$hpos->{strand}   = $self->strand();
		$hpos->{intspan} = $ps;
		$hpos->{utr} = $ps->diff($self->getSpanCoding);
		
		#$hpos->{name} .= "NC" if $hpos->{utr}->equals($hpos->{intspan});
		
		my $len = $hpos->{start}-$hpos->{end}+1;
		push( @$exons, $hpos );
	}
	if (scalar(@$exons) == 1){
		 $exons->[0]->{utr}->empty ;
		
	}
	my @temp = sort {$a->{start} <=> $b->{start}} @$exons;
	my $objs = $self->getProject()->flushObjects("introns",$exons);


	return $hreturn;	
}

sub getAllGenomicsParts {
	my $self = shift;
	my $objs =  $self->getExons();
	push(@$objs,@{$self->getIntrons});
	my @so = sort{$a->end <=> $b->end} @$objs;
	return \@so;
}

sub setExons {
	my ($self) = @_;
	my $span_genomic = $self->genomic_span();
	#my $span_coding = $self->span_coding();
	
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
		$hpos->{utr} = $ps->diff($self->getSpanCoding);
		$hpos->{name} .= "NC" if $hpos->{utr}->equals($hpos->{intspan});
		
		my $len = $hpos->{start}-$hpos->{end}+1;
		push( @$exons, $hpos );
	}
	#if (scalar(@$exons) == 1){
		
	#	 $exons->[0]->{utr}->empty ;
		
	#}
	my @temp = sort {$a->{start} <=> $b->{start}} @$exons;

	my $objs = $self->getProject()->flushObjects("exons",$exons);
	return $hreturn;
	my $h1 = $self->setIntrons();
	my %hash = (%$hreturn,%$h1);

	return \%hash;	
}

sub getGenomicPositionFromNomenclature {
	my ($self,$pos) = @_;
	my @coding_pos = $self->getSpanCoding()->as_array();
	@coding_pos = reverse(@coding_pos) if ($self->strand ==-1);
	return $coding_pos[$pos];
}

########################
# methods for coverage
########################

sub mean_exonic_coverage {
	my ($self,$patient) = @_;
	my $toto = $self->getGene()->get_coverage($patient)->coverage_intspan($self->getGenomicSpan);
	return int($toto->{mean});
}

sub mean_coding_coverage {
	my ($self,$patient) = @_;
#	warn $self->id;
	my $toto = $self->getGene()->get_coverage($patient)->coverage_intspan($self->getSpanCoding);
	return int($toto->{mean});
	
}





sub return_raw_coverage_obj{
	my ($self,$p) = @_;
	return $self->getGene()->get_coverage($p);
}


sub purge_patient{
	my ($self,$p) = @_;
	delete $self->{coverage_obj}->{$p->id} ;
}

sub getListCaptureDiag {
	my $self = shift;
	my @lTmp = split('_', $self->id());
	my $id = $lTmp[0];
	my $query = $self->getQuery();
	my $list = $query->listCaptureDiagnosticFromTranscriptId($id);
	return $list;
}

sub getIntspan {
	my ($self,$padding,$utr,$intronic) = @_;
	my $intspan  = $self->getSpanCoding;
	if ($utr == 1 or $intspan->is_empty ) {
		$intspan =  $self->genomic_span();
	}
	if ($intronic == 1){
		$intspan =  Set::IntSpan::Fast::XS->new($self->start."-".$self->end);
	}
	
	return $intspan if $padding == 0;
	my $ex_intspan = Set::IntSpan::Fast::XS->new();
	
	my $iter = $intspan->iterate_runs();
		while (my ( $from, $to ) = $iter->()) {
			$ex_intspan->add_range($from-$padding,$to+$padding);
		}
		warn$self->name unless  $ex_intspan->as_string;
	return $ex_intspan;
}
1;