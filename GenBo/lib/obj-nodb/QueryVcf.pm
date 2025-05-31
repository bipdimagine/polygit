package QueryVcf;
use strict;
use Moo;
use Data::Dumper;
use Carp;
use Clone 'clone';
use IPC::Open2;
use List::MoreUtils qw( uniq );
use Compress::Snappy;
use Storable qw/thaw freeze/;
use Scalar::Util qw(looks_like_number);
use Align::NW;
use Bio::DB::HTS::Tabix;
use Bio::DB::HTS::VCF;
use Scalar::Util qw(weaken);
use Time::HiRes qw ( time alarm sleep );
my %REV_IUB = (
	A    => 'A',
	T    => 'T',
	C    => 'C',
	G    => 'G',
	AC   => 'M',
	AG   => 'R',
	AT   => 'W',
	CG   => 'S',
	CT   => 'Y',
	'GT' => 'K',
	ACG  => 'V',
	ACT  => 'H',
	AGT  => 'D',
	CGT  => 'B',
	ACGT => 'N',
	N    => 'N'
);

has patient => (
	is     => 'ro',
	reader => 'getPatient',
	#required=> 1,
);

has method => (
	is       => 'ro',
	required => 1,
	required=> 1,
);

has file => (
	is       => 'rw',
	required => 1,

);

has parse_hgmd => (
	is      => 'rw',
	lazy    => 1,
	default => undef,
);

has vcf => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $vcf  = Vcf->new(
			file  => $self->file(),
			tabix => $self->buffer->config->software->tabix
		);
		$vcf->parse_header();
		return $vcf;
	}
);

has tabix => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self  = shift;
		my $tabix = Bio::DB::HTS::Tabix->new( filename => $self->file() );
		return $tabix;
	}
);

has buffer => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->getPatient->buffer();
	}
);

has project => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $project;
		if ( $self->noPatient() ) {
			$project = $self->buffer->newProject( -name => 'NGS2015_0794' );
		}
		else { return $self->getPatient->project(); }
	}
);

sub getThisChromosome {
	my ( $self, $chr_name ) = @_;
	return $self->project->getChromosome($chr_name);
}

has chromosomes => (
	is      => 'ro',
	lazy    => 1,
	default => sub {

		my ($self) = @_;
		my %chromosomes;

		foreach my $chr_id ( @{ $self->tabix->seqnames() } ) {
			my $chr_tmp = $chr_id;
			$chr_tmp =~ s/chr//;
			if (   $chr_tmp eq 'M'
				|| $chr_tmp eq 'MT'
				|| $chr_tmp eq 'X'
				|| $chr_tmp eq 'Y' )
			{
				$chromosomes{$chr_id} = undef;
				next;
			}

			next unless looks_like_number($chr_tmp);
			$chromosomes{$chr_id} = undef;
		}

		#	map{$chromosomes{$_} = undef } $self->tabix->getnames();
		return \%chromosomes;
	}
);

has isUcsc => (
	is      => 'ro',
	lazy    => 1,
	default => sub {

		my ($self) = @_;
		my @chrs = sort keys %{ $self->chromosomes() };

		#	warn "==> no variation or indel in file :".$self->file() unless @chrs;
		return 1 unless (@chrs);
		return $chrs[0] =~ /chr/;
	}
);

has noPatient => (
	is      => 'rw',
	default => undef,
);

has force_all_gt_he => (
	is      => 'rw',
	default => undef,
);

has ucsc_format => (
	is      => 'rw',
	lazy    => 1,
	default => sub {

		my ($self) = @_;
		my $chrs = $self->tabix->seqnames();
		return 1 unless @$chrs;
		return $self->getPatient->project->getChromosome( $chrs->[0] )
		  ->ucsc_name eq $chrs->[0];
	},

);

sub return_tabix_chromosome_name {
	my ( $self, $chr ) = @_;
	if ( $self->ucsc_format ) {
		return $chr->ucsc_name;
	}
	return $chr->name;
}

##### METHODS #####

sub parseVcfFile {
	die("\n\nERROR: no reason to use QueryVcf::parseVcfFile method !!!\n\n");
}

sub parseVcfFileForReference {
	my ( $self, $reference, $useFilter ) = @_;
	my $file = $self->file();

	#warn $file;
	#	warn $file if $file =~/mpileup/;
	die($file) unless -e $file;
	if ( $self->method() eq "melt" ) {
		return $self->parseVcfFileForReference_melt( $reference, $useFilter );
	}
	if ( $self->method() eq "manta" or $self->method() eq "dragen-sv" or  $self->method() eq "Sniffles2") {
		return $self->parseVcfFileForReference_manta( $reference, $useFilter );
	}
	
	if ( $file =~ /\.sam/ ) {
		my @res = `zgrep "#INFO=<ID" $file`;

		#unless (scalar @res){
		return $self->parseVcfFileForReference_gatk( $reference, $useFilter );

		#}
		#	return $self->parseVcfFileForReference_gatk($reference, $useFilter);

	}
	elsif ( $file =~ /.bcf.vcf/ ) {

		return $self->parseVcfFileForReference_samtools( $reference,
			$useFilter );
	}
	elsif ( $file =~ /.vcf/ ) {
		return $self->parseVcfFileForReference_gatk( $reference, $useFilter );
	}
	elsif ( $file =~ /.bcf/ ) {
		confess( $self->getPatient()->getProject->name() );
		return $self->parseVcfFileForReference_gatk( $reference, $useFilter );
	}
	elsif ( $file =~ /.gatk/ ) {
		return $self->parseVcfFileForReference_gatk( $reference, $useFilter );
	}
	elsif ( $file =~ /.casava/ ) {
		return $self->parseCasavaFile( $reference, $useFilter );
	}
	elsif ( $file =~ /.gff3/ ) {
		warn "WARNING: GFF3 file not supported yet (" . $file . ")";
		return {};
	}
	elsif ( $file =~ /.bed/ ) {
		return $self->parseBed( $reference, $useFilter );
	}
	elsif ( $file =~ /.lid/ ) {
		return $self->parseLid( $reference, $useFilter );
	}
	else { confess("\n\nERROR: no parser for this file ($file)\n\n") }
}

sub return_tabix_query {
	my ( $self, $chromosome, $start, $end ) = @_;
	my $res = $self->tabix->query_full( $chromosome->name, $start, $end );
	unless ($res) {
		$res = $self->tabix->query_full( $chromosome->ucsc, $start, $end );
	}

}

sub parseLid {
	my ( $self, $reference ) = @_;
	confess();
	#return {};

#	my $project    = $self->getPatient()->getProject();
#	my $file       = $self->file();
#	my $chromosome = $reference->getChromosome();
#	my $res        = $self->tabix->query_full(
#		$self->return_tabix_chromosome_name($chromosome),
#		$reference->start, $reference->end );
#
##my $res =  $self->tabix->query($chromosome->name.":".$reference->start."-".$reference->end);
##	my $res =  $self->tabix->query($chromosome->fasta_name);
#	return {} unless $res;
#
#	my @data;
#	my %hashRes;
#	my $patient_id = $self->getPatient()->id;
#	while ( my $line = $res->next ) {
#
##22	38540870	38540989	1/1	dup	ENSG00000184381-PLA2G6	2.66	1.1;0.8;0.9;1.2;0.7;0.8;1.3;1.0;0.7;1.0;1.0;1.3	117	45;30;36;48;26;41;63;49;25;34;39;41
#		my (
#			$chr_name, $start,  $end, $status, $type,
#			$gene,     $score1, $a,   $score2, $b
#		) = split( " ", $line );
#		if ( $project->isGenome ) {
#
#		#don't construct object > 10000 for genome project better use wisecondor
#			next if abs( $start - $end ) > 10_000;
#		}
#		my $ref_seq = $chromosome->getSequence( $start,     $end );
#		my $first   = $chromosome->getSequence( $start - 1, $start - 1 );
#		my $id;
#		my $sequence_id;
#		my $structType;
#
#		#warn $score1."-".$a."-".$score2;
#		if ( $type eq "del" ) {
#
#			$sequence_id = $first . $ref_seq . "_" . $first;
#			$id         = $chromosome->name . "_" . $start . "_" . $sequence_id;
#			$structType = 'del';
#			$hashRes{$structType}->{$id}->{'id'}              = $id;
#			$hashRes{$structType}->{$id}->{'vcf_id'}          = $id;
#			$hashRes{$structType}->{$id}->{'check_id'}        = $id;
#			$hashRes{$structType}->{$id}->{'isLargeDeletion'} = 1;
#			$hashRes{$structType}->{$id}->{'isLarge'}         = 1;
#			$hashRes{$structType}->{$id}->{'structuralType'}  = 'l_del';
#			$hashRes{$structType}->{$id}->{'structuralTypeObject'} =
#			  'large_deletions';
#			$hashRes{$structType}->{$id}->{'start'}      = $start;
#			$hashRes{$structType}->{$id}->{'end'}        = $end;
#			$hashRes{$structType}->{$id}->{'ref_allele'} = $ref_seq;
#			$hashRes{$structType}->{$id}->{'var_allele'} = "-";
#		}    #end deletion
#		else {
#			$sequence_id = $first . "_" . $first . $ref_seq;
#			$id         = $chromosome->name . "_" . $start . "_" . $sequence_id;
#			$structType = 'ins';
#			$hashRes{$structType}->{$id}->{'id'}                 = $id;
#			$hashRes{$structType}->{$id}->{'vcf_id'}             = $id;
#			$hashRes{$structType}->{$id}->{'check_id'}           = $id;
#			$hashRes{$structType}->{$id}->{'isLargeDuplication'} = 1;
#			$hashRes{$structType}->{$id}->{'isLarge'}            = 1;
#			$hashRes{$structType}->{$id}->{'structuralType'}     = 'l_dup';
#			$hashRes{$structType}->{$id}->{'structuralTypeObject'} =
#			  'large_duplications';
#			$hashRes{$structType}->{$id}->{'start'}         = $start;
#			$hashRes{$structType}->{$id}->{'end'}           = $start + 1;
#			$hashRes{$structType}->{$id}->{'allele_length'} = $end - $start;
#			$hashRes{$structType}->{$id}->{'ref_allele'}    = $ref_seq;
#			$hashRes{$structType}->{$id}->{'var_allele'}    = $ref_seq;
#			$hashRes{$structType}->{$id}->{'sequence'}      = $ref_seq;
#
#			#			$hashRes{$structType}->{$id}->{'var_allele'} = "-";
#		}
#		die() unless $sequence_id;
#		$hashRes{$structType}->{$id}->{'chromosomes_object'} =
#		  { $chromosome->id() => undef };
#
#		$hashRes{$structType}->{$id}->{'line_infos'} = {};
#		$hashRes{$structType}->{$id}->{'references_object'}->{ $reference->id }
#		  = undef;
#
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{method} =
#		  $self->method();
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{'ref_allele'} =
#		  $ref_seq;
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{'var_allele'} =
#		  "-";
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{dp} = '-';
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{nb_all_mut} = 0;
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{nb_all_ref} = 0;
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{score1} =
#		  $score1;
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{score2} =
#		  $score2;
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{score3} = 0;
#
##$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{score3} = int(($score1+$score2)*100);
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{info_vcf} = {};
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{score}    = 0;
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{method_calling}
#		  ->{ $self->method }->{nb_all_ref} = 0;
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{method_calling}
#		  ->{ $self->method }->{nb_all_mut} = 0;
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{method_calling}
#		  ->{ $self->method }->{dp} = '-';
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{method_calling}
#		  ->{ $self->method }->{score1} = $score1;
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{method_calling}
#		  ->{ $self->method }->{score2} = $score2;
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{method_calling}
#		  ->{ $self->method }->{score3} = 0;
#
##$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{method_calling}->{$self->method}->{score3} = int(($score1+$score2)*100);
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{method_calling}
#		  ->{ $self->method }->{info_vcf} = {};
#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{method_calling}
#		  ->{ $self->method }->{score} = 0;
#		if ( $status eq "1/1" ) {
#			$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{he} = 0;
#			$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{ho} = 1;
#			$hashRes{$structType}->{$id}->{annex}->{$patient_id}
#			  ->{method_calling}->{ $self->method }->{he} = 0;
#			$hashRes{$structType}->{$id}->{annex}->{$patient_id}
#			  ->{method_calling}->{ $self->method }->{ho} = 1;
#		}
#		else {
#			$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{he} = 1;
#			$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{ho} = 0;
#			$hashRes{$structType}->{$id}->{annex}->{$patient_id}
#			  ->{method_calling}->{ $self->method }->{he} = 1;
#			$hashRes{$structType}->{$id}->{annex}->{$patient_id}
#			  ->{method_calling}->{ $self->method }->{ho} = 0;
#		}
#
#		$hashRes{$structType}->{$id} =
#		  compress( freeze( $hashRes{$structType}->{$id} ) );
#
#	}
#	return \%hashRes;
}

sub parseBed {
	my ( $self, $reference ) = @_;
	my $project = $self->getPatient()->getProject();
	my $file    = $self->file();
	confess();
	my $chromosome = $reference->getChromosome();
	return {} if $chromosome->name eq "Y";
	my $res = $self->tabix->query_full(
		$self->return_tabix_chromosome_name($chromosome),
		$reference->start, $reference->end );
	return {} unless $res;
	my @data;
	my %hashRes;
	my $patient_id = $self->getPatient()->id;

	while ( my $line = $res->next ) {
		my ( $chr_name, $start, $end, $score1, $score2 ) = split( " ", $line );
		$score1 = 0 unless $score1;
		$score2 = 0 unless $score2;
		my $ref_seq = $chromosome->getSequence( $start,     $end );
		my $first   = $chromosome->getSequence( $start - 1, $start - 1 );
		my $sequence_id = $first . $ref_seq . "_" . $first;
		my $id          = $chromosome->name . "_" . $start . "_" . $sequence_id;

		my $structType = 'del';
		$hashRes{$structType}->{$id}->{'id'}              = $id;
		$hashRes{$structType}->{$id}->{'vcf_id'}          = $id;
		$hashRes{$structType}->{$id}->{'check_id'}        = $id;
		$hashRes{$structType}->{$id}->{'isLargeDeletion'} = 1;
		$hashRes{$structType}->{$id}->{'isLarge'}         = 1;
		$hashRes{$structType}->{$id}->{'structuralType'}  = 'del';
		$hashRes{$structType}->{$id}->{'structuralTypeObject'} =
		  'large_deletions';
		$hashRes{$structType}->{$id}->{'chromosomes_object'} =
		  { $chromosome->id() => undef };
		$hashRes{$structType}->{$id}->{'start'}      = $start;
		$hashRes{$structType}->{$id}->{'end'}        = $end;
		$hashRes{$structType}->{$id}->{'ref_allele'} = $ref_seq;
		$hashRes{$structType}->{$id}->{'var_allele'} = "-";
		$hashRes{$structType}->{$id}->{'line_infos'} = {};
		$hashRes{$structType}->{$id}->{'references_object'}->{ $reference->id }
		  = undef;
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{'ref_allele'} =
		  $ref_seq;
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{'var_allele'} =
		  "-";
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{dp}         = 0;
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{nb_all_mut} = 0;
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{nb_all_ref} = 0;

		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{score1} =
		  $score1;
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{score2} =
		  $score2;
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{score3} = 0;
		if ( $score2 > 0 ) {
			$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{score3} =
			  int( ( $score1 / $score2 ) * 100 );

		}
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{info_vcf} = {};

#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{complex_id} = $x->{COMPLEX};

		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{score} = 0;
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{he}    = 0;
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{ho}    = 1;
		$hashRes{$structType}->{$id} =
		  compress( freeze( $hashRes{$structType}->{$id} ) );

	}
	return \%hashRes;
}

sub getRandomInsertion {
	my ( $self, $row, $header ) = @_;

}

sub parseVCFLine {
	my ( $self, $line ) = @_;
	my ( $chr, $pos, $id, $ref, $alt, $qual, $filter, $info, @all ) = split( "\t", $line );
	confess() if scalar(@all) ne 2;
	my @alts = defined $alt ? split( ",", $alt ) : ();
	my %infos = defined $info ? map { my ( $k, $v ) = split '='; $k => $v } split ';', $info : ();
	my @a = split( ":", $all[0] );
	my @b = split( ":", $all[1] );
	my $gt;
	for ( my $i = 0 ; $i < @a ; $i++ ) {
		$gt->{ $a[$i] } = $b[$i];
	}
	$gt->{genotype}   = $self->genotype( $gt->{GT} );
	$gt->{he}         = 1;
	$gt->{ho}         = 0;
	$gt->{isref}      = 0;
	$gt->{is_cas_1_2} = 0;
		if ( $gt->{genotype}->[0] eq $gt->{genotype}->[1] ) {
		$gt->{he}    = 0;
		$gt->{ho}    = 1;
		$gt->{isref} = 1 if $gt->{genotype}->[0] == 0;

	}
	elsif ( $gt->{genotype}->[0] + $gt->{genotype}->[1] == 1 ) {
		$gt->{he}    = 0;
		$gt->{ho}    = 1;
		$gt->{isref} = 1 if $gt->{genotype}->[0] == 0;

	}
	elsif (($gt->{genotype}->[0] + $gt->{genotype}->[1] ) > 2) {
		$gt->{is_cas_1_2} = 1;
	}
	
	else {
		warn $line;
		warn Dumper  $gt->{genotype};
		
		die();
	}
	return {
		chr    => $chr,
		pos    => $pos,
		id     => $id,
		ref    => $ref,
		alt    => \@alts,
		qual   => $qual,
		filter => $filter,
		infos  => \%infos,
		gt     => $gt,
	};
}

sub genotype {
	my ( $self, $val ) = @_;
	my @t = split( "", $val );
	confess($val) if scalar(@t) ne 3;
	return [ $t[0], $t[-1] ] if $t[0] < $t[-1];
	return [ $t[-1], $t[0] ];    #if $t[0] < $t[-1];
}

sub parseBnd_alt {
	my ( $self, $alt ) = @_;

	#/([a-zA-Z])[\[\]](:\d)[\[\]]/
	my $ref;
	my $before_after;
	if ( $alt =~ /[\[\]]$/ ) {
		$alt =~ s/^([a-zA-Z]*)//;
		$before_after = 1;
		$ref          = $1;
	}
	else {
		$before_after = -1;
		$alt =~ s/([a-zA-Z]*)$//;
		$ref = $1;
	}

	$alt =~ s/^(.)//;
	my $fc = $1;
	$alt =~ s/(.)$//;
	my $bc = $1;
	my ( $chr, $pos ) = split( ":", $alt );
	return {
		ref          => $ref,
		first        => $fc,
		last         => $bc,
		chromosome   => $chr,
		pos          => $pos,
		ref_position => $before_after
	};
}

sub find_linked_bnd {
	my ( $self, $text, $chr, $start, $end ) = @_;
	my $v1         = Bio::DB::HTS::Tabix->new( filename => $self->file() );
	my $patient_id = $self->getPatient->id;
	my %hashRes;
	my $iter = $v1->query( $chr . ":" . $start . "-" . $end );

#my $iter = $v->query($chr->fasta_name.":".$reference->start."-".$reference->end);
	return undef unless $iter;
	my $xs;
	while ( my $row = $iter->next ) {
		next unless $row =~ /$text/;
		my $x       = $self->parseVCFLine($row);
		my $bnd_alt = $self->parseBnd_alt( $x->{alt}->[0] );
		$x->{bnd_alt} = $bnd_alt;
		push( @$xs, $x );
	}
	return $xs;
}

sub genericSVTransLoc {
	my ( $self, $x, $reference ) = @_;
	return;
	my $chr        = $reference->getChromosome();
	my $patient    = $self->getPatient();
	my $patient_id = $patient->id;
	my $text       = $x->{infos}->{MATEID};
	my $astart     = $x->{pos};
	my $aend       = $x->{bnd_alt}->{pos} + 2000;
	if ( $x->{bnd_alt}->{pos} < $x->{pos} ) {
		$astart = $x->{bnd_alt}->{pos} - 2000;
		$aend   = $x->{bnd_alt}->{pos} + 2000;
	}

	#$text = $x->{infos}->{EVENT} if exists $x->{infos}->{EVENT};

	my $hash;
	my $bnds = $self->find_linked_bnd( $text, $x->{bnd_alt}->{chromosome},
		$astart, $aend );

	#warn Dumper $x;
	return             unless $bnds;
	return             unless scalar(@$bnds) == 1;
	die( $self->file ) unless scalar(@$bnds) == 1;
	my $other = $bnds->[0];    #grep {$_->{ID} ne $x->{ID}} @$bnds;

	my $genbo_pos = $x->{pos}
	  ;    #+(length($bnds->[0]->{ref})*$bnds->[0]->{bnd_alt}->{ref_position});

	my $pos_end = $genbo_pos + 1
	  ; #$bnds->[-1]->{pos}+(length($bnds->[-1]->{ref})*$bnds->[-1]->{bnd_alt}->{ref_position});
	my $len = 1;
	my $ref = $x->{ref};
	my $id =
		$x->{chr} . "_"
	  . $x->{pos} . "_"
	  . $x->{ref} . "_bnd-"
	  . $other->{chr} . ":"
	  . $other->{pos};
	my $infos = $x->{infos};
	$hash->{'isSrPr'}               = 1;
	$hash->{'id'}                   = $id;
	$hash->{'isInversion'}          = 1;
	$hash->{'structuralType'}       = 'bnd';
	$hash->{'structuralTypeObject'} = 'boundaries';
	$hash->{'id'}                   = $id;
	$hash->{'vcf_id'} =
	  join( "_", $chr->name, $x->{pos}, $ref, $x->{alt}->[0] );
	$hash->{'isSV'}               = 1;
	$hash->{'chromosomes_object'} = { $chr->id => undef };
	$hash->{'start'}              = $genbo_pos;
	$hash->{'end'}                = $pos_end;
	$hash->{'length'}             = $len;
	$hash->{'ref_allele'}         = $ref;
	$hash->{'var_allele'}         = "++";
	$hash->{'line_infos'}         = "-";
	$hash->{'mate_pos'}           = $other->{pos};
	$hash->{'mate_chr'}           = $other->{chr};
	$hash->{'mate_id'} =
		$other->{chr} . "_"
	  . $other->{pos} . "_"
	  . $other->{ref} . "_bnd-"
	  . $x->{chr} . ":"
	  . $x->{pos};
	$hash->{'event_id'} = undef;
	$hash->{'event_id'} = $x->{infos}->{EVENT} if exists $x->{infos}->{EVENT};
	$hash->{'vcf_position'} = $x->{pos};
	$hash->{'cipos'}        = "0,0";
	$hash->{'cipos'} = $x->{infos}->{CIPOS} if exists $x->{infos}->{CIPOS};
	###OBJECTS
	$hash->{'references_object'}->{ $reference->id } = undef;
	### ANNEX
	$hash->{annex}->{$patient_id}->{Filter} = $bnds->[0]->{filter};

	$hash->{annex}->{$patient_id}->{'ref_allele'} = $ref;
	$hash->{annex}->{$patient_id}->{'var_allele'} = "inv";

	$self->add_DP_AD( $patient_id, $hash, $x );
	$self->add_SR_PR( $patient_id, $hash, $x );
	$hash->{annex}->{$patient_id}->{score}  = $x->{qual};
	$hash->{annex}->{$patient_id}->{he}     = $x->{gt}->{he};
	$hash->{annex}->{$patient_id}->{ho}     = $x->{gt}->{ho};
	$hash->{annex}->{$patient_id}->{method} = $self->method();
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }
	  ->{nb_all_other_mut} = 0;

	$hash->{annex}->{$patient_id}->{method} = $self->method();
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }
	  ->{nb_all_other_mut} = 0;
	i $hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }
	  ->{score} = $x->{qual};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }
	  ->{nb_all_mut} = $hash->{annex}->{$patient_id}->{nb_all_ref};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{he} =
	  $x->{gt}->{he};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{ho} =
	  $x->{gt}->{ho};
	return $hash;
}

sub genericSVInv {
	my ( $self, $x, $reference ) = @_;
	my $chr        = $reference->getChromosome();
	my $patient    = $self->getPatient();
	my $patient_id = $patient->id;
	my $text       = $x->{infos}->{MATEID};
	my $start      = $x->{pos};
	my $end        = $x->{bnd_alt}->{pos} + 2000;
	if ( $x->{bnd_alt}->{pos} < $x->{pos} ) {
		$start = $x->{bnd_alt}->{pos} - 2000;
		$end   = $x->{pos} + 2000;
	}
	$text = $x->{infos}->{EVENT} if exists $x->{infos}->{EVENT};

	my $hash;
	my $bnds =
	  $self->find_linked_bnd( $text, $x->{bnd_alt}->{chromosome}, $start,
		$end );
	$x = $bnds->[0];
	confess( scalar(@$bnds) ) if scalar(@$bnds) < 2 or scalar(@$bnds) > 4;

	my $genbo_pos = $bnds->[0]->{pos} +
	  ( length( $bnds->[0]->{ref} ) * $bnds->[0]->{bnd_alt}->{ref_position} );

	my $pos_end = $bnds->[-1]->{pos} +
	  ( length( $bnds->[-1]->{ref} ) * $bnds->[-1]->{bnd_alt}->{ref_position} );
	my $res = $chr->genesIntervalTree->fetch( $genbo_pos, $pos_end + 1 );
	return if scalar(@$res) > 2;
	my $var_allele = $chr->sequence( $genbo_pos, $pos_end );
	$var_allele = BioTools::complement_sequence($var_allele);

	my $len   = abs( $pos_end - $genbo_pos ) + 1;
	my $ref   = $bnds->[0]->{ref};
	my $id    = $chr->name . "_" . $genbo_pos . "_" . $ref . "_inv-" . $len;
	my $infos = $x->{infos};
	$hash->{'id'}                   = $id;
	$hash->{'isInversion'}          = 1;
	$hash->{'structuralType'}       = 'inv';
	$hash->{'structuralTypeObject'} = 'inversions';
	$hash->{'id'}                   = $id;
	$hash->{'vcf_id'} =
	  join( "_", $chr->name, $bnds->[0]->{pos}, $ref, $bnds->[0]->{alt}->[0] );
	$hash->{'isSV'}               = 1;
	$hash->{'chromosomes_object'} = { $chr->id => undef };
	$hash->{'start'}              = $genbo_pos;
	$hash->{'end'}                = $pos_end;

	$hash->{'length'}     = $len;
	$hash->{'ref_allele'} = $ref;

	#my $var_allele = $chr->sequence($start,$end);
	$hash->{'var_allele'}   = $var_allele;
	$hash->{'line_infos'}   = "-";
	$hash->{'vcf_position'} = $bnds->[0]->{pos};
	$hash->{'isSrPr'}       = 1;
	###OBJECTS
	$hash->{'references_object'}->{ $reference->id } = undef;
	### ANNEX
	$hash->{annex}->{$patient_id}->{Filter} = $bnds->[0]->{filter};
	$hash->{'cipos'}                        = "0,0";
	$hash->{'cipos'}                        = $bnds->[0]->{infos}->{CIPOS}
	  if exists $bnds->[0]->{infos}->{CIPOS};
	$hash->{annex}->{$patient_id}->{'ref_allele'} = $ref;
	$hash->{annex}->{$patient_id}->{'var_allele'} = "inv";

	$self->add_DP_AD( $patient_id, $hash, $x );
	$self->add_SR_PR( $patient_id, $hash, $x );
	$hash->{annex}->{$patient_id}->{score}  = $x->{qual};
	$hash->{annex}->{$patient_id}->{he}     = $x->{gt}->{he};
	$hash->{annex}->{$patient_id}->{ho}     = $x->{gt}->{ho};
	$hash->{annex}->{$patient_id}->{method} = $self->method();
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{nb_all_other_mut} = 0;

	$hash->{annex}->{$patient_id}->{method} = $self->method();
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }
	  ->{nb_all_other_mut} = 0;
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }
	  ->{nb_all_ref} = $hash->{annex}->{$patient_id}->{nb_all_mut};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{score}
	  = $x->{qual};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }
	  ->{nb_all_mut} = $hash->{annex}->{$patient_id}->{nb_all_ref};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{he} =
	  $x->{gt}->{he};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{ho} =
	  $x->{gt}->{ho};
	return $hash;
}

sub genericSVBnd {
	my ( $self, $x, $reference ) = @_;
	my $chr        = $reference->getChromosome();
	my $patient    = $self->getPatient();
	my $patient_id = $patient->id;
	my $hash;
	my $ref = $x->{ref};
	$ref = substr( $ref, 1 );
	my $alt = $x->{alt}->[0];

	#die($alt) if length($alt) > 1;

	my $bnd_alt = $self->parseBnd_alt($alt);
	$x->{bnd_alt} = $bnd_alt;

	if ( $chr->fasta_name eq $bnd_alt->{chromosome} ) {
		return $self->genericSVInv( $x, $reference );
	}
	else {
		return;
		return $self->genericSVTransLoc( $x, $reference );
	}
	confess();
	my $pos       = $x->{pos};
	my $genbo_pos = $pos + 1;

	my $structType = "deletion";

	my $id;

	my $len;
	my $var_allele = $alt;
	my $infos      = $x->{infos};
	$hash->{'isSrPr'} = 1;
	$hash->{'id'}     = $id;
	$hash->{'vcf_id'} = join( "_", $chr->name, $pos, $ref, $x->{alt}->[0] );
	$hash->{'structuralType'}       = "deletion";            #= $allele->{type};
	$hash->{'structuralTypeObject'} = 'deletions';
	$hash->{'isSV'}                 = 1;
	$hash->{'chromosomes_object'}   = { $chr->id => undef };
	$hash->{'start'}                = $genbo_pos;
	$hash->{'end'}                  = $x->{END};
	$hash->{'length'}               = abs( $x->{SVLEN} );
	$hash->{'ref_allele'}           = $ref;
	$hash->{'var_allele'}           = $var_allele;
	$hash->{'line_infos'}           = "-";
	$hash->{'vcf_position'}         = $pos;
	###OBJECTS
	$hash->{'references_object'}->{ $reference->id } = undef;
	### ANNEX
	$hash->{annex}->{$patient_id}->{Filter}       = "PASS";
	$hash->{annex}->{$patient_id}->{'ref_allele'} = $ref;
	$hash->{annex}->{$patient_id}->{'var_allele'} = $var_allele;
	$self->add_DP_AD( $patient_id, $hash, $x );
	$self->add_SR_PR( $patient_id, $hash, $x );
	$hash->{annex}->{$patient_id}->{score}  = $x->{qual};
	$hash->{annex}->{$patient_id}->{he}     = $x->{gt}->{he};
	$hash->{annex}->{$patient_id}->{ho}     = $x->{gt}->{ho};
	$hash->{annex}->{$patient_id}->{method} = $self->method();
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }
	  ->{nb_all_other_mut} = 0;

	$hash->{annex}->{$patient_id}->{method} = $self->method();
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{nb_all_other_mut} = 0;
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{nb_all_ref} = $hash->{annex}->{$patient_id}->{nb_all_mut};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{score} = $x->{qual};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{nb_all_mut} = $hash->{annex}->{$patient_id}->{nb_all_ref};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{he} = $x->{gt}->{he};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{ho} =	  $x->{gt}->{ho};
	return $hash;
}

sub genericSVDel {
	my ( $self, $x, $reference ) = @_;
	my $chr        = $reference->getChromosome();
	my $patient    = $self->getPatient();
	my $patient_id = $patient->id;
	my $hash;
	my $ref = $x->{ref};
	$ref = substr( $ref, 1 );
	my $alt = $x->{alt}->[0];
	my $pos = $x->{pos};
	my $genbo_pos = $pos + 1;
	my $id = $chr->name . "_" . $genbo_pos . "_del-" . abs( $x->{infos}->{SVLEN} );

	my $res = $chr->genesIntervalTree->fetch( $genbo_pos,
		$genbo_pos + abs( $x->{infos}->{SVLEN} ) );
	return if scalar(@$res) > 3;
	my $len;
	my $var_allele = $alt;
	my $infos      = $x->{infos};
	$hash->{'id'}                   = $id;
	$hash->{'isSrPr'}               = 1;
	#$hash->{'isLargeDeletion'}      = 1;
	$hash->{'structuralType'}       = 'deletion';
	$hash->{'structuralTypeObject'} = 'deletions';
	$hash->{'id'}                   = $id;
	$hash->{'vcf_id'} = join( "_", $chr->name, $pos, $ref, $x->{alt}->[0] );
	$hash->{'isSV'}   = 1;
	$hash->{'chromosomes_object'} = { $chr->id => undef };
	$hash->{'start'}              = $genbo_pos;
	$hash->{'end'}                = $x->{infos}->{END};
	$hash->{'length'}             = abs( $x->{infos}->{SVLEN} );
	$hash->{'ref_allele'}         = $ref;
	$hash->{'var_allele'}         = $var_allele;
	$hash->{'line_infos'}         = "-";
	$hash->{'vcf_position'}       = $pos;
	$hash->{'CUUCUUCU'}       = $pos;
	#$hash->{'isCnv'}              = 1;
	###OBJECTS
	$hash->{'references_object'}->{ $reference->id } = undef;
	### ANNEX
	$hash->{annex}->{$patient_id}->{Filter}       = "PASS";
	$hash->{annex}->{$patient_id}->{'ref_allele'} = $ref;
	$hash->{annex}->{$patient_id}->{'var_allele'} = $var_allele;
	$self->add_DP_AD( $patient_id, $hash, $x );
	$self->add_SR_PR( $patient_id, $hash, $x );
	$hash->{annex}->{$patient_id}->{score}  = $x->{qual};
	$hash->{annex}->{$patient_id}->{he}     = $x->{gt}->{he};
	$hash->{annex}->{$patient_id}->{ho}     = $x->{gt}->{ho};
	$hash->{annex}->{$patient_id}->{method} = $self->method();
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }
	  ->{nb_all_other_mut} = 0;

	$hash->{annex}->{$patient_id}->{method} = $self->method();
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }
	  ->{nb_all_other_mut} = 0;
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }
	  ->{nb_all_ref} = $hash->{annex}->{$patient_id}->{nb_all_mut};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{score}
	  = $x->{qual};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }
	  ->{nb_all_mut} = $hash->{annex}->{$patient_id}->{nb_all_ref};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{he} =
	  $x->{gt}->{he};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{ho} =
	  $x->{gt}->{ho};
	return $hash;
}

sub genericSVIns {
	my ( $self, $x, $reference ) = @_;
	my $chr        = $reference->getChromosome();
	my $patient    = $self->getPatient();
	my $patient_id = $patient->id;
	my $hash;
	my $ref       = $x->{ref};
	my $pos       = $x->{pos};
	my $genbo_pos = $pos + length($ref);
	my $id;
	my $alt = $x->{alt}->[0];
	my $len;
	my $var_allele;
	my $infos = $x->{infos};
	$hash->{'structuralType'}       = 'insertion';
	#$hash->{'structuralTypeObject'} = 'large_insertions';

	if ( $infos->{SVTYPE} eq "DUP" ) {
		$len = $infos->{SVLEN};
		$len = $infos->{DUPSVLEN} if exists $infos->{DUPSVLEN};
		return if $len > 2000;
		my $res =
		  $chr->genesIntervalTree->fetch( $genbo_pos, $genbo_pos + $len );
		return if @$res > 2;
	}

	if ( $infos->{SVTYPE} eq "DUP"
		or ( exists $infos->{DUPSVLEN} && $infos->{SVTYPE} eq "INS" ) )
	{
		$len = $infos->{SVLEN};
		$len = $infos->{DUPSVLEN} if exists $infos->{DUPSVLEN};

		#$genbo_pos += $len;
		$id = $chr->name . "_" . $genbo_pos . "_" . $ref . "_dup-" . $len;
		$var_allele = $chr->sequence( $genbo_pos, $genbo_pos + $len );
		$hash->{'isDup'} = 1;
		$hash->{var_allele} = $var_allele;
		#$hash->{'isLargeDuplication'} = 1;
		#$hash->{'isLarge'}              = 1;
		$hash->{'end'}                  = $genbo_pos + $len;
		$hash->{'length'}               = abs($len);
		$hash->{'structuralType'}       = 'insertion';
		$hash->{'structuralTypeObject'} = 'insertions';
		$hash->{'COUCOU'} = 'insertions';
		
		$hash->{'allele_length'}        = $len;

	}
	elsif ( $infos->{SVTYPE} eq "INS"
		&& ( exists $infos->{LEFT_SVINSSEQ} or exists $infos->{RIGHT_SVINSSEQ} )
	  )
	{
		#	return ;
		$len        = 1;
		$id         = $chr->name . "_" . $genbo_pos . "_" . $ref . "_ins-?";
		$var_allele = $infos->{LEFT_SVINSSEQ} . "+" . $infos->{RIGHT_SVINSSEQ};
		$hash->{'var_allele'}    = $var_allele;
		$hash->{'allele_length'} = "~" . length($var_allele);
	}
	elsif ( $infos->{SVTYPE} eq "INS" && $alt =~ /INS/ ) {
	
		$hash->{'CUICUI'} = 'insertions';
		if ( exists $infos->{SVINSSEQ} ) {
			$var_allele              = $infos->{SVINSSEQ};
			$hash->{'var_allele'}    = $var_allele;
			$hash->{'allele_length'} = "~" . length($var_allele);
		}
		else {
			confess( Dumper $infos) unless exists $infos->{LEFT_SVINSSEQ};
		}

	}
	elsif ( $infos->{SVTYPE} eq "INS" && $alt !~ /INS/ ) {
		$alt        = substr( $alt, length($ref) );
		$len        = length($alt);
		$ref        = substr( $ref, 0, 1 );
		$id         = $chr->name . "_" . $genbo_pos . "_" . $ref . "_" . $alt;
		$var_allele = $alt;

		$hash->{'var_allele'}    = $var_allele;
		$hash->{'allele_length'} = length($var_allele);
		
	}
	else {
		confess();
	}
		$hash->{'structuralTypeObject'} = 'insertions';
		$hash->{'structuralType'}       = 'insertion';
	$hash->{'start'}  = $genbo_pos;
	$hash->{'end'}    = $genbo_pos;
	$hash->{'id'}     = $id;
	$hash->{'vcf_id'} = join( "_", $chr->name, $pos, $ref, $x->{alt}->[0] );
	$hash->{'isSrPr'} = 1;
	$hash->{'isSV'}   = 1;
	
	$hash->{'chromosomes_object'} = { $chr->id => undef };

	$hash->{'ref_allele'}   = $ref;
	$hash->{'line_infos'}   = "-";
	$hash->{'vcf_position'} = $pos;

	###OBJECTS
	$hash->{'references_object'}->{ $reference->id } = undef;
	### ANNEX
	$hash->{annex}->{$patient_id}->{Filter}       = "PASS";
	$hash->{annex}->{$patient_id}->{'ref_allele'} = $ref;
	$hash->{annex}->{$patient_id}->{'var_allele'} = "";       #$var_allele;
	$self->add_DP_AD( $patient_id, $hash, $x );
	$self->add_SR_PR( $patient_id, $hash, $x );
	$hash->{annex}->{$patient_id}->{score}  = $x->{qual};
	$hash->{annex}->{$patient_id}->{he}     = $x->{gt}->{he};
	$hash->{annex}->{$patient_id}->{ho}     = $x->{gt}->{ho};
	$hash->{annex}->{$patient_id}->{method} = $self->method();
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }
	  ->{nb_all_other_mut} = 0;

	$hash->{annex}->{$patient_id}->{method} = $self->method();
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }
	  ->{nb_all_other_mut} = 0;
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }
	  ->{nb_all_ref} = $hash->{annex}->{$patient_id}->{nb_all_mut};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{score}
	  = $x->{qual};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }
	  ->{nb_all_mut} = $hash->{annex}->{$patient_id}->{nb_all_ref};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{he} =
	  $x->{gt}->{he};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{ho} =
	  $x->{gt}->{ho};
	return $hash;
}

sub add_DP_AD {
	my ( $self, $patient_id, $hash, $x ) = @_;
	my $a;
	$a .= $x->{gt}->{PR}       if exists $x->{gt}->{PR};
	$a .= "," . $x->{gt}->{SR} if exists $x->{gt}->{SR};
	my @aa = split( ",", $a );
	for ( my $i = 0 ; $i < @aa ; $i += 2 ) {
		$hash->{annex}->{$patient_id}->{nb_all_mut} = $aa[ $i + 1 ];
		$hash->{annex}->{$patient_id}->{nb_all_ref} = $aa[$i];
	}
		warn Dumper $hash unless $hash->{annex}->{$patient_id}->{nb_all_mut};
		warn Dumper $hash unless $hash->{annex}->{$patient_id}->{nb_all_ref};
		$hash->{annex}->{$patient_id}->{dp} =
	  $hash->{annex}->{$patient_id}->{nb_all_mut} +
	  $hash->{annex}->{$patient_id}->{nb_all_ref};

}

sub add_SR_PR {
	my ( $self, $patient_id, $hash, $x ) = @_;
	$hash->{annex}->{$patient_id}->{pr} = "-1,-1";
	$hash->{annex}->{$patient_id}->{sr} = "-1,-1";
	$hash->{annex}->{$patient_id}->{pr} = $x->{gt}->{PR}
	  if exists $x->{gt}->{PR};
	$hash->{annex}->{$patient_id}->{sr} = $x->{gt}->{SR}
	  if exists $x->{gt}->{SR};
}

sub parseVcfFileForReference_manta {
	my ( $self, $reference, $useFilter ) = @_;
	my $chr    = $reference->getChromosome();
	my $v      = Bio::DB::HTS::VCF->new( filename => $self->file() );
	my $v1     = Bio::DB::HTS::Tabix->new( filename => $self->file() );
	my $header = $v->header();
	confess() if $header->num_samples() ne 1;
#	die( $header->get_sample_names->[0] . " " . $self->getPatient->barcode." ".$self->file()." ".$header->num_samples())
#	  if $self->getPatient->name ne $header->get_sample_names->[0]
#	  and $header->get_sample_names->[0] ne $self->getPatient->barcode;

	#confess() if
	my $patient_id = $self->getPatient->id;
	my %hashRes;
	my $iter = $v1->query(
		$chr->fasta_name . ":" . $reference->start . "-" . $reference->end );

#my $iter = $v->query($chr->fasta_name.":".$reference->start."-".$reference->end);
	return {} unless $iter;
	while ( my $row = $iter->next ) {

		my $x = $self->parseVCFLine($row);
		confess( "alt " => Dumper $x) if scalar( @{ $x->{alt} } ) > 1;

		#	confess("ref =>".Dumper $x ) if scalar(@{$x->{ref}}) > 1;
		my $hash;
		if ( $x->{infos}->{SVTYPE} eq "INS" or $x->{infos}->{SVTYPE} eq "DUP" )
		{
			#next;
			$hash = $self->genericSVIns( $x, $chr );
			next unless $hash;

			#my ($self,$x,$patient,$chr,$reference) = @_;
		}
		elsif ( $x->{infos}->{SVTYPE} eq "DEL" ) {
			$hash = $self->genericSVDel( $x, $chr );
			next unless $hash;

			#my ($self,$x,$patient,$chr,$reference) = @_;
		}
		elsif ( $x->{infos}->{SVTYPE} eq "BND" ) {

			$hash = $self->genericSVBnd( $x, $chr );
			next unless $hash;

			#my ($self,$x,$patient,$chr,$reference) = @_;
		}
		else {
			confess( $x->{infos}->{SVTYPE} );
			next;
		}
		my $id = $hash->{id};
		$hash->{'isSrPr'} = 1;
		my $structType = $hash->{structuralType};
		warn Dumper $hash unless $id;
		next unless $id;
		$hashRes{$structType}->{$id} = compress( freeze($hash) );

	}

	return \%hashRes;
}

sub parseVcfFileForReference_melt {
	my ( $self, $reference, $useFilter ) = @_;
	my $chr    = $reference->getChromosome();
	my $v      = Bio::DB::HTS::VCF->new( filename => $self->file() );
	my $v1     = Bio::DB::HTS::Tabix->new( filename => $self->file() );
	my $header = $v->header();
	confess() if $header->num_samples() ne 1;

	#die() if $self->getPatient->name ne  $header->get_sample_names->[0];
	#confess() if
	my $patient_id = $self->getPatient->id;
	my %hashRes;
	my $iter = $v1->query(
		$chr->fasta_name . ":" . $reference->start . "-" . $reference->end );

#my $iter = $v->query($chr->fasta_name.":".$reference->start."-".$reference->end);
	return {} unless $iter;
	my $structType = "insertion";
	my $chr_name   = $chr->fasta_name;
	while ( my $row = $iter->next ) {
		my $x         = $self->parseVCFLine($row);
		my $ref       = $x->{ref};
		my $pos       = $x->{pos};
		my $genbo_pos = $pos + 1;
		confess() if scalar( @{ $x->{alt} } ) > 1;

		my $alt = $x->{alt}->[0];

		next if $alt =~ /DEL/;

		my $type_mei    = $x->{infos}->{SVTYPE};
		my $sequence_id = $ref . '_' . $type_mei;

		#my $type_mei = $type_mei_vcf;

		#$type_mei = "ALU" if $type_mei_vcf =~/ALU/i ;
		#$type_mei = "LINE" if $type_mei_vcf =~/LINE/i ;
		#$type_mei = "L1" if $type_mei_vcf =~/L1/i;

		my $id = $chr->name . "_" . $genbo_pos . "_" . $ref . "_" . $type_mei;

		$hashRes{$structType}->{$id}->{'id'} = $id;
		$hashRes{$structType}->{$id}->{'structuralType'} =
		  "insertion";    #= $allele->{type};
		$hashRes{$structType}->{$id}->{'structuralTypeObject'} = 'insertions';
		$hashRes{$structType}->{$id}->{'isMei'}                = 1;
		$hashRes{$structType}->{$id}->{'mei_type'}             = $type_mei;
		$hashRes{$structType}->{$id}->{'chromosomes_object'} =
		  { $chr->id => undef };
		$hashRes{$structType}->{$id}->{'start'} =
		  $genbo_pos;     # = $allele->{start};
		$hashRes{$structType}->{$id}->{'end'} = $genbo_pos;  # = $allele->{end};
		$hashRes{$structType}->{$id}->{'ref_allele'} =
		  $ref;    #; = $allele->{sequence_ref};
		$hashRes{$structType}->{$id}->{'var_allele'} =
		  "ALU";    #; = $allele->{sequence};
		$hashRes{$structType}->{$id}->{'line_infos'} =
		  "-";      #$allele->{vcf_parse};
		$hashRes{$structType}->{$id}->{'vcf_position'} =
		  $pos;     # = $allele->{vcf_parse}->{POS};
		$hashRes{$structType}->{$id}->{'vcf_id'} =
		  $chr_name . "_" . $pos . "_" . $ref . "_" . $type_mei;
		###OBJECTS
		$hashRes{$structType}->{$id}->{'references_object'}->{ $reference->id }
		  = undef;
		$hashRes{$structType}->{$id}->{'references_object'}->{ $reference->id }
		  = undef;

		### ANNEX

		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{Filter} = "PASS";
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{'ref_allele'} =
		  $ref;
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{'var_allele'} =
		  "ALU";
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{dp} =
		  $x->{gt}->{DP};    # = $allele->{dp}->{raw};
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{nb_all_mut} =
		  $x->{gt}->{AD};    # = $allele->{dp}->{alt};
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{nb_all_ref} =
		  abs( $x->{gt}->{DP} - $x->{gt}->{AD} );    # = $allele->{dp}->{ref};
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{score} =
		  $x->{qual};

		# = $allele->{score};
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{he} =
		  $x->{gt}->{he};
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{ho} =
		  $x->{gt}->{ho};
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{method} =
		  $self->method();    #$allele->{method};#  if exists $allele->{method};
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{method_calling}
		  ->{ $self->method }->{nb_all_other_mut} = 0;

		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{method_calling}
		  ->{ $self->method }->{nb_all_ref} =
		  abs( $x->{gt}->{DP} - $x->{gt}->{AD} );    #$allele->{dp}->{ref};
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{method_calling}
		  ->{ $self->method }->{score} = 0;          #$allele->{score};
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{method_calling}
		  ->{ $self->method }->{nb_all_mut} =
		  $x->{gt}->{AD};                            #$allele->{dp}->{alt};
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{method_calling}
		  ->{ $self->method }->{he} = $x->{gt}->{he};    #$allele->{he} ;
		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{method_calling}
		  ->{ $self->method }->{ho} = $x->{gt}->{ho};    #$allele->{ho}  ;
		$hashRes{$structType}->{$id} =
		  compress( freeze( $hashRes{$structType}->{$id} ) );

	}

	return \%hashRes;
}

sub parseVcfFileForReference_gatk {
	my ( $self, $reference, $useFilter ) = @_;
	my $hash_alleles;
	my $seqpilot;
	unless ($useFilter) { $useFilter = 0; }
	if ( $self->method() eq 'manta' ) {
				confess();
	}
	#	my $project = $self->getPatient()->getProject();
	return {} if $self->isUcsc() eq "empty";
	my $t = time;

	my $debug;
	my $he    = 0;
	my $ho    = 0;
	my $xtime = 0;
	my ( %hashRes, $chr, $idchr, $vcf );
	my $file = $self->file();

	if ($reference) {
		$chr   = $reference->getChromosome();
		$idchr = $chr->name();
		$idchr = $chr->ucsc_name() if $self->isUcsc();
		$vcf = Vcf->new(
			file   => $file,
			region => $idchr . ":" . $reference->start . "-" . $reference->end,
			tabix  => $self->buffer->software("tabix")
		);
	}
	else {
		$vcf =
		  Vcf->new( file => $file, tabix => $self->buffer->software("tabix") );
	}

	$vcf->parse_header();
	my $type;
	my ($find) =
	  grep { lc( $_->{value} ) =~ /seqpilot/ } @{ $vcf->{header}->{source} };
	$seqpilot = 1 if $find;

	#};
	#if ($@){
	#die();
	#return $self->parseVcfFileForReference_samtools($reference, $useFilter);
	#}
	$hash_alleles->{mnp}       = [];
	$hash_alleles->{snp}       = [];
	$hash_alleles->{insertion} = [];
	$hash_alleles->{deletion}  = [];
	$hash_alleles->{indel}     = [];
	while ( my $x = $vcf->next_data_hash() ) {
		
		if ( $self->method() eq 'manta' or $self->method() eq 'dragen-sv' ) {
			confess();

		}
		if ( $self->method() eq 'pbsv') {
			my $vcfRefAllele = $$x{'REF'};
			my $varAllele    = $$x{'ALT'};
			if ($varAllele =~ />/){
				confess();
			}
			
			confess();

		}
		

		if ( ( $useFilter == 1 ) and ( exists $$x{'FILTER'} ) ) {
			die();
			my ($find) = grep { uc($_) ne "PASS" } @{ $x->{'FILTER'} };
			next if $find;
		}
		my $chrom    = $$x{'CHROM'};
		my $varStart = $$x{'POS'};

		my $varEnd;
		my $vcfRefAllele = $$x{'REF'};
		my $varAllele    = $$x{'ALT'};
		my $alleles;
		my @ad;

		my $pat_name;
		unless ( $self->noPatient() ) {
			$pat_name = $self->getPatient->name;
		}

		my $gtypes;

# QueryOnlyVcf -Valeur aléatoire en 0/1 avec option force_all_gt_he pour parser toutes les lignes du vcf sans restriction de gt
		if ( $self->force_all_gt_he() ) {
			$gtypes->{GT} = '0/1';
			$gtypes->{GQ} = '99';
			$gtypes->{AD} = '59,42';
			$gtypes->{PL} = '961,0,1581';
			$gtypes->{DP} = '101';
		}
		elsif ( $file =~ /lofreq/ ) {
			if   ( $x->{INFO}->{AF} > 0.9 ) { $gtypes->{GT} = '1/1'; }
			else                            { $gtypes->{GT} = '0/1'; }
		}
		else {
			$gtypes = $x->{gtypes}->{$pat_name} unless ( $self->noPatient() );
		}

		unless ($gtypes) {
			my @agtypes = values %{ $x->{gtypes} };
			if ( scalar(@agtypes) > 1 ) {
				if ( $self->noPatient() ) {
					my $this_gt;
					foreach my $h (@agtypes) {
						if ( $h->{GT} ne '0/0' ) {
							$gtypes = $h;
							last;
						}
					}
					$gtypes = $agtypes[0] unless ($gtypes);
				}
				else {
					warn "\n\n";
					warn Dumper $x;
					warn "\n\n";
					warn Dumper @agtypes;
					warn "\n\n";
					warn 'Patient: ' . $pat_name;
					warn "\n\n";
					confess();
				}

			}
			else { $gtypes = $agtypes[0]; }
		}
		$gtypes->{GT} =~ s/\|/\//;
		confess($file) unless exists $gtypes->{GT};
		$gtypes->{GT} = "1/1" if $gtypes->{GT} eq "1";

		#$gtypes->{GT} = "1/1" if  $gtypes->{GT} eq "1";

		next if $gtypes->{GT} eq "./.";

		my @gt = split( "/", $gtypes->{GT} );
		if ( scalar(@gt) ne 2 ) { @gt = split( '|', $gtypes->{GT} ); }
		warn " SKIP VARIANT FOUND . " if $gtypes->{GT} eq ".";
		next                          if $gtypes->{GT} eq ".";
		next                          if $gtypes->{GT} eq "0/0";
		next                          if $gtypes->{GT} eq "./0";
		next                          if $gtypes->{GT} =~ /\./;

		if ( scalar(@gt) ne 2 and $file =~ /mutect2/ ) {
			my @lMyGt;
			$lMyGt[0] = $gt[0];

			for ( my $i = 1 ; $i <= scalar(@gt) ; $i++ ) {
				next unless $gt[$i];
				next if ( $gt[$i] eq "/" );

				next if ( $gt[$i] == 0 );
				push( @lMyGt, $gt[$i] );
			}
			@gt = @lMyGt;
			$gtypes->{GT} = join( '/', @gt );
		}

		if ( scalar(@gt) eq 1 ) {
			if ( $gt[0] =~ /[0-9]/ ) {
				push( @gt, $gt[0] );
				$gtypes->{'GT'} = join( '/', @gt );
			}
		}

		if ( scalar(@gt) ne 2 ) {
			warn Dumper $x;
			warn Dumper @gt;
			warn "\n\n";
			warn 'VCF FILE: ' . $self->file();
			die();
		}

		my %hgt;
		map { $hgt{$_}++ } @gt;

		my %alleles;
		my $index = 0;

		#############
		#find DP
		############
		my ( $pr, $sr );
		my $dp = 0;
		if ( exists $gtypes->{DP} ) {
			$dp = $gtypes->{DP};
		}
		elsif ( exists $x->{INFO}->{DP} ) {
			$dp = $x->{INFO}->{DP};
		}
		elsif ( exists $x->{INFO}->{TC} ) {
			$dp = $x->{INFO}->{TC};
		}

		if ( exists $gtypes->{AD} ) {
			@ad = split( ",", $gtypes->{AD} );
		}
		elsif ( exists $x->{INFO}->{DP4} ) {
			my @ad2 = split( ",", $x->{INFO}->{DP4} );
			my $z   = 0;
			for ( my $i = 0 ; $i < @ad2 ; $i += 2 ) {
				$ad[$z] = $ad2[$i] + $ad2[ $i + 1 ];
				$z++;
			}
		}
		elsif ( exists $gtypes->{AO} ) {
			if ( exists $gtypes->{RO} ) {
				$ad[0] = $gtypes->{RO};
			}
			else {
				$ad[0] = $x->{INFO}->{SRR} + $x->{INFO}->{SRF};
			}

			#die() unless $x->{INFO}->{SRR} <0;
			push( @ad, split( ",", $gtypes->{AO} ) );
		}
		elsif ( exists $gtypes->{'AF'} ) {
			$ad[0] = $dp;
			foreach my $p ( split( ",", $gtypes->{'AF'} ) ) {
				my $nb = int( $p * $dp );
				$ad[0] -= $nb;
				push( @ad, $nb );
			}

		}
		elsif ( $gtypes->{'GT'} eq '1/1' ) {
			$ad[0] = '0';
			if ( exists( $gtypes->{'DP'} ) ) { $ad[1] = $gtypes->{'DP'}; }
			elsif ( exists( $x->{'INFO'}->{'DP'} ) ) {
				$ad[1] = $x->{'INFO'}->{'DP'};
			}
			else { $ad[1] = '1'; }
		}
		elsif ( exists $x->{'INFO'}->{'AB'} ) {
			if ( $x->{'INFO'}->{'AB'} eq '0.00' ) {
				next;
			}
			else {
				my @AB    = split( ",", $x->{'INFO'}->{'AB'} );
				my $nbAlt = int( $gtypes->{'DP'} ) * $AB[0];
				$ad[0] = int( $gtypes->{'DP'} ) - int($nbAlt);
				$ad[1] = int($nbAlt);
			}
		}
		elsif ( exists $gtypes->{'NV'} ) {

			die() unless ( exists $gtypes->{'NR'} );
			my @NR = split( ",", $gtypes->{NR} );
			my @NV = split( ",", $gtypes->{NV} );
			@ad = ( $NR[0], @NV );

		}
		else {
			$ad[0] = '1';
			if ( exists( $gtypes->{'DP'} ) ) { $ad[1] = $gtypes->{'DP'}; }
			elsif ( exists( $x->{'INFO'}->{'DP'} ) ) {
				$ad[1] = $x->{'INFO'}->{'DP'};
			}
			else { $ad[1] = '1'; }
		}

		# Champ MANTA PR a cumuler avec SR pour le nb reads REF/ALT


		# correct bug in GATK allele with no reads
		my $test = 0;
		for ( my $i = 0 ; $i < @ad ; $i++ ) {
			next if $ad[$i] eq ".";
			$test += $ad[$i];
		}
		my @types_variation = ( "", "", "", "", "" );
		my @cigars          = ( "", "", "", "", "" );

		my @lTmp = split( "/", $gtypes->{'GT'} );
		unless ( scalar(@lTmp) == 2 ) {
			warn Dumper $gtypes;
			die;
		}

		unless ( $lTmp[0] =~ /^\d+?$/ ) {
			warn Dumper $gtypes;
			die;
		}
		unless ( $lTmp[1] =~ /^\d+?$/ ) {
			warn Dumper $gtypes;
			die;
		}

		my @gt_allele =
		  uniq( sort { $a <=> $b } split( "/", $gtypes->{'GT'} ) );

		my $max = $gt_allele[-1];
		for ( my $i = 0 ; $i <= $max + 1 ; $i++ ) {
			$cigars[$i]          = "";
			$types_variation[$i] = "";
		}

		my $nb_alleles = 0;
		my %tall;
		map { $tall{$_}++ if $_ ne '0' } split( "\/", $gtypes->{'GT'} );
		my $nb_all = scalar( keys %tall );

		if ( exists $x->{'INFO'}->{TYPE} ) {
			my @t = split( ",", $x->{'INFO'}->{TYPE} );
			my @c = split( ",", $x->{'INFO'}->{CIGAR} . "" )
			  if exists $x->{'INFO'}->{CIGAR};

			if ( scalar(@t) == scalar(@$varAllele) ) {
				@cigars          = @c if exists $x->{'INFO'}->{CIGAR};
				@types_variation = @t;
			}

		}
		my $set = "";
		$set = $x->{'INFO'}->{set} if exists $x->{'INFO'}->{set};

	  ALLELE: foreach my $num_allele (@gt_allele) {

			next if $num_allele eq 0;

			my $index = $num_allele - 1;

			my $vcfVarAllele;

			$vcfVarAllele = $varAllele->[$index];
			next if $vcfVarAllele eq "*";
			my $freebayes_type = $types_variation[$index];
			my $allele;
			if ( $self->method() eq 'HGMD' ) {
				$allele = $x;
				my @lHgmdCat = keys %{$x};
				$allele->{HGMD_CAT_PARSED} = \@lHgmdCat;
			}
			$allele->{ref_all_ori} = $vcfRefAllele;
			$allele->{var_all_ori} = $vcfVarAllele;
			$allele->{ref}         = $vcfRefAllele;
			$allele->{chromosome}  = $x->{'CHROM'};
			my $freebayes_cigar = $cigars[$index] . "";
			$allele->{cigar} = $freebayes_cigar;
			$allele->{start} = $x->{'POS'};
			$allele->{alt}   = $vcfVarAllele;

			$allele->{is_imprecise} = undef;
			my ( $type, $len, $ht );

				( $type, $len, $ht ) = $vcf->event_type( $vcfRefAllele, $vcfVarAllele );

			my $htype;
			warn "LEN ==>" . $len . " type : $type " if $debug;
			if ( $type eq 's' && length($vcfVarAllele) == 1 ) {
				$htype = "snp";
			}
			elsif ( $set eq 'freebayes' ) {
				$htype = "indel";
			}
			elsif ( $type eq "i" && $len >= 0 ) {
				$htype = "insertion";

			}
			elsif ( $type eq "i" && $len < 0 ) {
				$htype = "deletion";

			}
			else {
				$htype = "indel";
			}

			if ($set) {
				if ( lc($set) eq "intersection" ) {
					$allele->{method} = "inter";
				}
				else {
					$allele->{method} = substr( lc($set), 0, 4 );
				}
			}

			$allele->{vtype}     = $type;
			$allele->{ht}        = $ht;
			$allele->{len}       = $len;
			$allele->{index}     = $index;
			$allele->{dp}->{ref} = $ad[0];
			$allele->{dp}->{alt} = $ad[$num_allele];
			$allele->{gt_hash}   = $gtypes;
			$allele->{gt}        = $gtypes->{GT};

			$allele->{pos_vcf} = $x->{'POS'};
			$allele->{dp}->{raw} = $dp;

			$allele->{score} = $x->{'QUAL'};

			if ( scalar(@ad) == 0 ) {

				$allele->{dp}->{alt} = $dp;
				$allele->{dp}->{ref} = 0;
			}
			else {
				$allele->{dp}->{ref} = $ad[0];
				$allele->{dp}->{alt} = $ad[$num_allele];
			}
			$dp = 0 if $dp eq ".";
			if ( $dp == 0 ) {
				$allele->{dp}->{raw} =
				  $allele->{dp}->{ref} + $allele->{dp}->{alt};
			}

			if ( exists $x->{'FILTER'} ) {
				$allele->{filter} = $x->{'FILTER'}->[0];
			}
			$allele->{vcf_parse} = $x;
			##############
			#parse he /ho
			#############
			my $nb_alleles = 0;
			my @al         = sort { $a <=> $b } split( "/", $allele->{gt} );
			map { $nb_alleles++ if ( $_ ne 0 || $_ ne "." ) } split( "/", $allele->{gt} );
			if ( $al[0] eq $al[1] ) {
				$allele->{he} = 0;
				$allele->{ho} = 1;
				$ho++;
			}

			#					  #little trick for "1/2" view as ho
			#					  if ($al[0] ne"0"){
			#					   $allele->{he} = 0 ;
			#						$allele->{ho} = 1 ;
			#						$ho ++;
			#					  }
			else {
				# un seul allele donc homozygote
				$allele->{he} = 1;
				$allele->{ho} = 0;
				$he++;

			}

			if ( $allele->{he} == 1 and $al[0] ne 0 ) {
				$allele->{is_cas_1_2} = 1;
			}

			push( @{ $hash_alleles->{$htype} }, compress( freeze $allele) );
			$index++;
		}

	}
	$vcf->close();
	############
	# Complexes
	############

	while ( @{ $hash_alleles->{indel} } ) {
		my $debug;

		my $allele  = thaw( decompress( shift @{ $hash_alleles->{indel} } ) );
		my $allelec = clone($allele);

		$allelec->{sequence_id}  = $allele->{ref} . "_" . $allele->{alt};
		$allelec->{type}         = "complex";
		$allelec->{obj}          = "complex";
		$allelec->{end}          = $allelec->{start} + length( $allele->{ref} );
		$allelec->{sequence}     = $allele->{var_all_ori};
		$allelec->{iupac_allele} = "-";
		$allelec->{sequence_ref} = $allele->{ref};
		$allele->{ALT2}          = $allele->{alt};
		my $cigar = $allele->{cigar};

		#		$debug = 1 if  $allele->{pos_vcf} == 74315212;

		#		warn $cigar if $debug;

		if ( $cigar eq "" || $cigar =~ /^\dD/ ) {
			my $seq1 = "NNNNNNNNNN";

			unless ($reference) {
				my $chr_name = $allele->{chromosome};
				$chr_name =~ s/chr//;
				$chr = $self->getThisChromosome($chr_name);
			}

			#special case alignement start with a deletion
			my $len_add = 1;

#	my $first_base = 	 $chr->getSequence($allele->{start}-$len_add,$allele->{start}-$len_add);
#	warn $allele->{alt}."\n".$allele->{ref};
#	die();
			$allele->{alt} = "A" . $allele->{alt};
			$allele->{ref} = "A" . $allele->{ref};

			# $cigar;
			if ( length( $allele->{alt} ) eq length( $allele->{ref} ) ) {
				my @a = split( "", $allele->{alt} );
				my @b = split( "", $allele->{ref} );
				my $diff;
				for ( my $i = 0 ; $i < @a ; $i++ ) {
					if ( $a[$i] ne $b[$i] ) {
						$diff = 999 if $diff;
						$diff = $i;
					}
				}
				if ( $diff && $diff ne 999 ) {
					my $start = $diff;
					my $end   = length( $allele->{alt} ) - $start - 1;
					$cigar = ($start) . "M1X" . $end . "M";

				}
				$allele->{start} -= $len_add;
			}

			unless ($cigar) {
				$allele->{start} -= $len_add;
				my $t = time;
				$cigar = $self->create_cigar_sw( $seq1 . $allele->{ref} . $seq1,
					$seq1 . $allele->{alt} . $seq1, 1 );

				#warn "cigar ".$cigar." $len_add \n" if $debug;
				#warn Dumper parse_cigar($allele,$cigar);
				$xtime += abs( time - $t );
			}

		}
		else {

			$cigar = $allele->{cigar};

			#confess();
		}

		my $event;
		if ( $cigar =~ /^\dD/ ) {
			$allele->{start} -= 1;
			$event = parse_cigar( $allele, "1M" . $cigar );
		}
		else {
			$event = parse_cigar( $allele, $cigar );
		}

		my $complex_id =
			$allele->{chromosome} . "_"
		  . $allele->{start} . "_"
		  . $allele->{ref} . "_"
		  . $allele->{alt};

		foreach my $e (@$event) {
			my $new_allele = clone($allele);
			$new_allele->{start} = $e->{start};
			$new_allele->{ref}   = $e->{ref};
			$new_allele->{alt}   = $e->{alt};
			$new_allele->{cigar} = $e->{cigar};
			my ( $type, $len, $ht ) = $vcf->event_type( $e->{ref}, $e->{alt} );
			$new_allele->{ht}  = $ht;
			$new_allele->{len} = $len;

			$new_allele->{complex_id} = $complex_id;
			my $type1 = $e->{type};
			push( @{ $hash_alleles->{$type1} }, compress( freeze $new_allele) );
		}

	}
	delete $hash_alleles->{indel};

	############
	# mnp
	############
	while ( @{ $hash_alleles->{mnp} } ) {
		my $allele = thaw( decompress( shift @{ $hash_alleles->{mnp} } ) );

	  #foreach my $allele ( map{thaw(decompress($_))} @{$hash_alleles->{mnp}} ){
		$allele->{type} = 'mnp';
		$allele->{obj}  = 'mnps';
		my $addlen = length( $allele->{ref} ) - 1;
		confess() if length( $allele->{alt} ) ne length( $allele->{ref} );
		my $complex_id =
			$allele->{chromosome} . "_"
		  . $allele->{start} . "_"
		  . $allele->{ref} . "_"
		  . $allele->{alt};
		my @ref_bases = split( "", $allele->{ref} );
		my @alt_bases = split( "", $allele->{alt} );
		my $find;

		for ( my $i = 0 ; $i < @alt_bases ; $i++ ) {

			my $new_allele = clone($allele);
			next if $ref_bases[$i] eq $alt_bases[$i];
			$find++;

			$new_allele->{start}      = $allele->{start} + $i;
			$new_allele->{ht}         = $alt_bases[$i];
			$new_allele->{len}        = 1;
			$new_allele->{complex_id} = $complex_id;
			$new_allele->{ref}        = $ref_bases[$i];
			$new_allele->{alt}        = $alt_bases[$i];
			push( @{ $hash_alleles->{snp} }, compress( freeze $new_allele) );

		}

	}

	#	delete $hash_alleles->{mnp};
	#	 $hash_alleles->{mnp} =[];
	############
	# Deletions
	############
#	chr17	29587480	.	TC	G	1980.8	PASS	.	GT:GQ:DP:FDP:RO:FRO:AO:FAO:AF:SAR:SAF:SRF:SRR:FSAR:FSAF:FSRF:FSRR	0/1:99:380:377:180:183:193:194:0.514589:96:97:101:79:97
#:97:102:81
	my $temp_all;
	while ( @{ $hash_alleles->{deletion} } ) {
		my $allele = thaw( decompress( shift @{ $hash_alleles->{deletion} } ) );

#	foreach my $allele (map{thaw(decompress($_))} shift @{$hash_alleles->{deletion}} ){
		$allele->{type} = 'del';
		$allele->{obj}  = 'deletions';
		if ( $self->method() eq 'manta' or $self->method() eq 'dragen-sv' or $self->method() eq 'Sniffles2' ) {
			$allele->{method} = 'manta';
			
			if ( abs( $allele->{len} ) >= 50 ) {
				$allele->{type} = 'l_del';
				$allele->{obj}  = 'deletions';
			}
		}
		die() if $allele->{len} == 0;
		$allele->{len} = abs( $allele->{len} );
		if ( $self->method() eq 'manta' or $self->method() eq 'dragen-sv' ) {
			$allele->{sequence_ref} = $allele->{'ref'};
			$allele->{sequence_id}  = $allele->{'ref'} . '_' . $allele->{'alt'};
			$allele->{sequence}     = "-";
			$allele->{iupac_allele} = "-";
			$allele->{PARSE}        = "tutu";
			$allele->{end}          = $allele->{vcf_parse}->{INFO}->{END};
			$allele->{start}        = $allele->{start} + 1;
			push( @$temp_all, compress( freeze $allele) );
			next;
		}
		my $first = substr( $allele->{alt}, 0, 1 );
		$allele->{sequence_id}  = $first . $allele->{ht} . "_" . $first;
		$allele->{sequence}     = "-";
		$allele->{iupac_allele} = "-";
		my @lTmpRef = split( '', $allele->{ref} );
		$allele->{sequence_ref} = $allele->{ht};
		my $dpos = index( $allele->{ref}, $allele->{ht} );

		if ( $dpos == -1 ) {
			warn "\n\n";
			warn Dumper $allele;
			warn "\n\n";
			warn $allele->{ref};
			warn "\n\n";
			warn $allele->{ht};
			confess();
		}
		$allele->{start} = $allele->{start} + $dpos;
		$allele->{end}   = $allele->{start} + $allele->{len} - 1;
		confess() if $allele->{end} < $allele->{start};
		$allele->{PARSE} = "tutu";
		push( @$temp_all, compress( freeze $allele) );

	}
	$hash_alleles->{deletion} = $temp_all;

	$temp_all = [];
	############
	# insertions
	############
	my @iall;
	while ( @{ $hash_alleles->{insertion} } ) {
		my $allele =
		  thaw( decompress( shift @{ $hash_alleles->{insertion} } ) );

		# insertion
		$allele->{type} = 'ins';
		$allele->{obj}  = 'insertions';
		if ( $self->method() eq 'manta' or $self->method() eq 'dragen-sv' ) {
			$allele->{method} = $self->method();
			if ( abs( $allele->{len} ) >= 50 ) {
				$allele->{type} = 'insertion';
				$allele->{obj}  = 'insertions';
			}
		}
		my $dpos = index( $allele->{alt}, $allele->{ht} );

		#						$dpos = 0 if ($self->method() eq 'manta');
		if ( $dpos == 0
			and
			( $self->method() eq 'manta' or $self->method() eq 'dragen-sv' ) )
		{
			unless ($reference) {
				my $chr_name = $allele->{chromosome};
				$chr_name =~ s/chr//;
				$chr = $self->getThisChromosome($chr_name);
			}
			my $first_base =
			  $chr->getSequence( $allele->{start} - 1, $allele->{start} - 1 );
			$allele->{start} -= 1;
			$allele->{ref} = $first_base . $allele->{ref};
			$allele->{alt} = $first_base . $allele->{alt};
			$dpos = index( $allele->{alt}, $allele->{ht} );
		}
		my $first = substr( $allele->{ref}, 0, 1 );

		$allele->{sequence_id}  = $first . "_" . $first . $allele->{ht};
		$allele->{sequence_ref} = "-";
		$allele->{sequence}     = $allele->{ht};

		if ( $dpos == -1 ) {

			#cas de l'insertion au milieu de la sequence
			# exemple : TATTTTA ref :TTTTTA
			warn "\n\n";
			warn Dumper $allele;
			warn "\n\n";
			warn 'ALT: ' . $allele->{alt};
			warn 'HT: ' . $allele->{ht};
			confess() if $dpos == -1;
		}

		else {
			$allele->{start} = $allele->{start} + $dpos;
		}
		$allele->{end}   = $allele->{start};
		$allele->{PARSE} = "tutu";

		# ID: vcf position
		# ID: start position
		push( @$temp_all, compress( freeze $allele) );

	}
	$hash_alleles->{insertion} = $temp_all;

	$temp_all = [];
	############
	# snps
	############
	while ( @{ $hash_alleles->{snp} } ) {
		my $allele = thaw( decompress( shift @{ $hash_alleles->{snp} } ) );

#	foreach my $allele (map{thaw(decompress($_))} shift @{$hash_alleles->{snp}} ){
# SNP
		$allele->{type} = 'snp';
		$allele->{obj}  = 'variations';
		my $addlen = length( $allele->{ref} ) - 1;
		if ( length( $allele->{ref} ) > 1 ) {

			#warn $vcfRefAllele;
			$allele->{sequence_ref} = substr $allele->{ref}, -1;
		}
		else {
			$allele->{sequence_ref} = $allele->{ref};
		}
		$allele->{sequence_id} = $allele->{sequence_ref} . "_" . $allele->{ht};

		$allele->{sequence} = $allele->{ht};

		$allele->{start} += $addlen;

		$allele->{end} = $allele->{start};
		push( @$temp_all, compress( freeze $allele) );
	}
	$hash_alleles->{snp} = $temp_all;

	$temp_all = [];
	############
	# mnp
	############
	my @mall;

#	foreach my $allele ( map{thaw(decompress($_))} shift @{$hash_alleles->{mnp}} ){
#					confess(Dumper($allele));
#						$allele->{type} = 'mnp';
#						$allele->{obj} = 'mnps';
#						my $addlen = length($allele->{ref}) -1;
#						confess() if  length($allele->{alt}) ne  length( $allele->{ref});
#						$allele->{sequence_ref}  = $allele->{ref};
#						$allele->{sequence_id} = $allele->{sequence_ref}."_".$allele->{ht};
#						if (length ($allele->{sequence_ref}) >2){
#						my $dpos = index($allele->{alt},$allele->{ht});
#						confess(Dumper $allele->{vcf_parse}) if $dpos == -1;
#
#						$allele->{start} = $allele->{start} +$dpos;
#						$allele->{end} = $allele->{start} + $allele->{len} - 1;
#						$allele->{sequence_ref} = substr($allele->{sequence_ref},$dpos,$allele->{len});
#						}
#						$allele->{sequence} =  $allele->{ht};
#						my @seq_alt = split('',$allele->{sequence});
#						my @seq_ref = split('',$allele->{sequence_ref});
#						$allele->{iupac_allele} ='';
#						for (my $i=0;$i< @seq_alt;$i++){
#							eval {
#							$allele->{iupac_allele} .= $reference->getProject()->biotools->getIUPAC([$seq_ref[$i] ,$seq_alt[$i]]);
#						};
#						}
#						$allele->{PARSE} = "tutu";
#
#						#$allele->{start} =  $$x{'POS'};
#						$allele->{end} = $allele->{start}+$addlen;
#					push(@mall,compress(freeze($allele)));
#
#}

	#construct hash for GenBo
	my $patient_id;
	unless ( $self->noPatient() ) { $patient_id = $self->getPatient->{id}; }

	#warn $xtime;
	#die();

	#foreach my $t (keys %{$hash_alleles}) {
	#	warn "$t: ".scalar (@{$hash_alleles->{$t}});
	#}
	#die;
	foreach my $type ( keys %{$hash_alleles} ) {
		die() if $type eq "indel";
		next  if $type eq "indel";
		next unless $hash_alleles->{$type};

		while ( @{ $hash_alleles->{$type} } ) {
			my $allele =
			  thaw( decompress( shift @{ $hash_alleles->{$type} } ) );

 #foreach my $allele (map{thaw(decompress($_))} shift@{$hash_alleles->{$type}}){
			my $structType = $allele->{type};

			# ID: vcf position
			#my $id = $chr->name."_".$$x{'POS'}."_".$allele->{sequence_id};
			# ID: start position
			my ( $chr_name, $chr_id );
			if ($reference) {
				$chr_name = $chr->name;
				$chr_id   = $chr->id();
			}
			else {
				$chr_name = $allele->{chromosome};
				$chr_name =~ s/chr//;
				$chr_id = $chr_name;
			}
			my $id =
			  $chr_name . "_" . $allele->{start} . "_" . $allele->{sequence_id};
			die()                    unless exists $allele->{start};
			die()                    unless exists $allele->{sequence_id};
			die()                    unless exists $allele->{ref_all_ori};
			die()                    unless exists $allele->{var_all_ori};
			die()                    unless exists $allele->{pos_vcf};
			die()                    unless exists $allele->{chromosome};
			die()                    unless exists $allele->{type};
			die()                    unless exists $allele->{obj};
			confess( Dumper $allele) unless exists $allele->{end};
			die()                    unless exists $allele->{sequence_ref};
			die()                    unless exists $allele->{sequence};
			die()                    unless exists $allele->{vcf_parse};
			die()                    unless exists $allele->{dp}->{raw};
			die()                    unless exists $allele->{dp}->{ref};
			die()                    unless exists $allele->{dp}->{alt};
			die()                    unless exists $allele->{he};
			die()                    unless exists $allele->{ho};
			die()                    unless exists $allele->{score};

			next if exists $hashRes{$structType}->{$id};

			if ( $self->method() eq 'HGMD' ) {
				foreach my $cat ( @{ $allele->{HGMD_CAT_PARSED} } ) {
					$hashRes{$structType}->{$id}->{$cat} = $allele->{$cat};
				}
			}

			$hashRes{$structType}->{$id}->{'id'} = $id;
			$hashRes{$structType}->{$id}->{'vcf_id'} =
				$allele->{chromosome} . "_"
			  . $allele->{pos_vcf} . "_"
			  . $allele->{sequence_id};
			$hashRes{$structType}->{$id}->{'check_id'} =
				$allele->{chromosome} . "_"
			  . $allele->{pos_vcf} . "_"
			  . $allele->{ref_all_ori} . "_"
			  . $allele->{var_all_ori};

			$hashRes{$structType}->{$id}->{'structuralType'} = $allele->{type};
			$hashRes{$structType}->{$id}->{'structuralTypeObject'} =
			  $allele->{obj};
			$hashRes{$structType}->{$id}->{'chromosomes_object'} =
			  { $chr_id => undef };
			$hashRes{$structType}->{$id}->{'start'} = $allele->{start};
			$hashRes{$structType}->{$id}->{'end'}   = $allele->{end};
			$hashRes{$structType}->{$id}->{'ref_allele'} =
			  $allele->{sequence_ref};
			$hashRes{$structType}->{$id}->{'var_allele'} = $allele->{sequence};
			$hashRes{$structType}->{$id}->{'line_infos'} =
			  "-";    #$allele->{vcf_parse};
			$hashRes{$structType}->{$id}->{'vcf_position'} =
			  $allele->{vcf_parse}->{POS};

#$hashRes{$structType}->{$id}->{'vcf_sequence'} =  $allele->{vcf_parse}->{ALT}->[0]."";
			$hashRes{$structType}->{$id}->{'is_imprecise'} =
			  $allele->{is_imprecise};

			#warn Dumper $hashRes{$structType}->{$id}->{'line_infos'};
			if ($reference) {
				$hashRes{$structType}->{$id}->{'references_object'}
				  ->{ $reference->id } = undef;
			}
			else {
				$hashRes{$structType}->{$id}->{'references_object'}
				  ->{'no_reference'} = undef;
			}
			unless ( $self->noPatient() ) {
				if ( exists $allele->{filter} ) {
					$hashRes{$structType}->{$id}->{annex}->{$patient_id}
					  ->{Filter} = $allele->{filter};
				}

				#die();
				my $debug;
				$hashRes{$structType}->{$id}->{annex}->{$patient_id}
				  ->{'ref_allele'} = $allele->{sequence_ref};

#my $debug;
#$debug = 1 if $hashRes{$structType}->{$id}->{annex}->{$patient_id}->{'ref_allele'} eq "CTGTGAAAATGCTTTGTGATGTGTGGATTCATCTCACAGAATGAAACCTGTGTTTTGATTCACTAGATTGGAAACACTCTTTTTGTGAATCTAGGAAGGGACATTTCTGAGCCCAATGTGGCCTATAGCAAGAAACCAAATATCCCACGATAAAAAGTGGAAGCAAGCTATCTTTGAAAATGCTTTGTCATACGTTAATTCATCTTACAGAATGAAACTTGTGTTTTGATTCACCAGGTAAAAAACAATTTTTTTTGTAGACTCTAAATAGGGACATTTCTAAGCCTATTGAAGCCTATAGTGAAGAACCAAATATCACGTGATGAAACCTCAAAACAAACTATCTGTGAAAATGCTATGTGATGTGTGGATTCAATTCACTGAAGGGAACCTGTGTTTTAATTCACCAGGTTAGAAACACTCTTTCTGTAGAATCTACAAAATGACATTCCTGAGCTCATTGAAGCCTACAGTGAGAAACTATATATCCCTTGATGAAACTAGAAACCAACTACCTGTGAAAATGTTTTGTGATGGGTGGATTCATCTCACAGAATGGAATCTGTGTTTTGATTTACCAGGTAGGAAACACTGTTCTTGCGGAATCTACTAAGGGACATTTCATAGCCCATTGAGGCCTATAGTGAAAACCAAATATCCCGCGATACAAACTAGAAACAAGCTATCTGTGAAAAATGCTCTGTGCCGTGTGGATTCACCTCACAGGGCTAAACCTATGTTTTCATTTACCAGGTTAGAAACACTTTTTTTTTATATATAGAATCTACAAATGGACATTTCTGAGCTCATTGAAGTCCATAGTGAAAAACTGAATATCCCACGATAAAAACTAAAACAAGCTATCTAAAAAAATGCTTTGTGACATGTGAATTCTTCTCACAGAATGGAAACTGTGTTTTTATTCACCAGGTTGGAAACATTCTTTTTGTAGTATATACAAAGGGACATTCCTGTTCCCAGTGAAACCTATATTATAAAACCAAATATCCTGTGATTTAAACAGAAACAAGCTATGTGTAAAACTGCTTTGTGATGTGTGAATTCATCTCACAGAATAGAACCTGTGTTTTGATTCACCTGGTTGGAAACACTCCTTCTGTACTATATGCGAAGGAACATTTCTGAGTCCATTGAAGCCCATAGTGAAAAACCGATTACCTGTGATAAAAACTAGAAACAAGCTATCTGTGAAAATACTTTGTGATGTGTGAATTTATCTCACAGAATGAAACCTTTGTTTTCACTCACAAGGTTGGAAACATTCTGTTTGTAGTATATACGAAGAAACATTCCTGTTCCCGTTGAAGCCTATAATGTAAAACCGATTATGCTGTGATAAAAACTAGAAACAATCAATCTGTCAAAATGCTTTGTTCTGTGTGGATTCATCTCACAATGGAACCAGTGTTTTGATTCATCAGGTAGGAAACTATCTTTTTTGTAGAATCTATGAAGGGACATTTCAGAGCCCATTGAAGCCTACAGTTATAAATCCAATATCCTGCGATAAAAACTGGAAACAAACTATCTGTGAAAATGCTTTGTGATGTGTGGATTCATCTCCATGATCAAACCTGTGTTGTGATTTACCAGATTAGAAAAGCTTTTCTTTTAAGATCTACAAAGGGACATTTCTAAGCCCATTGAGTCCTATAGTGAAAAACATAATATACCATGATAAAAACTAGAAACAAGCTATAGGTGAAAATACTTTGCACTGTGTGGATTCATCTCACAGAATGGAAACTGTGTTTTGATTCACCAAGTTGGAAACAATCTTTTTGTAGAATCTACAGAGGGACATTTCTGAGCCCATTATGGCCTATAGTGAAAAACTGAATATCCCACGCTAAAAACTAGAAAAAAAACTATCTGTGAAAATGCTTTGTGACAGGTGGATTCATCTAATAGAATGGAAACTGGTTTTGATTCACCAGGTTGGCAACACTCCTTTTGTAGAATTTATGAAGGGACATTTCTGAGCCCATTGAAGCCTATATTGAAAACCAAATATCCCCCAGTAAAAACTTGAAACAAGCTATCTGTGAAAATGCTTTGTGATGTGTGGATTCATATCAGAGAATGGAACCTGTGTTTTCATTCACCAGGTTGGAAACGATCTTTTTTGTGGAATCTAGGAAGGGACATTTCTGAGCCCAATGAGGCCTAAAGCCAAAATCCAAATATCCCATGAGAAAAAGTAGAAACAAGTT";
#warn Dumper  $allele if $debug;
#die() if $debug;
				$hashRes{$structType}->{$id}->{annex}->{$patient_id}
				  ->{'var_allele'} = $allele->{sequence};
				$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{dp} =
				  $allele->{dp}->{raw};
				$hashRes{$structType}->{$id}->{annex}->{$patient_id}
				  ->{nb_all_mut} = $allele->{dp}->{alt};
				$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{nb_all_ref} = $allele->{dp}->{ref};
				$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{pr} = $allele->{pr}	if ( exists $allele->{pr} and $allele->{pr} );

				$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{sr} =$allele->{sr} if ( exists $allele->{sr} and $allele->{sr} );

				$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{set} = $allele->{vcf_parse}->{INFO}->{set} if exists $allele->{vcf_parse}->{INFO}->{set};

#$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{info_vcf} = $allele->{vcf_parse}->{INFO};#compress(freeze($allele->{vcf_parse}->{INFO}));

#		$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{complex_id} = $x->{COMPLEX};

				$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{score} =
				  $allele->{score};
				$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{he} =
				  $allele->{he};
				$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{ho} =
				  $allele->{ho};
				$hashRes{$structType}->{$id}->{annex}->{$patient_id}
				  ->{is_cas_1_2} = 1
				  if ( exists $allele->{is_cas_1_2}
					and $allele->{is_cas_1_2} == 1 );

				#delete $hashRes{$structType}->{$id}->{annex};
				#	warn $allele->{method} if exists $allele->{method};
				$hashRes{$structType}->{$id}->{annex}->{$patient_id}->{method}
				  = $allele->{method}
				  if exists $allele->{method};
				$hashRes{$structType}->{$id}->{annex}->{$patient_id}
				  ->{method_calling}->{ $self->method }->{nb_all_other_mut} = 0;
				if ( exists $allele->{is_cas_1_2}
					and $allele->{is_cas_1_2} == 1 )
				{
					$allele->{dp}->{alt} = 0 unless $allele->{dp}->{alt};
					$hashRes{$structType}->{$id}->{annex}->{$patient_id}
					  ->{method_calling}->{ $self->method }->{nb_all_other_mut}
					  = $allele->{dp}->{raw};
					$hashRes{$structType}->{$id}->{annex}->{$patient_id}
					  ->{method_calling}->{ $self->method }->{nb_all_other_mut}
					  -= $allele->{dp}->{alt};
					$hashRes{$structType}->{$id}->{annex}->{$patient_id}
					  ->{method_calling}->{ $self->method }->{nb_all_other_mut}
					  -= $allele->{dp}->{ref};
					$hashRes{$structType}->{$id}->{annex}->{$patient_id}
					  ->{nb_all_other_mut} =
					  $hashRes{$structType}->{$id}->{annex}->{$patient_id}
					  ->{method_calling}->{ $self->method }->{nb_all_other_mut};
				}
				$hashRes{$structType}->{$id}->{annex}->{$patient_id}
				  ->{method_calling}->{ $self->method }->{nb_all_ref} =
				  $allele->{dp}->{ref};
				$hashRes{$structType}->{$id}->{annex}->{$patient_id}
				  ->{method_calling}->{ $self->method }->{score} =
				  $allele->{score};
				$hashRes{$structType}->{$id}->{annex}->{$patient_id}
				  ->{method_calling}->{ $self->method }->{nb_all_mut} =
				  $allele->{dp}->{alt};
				$hashRes{$structType}->{$id}->{annex}->{$patient_id}
				  ->{method_calling}->{ $self->method }->{he} = $allele->{he};
				$hashRes{$structType}->{$id}->{annex}->{$patient_id}
				  ->{method_calling}->{ $self->method }->{ho} = $allele->{ho};

			}
			else {
				$hashRes{$structType}->{$id}->{annex} = undef;

				#$hashRes{$structType}->{$id}->{line_infos} = undef;
			}
			$allele->{dp}->{alt} = "." unless exists $allele->{dp}->{alt};
			$allele->{dp}->{ref} = "." unless exists $allele->{dp}->{ref};

#					if ($self->parse_hgmd() or $reference->parse_hgmd()) {
#						$hashRes{$structType}->{$id}->{hgmd_infos} = $allele->{vcf_parse};
#						$hashRes{$structType}->{$id}->{hgmd_id} = $allele->{vcf_parse}->{ID};
#						$hashRes{$structType}->{$id}->{isDM} = undef;
#						if ($allele->{vcf_parse}->{INFO}->{CLASS} eq 'DM') { $hashRes{$structType}->{$id}->{isDM} = 1; }
#					}

			$hashRes{$structType}->{$id} =
			  compress( freeze( $hashRes{$structType}->{$id} ) );

#					my $varObj = $project->flushObject($allele->{obj}, $hashRes{$structType}->{$id});
#					$allele = undef;
#					delete $hashRes{$structType}->{$id};

			#die()  if $allele->{chromosome} eq "chr17";

		}
	}

	$self->buffer->closeNeedlemanWunsch();

	return \%hashRes;

}

sub ucsc2ensembl {
	my ($chr) = @_;
	if ( $chr =~ /chr/ ) {
		$chr =~ s/chr//;
		$chr = 'MT' if $chr eq 'M';
	}
	return $chr;
}

sub create_cigar_sw2 {
	my ( $self, $ref1, $alt1, $needle, $debug ) = @_;
	my @opt = ( 'maxhits' => 1 );

	#push(@opt,'nogaps'=>'1') if (length($ref1) eq length ($alt1));
	my $gap  = -3;
	my $type = 0;
	if ( length($alt1) < length($ref1) * 0.75 ) {
		$needle = "needle";
		$gap    = -2;
		$type   = 1;
	}
	my $sw;

	#my $sw2;
	unless ($needle) {
		@opt = ( 'maxhits' => 1, 'gapopen' => -2 );
		$sw  = $self->buffer->getNeedlemanWunsch( \@opt, 2 );

		#$sw2 = new SmithWaterman($self->buffer,@opt);#,'gapopen'=>-3
	}
	else {
		@opt = ( 'gapopen' => $gap );
		#
		$sw = $self->buffer->getNeedlemanWunsch( \@opt, $type );

		#$sw2 = new NeedlemanWunsch($self->buffer,@opt);
	}
	warn "gap " if $debug;
	warn Dumper $sw if $debug;

	#my $sw = new NeedlemanWunsch($self->buffer,@opt);#,'gapopen'=>-3
	eval { $sw->do_alignment( $ref1, $alt1 ); };
	if ($@) {

		$self->buffer->closeNeedlemanWunsch();
		$sw = $self->buffer->getNeedlemanWunsch( \@opt, $type );
		$sw->do_alignment( $ref1, $alt1 );
	}

	#$sw2->do_alignment($ref1,$alt1);
	#warn Dumper $sw if $debug;
	my $string;
	my $hit;
	my $startN = 0;
	my $endN   = 0;
	if ( defined( $hit = $sw->get_next_hit() ) ) {

		#my $hit2 = $sw2->get_next_hit();
		#	die() if $hit2->{'align1'} ne $hit->{'align1'};
		if ($needle) {
			$hit->{pos1} = 0;
			$hit->{pos2} = 0;
		}
		if ( $hit->{pos1} > $hit->{pos2} ) {
			die() if $hit->{pos2} ne 10;
			my $string = substr( $alt1, 0, 11 );
			$string .= "-" x ( $hit->{pos1} - 10 );
			$hit->{'align2'} = $string . substr( $hit->{'align2'}, 1 );
			$string .= $hit->{'align2'};
			my $string2 = substr( $ref1, 0, $hit->{pos1} );
			$hit->{'align1'} = $string2 . $hit->{'align1'};

			#	$hit->{'align2'}
			#	warn Dumper $hit;
			#warn "coucou";
			#die();
			#$sw = new NeedlemanWunsch($self->buffer);
			#$sw->do_alignment($ref1,$alt1);
			#$hit = $sw->get_next_hit();
			#warn Dumper $hit;
			#	die();

		}
		if ( $hit->{pos2} > $hit->{pos1} ) {

			confess( Dumper($hit) ) if $needle;
			return $self->create_cigar_sw( $ref1, $alt1, "needle" );
		}
		my @s1 = split( "", $hit->{'align1'} );
		my @s2 = split( "", $hit->{'align2'} );
		my $M  = 0;
		my $I  = 0;
		my $D  = 0;

		for ( my $i = 0 ; $i < @s1 ; $i++ ) {
			if ( $startN == 0 && $s1[$i] ne 'N' ) {
				$startN = $i;
			}
			if ( $startN > 0 && $s1[$i] eq 'N' && $endN == 0 ) {
				$endN = $i;
			}
			if ( $s1[$i] eq $s2[$i] ) {
				$string .= "M";
			}
			elsif ( $s1[$i] eq "-" ) {
				$string .= "I";
			}
			elsif ( $s2[$i] eq "-" ) {
				$string .= "D";
			}
			else {
				$string .= "X";
			}

		}
	}
	else {
		confess();
	}

	my $l            = ( abs( $startN - $endN ) );
	my $tt           = substr( $hit->{'align1'}, $startN, $l );
	my $tt2          = substr( $hit->{'align2'}, $startN, $l );
	my $cigar_string = substr( $string, $startN, $l );

	my @res = split( "", $cigar_string );

	#warn Dumper @res;
	my $cigar        = "";
	my $current_code = $res[0];
	my $code_length  = 1;
	for ( my $i = 1 ; $i < @res ; $i++ ) {
		if ( $res[$i] eq $current_code ) {
			$code_length++;
		}
		else {

			$cigar .= $code_length . $current_code;
			$current_code = $res[$i];
			$code_length  = 1;
		}
	}
	$cigar .= $code_length . $current_code;

	warn $cigar if $debug;

	#	$sw2->destructor();
	#	$sw2 = undef;
	warn $cigar if $debug;

	#die() if $debug;
	return $cigar;
}

sub getNeedleAlign {
	my ( $self, $ref, $alt ) = @_;
	unless ( exists $self->{nw} ) {
		my $p = {
			match      => 5,
			mismatch   => -2,
			gap_open   => -1,
			gap_extend => -1
		};

		$self->{nw} = Align::NW->new( $ref, $alt, $p );
	}
	my $gap = -3;
	if ( length($alt) < length($ref) * 0.75 ) {

		$gap = -2;
	}

	return $self->{nw}->align( $ref, $alt, $gap );
}

sub create_cigar_sw {
	my ( $self, $ref1, $alt1, $needle, $debug ) = @_;

	#warn $ref1."\n".$alt1 if length($ref1) ne length($alt1);;

	#warn Dumper $self->getNeedleAlign($ref1,$alt1);
	#die();

	#push(@opt,'nogaps'=>'1') if (length($ref1) eq length ($alt1));
	return $self->create_cigar_sw2( $ref1, $alt1, $needle, $debug );
	if ( length($alt1) < length($ref1) * 0.75 ) {
		return $self->create_cigar_sw2( $ref1, $alt1, $needle, $debug );

	}
	my $hit    = $self->getNeedleAlign( $ref1, $alt1 );
	my $startN = 0;
	my $endN   = 0;

	my @s1 = split( "", $hit->{'align1'} );
	my @s2 = split( "", $hit->{'align2'} );

	#warn Dumper $hit unless @s2;
	die() unless @s2;
	my $M = 0;
	my $I = 0;
	my $D = 0;
	my $string;
	for ( my $i = 0 ; $i < @s1 ; $i++ ) {
		if ( $startN == 0 && $s1[$i] ne 'N' ) {
			$startN = $i;
		}
		if ( $startN > 0 && $s1[$i] eq 'N' && $endN == 0 ) {
			$endN = $i;
		}
		if ( $s1[$i] eq $s2[$i] ) {
			$string .= "M";
		}
		elsif ( $s1[$i] eq "-" ) {
			$string .= "I";
		}
		elsif ( $s2[$i] eq "-" ) {
			$string .= "D";
		}
		else {
			$string .= "X";
		}

	}

	my $l            = ( abs( $startN - $endN ) );
	my $tt           = substr( $hit->{'align1'}, $startN, $l );
	my $tt2          = substr( $hit->{'align2'}, $startN, $l );
	my $cigar_string = substr( $string, $startN, $l );

	my @res          = split( "", $cigar_string );
	my $cigar        = "";
	my $current_code = $res[0];
	my $code_length  = 1;
	for ( my $i = 1 ; $i < @res ; $i++ ) {
		if ( $res[$i] eq $current_code ) {
			$code_length++;
		}
		else {

			$cigar .= $code_length . $current_code;
			$current_code = $res[$i];
			$code_length  = 1;
		}
	}
	$cigar .= $code_length . $current_code;

	# my $c2 = $self->create_cigar_sw2 ($ref1,$alt1,$needle,$debug);
	#warn Dumper $hit if $c2 ne $cigar;
	#$self->create_cigar_sw2 ($ref1,$alt1,$needle,1) if $c2 ne $cigar;
	#warn $c2  if $c2 ne $cigar;
	#die() if $c2 ne $cigar;
	return $cigar;
}

sub parse_cigar {
	my ( $allele, $cigar, $debug ) = @_;
	my $alt   = $allele->{alt};
	my $ref   = $allele->{ref};
	my $start = $allele->{start};

	#	warn $cigar;
	my $res     = parse_cigar_line($cigar);
	my $num_alt = 0;
	my $pos_ref = 0;
	my $positions;
	my $pos_alt = 0;
	my $len_ref = 0;
	my $len_alt = 0;
	for ( my $i = 0 ; $i < @$res ; $i++ ) {
		my $c = $res->[$i];
		$c->{ref_start} = $pos_ref;
		$c->{alt_start} = $pos_alt;
		if ( $c->{code} ne 'I' ) {
			$pos_ref += $c->{len};
			$len_ref += $c->{len};
		}
		if ( $c->{code} ne 'D' ) {
			$pos_alt += $c->{len};
			$len_alt += $c->{len};
		}

		$c->{alt_end} = $pos_alt;
		$c->{ref_end} = $pos_ref;

		if ( $c->{code} ne 'M' ) {

			push( @$positions, $i );
			$num_alt++;
		}

		#$num_alt ++ if $c->{code} ne 'M';
	}

	warn " $len_ref => $ref  $alt $cigar" if $len_ref ne length($ref);
	confess( Dumper $allele)              if $len_ref ne length($ref);
	
	die()                                 if $len_ref ne length($ref);
	confess( Dumper $allele)              if $len_ref ne length($ref);

	if ( $len_alt ne length($alt) ) {
		confess();

	}

	confess( Dumper $allele) if $len_alt ne length($alt);

	my $cut_event;
	foreach my $i (@$positions) {
		my $c = $res->[$i];
		my $ref1;
		my $alt1;
		my $start1;
		my $end1;
		my $type1;
		my $cigar1;
		my $h;

		if ( $c->{code} eq 'D' ) {

			#my $s = $c->{ref_start}-1;
			$alt1 = substr( $ref, $c->{ref_start} - 1, 1 );
			$ref1 = substr( $ref, $c->{ref_start} - 1, $c->{len} + 1 );
			$start1 = $start + $c->{ref_start} - 1;
			$type1  = "deletion";
			$cigar1 = $c->{len} . $c->{code};

			#warn "deletion $ref/$alt ". $ref1."/".$alt1;

		}
		elsif ( $c->{code} eq 'X' ) {
			$ref1 = substr( $ref, $c->{ref_start}, $c->{len} );
			$alt1 = substr( $alt, $c->{alt_start}, $c->{len} );
			$start1 = $start + $c->{ref_start};
			$type1  = "snp";
			$type1  = "mnp" if ( $c->{len} > 1 );
			$cigar1 = $c->{len} . $c->{code};

			#	warn "snp $ref/$alt ". $ref1."/".$alt1;
		}
		elsif ( $c->{code} eq 'I' ) {
			$ref1 = substr( $ref, $c->{ref_start} - 1, 1 );
			$alt1 = substr( $alt, $c->{alt_start} - 1, $c->{len} + 1 );
			$start1 = $start + $c->{ref_start} - 1;
			$type1  = "insertion";
			$cigar1 = $c->{len} . $c->{code};

			#	warn "insertion $ref/$alt ". $ref1."/".$alt1;

		}
		else {
			confess( $c->{code} );
		}

		$h->{ref}   = $ref1;
		$h->{alt}   = $alt1;
		$h->{start} = $start1;
		$h->{type}  = $type1;
		$h->{cigar} = $cigar1;
		push( @{$cut_event}, $h );
	}

	return $cut_event;

	#exit();
	#die();

}

sub parse_cigar_line {
	my ($cigar) = @_;
	my $res;
	while ( $cigar =~ m/(\d+)([A-Z])/g ) {

		my $c;
		$c->{len}  = $1;
		$c->{code} = $2;
		push( @$res, $c );

	}
	return $res;

}

sub parseVcfFileForReference_samtools {
	my ( $self, $reference, $useFilter ) = @_;
	my $file   = $self->file();
	my $chrRef = $reference->getChromosome()->name();
	my $chr    = $reference->getChromosome();
	my $idchr  = $chr->name();
	$idchr = $chr->ucsc_name() if $self->isUcsc();
	my $tabix  = $self->buffer->config->{software}->{tabix};
	my $region = $idchr . ":" . $reference->start . "-" . $reference->end;
	my $cmd    = "$tabix $file $region | ";
	my %hashRes;
	my %skip;
	warn $cmd;
	open( FILE, $cmd ) || confess("can t open $file");
	my %verif_type;
	my $nb        = 0;
	my $line_skip = 0;

	while (<FILE>) {
		next if $_ =~ /^#/;
		chomp();
		my $line_file = $_;
		warn $line_file;
		my @data = split( " ", $line_file );

		my $var;

		my @test = split( ",", $data[4] );
		$data[4] = $test[0];

		my ( $structType, $structTypeObj );
		unless ( $line_file =~ /INDEL/ ) {
			warn "var";
			$var = $self->variationVcfFile( \@data );
			warn Dumper $var;

			#die();
		}
		else {
			if ( length( $data[3] ) > length( $data[4] ) ) {
				$structType    = 'del';
				$structTypeObj = 'deletions';
				my $seqshort = $data[4];
				my $seqlong  = $data[3];
				my $index    = index( $seqlong, $seqshort );
				if ( $index == -1 ) { $line_skip++; }
				$var = $self->deletionVcfFile( \@data );
			}
			elsif ( length( $data[3] ) == length( $data[4] ) ) {
				my $seqshort = $data[4];
				my $seqlong  = $data[3];
				warn $seqlong . " " . $seqshort;
				die("\n\nERROR: can't have equal length...\n\n");
			}
			else {
				$structType    = 'ins';
				$structTypeObj = 'insertions';
				my $seqshort = $data[3];
				my $seqlong  = $data[4];
				my $index    = index( $seqlong, $seqshort );
				if ( $index == -1 ) {
					$line_skip++;
					next;
				}
				$var = $self->insertionVcfFile( \@data );
			}
		}
		my $pat = $self->getPatient->id;

		#my $id = $chr->name."_".$data[1]."_".$var->{sequence_id};
		my $id = $chr->name . "_" . $var->{start} . "_" . $var->{sequence_id};
		warn $structTypeObj;
		$structType = $var->{structType};
		$hashRes{$structType}->{$id}->{'id'} = $id;
		$hashRes{$structType}->{$id}->{'vcf_id'} =
		  $data[0] . "_" . $data[1] . "_" . $data[3] . "_" . $data[4];
		$hashRes{$structType}->{$id}->{'structuralType'} = $var->{structType};
		$hashRes{$structType}->{$id}->{'structuralTypeObject'} =
		  $var->{structTypeObj};
		$hashRes{$structType}->{$id}->{'chromosomes_object'} =
		  { $chr->id() => undef };
		$hashRes{$structType}->{$id}->{'start'}      = $var->{start};
		$hashRes{$structType}->{$id}->{'end'}        = $var->{end};
		$hashRes{$structType}->{$id}->{'ref_allele'} = $var->{ref_allele};
		$hashRes{$structType}->{$id}->{'var_allele'} = $var->{sequence};
		$hashRes{$structType}->{$id}->{'line_infos'}->{$pat} = $line_file;
		$hashRes{$structType}->{$id}->{'references_object'}
		  ->{ $reference->id() } = undef;
		$hashRes{$structType}->{$id}->{'annex'}->{$pat}->{'ref_allele'} =
		  $var->{ref_allele};
		$hashRes{$structType}->{$id}->{'annex'}->{$pat}->{'var_allele'} =
		  $var->{mut_allele};
		$hashRes{$structType}->{$id}->{annex}->{$pat}->{he} =
		  $var->{heterozygote};
		$hashRes{$structType}->{$id}->{annex}->{$pat}->{ho} =
		  $var->{homozygote};
		$hashRes{$structType}->{$id}->{annex}->{$pat}->{score} = $var->{score};

		$hashRes{$structType}->{$id}->{annex}->{$pat}->{dp} = $var->{score2};
		$hashRes{$structType}->{$id}->{annex}->{$pat}->{nb_all_ref} =
		  $var->{score3};
		$hashRes{$structType}->{$id}->{annex}->{$pat}->{nb_all_mut} =
		  $var->{score4};
		$hashRes{$structType}->{$id} =
		  compress( freeze( $hashRes{$structType}->{$id} ) );
		$nb++;

	}

	#warn "Chromosome ".$chrRef." -> Nb var: $nb | Skip : $line_skip\n";
	close FILE;
	return \%hashRes;
}

sub variationVcfFile {
	my ( $self, $data ) = @_;
	my %var;
	$var{structType}    = 'snp';
	$var{structTypeObj} = 'variations';
	$var{start}         = $data->[1];
	$var{end}           = $data->[1];
	$var{position}      = $data->[1];

	die() unless $var{position};
	$var{type}        = "variation";
	$var{ref_allele}  = uc( $data->[3] );
	$var{mut_allele}  = uc( $data->[4] );
	$var{sequence}    = uc( $data->[4] );
	$var{sequence_id} = $var{ref_allele} . "_" . $var{sequence};
	$var{score}       = $var{score2} = $var{score3} = 0;
	getCommonData( $data, \%var );

	if ( $var{heterozygote} == 1 ) {
		$var{mut_allele} = $self->getPatient()->getProject()
		  ->biotools->getIUPAC( [ $var{ref_allele}, $var{mut_allele} ] );
	}
	return ( \%var );
}

sub insertionVcfFile {
	my ( $self, $data ) = @_;
	my %var;
	my $length_insertion = length( $data->[4] ) - length( $data->[3] );
	confess("\n\nERROR: length deletion < 0...\n\n") if $length_insertion < 0;
	my $sequence = uc( substr( $data->[4], 1, $length_insertion ) );
	my $start    = $data->[1] + 1;
	my $end      = $start + length($sequence);
	$var{start}      = $start;
	$var{end}        = $start;
	$var{mut_allele} = $sequence;
	$var{sequence}   = $sequence;
	$var{ref_allele} = "-";
	my $first = substr( $data->[3], 0, 1 );
	$var{sequence_id}   = $first . "_" . $first . $var{mut_allele};
	$var{structType}    = 'ins';
	$var{structTypeObj} = 'insertions';
	getCommonData( $data, \%var );
	return \%var;
}

sub deletionVcfFile {
	my ( $self, $data ) = @_;
	my %var;
	my $length_deletion = length( $data->[3] ) - length( $data->[4] );
	confess("\n\nERROR: length deletion < 0...\n\n") if $length_deletion < 0;
	my $sequence = uc( substr( $data->[3], 1, $length_deletion ) );
	my $start    = $data->[1] + length( $data->[4] );
	my $end      = $data->[1] + length($sequence);
	$var{old_start} = $start;
	$var{old_end}   = $end;

	if ( $start > $end ) {
		$var{old_start} = $end;
		$var{old_end}   = $start;
	}
	$var{start}      = $data->[1] + 1;
	$var{end}        = $var{start} + length($sequence) - 1;
	$var{mut_allele} = "-";                                   #join(@seq);
	$var{sequence}   = "-";
	$var{ref_allele} = $sequence;
	my $first = substr( $data->[3], 0, 1 );

	$var{sequence_id}   = $first . $var{ref_allele} . "_" . $first;
	$var{score}         = $var{score2} = $var{score3} = 0;
	$var{structType}    = 'del';
	$var{structTypeObj} = 'deletions';
	getCommonData( $data, \%var );
	return \%var;
}

sub getCommonData {
	my ( $data, $var ) = @_;
	$var->{score}  += parse_score($data);
	$var->{score2} += parse_dp($data);
	my $homo = parse_homozygote($data);
	if ( $homo == 1 ) {
		$var->{heterozygote} = 0;
		$var->{homozygote}   = 1;
	}
	else {
		$var->{heterozygote} = 1;
		$var->{homozygote}   = 0;
	}
	( $var->{score3}, $var->{score4} ) = parse_dp_by_allele($data);
	warn $data if $var->{score2} == 0;
	die() if $var->{score2} == 0;
}

sub parse_score {
	my ($data) = @_;
	return $data->[5];
}

sub parse_dp {
	my ($data) = @_;
	my $string = $data->[7];
	my %info;
	foreach my $tt ( split( ";", $string ) ) {
		my ( $scn, $val ) = split( "=", $tt );
		$info{$scn} = $val;
	}
	return $info{DP};
}

sub parse_dp_by_allele {
	my ($data) = @_;
	my $string = $data->[7];
	my %info;
	foreach my $tt ( split( ";", $string ) ) {
		my ( $scn, $val ) = split( "=", $tt );
		$info{$scn} = $val;
	}
	my @dp = split( ",", $info{DP4} );
	return ( ( $dp[0] + $dp[1] ), ( $dp[2] + $dp[3] ) );
}

sub parse_homozygote {
	my ($data) = @_;
	my @ids  = split( ":", $data->[8] );
	my @vals = split( ":", $data->[9] );
	for ( my $i = 0 ; $i < @ids ; $i++ ) {
		if ( $ids[$i] eq "GT" ) {
			if    ( $vals[$i] eq "1/1" ) { return 1; }
			elsif ( $vals[$i] eq "0/1" ) { return -1; }
			elsif ( $vals[$i] eq "1/0" ) { return -1; }
			elsif ( $vals[$i] eq "0/0" ) {
				warn Dumper $data;
				confess(
					"\n\nERROR: parse_homozygote error -> '0/0' case...\n\n");
			}
			else { confess( "n\nERROR: he/ho : " . Dumper $data . "\n\n" ); }
		}
	}
	confess( join( ";", @$data ) );
}

sub parseCasavaFile {
	my ( $self, $reference, $useFilter ) = @_;
	confess();
	my $file   = $self->file();
	my $chr    = $reference->getChromosome();
	my $chrRef = $reference->getChromosome()->name();
	confess("file not indexed $file") unless -e $file . ".tbi";

	#my $tabix = $self->buffer->config->{software}->{tabix};
	my $idchr = $chr->name();
	$idchr = $chr->ucsc_name() if $self->isUcsc();

	my $tabix = new Tabix( -data => $file );

#my ($idchr) = grep {$_ eq $chr->ucsc_name || $_ eq $chr->name } $tabix->getnames();
	my ($find) = grep { $_ eq $idchr } $tabix->getnames();
	warn "can't find $idchr $file " unless $find;
	return {} unless $find;
	my $region = $idchr . ":" . $reference->start . "-" . $reference->end;
	my $res = $tabix->query_full( $idchr, $reference->start, $reference->end );
	my %verif_type;
	my $nb;
	my $line_skip;
	my $nb_variation;
	my %hashRes;

	while ( my $line_file = $res->next ) {
		next if $line_file =~ /^#/;
		chomp($line_file);
		my @data = split( " ", $line_file );

		my $pat = $self->getPatient->id;
		my $pos = $data[1];

		if ( $line_file =~ /snp_/i ) {

			my $ref = uc( $data[10] );

			my @alleles = split( "", uc( $data[6] ) );

			my $homo = 0;
			my $he   = 1;
			my $mut_all;
			my $mut_pat = $REV_IUB{ join( "", sort @alleles ) };
			if ( scalar(@alleles) == 1 ) {
				$homo    = 1;
				$he      = 0;
				$mut_all = uc( $data[6] );
			}
			else {
				$mut_all = uc( $data[6] );
				$mut_all =~ s/$ref//;

				if ( length($mut_all) > 1 ) {
					$mut_all = chop($mut_all);
				}

			}
			die() if scalar(@alleles) > 2;
			my $id = $chr->name . "_" . $pos . "_" . $ref . "_" . $mut_all;

			my $structType    = 'snp';
			my $structTypeObj = 'variations';

			$hashRes{$structType}->{$id}->{'id'}     = $id;
			$hashRes{$structType}->{$id}->{'vcf_id'} = $id;
			$hashRes{$structType}->{$id}->{'chromosomes_object'} =
			  { $chr->id() => undef };
			$hashRes{$structType}->{$id}->{'structuralType'} = $structType;
			$hashRes{$structType}->{$id}->{'structuralTypeObject'} =
			  $structTypeObj;

			$hashRes{$structType}->{$id}->{'start'} = $pos;
			$hashRes{$structType}->{$id}->{'end'}   = $pos;

			$hashRes{$structType}->{$id}->{'ref_allele'}         = $ref;
			$hashRes{$structType}->{$id}->{'var_allele'}         = $mut_all;
			$hashRes{$structType}->{$id}->{'line_infos'}->{$pat} = $line_file;
			$hashRes{$structType}->{$id}->{'references_object'}
			  ->{ $reference->id() } = undef;
			$hashRes{$structType}->{$id}->{'annex'}->{$pat}->{'ref_allele'} =
			  $ref;
			$hashRes{$structType}->{$id}->{'annex'}->{$pat}->{'var_allele'} =
			  $mut_pat;
			$hashRes{$structType}->{$id}->{annex}->{$pat}->{he} = $he;
			$hashRes{$structType}->{$id}->{annex}->{$pat}->{ho} = $homo;
			my @sc = split( ":", $data[9] );

			my %DP = (
				A => 2,
				C => 3,
				G => 4,
				T => 5
			);
			$hashRes{$structType}->{$id}->{annex}->{$pat}->{score} = $sc[-1];
			$hashRes{$structType}->{$id}->{annex}->{$pat}->{dp}    = $data[8];
			$hashRes{$structType}->{$id}->{annex}->{$pat}->{nb_all_ref} =
			  $data[ $DP{$ref} ];
			$hashRes{$structType}->{$id}->{annex}->{$pat}->{nb_all_mut} =
			  $data[ $DP{$mut_all} ];

			$hashRes{$structType}->{$id} =
			  compress( freeze( $hashRes{$structType}->{$id} ) );

		}
		else {
			my $first_char = chop( $data[3] );
			my ( $ref, $mut_all ) = split( "/", $data[4] );

			my $homo = 0;
			my $he   = 1;
			if ( uc( $data[9] ) eq "HOM" ) {
				$homo = 1;
				$he   = 0;
			}
			my $score = $data[6];
			my $dp    = 0;
			$dp += $data[10];
			my $nb_ref = $data[11];
			my $nb_mut = $data[12];
			my $id;
			my $structType;
			my $structTypeObj;

			if ( $data[2] =~ /D/ && $data[2] !~ /I/ ) {

				#deletion
				$structType    = 'del';
				$structTypeObj = 'deletions';
				$ref           = $first_char . $ref;
				$mut_all       = $first_char;
				$id =
				  $chr->name . "_" . ( $pos - 1 ) . "_" . $ref . "_" . $mut_all;
			}
			elsif ( $data[2] =~ /I/ && $data[2] !~ /D/ ) {
				$structType    = 'ins';
				$structTypeObj = 'insertions';
				$mut_all       = $first_char . $mut_all;
				$ref           = $first_char;
				$id =
				  $chr->name . "_" . ( $pos - 1 ) . "_" . $ref . "_" . $mut_all;
			}

			else { next; }
			$hashRes{$structType}->{$id}->{'id'}     = $id;
			$hashRes{$structType}->{$id}->{'vcf_id'} = $id;
			$hashRes{$structType}->{$id}->{'chromosomes_object'} =
			  { $chr->id() => undef };
			$hashRes{$structType}->{$id}->{'structuralType'} = $structType;
			$hashRes{$structType}->{$id}->{'structuralTypeObject'} =
			  $structTypeObj;
			$hashRes{$structType}->{$id}->{'start'}              = $pos;
			$hashRes{$structType}->{$id}->{'end'}                = $pos;
			$hashRes{$structType}->{$id}->{'ref_allele'}         = $ref;
			$hashRes{$structType}->{$id}->{'var_allele'}         = $mut_all;
			$hashRes{$structType}->{$id}->{'line_infos'}->{$pat} = $line_file;
			$hashRes{$structType}->{$id}->{'references_object'}
			  ->{ $reference->id() } = undef;
			$hashRes{$structType}->{$id}->{'annex'}->{$pat}->{'ref_allele'} =
			  $ref;
			$hashRes{$structType}->{$id}->{'annex'}->{$pat}->{'var_allele'} =
			  $mut_all;
			$hashRes{$structType}->{$id}->{annex}->{$pat}->{he} = $he;
			$hashRes{$structType}->{$id}->{annex}->{$pat}->{ho} = $homo;
			$hashRes{$structType}->{$id}->{annex}->{$pat}->{method_calling}
			  ->{ $self->method }->{nb_all_ref} = $nb_ref;
			$hashRes{$structType}->{$id}->{annex}->{$pat}->{method_calling}
			  ->{ $self->method }->{nb_all_mut} = $nb_mut;
			$hashRes{$structType}->{$id}->{annex}->{$pat}->{score} = $score;
			$hashRes{$structType}->{$id}->{annex}->{$pat}->{dp}    = $dp;
			$hashRes{$structType}->{$id}->{annex}->{$pat}->{nb_all_ref} =
			  $nb_ref;
			$hashRes{$structType}->{$id}->{annex}->{$pat}->{nb_all_mut} =
			  $nb_mut;
			$hashRes{$structType}->{$id} =
			  compress( freeze( $hashRes{$structType}->{$id} ) );
		}

	}
	$tabix->close();

	#close FILE;
	return \%hashRes;

}

### TODO: Regarder ce cas et voir ce quil donne
### INS [chr13-49034022] REF: ATTC - VAR: ATTT

# CAS TORDU:

# chr13	49034022	rs10648129	ATTC	ATTTC,A	228.59	PASS	AC=3,1;AF=0.750,0.250;AN=4;DB;DP=13;FS=0.000;HaplotypeScore=106.6404;MLEAC=3,1;MLEAF=0.750,0.250;MQ=55.58;MQ0=0;QD=1.06;VQSLOD=1.95;culprit=HaplotypeScore	GT:AD:DP:GQ:PL	1/2:0,8,5:13:93:4392,418,93,442,0,156

# chr5	93388777	.	TATAG	T,TATAGATAG	1745.19	.	AC=1,1;AF=0.500,0.500;AN=2;DP=31;FS=0.000;HaplotypeScore=4.4324;MLEAC=1,1;MLEAF=0.500,0.500;MQ=62.13;MQ0=0;QD=49.86;RPA=10,9,11;RU=ATAG;STR	GT:AD:DP:GQ:PL	1/2:0,22,9:31:99:1767,576,880,1281,0,1961

1;
