package QueryDragenSv;
use strict;
use Moo;
use Data::Dumper;
use Carp;
use Compress::Snappy;
use Storable qw/thaw freeze/;
use Align::NW;
use Bio::DB::HTS::Tabix;
extends "QueryVcf";




sub parseVcfFileForReference {
	my ( $self, $reference, ) = @_;
	my $chr    = $reference->getChromosome();
	my $v      = Bio::DB::HTS::VCF->new( filename => $self->file() );
	my $v1     = Bio::DB::HTS::Tabix->new( filename => $self->file() );
	warn $self->file();
	my $header = $v->header();
	warn Dumper $header->get_sample_names();
	confess() if $header->num_samples() ne 1;
	
	#die( $header->get_sample_names->[0] . " " . $self->getPatient->barcode )
	 # if $self->getPatient->name ne $header->get_sample_names->[0]
	 # and $header->get_sample_names->[0] ne $self->getPatient->barcode;

	my $patient_id = $self->getPatient->id;
	my %hashRes;
	my $iter = $v1->query($chr->fasta_name . ":" . $reference->start . "-" . $reference->end );
	return {} unless $iter;
	while ( my $row = $iter->next ) {

		my $x = $self->parseVCFLine($row);
		confess( "alt " => Dumper $x) if scalar( @{ $x->{alt} } ) > 1;
		my $hash;
		if ( $x->{infos}->{SVTYPE} eq "INS")
		{
			$hash = $self->DragenIns( $x, $chr );
			next unless $hash;

			#my ($self,$x,$patient,$chr,$reference) = @_;
		}
		elsif ( $x->{infos}->{SVTYPE} eq "DUP" ) {
			$hash = $self->DragenDup( $x, $chr );
			next unless $hash;

			#my ($self,$x,$patient,$chr,$reference) = @_;
		}
		elsif ( $x->{infos}->{SVTYPE} eq "DEL" ) {
			$hash = $self->genericSVDel( $x, $chr );
			next unless $hash;

			#my ($self,$x,$patient,$chr,$reference) = @_;
		}
		elsif ( $x->{infos}->{SVTYPE} eq "INV" ) {
			$hash = $self->genericSVInv( $x, $chr );
			next unless $hash;

			#my ($self,$x,$patient,$chr,$reference) = @_;
		}
		elsif ( $x->{infos}->{SVTYPE} eq "BND" ) {
			next ;
			#my ($self,$x,$patient,$chr,$reference) = @_;
		}
		else {
			confess( $x->{infos}->{SVTYPE} );
			next;
		}
		my $id = $hash->{id};
		$hash->{'isSrPr'} = 1;
		my $structType = $hash->{structuralType};
		next unless $id;
		$hashRes{$structType}->{$id} = compress( freeze($hash) );

	}
	return \%hashRes;
}

sub DragenDup {
	my ( $self, $x, $reference ) = @_;
	my $chr        = $reference->getChromosome();
	my $patient    = $self->getPatient();
	my $patient_id = $patient->id;
	my $hash;
	my $ref       = $x->{ref};
	$ref        = substr( $ref, 0, 1 );
	my $pos       = $x->{pos};
	my $genbo_pos = $pos + length($ref);
	my $id;
	my $alt = $x->{alt}->[0];
	my $len  = $x->{infos}->{SVLEN};
	die() unless $len;
	my $var_allele;
	my $infos = $x->{infos};
	$len = $infos->{SVLEN};
	$len = $infos->{DUPSVLEN} if exists $infos->{DUPSVLEN};
 	my $res = $chr->genesIntervalTree->fetch( $genbo_pos, $genbo_pos + $len );
 	
		$id = $chr->name . "_" . $genbo_pos . "_" . $ref . "_dup-" . $len;
		$var_allele = $chr->sequence( $genbo_pos, $genbo_pos + $len );
		$hash->{'isDup'} = 1;

		#$hash->{'isLargeDuplication'} = 1;
		$hash->{'isLarge'}              = 1;
		$hash->{'end'}                  = $genbo_pos + $len;
		$hash->{'length'}               = abs($len);
		$hash->{'structuralType'}       = 'l_dup';
		$hash->{'structuralTypeObject'} = 'large_duplications';
		$hash->{'allele_length'}        = $len;
	
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
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{nb_all_other_mut} = 0;
	$hash->{annex}->{$patient_id}->{method} = $self->method();
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method } ->{nb_all_other_mut} = 0;
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{nb_all_ref} = $hash->{annex}->{$patient_id}->{nb_all_mut};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{score} = $x->{qual};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{nb_all_mut} = $hash->{annex}->{$patient_id}->{nb_all_ref};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{he} = $x->{gt}->{he};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{ho} = $x->{gt}->{ho};
	return $hash;
}

sub DragenIns {
	my ( $self, $x, $reference ) = @_;
	my $chr        = $reference->getChromosome();
	my $patient    = $self->getPatient();
	my $patient_id = $patient->id;
	my $hash;
	my $ref       = $x->{ref};
	$ref        = substr( $ref, 0, 1 );
	my $pos       = $x->{pos};
	my $genbo_pos = $pos + length($ref);
	my $id;
	my $alt = $x->{alt}->[0];
	my $len = 1;
	
	my $var_allele;
	my $infos = $x->{infos};
	$hash->{'structuralType'}       = 'ins';
	$hash->{'structuralTypeObject'} = 'insertions';
	if ($alt ne "<INS>"){
		$alt        = substr( $alt, length($ref) );
		$len        = length($alt);
		$id         = $chr->name . "_" . $genbo_pos . "_" . $ref . "_" . $alt;
		$var_allele = $alt;
		$hash->{'var_allele'}    = $var_allele;
		$hash->{'allele_length'} = length($var_allele);
	}
	elsif ( exists $infos->{DUPSVLEN} ){
		$var_allele = $chr->sequence( $genbo_pos, $genbo_pos + $infos->{DUPSVLEN});
		
		 if   ($infos->{DUPSVLEN} < $infos->{SVLEN} ) { 
				$var_allele .= $infos->{DUPSVINSSEQ};
			
	 	}
	 	elsif  ($infos->{DUPSVLEN} == $infos->{SVLEN} ) {
	 		$hash->{'structuralType'}       = 'l_dup';
			$hash->{'structuralTypeObject'} = 'large_duplications';
	 	}
	 	else {
	 		confess();
	 	}
	 	$id         = $chr->name . "_" . $genbo_pos . "_" . $ref . "_" . $alt;
		$hash->{'var_allele'}    = $var_allele;
		$hash->{'allele_length'} = $infos->{SVLEN};
		
	}
	elsif ( exists $infos->{SVINSSEQ} ) {
			$var_allele   = $infos->{SVINSSEQ};
			$alt = $infos->{SVINSSEQ};
			$id         = $chr->name . "_" . $genbo_pos . "_" . $ref . "_" . $alt;
			$hash->{'var_allele'}    = $var_allele;
			$hash->{'allele_length'} = length($var_allele);
	}
	elsif ( exists $infos->{LEFT_SVINSSEQ} or exists $infos->{RIGHT_SVINSSEQ} ){
		#	return ;
		$hash->{'structuralType'}       = 'l_ins';
		$hash->{'structuralTypeObject'} = 'large_insertions';
		$infos->{LEFT_SVINSSEQ} = "" unless $infos->{LEFT_SVINSSEQ};
		$infos->{RIGHT_SVINSSEQ} = ""  unless $infos->{RIGHT_SVINSSEQ};
		$var_allele = $infos->{LEFT_SVINSSEQ} . "+" . $infos->{RIGHT_SVINSSEQ};
		$id         = $chr->name . "_" . $genbo_pos . "_" . $ref . "_ins-?";
	
		$hash->{'var_allele'}    = $var_allele;
		$hash->{'allele_length'} = "~" . length($var_allele);
	}
	else {confess();}

	
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
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{nb_all_other_mut} = 0;
	$hash->{annex}->{$patient_id}->{method} = $self->method();
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method } ->{nb_all_other_mut} = 0;
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{nb_all_ref} = $hash->{annex}->{$patient_id}->{nb_all_mut};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{score} = $x->{qual};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{nb_all_mut} = $hash->{annex}->{$patient_id}->{nb_all_ref};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{he} = $x->{gt}->{he};
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{ho} = $x->{gt}->{ho};
	return $hash;
}
1;
