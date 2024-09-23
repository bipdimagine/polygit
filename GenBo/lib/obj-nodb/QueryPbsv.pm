package QueryPbsv;
use strict;
use Moo;
use Data::Dumper;
use Carp;
use Clone 'clone';
use IPC::Open2;
use List::MoreUtils qw( uniq );
use Compress::Snappy;
use Storable qw/thaw freeze/;
use Align::NW;
use Bio::DB::HTS::Tabix;
use Bio::DB::HTS::VCF;
extends "QueryVcf";

sub parseVcfFileForReference {
	my ( $self, $reference, ) = @_;
	my $chr    = $reference->getChromosome();
	my $v      = Bio::DB::HTS::VCF->new( filename => $self->file() );
	my $v1     = Bio::DB::HTS::Tabix->new( filename => $self->file() );
	my $header = $v->header();
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
		if ( $x->{infos}->{SVTYPE} eq "INS" or $x->{infos}->{SVTYPE} eq "DUP" )
		{
			$hash = $self->PbsvIns( $x, $chr );
			next unless $hash;

			#my ($self,$x,$patient,$chr,$reference) = @_;
		}
		elsif ( $x->{infos}->{SVTYPE} eq "DEL" ) {
			$hash = $self->PbsvDel( $x, $chr );
			next unless $hash;

			#my ($self,$x,$patient,$chr,$reference) = @_;
		}
		elsif ( $x->{infos}->{SVTYPE} eq "INV" ) {

			$hash = $self->PbsvINV( $x, $chr );
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
		warn Dumper $hash unless $id;
		next unless $id;
		$hashRes{$structType}->{$id} = compress( freeze($hash) );

	}

	return \%hashRes;
	
}

sub PbsvINV {
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
	my $pos_end = $genbo_pos + abs( $x->{infos}->{SVLEN});
	
	my $res = $chr->genesIntervalTree->fetch( $genbo_pos,$genbo_pos + abs( $x->{infos}->{SVLEN} ) );
	return if scalar(@$res) > 2;
	my $var_allele = $chr->sequence( $genbo_pos, $pos_end );
	$var_allele = BioTools::complement_sequence($var_allele);

	my $len   = abs( $pos_end - $genbo_pos ) + 1;
	my $id    = $chr->name . "_" . $genbo_pos . "_" . $ref . "_inv-" . $len;
	my $infos = $x->{infos};
	$hash->{'id'}                   = $id;
	$hash->{'isInversion'}          = 1;
	$hash->{'structuralType'}       = 'inv';
	$hash->{'structuralTypeObject'} = 'inversions';
	$hash->{'id'}                   = $id;
	$hash->{'vcf_id'} = join( "_", $chr->name,$x->{ref}, $alt );
	
	$hash->{'isSV'}               = 1;
	$hash->{'chromosomes_object'} = { $chr->id => undef };
	$hash->{'start'}              = $genbo_pos;
	$hash->{'end'}                = $pos_end;

	$hash->{'length'}     = $len;
	$hash->{'ref_allele'} = $ref;

	#my $var_allele = $chr->sequence($start,$end);
	$hash->{'var_allele'}   = $var_allele;
	$hash->{'line_infos'}   = "-";
	$hash->{'vcf_position'} = $pos;
	$hash->{'references_object'}->{ $reference->id } = undef;
	### ANNEX
	$hash->{annex}->{$patient_id}->{Filter} =$x->{filter};
	$hash->{annex}->{$patient_id}->{'ref_allele'} = $ref;
	$hash->{annex}->{$patient_id}->{'var_allele'} = "inv";

	$self->add_DP_AD( $patient_id, $hash, $x );
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

sub  PbsvDel {
	
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
	return if  $x->{infos}->{SVLEN} > 2000;
	
	my $res = $chr->genesIntervalTree->fetch( $genbo_pos,
		$genbo_pos + abs( $x->{infos}->{SVLEN} ) );
	return if scalar(@$res) > 2;
	my $len;
	my $var_allele = $alt;
	my $infos      = $x->{infos};
	$hash->{'id'}                   = $id;
	$hash->{'isSrPr'}               = 1;
	$hash->{'isLargeDeletion'}      = 1;
	$hash->{'structuralType'}       = 'del';
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
	$hash->{'isCnv'}              = 1;
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
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{nb_all_other_mut} = 0;

	$hash->{annex}->{$patient_id}->{method} = $self->method();
	$hash->{annex}->{$patient_id}->{method_calling}->{ $self->method }->{nb_all_other_mut} = 0;
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
sub PbsvIns {
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
	my $len  = $x->{infos}->{SVLEN};
	die() unless $len;
	
	my $var_allele;
	my $infos = $x->{infos};
	$hash->{'structuralType'}       = 'ins';
	$hash->{'structuralTypeObject'} = 'insertions';
	
	if ( $infos->{SVTYPE} eq "DUP" ) {
		my $res = $chr->genesIntervalTree->fetch( $genbo_pos, $genbo_pos + $len );
		return if @$res > 2;
		$id = $chr->name . "_" . $genbo_pos . "_" . $ref . "_dup-" . $len;
		$var_allele = $chr->sequence( $genbo_pos, $genbo_pos + $len );
		$hash->{'isDup'} = 1;
		$hash->{'isLarge'}              = 1;
		$hash->{'end'}                  = $genbo_pos + $len;
		$hash->{'length'}               = abs($len);
		$hash->{'structuralType'}       = 'l_dup';
		$hash->{'structuralTypeObject'} = 'large_duplications';
		$hash->{'allele_length'}        = $len;
	}
	elsif ( $infos->{SVTYPE} eq "INS")
	{
		$len        = 1;
		$alt        = substr( $alt, length($ref) );
		$len        = length($alt);
		$ref        = substr( $ref, 0, 1 );
		if ($alt ne "<INS>"){
			$var_allele        = substr( $alt, length($ref) );
			$len        = length($var_allele);
			$id         = $chr->name . "_" . $genbo_pos . "_" . $ref . "_" . $var_allele;
			$hash->{'allele_length'} = length($var_allele);
			
		}
		else {
			die Dumper $x;
		}
	}
		$hash->{'var_allele'}    = $var_allele;
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



sub add_DP_AD {
	my ( $self, $patient_id, $hash, $x ) = @_;
	my @ad = split( ",",  $x->{gt}->{AD} );
	
	$hash->{annex}->{$patient_id}->{nb_all_mut} = $ad[1];
	$hash->{annex}->{$patient_id}->{nb_all_ref} = $ad[0];
	$hash->{annex}->{$patient_id}->{dp} =  $x->{gt}->{DP};
}

1;