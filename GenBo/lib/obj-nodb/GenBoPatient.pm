package GenBoPatient;

use strict;
use Moo;
use Data::Dumper;
use Config::Std;
use GenBoCapture;
use GenBoPcrMultiplex;
use QueryVcf;
use QueryJunctionFile;
use File::Util;
use Bio::DB::Sam;
use Storable;
use validationQuery;
use Compress::Snappy;
use GenBoBinaryFile;
use Storable qw/thaw freeze/;
use List::Util qw(max sum);
use JSON::XS;
use Statistics::Descriptive;
use GenBoNoSqlLmdbCache;
use Carp;
use List::MoreUtils qw(firstidx );
extends "GenBo";

has patient_id => (
	is => 'ro',

	#isa		=> 'Str',
	reader   => 'getPatientId',
	required => 1,
);
has bar_code => (
	is     => 'ro',
	reader => 'barcode',

	#isa		=> 'Str',
);

has capture_id => (
	is => 'rw',

	#isa		=> 'Str',
	reader   => 'getCaptureId',
	required => 1,
);

has run_id => (
	is => 'ro',

	#isa		=> 'Str',
	reader   => 'getRunId',
	required => 1,
);

has genbo_id => (
	is => 'ro',

	#isa		=> 'Str',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->id();
	}
);

has project_id => (
	is => 'ro',

	#isa		=> 'Str',
	reader => 'getProjectId',
);

has control => (
	is     => 'ro',
	reader => 'is_control',
);

has identity_vigilance_vcf => ( is => 'ro', );

has identity_vigilance => ( is => 'ro', );

has captureFiles => (
	is => 'ro',

	#isa		=> 'ArrayRef[Str]',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my @lFiles;
		my $lCaptureObj = $self->getCaptures();
		foreach my $obj (@$lCaptureObj) {
			push( @lFiles, $obj->gzFileName() );
		}
		return \@lFiles;
	},
);

has transcriptsCoverageFile => (
	is => 'ro',

	#isa		=> 'ArrayRef[Str]',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $f =
			$self->getProject()->getCoverageDir()
		  . $self->name
		  . ".transcript.cov.kct";
		return $f;
	},
);

has kyoto_polydiag_cache => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $f =
			$self->project->getDiagCacheDir() . "/"
		  . $self->name
		  . ".polydiag.kct";
		return $f;
	},
);


has isGenome => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		foreach my $c ( @{ $self->getCaptures } ) {
			return 1 if lc( $c->analyse ) =~ /genome/;
			return 1 if lc( $c->analyse ) =~ /all_exons/;
		}
		return undef;
	},
);

sub transcriptsCoverageLite {
	my ( $self, $mode, $reopen ) = @_;
	my $project = $self->project;
	return $project->transcriptsCoverageLite($mode);

}

sub transcriptsCoverage {
	my ( $self, $mode, $reopen ) = @_;
	$mode = '' unless $mode;
	confess();

	#warn "KYOTO ----";
	if ($reopen) {
		warn " ****** REOPEN COVERAGE with mode $mode ***** ";
		$self->{transcriptsCoverage}->close()
		  if exists $self->{transcriptsCoverage};
		delete $self->{transcriptsCoverage};
	}
	return $self->{transcriptsCoverage} if exists $self->{transcriptsCoverage};
	my $f   = $self->transcriptsCoverageFile();
	my $db1 = new KyotoCabinet::DB;
	my $kmode;
	if ( $mode eq 'c' ) {
		if (
			!$db1->open(
				$f,
				$db1->OWRITER | $db1->ONOLOCK | $db1->OCREATE | $db1->OTRUNCATE
			)
		  )
		{

			printf STDERR ( "open error: %s %s\n", $db1->error, $f );
			confess();
		}
		system("chmod a+w $f");

	}
	elsif ( $mode eq 'w' ) {
		if ( !$db1->open( $f, $db1->OWRITER | $db1->ONOLOCK ) ) {

			printf STDERR ( "open error: %s %s\n", $db1->error, $f );
			confess();
		}

	}
	else {
		return unless -e $f;
		if ( !$db1->open( $f, $db1->OREADER | $db1->ONOLOCK ) ) {
			printf STDERR ( "open error: %s %s\n", $db1->error, $f );
			confess();
		}
	}

	$self->{transcriptsCoverage} = $db1;
	return $self->{transcriptsCoverage};
}

sub getChromosomeCoverage {
	my ( $self, $chr ) = @_;
	confess();
	return $self->{chr_coverage}->{ $chr->id }
	  if exists $self->{chr_coverage}->{ $chr->id };
	my $f =
		$self->getProject()->getCoverageDir() . "/"
	  . $self->name . "/"
	  . $chr->name
	  . ".freeze";
	unless ( -e $f ) {
		$self->{chr_coverage}->{ $chr->id } = undef;
		return;
	}
	my $array = retrieve($f);
	$self->{chr_coverage}->{ $chr->id } = $array;
	return $array;

}

has tabix_coverage => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $coverage_file;
		$coverage_file = $self->getCoverageFile();
		return unless -e $coverage_file;
		return Bio::DB::HTS::Tabix->new( filename => $coverage_file );

	}
);

#has bio_db_sam => (
#	is		=> 'ro',
#	lazy	=> 1,
#	default => sub {
#		my $self = shift;
#		my $hSam;
#		my $fasta_file = $self->project->getGenomeFasta();
#		foreach my $bam_file (@{$self->getBamFiles()}) {
#			$hSam->{$bam_file} = Bio::DB::Sam->new(-bam=>$bam_file, -fasta=>$fasta_file);
#		}
#		return $hSam;
#	},
#);

has bio_db_sam => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $bam  = $self->getBamFile();
		return unless -e $bam;
		my $sam = Bio::DB::Sam->new( -bam => $bam );
		return $sam;
	}
);

has is_ucsc_nomenclature => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $sam  = $self->bio_db_sam;
		my $find = grep { $_ =~ /^chr/ } $sam->seq_ids;
		return $find;
	}
);

#sub get_raw_coverage {
#	my ($self,$chr,$start,$end) = @_;
#	my $sam = $self->bio_db_sam();
#
#	my ($coverage) = $sam->features(-type=>'coverage',-seqid=>$chr,-start=>$start,-end=>$end);
#	my $array_intspan =  Array::IntSpan->new();
#	my $sum =0;
#	my $nb;
#	for (my $i=0;$i< @{$coverage->coverage};$i++){
#
#		#warn $coverage->coverage->[$i];
#			my $pos = $i+$coverage->start;
#			$array_intspan->set_range($pos,$pos,$coverage->coverage->[$i]);
#			$sum += $coverage->coverage->[$i];
#			$nb++;
#
#	}
#	 $array_intspan->consolidate(1,$coverage->end);
#	 my $mean =0;
#	 $mean = $sum/$nb if $sum >0;
#	 my $res;
#	 $res->{array_intspan} = $array_intspan;
#	 $res->{sum} = $sum;
#	 $res->{mean} = $mean;
#	  $res->{data} = $coverage->coverage;
#
#	return $res;
#}

has coverage => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;

		my $item;
		eval {
			my $tabix = $self->tabix_coverage;
			return unless $tabix;
			my $res   = $tabix->query_full( "mean_all" ) ;
			my @data;

			while ( my $line = $res->next ) {
				my ( $a, $b, $c ) = split( " ", $line );
				if ( $b == 99 ) {
					$b = "mean";

				}
				else { $b .= "x"; $c *= 100; }
				$item->{"$b"} = int( $c * 10 ) / 10;

			}
		};
		return $item;
	},
);

has sequencesDir => (
	is      => 'ro',
	lazy    => 1,
	reader  => 'getSequencesDirectory',
	default => sub {
		my $self        = shift;
		my $run         = $self->getRun();
		my $machine     = $run->machine;
		my $run_name    = $run->plateform_run_name();
		my $constructor = $run->machine_constructor();
		my $plateform   = $run->plateform();
		my $path        = $self->buffer()->getDataDirectory("sequences");

		my $seq_dir =
			$path . "/"
		  . $constructor . "/"
		  . $machine . "/"
		  . $plateform
		  . "/$run_name/";
		my $seq_dir2 =
			$path . "/"
		  . $constructor . "/"
		  . $machine . "/"
		  . $plateform
		  . "/$run_name.saved/";
		return $seq_dir2 if -e $seq_dir2;
		return $seq_dir;

	},
);

sub getCaptureFile {
	my $self = shift;
	return $self->captureFiles()->[0]
	  if scalar( @{ $self->captureFiles() } ) == 1;
	confess();
}

has captureBedFiles => (
	is => 'ro',

	#isa		=> 'ArrayRef[Str]',
	#	reader	=> 'getCaptureBedFiles',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my @lFiles;
		my $lCaptureObj = $self->getCaptures();
		foreach my $obj (@$lCaptureObj) {
			push( @lFiles, $obj->bedFileName() );
		}
		return \@lFiles;
	},
);

sub getCaptureBedFile {
	my $self = shift;
	return $self->captureBedFiles->[0]
	  if scalar( @{ $self->captureBedFiles } ) == 1;
	confess();
}

has captureIntSpan => (
	is => 'ro',

	#isa		=> 'Any',
	reader  => 'getCapturesIntSpan',
	lazy    => 1,
	default => sub {
		my $self       = shift;
		my $captureObj = $self->getGenBoCapture();
		return $captureObj->getIntSpan();
	},
);

has 'sex' => (
	is => 'rw',

	#isa 	=> 'Int',
	default => -1,
);

has 'coverage_SRY' => (
	is => 'rw',

	#isa 	=> 'Int',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my ($find) =
		  grep { $_->name eq 'Y' } @{ $self->project->getChromosomes };
		return -1 unless $find;
	
		my $sambamba = $self->buffer->software("sambamba");
		my $intspan_capture =
		  $self->project->getChromosome('Y')->getCapturesGenomicSpan();
		my $intSpan = Set::IntSpan::Fast::XS->new(
			"6736370-6736371,6736442-6736443,6737844-6737845,2654896-2655740");

		my $zint = $intspan_capture->intersection($intSpan);

		return -1 if $zint->is_empty;
		my $iter = $zint->iterate_runs();
		my @tt;
		my $chr_name = $self->project->getChromosome('Y')->fasta_name();
		my $bam      = $self->getBamFile();
		my $max_mean = -50;
		while ( my ( $from, $to ) = $iter->() ) {
			my ($res) =
`$sambamba depth region $bam -L $chr_name:$from-$to 2>/dev/null | grep $chr_name | cut -f 5`;
			$res = 0 unless $res;
			chomp($res);

			#$res = 0 unless $res;

			$max_mean = $res if $res > $max_mean;
			last if $max_mean > 20;
		}
		return int( $max_mean * 10 ) / 10;


	},
);

has 'compute_sex' => (
	is => 'rw',

	#isa 	=> 'Int',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $covm = $self->coverage_SRY();
		return 1 if $covm > 30;
		return -1 if $covm == -1;
		my $covh = $self->coverage();
		$covh->{mean} += 0;
		return 2 if $covm < 5;
		return 1 if $covm > $covh->{mean} * 0.10;
		return 2 if $covm < $covh->{mean} * 0.05;
		if ( $covh->{mean} > 300 ) {
			return 1 if $covm > $covh->{mean} * 0.8;
		}

		#return 1 if $covm > 70;
		return (2);

	},
);

has somatic_group => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return unless ( $self->project->isSomatic() );
		return unless $self->project->somatic_details();
		foreach my $group_name ( keys %{ $self->project->somatic_details() } ) {
			return $group_name
			  if (
				exists $self->project->somatic_details->{$group_name}
				->{ $self->name() } );
		}
		return;
	},
);

# tissue du patient (somatic model)
has tissue => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return unless ( $self->project->isSomatic() );

		return $self->project->somatic_details->{ $self->somatic_group() }
		  ->{ $self->name() }->{tissue};
	},
);

has isSomatic => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return 1 if ( $self->tissue() eq 'T' );
		return;
	}
);

has isGerminal => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return 1 if ( $self->tissue() eq 'C' );
		return;
	}
);

has status => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->project->pedigree_details->{ $self->family() }
		  ->{ $self->name }->{status};
	}
);
has vstatus => (
	is      => 'ro',
	
);
has isHealthy => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return 1 if ( $self->status == 1 );
		return;
	}
);

has isIll => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return 1 unless ( $self->status() );
		return 1 if ( $self->status == 2 );
		return;
	}
);

has pedigreeLine => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $line = $self->family() . "\t" . $self->name();
		if ( $self->isChild() ) {
			my $fam    = $self->getFamily();
			my $father = $self->name() . '_father';
			$father = $fam->father() if ( $fam->father() );
			my $mother = $self->name() . '_mother';
			$mother = $fam->mother() if ( $fam->mother() );
			$line .= "\t" . $father . "\t" . $mother;
		}
		else { $line .= "\t0\t0"; }
		$line .= "\t" . $self->sex() . "\t" . $self->status();
		return $line;
	},
);

has family => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		foreach my $fam_name ( keys %{ $self->project->somatic_details() } ) {
			return $fam_name
			  if (
				exists $self->project->pedigree_details->{$fam_name}
				->{ $self->name() } );
		}
		return;
	},
);

has isChild => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return 1 unless $self->project->isFamilial;
		return 1 unless ( $self->getFamily() );
		return 1 if ( exists $self->getFamily->children->{ $self->name() } );
		return 0;
	},
);

has isKid => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->isChild;

	},
);

has callingMethods => (
	is      => 'rw',
	reader  => 'getCallingMethods',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my @lMethods;
		my $query    = $self->getProject()->buffer->getQuery();
		my $lMethods = $query->getCallingMethods(
			patient_name => $self->name() . "",
			project_id   => $self->getProject->id
		);
		my %order_methods = (
			ion_merge        => 1,
			unifiedgenotyper => 2,
			mpileup          => 3,
			gatk             => 4,
			dibayes          => 5,
			lifinder         => 6,
		);
		foreach my $method_name (@$lMethods) {

			$order_methods{$method_name} = 999
			  unless exists $order_methods{$method_name};
		}

		foreach my $methName (@$lMethods) {
			if (   ( $self->getProject()->getVariationsDir($methName) )
				or ( $self->getProject()->getIndelsDir($methName) )
				or ( $self->getProject()->getLargeIndelsDir($methName) ) )
			{
				push( @lMethods, $methName );
			}
		}
		my @sort_methods =
		  sort { $order_methods{$a} <=> $order_methods{$b} } @lMethods;
		return \@sort_methods;
	},
);

has callingSVMethods => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my @lMethods;
		my $query    = $self->getProject()->buffer->getQuery();
		my $lMethods = $query->getCallingSVMethods(
			patient_name => $self->name() . "",
			project_id   => $self->getProject->id
		);
		my @final;
		foreach my $methName (@$lMethods) {
			if ( $self->getProject()->getVariationsDir($methName) ) {
				push( @final, $methName );
			}
		}

		return \@final;
	},
);

has alignmentMethods => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self  = shift;
		my $query = $self->getProject()->buffer->getQuery();
		return $query->getAlignmentMethods(
			patient_name => $self->name() . "",
			project_id   => $self->getProject->id
		);
	},

);

sub alignmentMethod {
	my $self    = shift;
	my $methods = $self->alignmentMethods();
	confess("more than one alignmeent method use default mapread")
	  if scalar(@$methods) > 1;
	return $methods->[0];
}



sub setRuns {
	my $self = shift;
	$self->getProject->getRuns();
	my $hids;
	$hids->{ $self->getRunId } = undef;
	return $hids;
}

sub setCaptures {
	my $self = shift;
	my $hashObjIds;
	my $args;
	$args->{id} = $self->getCaptureId();

	#$args->{name} = "capture_" . $self->getCaptureId();

	my $obj = $self->getProject()->flushObject( 'captures', $args );
	$hashObjIds->{ $self->getCaptureId() } = $self->getCaptureId();
	return $hashObjIds;
}

sub setPrimers {
	my ($self) = @_;

	#	warn "coucou";
	#	die();
	my %hash;
	foreach my $c ( $self->getCaptures ) {
		my $primers = $c->getPrimers();
		map { $hash{ $_->id }++ } @$primers;
	}

	return \%hash;

}

has deletions_object => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $hRes = $self->setDeletions();
		unless ($hRes) { $hRes->{none} = 'none'; }
		return $hRes;
	}
);

sub setVariants {
	my ( $self, $typeVar ) = @_;
	my @objs;
	foreach my $ref ( @{ $self->getProject()->getReferences() } ) {
		foreach my $o ( @{ $self->setVariantsForReference( $ref, $typeVar ) } )
		{
			$o->{references_object}->{ $self->id } = undef;
			$self->{ $o->type_object }->{ $o->id } = undef;
		}
	}
}

sub setValidatedVariants {
	my ($self) = @_;
	$self->{validated_variants} = {};

	return $self->{validated_variants} if exists $self->{validated_variants};
	my $vs = {};
	if ( $self->getCapture->validation_db() ) {

#my $vquery = validationQuery->new(dbh=>$self->buffer->dbh,capture_name=>$self->getCapture->validation_db());
		$vs =
		  $self->project->validations_query->get_variations_in_validation_table(
			sample_name  => $self->name,
			project_name => $self->getProject->name
		  );
	}
	my $tt;
	foreach my $v ( values %$vs ) {

		my ( $chr_name, $start, $ref, $alt ) = split( '_', $v->{vcfid} );
		$alt = "" unless $alt;
		my $hash;
		my $type;
		my $strucType;
		$hash->{id}         = $v->{polyid};
		$hash->{annex}      = undef;
		$hash->{line_infos} = "";
		$hash->{chr_name}   = $chr_name;
		$hash->{start}      = $start;

		$hash->{end}               = $hash->{start};
		$hash->{strand}            = 1;
		$hash->{vcf_id}            = $v->{vcfid};
		$hash->{validation_method} = $v->{method};
		$hash->{validation_sanger} = $v->{validation_sanger};
		$hash->{validation_ngs}    = $v->{validation};
		$hash->{user_name}         = $v->{user_name};
		if ( length($ref) == length($alt) ) {
			$type               = "variations";
			$strucType          = "snp";
			$hash->{ref_allele} = $ref;
			$hash->{var_allele} = $alt;
		}
		elsif ( length($ref) > length($alt) ) {
			$type      = "deletions";
			$strucType = "del";
			my $t    = 0;
			my @lRef = split( "", $ref );
			while ( $t < length($alt) ) {
				shift(@lRef);
				$t++;
			}
			$hash->{ref_allele} = join( "", @lRef );
			$hash->{var_allele} = "-";
			$hash->{end} = $hash->{start} + length( $hash->{ref_allele} ) - 1;
		}
		elsif ( length($ref) < length($alt) ) {
			$type      = "insertions";
			$strucType = "ins";
			my $t    = 0;
			my @lAlt = split( "", $alt );
			while ( $t < length($ref) ) {
				shift(@lAlt);
				$t++;
			}
			$hash->{ref_allele} = $lAlt[0];
			$hash->{var_allele} = join( "", @lAlt );
			$hash->{end}        = $start + length($alt) - 1;
		}
		$hash->{structuralTypeObject}                     = $type;
		$hash->{structuralType}                           = $strucType;
		$tt->{$strucType}->{ $hash->{id} }                = $hash;
		$self->{hash_validated_variants}->{ $hash->{id} } = $hash;
	}

	$self->{validated_variants} = $tt;
}

sub getListValidatedVariants {
	my ($self) = @_;
	$self->setValidatedVariants();
	return $self->{hash_validated_variants};
}

has infos_validation => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self  = shift;
		my $hash  = $self->getValidatedVariants();
		my $hash2 = {};
		foreach my $h ( values %$hash ) {
			foreach my $k ( keys %$h ) {
				$hash2->{$k} = $h->{$k};
			}
		}
		return $hash2;
	}
);

sub getValidatedVariants {
	my ($self) = @_;

	return $self->{validated_variants} if exists $self->{validated_variants};
	$self->setValidatedVariants();
	return $self->{validated_variants};

}

sub _addInfosValidation {
	my ( $self, $variants, $t, $reference ) = @_;
	$t = "ins" if $t eq "insertions";
	$t = "del" if $t eq "deletions";
	$t = "snp" if $t eq "variations";
	return unless $self->{validated_variants};
	my @types = ( "ins", "del", "snp" );
	foreach my $t (@types) {

		my $vs = $self->{validated_variants};
		my $cc;
		foreach my $v ( values %{ $variants->{$t} } ) {
			my $id = $v->{id};
			if ( exists $vs->{$t}->{ $v->{id} } ) {
				my $v1 = $vs->{$t}->{ $v->{id} };
				$v->{validation_method}           = $v1->{validation_method};
				$v->{validation_sanger}           = $v1->{validation_sanger};
				$v->{validation_ngs}              = $v1->{validation_ngs};
				$self->{validation_method}->{$id} = $v1->{validation_method};
				$self->{validation_sanger}->{$id} = $v1->{validation_sanger};
				$self->{validation_ngs}->{$id}    = $v1->{validation_ngs};
				delete $vs->{$t}->{ $v->{id} };
			}
			else {

				#$self->{validation_method}->{$id}  ="none";
				#$self->{validation_sanger} ->{$id} =0;
				#$self->{validation_ngs}->{$id}  =0;
			}

		}

	}
}

my $already_parse;

sub setVariantsForReference {
	my ( $self, $reference, $typeVar, $cursor ) = @_;
	$self->setValidatedVariants() if $self->project->isDiagnostic();
	my @objs;
	my $hashRes;
	my $z = time;
	foreach my $this_typeVar ( keys %{ $self->callingFiles() } ) {
		my $hfiles = $self->callingFiles->{$this_typeVar};
		foreach my $method ( keys %{$hfiles} ) {
			#next unless $method eq 'manta';
			#next unless $method eq 'haplotypecaller4';
			my $vcfFile = $hfiles->{$method};
			next if exists $already_parse->{ $reference->name }->{$vcfFile};

			#	warn "coucou ".$reference->name;
			$already_parse->{ $reference->name }->{$vcfFile}++;
			$self->{queryVcf}->{$vcfFile} = $self->getQueryVcf( $vcfFile, $method ) unless exists $self->{queryVcf}->{$vcfFile};
			my $queryVcf = $self->{queryVcf}->{$vcfFile};
			my $z        = $queryVcf->parseVcfFileForReference($reference);
			foreach my $type ( keys %$z ) {
				foreach my $id ( keys %{ $z->{$type} } ) {
					if ( exists $hashRes->{$type}->{$id} ) {
						my $h1 = thaw( decompress( $hashRes->{$type}->{$id} ) );
						my $h2 = thaw( decompress( $z->{$type}->{$id} ) );

						$h1->{annex}->{ $self->id }->{method_calling}->{$method}= $h2->{annex}->{ $self->id }->{method_calling}->{$method};
						$hashRes->{$type}->{$id} = compress( freeze $h1);

					}
					else {
						$hashRes->{$type}->{$id} = $z->{$type}->{$id};
					}

				}
			}
		}
	}

	my $o = [];
	#warn "\t end parsing :".abs(time-$z)." ".$self->name;
	$o = $self->myflushobjects2( $hashRes, $cursor );
	
	return $o;
}

has tempArray => (
	is => 'rw',

	#isa		=> 'ArrayRef[Str]',
	lazy    => 1,
	default => sub {
		[],;
	},
);

#sub setVariantsForReference2 {
#	my ( $self, $reference, $typeVar ) = @_;
#	$self->setValidatedVariants() if $self->project->isDiagnostic();
#	my @objs;
#	my $hashRes;
#	my $z = time;
#	foreach my $this_typeVar ( keys %{ $self->callingFiles() } ) {
#		my $hfiles = $self->callingFiles->{$this_typeVar};
#		foreach my $method ( keys %{$hfiles} ) {
#			my $vcfFile = $hfiles->{$method};
#			next if exists $already_parse->{ $reference->name }->{$vcfFile};
#
#			#	warn "coucou ".$reference->name;
#			$already_parse->{ $reference->name }->{$vcfFile}++;
#			$self->{queryVcf}->{$vcfFile} =
#			  $self->getQueryVcf( $vcfFile, $method )
#			  unless exists $self->{queryVcf}->{$vcfFile};
#			my $queryVcf = $self->{queryVcf}->{$vcfFile};
#			my $z        = $queryVcf->parseVcfFileForReference($reference);
#			foreach my $type ( keys %$z ) {
#				foreach my $id ( keys %{ $z->{$type} } ) {
#					if ( exists $hashRes->{$type}->{$id} ) {
#						my $h1 = thaw( decompress( $hashRes->{$type}->{$id} ) );
#						my $h2 = thaw( decompress( $z->{$type}->{$id} ) );
#
#						$h1->{annex}->{ $self->id }->{method_calling}->{$method}
#						  = $h2->{annex}->{ $self->id }->{method_calling}
#						  ->{$method};
#						$hashRes->{$type}->{$id} = compress( freeze $h1);
#
#					}
#					else {
#						$hashRes->{$type}->{$id} = $z->{$type}->{$id};
#					}
#
#				}
#			}
#		}
#	}
#	my $o = [];
#
#	warn "\t end parsing :" . abs( time - $z ) . " " . $self->name;
#	$o = $self->myflushobjects2($hashRes);
#	return $o;
#}

sub _newDeletions {
	my ( $self, $pos, $ref, $alt, $len, $chr ) = @_;

}

sub returnListVarObject {
	my ( $self, $type, $hchrs ) = @_;
	my @lVarToReturn = ();
	foreach my $var ( @{ $self->getObjects($type) } ) {
		foreach my $hchr (@$hchrs) {
			next if $var->getChromosome->id ne $hchr->{obj}->id;
			if ( $var->getGenomicSpan()->intersection( $hchr->{intspan} ) ) {
				push( @lVarToReturn, $var );
			}
		}
	}
	return \@lVarToReturn;
}

has vcfParsed => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $lChr = $self->getProject()->getChromosomes();
		my %hashVcf;
		foreach my $chrObj (@$lChr) {
			my $lVcfSnp = $self->getVariationsFiles();
			foreach my $vcfFile (@$lVcfSnp) {
				$hashVcf{$vcfFile}->{ $chrObj->id } =
				  Set::IntSpan::Fast::XS->new();
			}
			my $lVcfIndel = $self->getIndelsFiles();
			foreach my $vcfFile (@$lVcfIndel) {
				$hashVcf{$vcfFile}->{ $chrObj->id } =
				  Set::IntSpan::Fast::XS->new();
			}
		}
		return \%hashVcf;
	}
);

sub _completHashesForOneChromosome {
	my ( $self, $vcfFile, $hchr ) = @_;
	confess() unless $hchr->{obj};
	my $intSpan_alreadyParsed =
	  $self->vcfParsed()->{$vcfFile}->{ $hchr->{obj}->id };
	my $intSpan_posToParse = $hchr->{intspan}->diff($intSpan_alreadyParsed);
	$self->vcfParsed()->{$vcfFile}->{ $hchr->{obj}->id } =
	  $hchr->{intspan}->union($intSpan_alreadyParsed);
	return $intSpan_posToParse;
}

my %test;

sub myflushobjects2 {
	my ( $self, $hashRes, $cursor ) = @_;
	my @objs;
	my $z = time;
	my $x = 0;
	if ( scalar( keys %$hashRes ) > 0 ) {
		foreach my $structType ( keys %$hashRes ) {
			while ( my ( $keyId, $valHash ) = each( %{ $hashRes->{$structType} } ) )
			{
				if ($cursor) {
					push( @{ $self->tempArray }, $valHash );
					next;
				}

				$x++;
				$valHash = thaw( decompress($valHash) );
				my $varObj = $self->getProject()->flushObject( $valHash->{structuralTypeObject}, $valHash );
				#warn $varObj;
				$varObj->patients_object()->{ $self->id() } = undef;

				#here add all method for this variant in method calling hash

				$varObj->annex()->{ $self->id } =  $valHash->{annex}->{ $self->id };
				
				$valHash = undef;
				delete $hashRes->{$structType}->{$keyId};
				push( @objs, $varObj );
			}
		}
		$hashRes = {};
	}
	my $t = abs( $z - time );
	
#	warn $t." nb objs => ".scalar(@objs)." ".$self->name." construct ".$x if scalar(@objs) >0 ;
	
	#else { warn 'No object...'; }
	return \@objs;
}

sub setStructuralVariants {
	my ($self) = @_;
	foreach my $method ( @{ $self->callingSVMethods } ) {
		foreach my $chr ( @{ $self->project->getChromosomes } ) {
			my $ao = $self->setStructuralVariantsForReference( $chr, $method );

		}
	}

}

sub setStructuralVariantsForReference {
	my ( $self, $chr, $caller ) = @_;
	my $file  = $self->getSVFile($caller);
	my $query = QueryVcfSV->new( patient => $self );
	my $array = $query->parseVcfFile( $self, $chr, $file, $caller );
	my $ao;
	foreach my $hash (@$array) {
		my $obj = $self->project->flushObjectSV($hash);
		$obj->add_event( $self, $caller, $hash );
		$obj->{patients_object}->{ $self->id() }   = undef;
		$chr->{ $obj->type_object }->{ $obj->id }  = undef;
		$self->{ $obj->type_object }->{ $obj->id } = undef;

		push( @$ao, $obj );

	}
	return $ao;

}

sub getStructuralVariations {
	my ( $self, $chrName, $start, $end ) = @_;
	my $lRes = $self->getVariations( $chrName, $start, $end );
	push( @$lRes, @{ $self->getInsertions( $chrName, $start, $end ) } );
	push( @$lRes, @{ $self->getDeletions( $chrName, $start, $end ) } );
	push( @$lRes, @{ $self->getMnps( $chrName, $start, $end ) } );
	push( @$lRes, @{ $self->getLargeDeletions( $chrName, $start, $end ) } );
	push( @$lRes, @{ $self->getLargeDuplications( $chrName, $start, $end ) } );
	push( @$lRes, @{ $self->getInversions( $chrName, $start, $end ) } );
	push( @$lRes, @{ $self->getBoundaries( $chrName, $start, $end ) } );
	return $lRes;
}

sub getQueryVcf {
	my ( $self, $fileName, $method ) = @_;
	my %args;
	$args{patient} = $self;
	$args{file}    = $fileName;
	$args{method}  = $method;
	my $queryVcf = QueryVcf->new( \%args );
	return $queryVcf;
}

sub isMale {
	my $self = shift;
	return 1 if $self->sex == 1;
}

sub isFemale {
	my $self = shift;
	return 1 if $self->sex == 2;
}

sub isChildren {
	my $self = shift;
	confess();
	return 1 if ( exists $self->getFamily->children->{ $self->name() } );
}

sub isParent {
	my $self = shift;
	return !( $self->isChild );
}

sub isMother {
	my $self = shift;
	return if $self->isChild();
	return if $self->isMale();
	return 1;
}

sub isFather {
	my $self = shift;
	return if $self->isChild();
	return if $self->isFemale();
	return 1;
}

has callingFiles => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		my $methods = $self->getCallingMethods();
		my $files   = {};
		foreach my $method_name (@$methods) {
			my $file = $self->_getCallingFileWithMethodName( $method_name, "variations" );
			$files->{variations}->{$method_name} = $file if $file;

		}
		foreach my $method_name (@$methods) {
			my $file2 = $self->_getCallingFileWithMethodName( $method_name, "indels" );
			$files->{indels}->{$method_name} = $file2 if $file2;
		}
		my $file3 = $self->getLargeIndelsFile();
		$files->{large_indels}->{'lifinder'} = $file3 if ( -e $file3 );

#		die();
#		foreach my $method_name (@$methods) {
#			my $file2 = $self->_getCallingFileWithMethodName($method_name,"large_indels");
#			$files->{large_indels}->{$method_name} = $file2 if $file2;
#		}

		return $files;
	},

);



has callingSVFiles => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		my $methods = $self->callingSVMethods();
		my $files   = {};
		foreach my $method_name (@$methods) {
			my $file = $self->_getCallingSVFileWithMethodName( $method_name,"variations" );
			$files->{sv}->{$method_name} = $file if $file;
		}
		return $files;
	},
);

has bamUrl => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		my $methods = $self->alignmentMethods();
		die() if ( scalar( @$methods > 1 ) );
		my $method_name = $methods->[0];
		my $bam_dir     = $self->getProject->getAlignmentUrl($method_name);
		my $bamf        = $self->getBamFileName();
		my (@t) = split( "/", $bamf );

		#warn $bam_dir;
		return $bam_dir . "/" . $t[-1];
	},
);

has hasBamFile => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		my ($no) =  grep {$_ eq "no_align"} @{$self->alignmentMethods()};
		return undef if $no;
		return 1;
	},
);

has bamFiles => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		#warn $self->name."----";
		my $methods = $self->alignmentMethods();
		my $files   = {};
		foreach my $method_name (@$methods) {
			my $bam_dir = $self->getProject->getAlignmentDir($method_name);
			next unless -e $bam_dir;

			my $bam_file = $self->_getFileByExtention( $bam_dir, "align" );
			next unless $bam_file;

			$files->{$method_name} = $bam_dir . "/" . $bam_file;
		}
		return $files;
	},

);

sub get_fastq {
	my ($self) = @_;
	require "fastq_util.pm";
	my $name=$self->name();
	my $dir = $self->getSequencesDirectory();
	my @names;
	my $couple;
	push(@names,$self->name);
	push(@names,$self->barcode) if length($self->barcode)>1 ;
	unless (exists $self->buffer->{cached_dir}->{$dir}){
		$self->buffer->{cached_dir}->{$dir}= [];
		opendir(DIR,$dir);
		my @allFiles= readdir(DIR);

		$self->buffer->{cached_dir}->{$dir} = \@allFiles;
	}
	NAME: foreach my $name (@names) {
	my @pattern = ("^".$name."_[ATGC][ATGC][ATGC]","^".$name."_S[1-9]+","^".$name."_","$name");
	foreach my $find (@pattern){
		my (@titi) = grep { /$find/} grep { /fastq/} @{$self->buffer->{cached_dir}->{$dir}} ;
	
		if (@titi) {
			$couple = fastq_util::find_paired_files(\@titi,$dir);
			last NAME if ($couple);
		}
	}
	}
	warn "NO fastq file for : -".$self->name()."- ".$self->barcode." ".$dir unless $couple;	
	foreach my $cp (@$couple) {
		$cp->{R1} = $self->getSequencesDirectory()."/".$cp->{R1};
		$cp->{R2} = $self->getSequencesDirectory()."/".$cp->{R2};
		delete $cp->{dif};
		delete $cp->{pos};
	}
	return $couple if scalar(@$couple)>0;
	die();	
	

}




has trackingFile => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self         = shift;
		my $tracking_dir = $self->getProject->getPipelineTrackingDir();
		return
			$self->getProject->getPipelineTrackingDir() . "/"
		  . $self->name . ".json";
	},
);

sub getTracking {
	my ($self) = @_;
	my $file = $self->trackingFile;
	warn $file;
	return {} unless -e $file;
	return {} if -z $file;
	open( JSON, $file );
	my $desc = decode_json <JSON>;
	close(JSON);
	return $desc;
}

sub getMetricsFile {
	my $self = shift;
	return $self->project->getMetricsDir() . "/" . $self->name() . ".stats";
}

sub getCaptureMetricsFile {
	my ( $self, $headerHSmetrics ) = @_;
	return $self->setCaptureMetricsFile($headerHSmetrics);
}

sub setCaptureMetricsFile {
	my ( $self, $headerHSmetrics ) = @_;
	my $bedFileName = $self->getCaptureBedFile();

	#warn  $bedFileName ;
	my $metricsFileName = $bedFileName;
	$metricsFileName =~ s/\.bed/\.metrics/;
	my $intervalFileName = $metricsFileName . ".intervals";
	return $intervalFileName if -e $intervalFileName;

	my $cmd =
"cat  $bedFileName |awk  -F \"\\t\" 'BEGIN {OFS=\"\\t\"} {print \$1,\$2,\$3,\"+\",(\$3-\$2)+1}' > $metricsFileName";
	`$cmd`;
	my $cmd2 = "cat $headerHSmetrics $metricsFileName > $intervalFileName";
	`$cmd2`;
	return $intervalFileName;

}

sub getBamFiles {
	my ($self) = @_;
	my (@t)    = values %{ $self->bamFiles() };
	return \@t;
}

sub getRecalFile {
	my ($self) = @_;
	my $methods = $self->alignmentMethods();
	return if scalar(@$methods) > 1;
	my $recal =
		$self->project->getRecalDir( $methods->[0] ) . "/"
	  . $self->name
	  . ".recal.table";
	warn( $recal . " is empty" ) if ( -z $recal );
	return
		$self->project->getRecalDir( $methods->[0] ) . "/"
	  . $self->name
	  . ".recal.table";
}



sub getBamFileName {
	my ( $self, $method_name ) = @_;
	my $bam_dir;
	if($method_name){
		 $bam_dir = $self->getProject->getAlignmentDir( $method_name );
		 
	}
	else {
	my $methods = $self->alignmentMethods();
	die( $self->project->name." ".Dumper($methods)) if scalar(@$methods) > 1;
	 $bam_dir = $self->getProject->getAlignmentDir( $methods->[0] );
	}
	die() unless $bam_dir;
	my $bam     = $bam_dir . "/" . $self->name . ".bam";
	return $bam;
}
sub getCramFileName {
	my ( $self, $method_name,$version ) = @_;
	my $bam_dir;
	if($method_name){
		 $bam_dir = $self->getProject->getAlignmentDir( $method_name,$version );
		 
	}
	else {
	my $methods = $self->alignmentMethods();
	die( $self->project->name." ".Dumper($methods)) if scalar(@$methods) > 1;
	 $bam_dir = $self->getProject->getAlignmentDir( $methods->[0],$version );
	}
	die() unless $bam_dir;
	my $bam     = $bam_dir . "/" . $self->name . ".cram";
	return $bam;
}
sub getBamFile {
	my ( $self, $method_name, $nodie ) = @_;
	
	unless ($method_name) {
		my $files = $self->getBamFiles();
		return $files->[0] if scalar(@$files) == 1;
		if ($nodie) {
			warn "NO BAM FILES " . $self->name;
			return;
		}

	  		confess($self->getBamFileName." ".$self->name." :: "." \n:: ".Dumper  $self->alignmentMethods());
	}
	 
	confess("ERROR: no bam file with $method_name method name. Exit. "
		  . $self->name." "
		  . $self->project->name
		  . "\n\n" )
	  unless exists $self->bamFiles()->{$method_name};
		return $self->bamFiles()->{$method_name};
	#return  $self->{files}->{alignment}->{alignment}; # Pas rempli cette table
}

sub getSVFiles {
	my $self      = shift;
	my @lVcfFiles = values %{ $self->callingSVFiles()->{sv} };
	warn( " warn I was unable to find variation vcf file for :  "
		  . $self->name() )
	  unless scalar @lVcfFiles;
	return \@lVcfFiles;
}

sub getSVFile {
	my ( $self, $method ) = @_;
	if ($method) {

		#	warn Dumper $self->callingSVFiles();
		confess( "can t find vcf for $method " . $self->name )
		  unless exists( $self->callingSVFiles()->{sv}->{$method} );
		return $self->callingSVFiles()->{sv}->{$method};
	}
	my @all = values %{ $self->callingSVFiles()->{sv} };
	return "" if scalar(@all) eq 0;

	confess($self->name
		  . "you have exactly "
		  . scalar(@all)
		  . " methods defined on your project "
		  . $self->getProject->name() )
	  if scalar(@all) ne 1;
	return $all[0];
}

sub getAnnotSVFileName {
	my ( $self, $method ) = @_;
	my $annot_dir = $self->getProject->getAnnotSVDir($method);
	my $annot_file =
	  $annot_dir . "/" . $method . "_SV_" . $self->name . ".annotated.tsv";
	return $annot_file;
}

sub getBestOneFileName {
	my $self = shift;
	my $dir  = $self->project->getCNVDir();
	return $dir . "/" . $self->name . ".bestCNV";
}

sub getBestOne {
	my $self = shift;
	my $file = $self->getBestOneFileName();
	return [] unless -e $file;
	open( JSON, $file );
	my $best = decode_json <JSON>;
	close(JSON);
	return $best->{items};
}

sub getVariationsFiles {
	my $self      = shift;
	my @lVcfFiles = values %{ $self->callingFiles()->{variations} };
	warn( " warn I was unable to find variation vcf file for :  "
		  . $self->name() )
	  unless scalar @lVcfFiles;
	return \@lVcfFiles;
}

sub getVariationsFileName {
	my ( $self, $method, $nodie ) = @_;
	my $dir = $self->project->getVariationsDir($method);
	return $dir . "/" . $self->name . ".vcf.gz";
}

sub vcfFileName {
	my ( $self, $method, $nodie ) = @_;
	return $self->getVariationsFileName($method);
}

sub getVariationsFile {
	my ( $self, $method, $nodie ) = @_;
	if ($method) {
		if ($nodie) {
			return ""
			  unless exists( $self->callingFiles()->{variations}->{$method} );
		}
		confess( "can t find vcf for $method" . Dumper $self->callingFiles() )
		  unless exists( $self->callingFiles()->{variations}->{$method} );
		return $self->callingFiles()->{variations}->{$method};
	}
	my @all = values %{ $self->callingFiles()->{variations} };
	return "" if scalar(@all) eq 0;

	confess($self->name
		  . "you have exactly "
		  . scalar(@all)
		  . " !$method! methods defined on your project "
		  . $self->getProject->name() )
	  if scalar(@all) ne 1;
	return $all[0];
}

sub getIndelsFiles {
	my $self      = shift;
	my @lVcfFiles = values %{ $self->callingFiles()->{indels} };

#warn ("no indels file ".$self->name." ".$self->getProject->name) unless scalar @lVcfFiles;
	return \@lVcfFiles;
}

sub getIndelsFile {
	my ( $self, $method, $nodie ) = @_;
	if ($method) {
		confess("can t find vcf for $method")
		  unless exists( $self->callingFiles()->{indels}->{$method} );
		return $self->callingFiles()->{indels}->{$method};
	}
	my @all = values %{ $self->callingFiles()->{indels} };
	return "" if scalar(@all) eq 0;
	confess("you have exactly "
		  . scalar(@all)
		  . " methods defined on your project " )
	  if scalar(@all) ne 1;
	return $all[0];
}

sub getCnvsFiles {
	my $self = shift;
	my $file =
		$self->project->getLargeIndelsDir("lifinder") . "/"
	  . $self->name
	  . ".bed.gz";
	return [$file] if -e $file;
	return [];

}

sub getLargeIndelsFiles {
	my $self = shift;
	my $file =
		$self->project->getLargeIndelsDir("lifinder") . "/"
	  . $self->name
	  . ".bed.gz";
	return [$file] if -e $file;
	return [];

}
sub targetGCFile {
	my ( $self ) = @_;
	my $file = $self->project->getTargetCountDir()."/".$self->name.".target.counts.gc-corrected.gz";
	return $file;
}
sub targetFile {
	my ( $self ) = @_;
	my $file = $self->project->getTargetCountDir()."/".$self->name.".target.counts.gz";
	return $file;
}
sub gvcfFileName {
	my ( $self, $method ) = @_;
	$method = "haplotypecaller4" unless $method;
	my $dir = $self->project->getGvcfDir($method);
	my $f   = $self->name . ".g.vcf.gz";
	return $dir . "/" . $f;

}

sub getGvcfFile {
	my ( $self, $method, $nodie ) = @_;
	$method = "haplotypecaller4" unless $method;
	my $f = $self->gvcfFileName($method);
	return $f if -e $f;
	return undef;

}

sub getCnvsFile {
	my ( $self, $method, $nodie ) = @_;
	my $file =
		$self->project->getLargeIndelsDir("lifinder") . "/"
	  . $self->name
	  . ".bed.gz";
	return $file;
}

sub getLargeIndelsFile {
	my ( $self, $method, $nodie ) = @_;
	my $file =
		$self->project->getLargeIndelsDir("lifinder") . "/"
	  . $self->name
	  . ".bed.gz";
	return $file;
}

sub _getCallingFileWithMethodName {
	my ( $self, $method_name, $type ) = @_;
	confess() unless $method_name;
	my $project = $self->getProject();
	my $dir;
	if ( $type eq 'variations' ) {
		$dir = $project->getVariationsDir($method_name);
	}
	elsif ( $type eq 'large_indels' ) {
		$dir = $project->getLargeIndelsDir($method_name);
	}
	else { $dir = $project->getIndelsDir($method_name); }

	unless ( -e $dir ) {
		warn "\n\nWARNING: no directory named $method_name found for "
		  . $self->name()
		  . " patient !!\n\n";
		return;
	}
	my $file = $self->_getFileByExtention( $dir, "variation" );

	return unless $file;

	confess( "no vcf file for " . $self->name() . " method is " . $method_name )
	  unless $file;
	return $dir . "/" . $file;
}

sub _getCallingSVFileWithMethodName {
	my ( $self, $method_name, $type ) = @_;
	confess() unless $method_name;
	my $project = $self->getProject();
	my $dir;
	$dir = $project->getVariationsDir($method_name);

	unless ( -e $dir ) {
		warn "\n\nWARNING: no directory named $method_name found for "
		  . $self->name()
		  . " patient !!\n\n";
		return;
	}
	my $SV_vcf = $self->name() . ".vcf.gz";
	my $SV_bed = $self->name() . "_aberrations.bed.gz";

	my $file;
	$file = $SV_vcf if -e $dir . "/" . $SV_vcf;
	$file = $SV_bed if -e $dir . "/" . $SV_bed;

	confess("no vcf or bed file for "
		  . $self->name()
		  . " method is "
		  . $method_name )
	  unless $file;
	return $dir . "/" . $file;
}

sub _getFileByExtention {
	my ( $self, $dir, $type ) = @_;
	my $name = $self->name();
	confess() unless $dir;
	my $exts = {
		"variation" => [
			"vcf.gz", "gz",  "vcf", "gff3", "sam", "txt",
			"casava", "tab", "casava2"
		],
		"align" => ["bam","cram"]
	};

	#my $f = File::Util->new();
	foreach my $ext ( @{ $exts->{$type} } ) {
		my $theoric = $self->name() . ".$ext";
		return $theoric if -e $dir . "/" . $theoric;
		my @files = grep { /^$name\./ }
		  grep { /\.$ext$/ } @{ $self->buffer->getListFiles($dir) };
		if ( scalar(@files) > 1 ) {
			my @files1 = grep { $_ !~ /gatk/ } @files;
			return $files1[0];

		}
		return $files[0] if scalar(@files);
	}
	return;
}

sub getCoverageFileName {
	my $self = shift;
	my $dir  = $self->getProject->getCoverageDir();
	return $dir . $self->name() . ".cov.gz";
}

sub getCoverageFile {
	my $self = shift;
	my $dir  = $self->getProject->getCoverageDir();
	return $dir . $self->name() . ".cov.gz";
}

sub fileNoSqlDepth {
	my ($self) = @_;
#	warn $self->NoSqlDepthDir;
	system( "mkdir " . $self->NoSqlDepthDir ) unless -d $self->NoSqlDepthDir;
	return $self->NoSqlDepthDir . "/" . $self->name . ".depth.lmdb";
}

sub fileWiseCondor {
	my ($self) = @_;
	system( "mkdir " . $self->getProject->getCoverageDir() . "/wisecondor/" )
	  unless -d $self->getProject->getCoverageDir() . "/wisecondor/";
	return
		$self->getProject->getCoverageDir()
	  . "/wisecondor/"
	  . $self->name . ".npz";
}
sub rawDataWiseCondor {
 my ($self) = @_;	
 my $dir = $self->project->getVariationsDir("wisecondor");
 return $dir."/".$self->name."_bins.bed.gz";
}
has tabix_wisecondor => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $coverage_file;
		$coverage_file = $self->rawDataWiseCondor();
		return unless -e $coverage_file;
		
		return  Bio::DB::HTS::Tabix->new( filename => $coverage_file );;

	}
);


has vntyperTsv => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $file = $self->project->getVariationsDir("vntyper")."/muc1/".$self->name."_Final_result.tsv";
	}
);
has isVntyperPositif => (
is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return 0 unless -e  $self->project->getVariationsDir("vntyper")."/muc1/".$self->name."_Final_result.tsv";
		
		my $file = $self->project->getVariationsDir("vntyper")."/muc1/".$self->name."_Final_result.tsv";
		my @lines = `tail -n +4 $file`;
		#warn Dumper @lines;
		chomp(@lines);
		#confess() if scalar(@lines) > 1;
		 return undef unless $lines[0];    
		 return 1;
	}
);

sub kestrel {
	my ($self) = @_;
	return $self->{kestrel} if exists $self->{kestrel};
	$self->vntyperResults();
	
	if (@{$self->{kestrel}->[0]} > 3){
		my $l = length ($self->{kestrel}->[0]->[7]);
		my $pos = $self->{kestrel}->[0]->[4];
		my $position = ($l - $pos);
		my $reverse = BioTools::complement_sequence($self->{kestrel}->[0]->[7]);
		my $ref = $self->{kestrel}->[0]->[5] ;
		my $alt = $self->{kestrel}->[0]->[6];
		my $left = substr($reverse, 0, $position+1);
		my $right = substr($reverse, $position+1);
		if (length($alt)> length($ref)){
			#insertion 
			 my $ralt = BioTools::complement_sequence(substr($alt, 1));
		my $ins = qq{<span style="color:red">[<span style="text-emphasis: double-circle red; ">$ralt</span/>]</span/>};
			#  my $ins = qq{<span style="color:blue;text-emphasis: double-circle blue; ">$ralt</span/>};
			$self->{kestrel}->[0]->[7] =$left.$ins.$right;
		}
		elsif (length($alt)< length($ref)){
			my $ralt = BioTools::complement_sequence(substr($alt, 1));
			 my $ins = qq{<span style="color:red">[$ralt]</span/>};
			my $right = substr($reverse, $position+length($ralt));
			$self->{kestrel}->[0]->[7] =$left.$ins.$right;
		}
		elsif (length($alt) == length($ref)){
			my $ralt = BioTools::complement_sequence($alt);
			 my $ins = qq{<span style="color:red">[$ref/$ralt]</span/>};
			  #my $ins = qq{<span style="color:red;text-decoration=underline overline">[$ref/$ralt]</span/>};
			 $left = substr($reverse, 0, $position);
			my $right = substr($reverse, $position+2);
			$self->{kestrel}->[0]->[7] =$left.$ins.$right;
		}
		
		#my @ins = split("",$self->{kestrel}->[0]->[6]); 
		#shift(@ins);
		#my $t = join("",@ins);
		#$self->{kestrel}->[0]->[5] = BioTools::complement_sequence($self->{kestrel}->[0]->[5]);
		#$self->{kestrel}->[0]->[6] = BioTools::complement_sequence($self->{kestrel}->[0]->[6]);
		#$self->{kestrel}->[0]->[7] = $self->{kestrel}->[0]->[7]."<BR>".$left."[".$position.$self->{kestrel}->[0]->[5]."/".$self->{kestrel}->[0]->[6]."]".$right."<BR>".$self->{kestrel}->[0]->[7];
		#$self->{kestrel}->[0]->[7] = $left.$right."<BR>".$self->{kestrel}->[0]->[7];
		
	}
	$self->{kestrel}->[0]->[0] .= "XXX";
	
	return $self->{kestrel};
}
sub adVNTR {
	my ($self) = @_;
	return $self->{adVNTR} if exists $self->{adVNTR};
	$self->vntyperResults();
	return $self->{adVNTR};
}
	
sub vntyperResults {
		my $self = shift;
		my $file = $self->project->getVariationsDir("vntyper")."/muc1/".$self->name."_Final_result.tsv";
		unless( -e $file){
			$self->{kestrel} =[];
			$self->{adVNTR} =[];
			return;
 		}
 		
		
			my $date = POSIX::strftime( 
             "%d/%m/%y", 
             localtime( 
               		(stat $file )[10]
                 )
             );
		my @t;
		my @lines = `cat $file`;
		chomp(@lines);
		my $kestrel = firstidx { $_ =~ /Kestrel_Result/ } @lines;
		my $limit2 = $kestrel - 1;
		
		my $advntr = firstidx { $_ =~ /adVNTR_Result/ } @lines;
		my $limit = $advntr ;
		$limit = scalar(@lines) if $advntr < $kestrel;
	
		for (my $i = $kestrel +2 ; $i<$limit;$i++){
			push(@t,[$date,"kestrel",split(" ",$lines[$i])] ) ;
		}
		push(@t,[$date]) unless @t; 
		$self->{header_adVNTR} = [];
		$self->{header_adVNTR} = ["date",split(" ",$lines[$advntr+1])] if $advntr > 0;
		$self->{header_kestrel} = ["date",split(" ",$lines[$kestrel+1])]; 
		$self->{kestrel} =\@t;
		$limit2 = scalar(@lines) if $advntr >= $kestrel ;
		warn $advntr." ".$limit2;
		my @t2;
		
		for (my $i = $advntr+2 ; $i<$limit2;$i++) {
			push(@t2,[$date,"adVNTR",split(" ",$lines[$i])] ) ;
		}
		$date = "-" if $advntr == -1;
		push(@t2,[$date,"adVNTR"]) unless @t2; 
		$self->{adVNTR} =\@t2;
		return;
		#confess() if scalar(@lines) > 1;

      #      
	
		
		#return \@t;
	}


sub isNoSqlDepth {
	my ($self) = @_;
	return -e $self->fileNoSqlDepth;
}

has NoSqlDepthDir => (
	is => 'rw',

	#isa		=> 'ArrayRef[Str]',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->getProject->getCoverageDir() . "/lmdb_depth";
		return $dir;
	},
);

sub getNoSqlDepth {
	my ( $self, $mode ) = @_;
	$mode = "r" unless $mode;
	my $buffer = $self->buffer;
	return $buffer->{lmdb_hash}->{depth}->{ $self->name }
	  if exists $buffer->{lmdb_hash}->{depth}->{ $self->name };
	my $dir = $self->NoSqlDepthDir;
	
	unless ( -e $dir ) {
		system("mkdir $dir;chmod a+rwx $dir");
	}
	$buffer->{lmdb_hash}->{depth}->{ $self->name } = GenBoBinaryFile->new(
		name => $self->name . ".depth.lmdb",
		dir  => $dir,
		mode => $mode
	);
	return $buffer->{lmdb_hash}->{depth}->{ $self->name };
}

sub getTranscriptsCoverageDepth {
	my ( $self, $mode ) = @_;
	my $buffer = $self->buffer;
	$mode = "r" unless $mode;
	if ( $mode eq "d" ) {
		$buffer->{lmdb_hash}->{transcripts}->{ $self->name }->close
		  if exists $buffer->{lmdb_hash}->{transcripts}->{ $self->name };
		delete $buffer->{lmdb_hash}->{transcripts}->{ $self->name };
		return;
	}

	return $buffer->{lmdb_hash}->{transcripts}->{ $self->name }
	  if exists $buffer->{lmdb_hash}->{transcripts}->{ $self->name };

	my $dir = $self->project->transcriptsCoverageDir();
	$buffer->{lmdb_hash}->{transcripts}->{ $self->name } = GenBoNoSqlLmdb->new(
		name        => $self->name . ".transcripts",
		dir         => $dir,
		mode        => $mode,
		is_compress => 1,
		vmtouch     => $self->buffer->vmtouch
	);
	return $buffer->{lmdb_hash}->{transcripts}->{ $self->name };
}

sub getGenesDude {
	my ( $self, $mode ) = @_;
	$mode = "r" unless $mode;
	my $buffer = $self->buffer;
	if ( $mode eq "d" ) {
		$buffer->{lmdb_hash}->{genes_dude}->{ $self->name }->close
		  if exists $buffer->{lmdb_hash}->{genes_dude}->{ $self->name };
		delete $buffer->{lmdb_hash}->{genes_dude}->{ $self->name };
		return;
	}

	return $buffer->{lmdb_hash}->{genes_dude}->{ $self->name }
	  if exists $buffer->{lmdb_hash}->{genes_dude}->{ $self->name };

	my $dir = $self->project->transcriptsDudeDir();
	$buffer->{lmdb_hash}->{genes_dude}->{ $self->name } = GenBoNoSqlLmdb->new(
		name        => $self->name . ".dude.genes",
		dir         => $dir,
		mode        => $mode,
		is_compress => 1,
		vmtouch     => $self->buffer->vmtouch
	);
	return $buffer->{lmdb_hash}->{genes_dude}->{ $self->name };
}

sub getTranscriptsDude {
	my ( $self, $mode ) = @_;
	$mode = "r" unless $mode;
	my $buffer = $self->buffer;
	if ( $mode eq "d" ) {
		$buffer->{lmdb_hash}->{transcripts_dude}->{ $self->name }->close
		  if exists $buffer->{lmdb_hash}->{transcripts_dude}->{ $self->name };
		delete $buffer->{lmdb_hash}->{transcripts_dude}->{ $self->name };
		return;
	}

	return $buffer->{lmdb_hash}->{transcripts_dude}->{ $self->name }
	  if exists $buffer->{lmdb_hash}->{transcripts_dude}->{ $self->name };

	my $dir = $self->project->transcriptsDudeDir();
	$buffer->{lmdb_hash}->{transcripts_dude}->{ $self->name } =
	  GenBoNoSqlLmdb->new(
		name        => $self->name . ".dude.transcripts",
		dir         => $dir,
		mode        => $mode,
		is_compress => 1,
		vmtouch     => $self->buffer->vmtouch
	  );
	return $buffer->{lmdb_hash}->{transcripts_dude}->{ $self->name };
}

sub depth {
	my ( $self, $chr_name, $start, $end ) = @_;
	my $chr   = $self->project->getChromosome($chr_name);
	my $array = $self->getNoSqlDepth->getDepth( $chr->name, $start, $end );
	return $array;
}
sub normalize_depth {
	my ( $self, $chr_name, $start, $end ) = @_;
	my $chr   = $self->project->getChromosome($chr_name);
	my @array = map{int($_/$self->normalized_reads *100)/100} @{$self->getNoSqlDepth->getDepth( $chr->name, $start, $end )};
	return \@array;
}
sub meanDepth {
	my ( $self, $chr, $start, $end ) = @_;
	$end = $start +1 if $start == $end;
	return $self->getNoSqlDepth->getMean( $chr, $start, $end );
}

sub minDepth {
	my ( $self, $chr_name, $start, $end ) = @_;
	my @lcov = sort {$a <=> $b} @{$self->depth($chr_name, $start, $end)};
	return $lcov[0];
}

sub maxDepth {
	my ( $self, $chr_name, $start, $end ) = @_;
	my @lcov = sort {$a <=> $b} @{$self->depth($chr_name, $start, $end)};
	return $lcov[-1];
}

sub depthIntspan {
	my ( $self, $chr, $intspan ) = @_;
	my $array = $self->getNoSqlDepth->getDepthForIntspan( $chr, $intspan );
	return $array;
}

sub statistics_multiplex {
	my ( $self, $multiplex ) = @_;

	return $self->{statistics_multiplex}->{$multiplex}
	  if exists $self->{statistics_multiplex}->{$multiplex};
	my $debug;
	my $primers =
	  $self->getProject->getCapture->getPrimersByMultiplex($multiplex);
	my $score = 0;
	my %genes;
	my $nb = 0;
	$self->{statistics_multiplex}->{$multiplex} =
	  new Statistics::Descriptive::Sparse;
	my @data;

	foreach my $primer (@$primers) {
		next if $primer->getChromosome()->name eq "X";
		push( @data, $primer->cnv_score($self) );

	}
	$self->{statistics_multiplex}->{$multiplex}->add_data(@data);

#		 sub outlier_filter2 { return 1; }
#		 $self->{statistics_multiplex}->{$patient->id} ->set_outlier_filter( \&outlier_filter2 );
#		my @data = $self->{statistics_multiplex}->{$patient->id} ->get_data_without_outliers();
#		$self->{statistics_multiplex}->{$patient->id} ->clear();
#		$self->{statistics_multiplex}->{$patient->id} ->add_data(@data);
#		 $self->{statistics_multiplex}->{$patient->id} ->set_outlier_filter( \&outlier_filter2 );
#		 @data = $self->{statistics_multiplex}->{$patient->id} ->get_data_without_outliers();
#		$self->{statistics_multiplex}->{$patient->id} ->clear();
#		$self->{statistics_multiplex}->{$patient->id} ->add_data(@data);

	return $self->{statistics_multiplex}->{$multiplex};
}

sub get_data_primers {
	my ( $self, $primer_id, $type ) = @_;

	my $db1 = $self->transcriptsCoverage();
	if ($db1) {
		my $count = $db1->get( $primer_id . "$type" );
		if ( defined $count ) {

			#$self->{zscore}->{$patient->id} = $count;
			return $count;
		}
	}
	return undef;
}

sub getCoverage {
	my ( $self, $chr, $start, $end ) = @_;
	my $tabix = $self->tabix_coverage();
	my $len   = abs( $start - $end ) + 1;
	my @v     = ( (0) x $len );
	my $res   = $tabix->query_full( $chr, $start - 1, $end ) if $start;
	my @data;

	while ( my $line = $res->next ) {

		my ( $a, $p, $c ) = split( " ", $line );

		confess() if $a ne $chr;
		my $pv = $p - $start;
		$v[$pv] = $c;
	}

	my $gc = GenBoCoverage->new( start => $start, end => $end, array => \@v );
	return ($gc);
}

sub count {
	my ( $self, $primer_id, $capture ) = @_;
	$self->load_cached_statistics();
	return $self->{statistics}->{$primer_id}->{corrected_mean}
	  if exists $self->{statistics}->{$primer_id}->{corrected_mean};
	if ( ref( $self->getProject() ) eq 'GenBoProjectCache' ) {
		warn $primer_id . " " . $self->name;
		warn Dumper $self->{statistics}->{$primer_id};
		confess();
	}
	unless ($capture) {
		$capture = $self->getCapture();
	}
	my ( $chr_name, $start, $end ) = split( "_", $primer_id );

	$chr_name =~ s/primer//;

	my $chr = $self->project->getChromosome($chr_name);
	my $len = $capture->primer_size($primer_id);

	$self->set_statistics_coverage( $chr, $primer_id, $start,
		( $start + $len ) + 1 );

	#my $tabix = $self->tabix_coverage();

	confess() unless exists $self->{statistics}->{$primer_id}->{corrected_mean};
	return $self->{statistics}->{$primer_id}->{corrected_mean};

}

sub load_cached_statistics {
	my ($self) = @_;
	return if exists $self->{cached};
	my $no   = $self->project->noSqlCoverage();
	my $hash = $no->get_like( $self->name, "primer%" );

	foreach my $key ( keys %$hash ) {
		my ( $chr_name, $start, $end ) = split( "_", $key );
		$self->{statistics}->{$key} = $hash->{$key};
		$self->{statistics}->{$key}->{corrected_mean} =
		  $self->{statistics}->{$key}->{mean};

		if ( $chr_name =~ /X/ ) {
			my $chr = $self->project->getChromosome("X");
			unless ( $chr->isPseudoAutosomal( $start, $end ) ) {
				$self->{statistics}->{$key}->{corrected_mean} *= 2
				  if $self->isMale();
			}

		 #$self->{statistics}->{$key}->{corrected_mean} *= 2 if $self->isMale();
		}

#$self->{statistics}->{$id}->{corrected_mean} = $self->{statistics}->{$id}->{mean};

	}

	$self->{cached} = 1;
	return $self->{statistics};
}

sub set_statistics_coverage {
	my ( $self, $chr, $id, $start, $end ) = @_;

	my $span = Set::IntSpan::Fast::XS->new( $start . "-" . $end );
	
	my $gc;
	my $inter = $chr->intspan_duplicate_region->intersection($span);
	if ( $inter->is_empty ) {
		$gc = GenBoCoverageTabix->new(
			start      => $start,
			end        => $end,
			chromosome => $chr,
			patient    => $self
		);

	}
	else {
#$gc = GenBoCoverageTabix->new(start=>$start,end=>$end,chromosome=>$chr,patient=>$self);
		$gc = GenBoCoverageSamtools->new(
			start      => $start,
			end        => $end,
			chromosome => $chr,
			patient    => $self,
			raw        => 0
		);
	}
	eval {
		$self->{statistics}->{$id}->{mean} = $gc->mean;
		$self->{statistics}->{$id}->{min}  = $gc->minimum;
		$self->{statistics}->{$id}->{sum}  = $gc->sum;
	};
	if ($@) {
		return;
	}
	$self->{statistics}->{$id}->{corrected_mean} =
	  $self->{statistics}->{$id}->{mean};

	if ( $chr->name eq "X" ) {
		unless ( $chr->isPseudoAutosomal( $start, $end ) ) {
			$self->{statistics}->{$id}->{corrected_mean} *= 2
			  if $self->isMale();
		}

	}

  #warn "coucou" if $chr->name eq "X" && ($chr->isPseudoAutosomal($start,$end));

	return $self->{statistics}->{$id};

}

sub set_statistics_coverage_obj {
	my ( $self, $obj ) = @_;
	confess();
	my $id = $obj->id;
	unless ( $obj->getGenes ) {
		$self->set_statistics_coverage( $obj->getChromosome, $obj->id,
			$obj->start, $obj->end );
		return;
	}
	my $res = $obj->getGenes->[0]->get_coverage($self)
	  ->coverage( $obj->start, $obj->end );

#my $gc = GenBoCoverageTabix->new(start=>$start,end=>$end,chromosome=>$chr,patient=>$self);
	$self->{statistics}->{$id}->{mean} = $res->{mean};
	$self->{statistics}->{$id}->{min}  = $res->{min};
	$self->{statistics}->{$id}->{sum}  = $res->{sum};
	$self->{statistics}->{$id}->{corrected_mean} =
	  $self->{statistics}->{$id}->{mean};
	if ( $obj->getChromosome->name eq "X" )
	{    #&& !($obj->getChromosome->isPseudoAutosomal($obj->start,$obj->end))){
		$self->{statistics}->{$id}->{corrected_mean} *= 2 if $self->isMale();

	}
	return 1;

}

sub minimum {
	my ( $self, $obj, $capture ) = @_;
	my $id = $obj->id;
	return $self->{statistics}->{$id}->{min}
	  if exists $self->{statistics}->{$id}->{min};
	$self->set_statistics_coverage( $obj->getChromosome, $obj->id, $obj->start,
		$obj->end );
	die() unless exists $self->{statistics}->{$id}->{min};
	return $self->{statistics}->{$id}->{min};
}

sub mean {
	my ( $self, $obj, $capture ) = @_;
	my $id = $obj->id;
	return $self->{statistics}->{$id}->{mean}
	  if exists $self->{statistics}->{$id}->{mean};

	#	warn $obj->start." ".$obj->end;
	$self->set_statistics_coverage( $obj->getChromosome, $obj->id, $obj->start,
		$obj->end );
	die($self->name." ".$obj->getChromosome->name) unless exists $self->{statistics}->{$id}->{mean};
	return $self->{statistics}->{$id}->{mean};
}

sub set_count_multiplex {
	my ( $self, $multiplex, $value ) = @_;
	$self->{count_multiplex}->{$multiplex} = $value;
}

sub count_multiplex {
	my ( $self, $multiplex, $gene_id ) = @_;

	return $self->{count_multiplex}->{$multiplex}
	  if exists $self->{count_multiplex}->{$multiplex};

	my $debug;

#	my $primers = $self->getProject->getCapture->getListPrimers($multiplex,$gene_id);

	my $primers = $self->getCapture->getListPrimers($multiplex);
	my @ids     = keys %$primers;

	my @data;
	foreach my $primer_id ( keys %$primers ) {
		$self->{count_multiplex}->{$multiplex} += $self->count($primer_id);
	}

	#$self->{count_multiplex}->{$multiplex} ->{$gene_id} = sum(@data);
	return $self->{count_multiplex}->{$multiplex};

}

sub checkParseVcf {
	my $self      = shift;
	my $nbVar_vcf = 0;
	foreach my $fileName ( @{ $self->getVariationsFiles } ) {
		my $hIds = $self->parseVcfFile($fileName);
		$nbVar_vcf += scalar keys(%$hIds);
	}
	my $nbInd_vcf = 0;
	foreach my $fileName ( @{ $self->getIndelsFiles } ) {
		my $hIds = $self->parseVcfFile($fileName);
		$nbInd_vcf += scalar keys(%$hIds);
	}
	my $nbVar_obj = scalar( @{ $self->getVariations() } );
	my $nbInd_obj = scalar( @{ $self->getIndels() } );
	my $hRes;
	if ( $nbVar_obj != $nbVar_vcf ) {
		my $msg =
		  "\n\nERROR: $nbVar_obj var objects for $nbVar_vcf var in VCF...\n";
		$msg .=
"To see missing ID(s), please use debug mode for GenboPatient::checkParseVcf. Die...\n\n";
		warn $msg;
		return;
	}
	elsif ( $nbInd_obj != $nbInd_vcf ) {
		warn "Nb var obj: $nbVar_obj  -  Nb var vcf: $nbVar_vcf\n";
		warn
		  "\n\nWARN: $nbInd_obj ind objects for $nbInd_vcf ind in VCF...\n\n";
		return;
	}
	warn "[SNPs] Obj: $nbVar_obj  -  Vcf: $nbVar_vcf\n";
	warn "[INDELs] Obj: $nbInd_obj  -  Vcf: $nbInd_vcf\n";
	warn "QueryVcf parsing seems ok !\n";
	return 1;
}

sub parseVcfFile {
	my ( $self, $vcfFileName ) = @_;
	my $hIds;
	open( FILE, "zcat $vcfFileName |" );
	while (<FILE>) {
		next if ( $_ =~ /#/ );
		chomp($_);
		my $line     = $_;
		my @lTmp     = split( "\t", $line );
		my $ref_all  = $lTmp[3];
		my $var_all  = $lTmp[4];
		my @lVarAll  = split( ',', $var_all );
		my @lHeaders = split( ':', $lTmp[8] );
		my ( $iGt, $iAd );
		my $i = 0;

		foreach my $elmt (@lHeaders) {
			$iGt = $i if ( $elmt eq 'GT' );
			$iAd = $i if ( $elmt eq 'AD' );
			$i++;
		}
		my @lRes = split( ':', $lTmp[9] );
		my @lGt  = split( '/', $lRes[$iGt] );
		my $max_ad;
		if ($iAd) {
			foreach my $ad ( split( ',', $lRes[$iAd] ) ) {
				$max_ad += int($ad);
			}
		}
		next if ( $max_ad == 0 );
		foreach my $gt (@lGt) {
			next if ( $gt eq '0' );
			my $id;
			my $var_this = $lVarAll[ int($gt) - 1 ];
			if (    length($ref_all) == length($var_this)
				and length($ref_all) == 1 )
			{
				$lTmp[0] =~ s/chr//;
				$id =
				  $lTmp[0] . '_' . $lTmp[1] . '_' . $ref_all . '_' . $var_this;
			}
			else {
				$id =
				  $lTmp[0] . '_' . $lTmp[1] . '_' . $ref_all . '_' . $var_this;
			}
			$hIds->{$id} =
			  $lTmp[0] . " - " . $lTmp[1] . " - " . $ref_all . " - " . $var_all;
		}
	}
	close(FILE);
	return $hIds;
}

sub is_multiplex_ok {
	my ( $self, $multi ) = @_;
	return 1;
	return $self->{is_multiplex_ok}->{$multi}
	  if exists $self->{is_multiplex_ok}->{$multi};
	my $debug;
	my $primers = $self->getProject->getCapture->getPrimersByMultiplex($multi);
	my $score   = 0;
	my %genes;
	my $nb = 0;

	foreach my $primer (@$primers) {
		next if $primer->getChromosome()->name eq "X";
		$nb++;
		if ( abs( $primer->zscore($self) ) > 2.5 ) {

#if ($primer->cnv_score($patient) > 1.5 || $primer->cnv_score($patient) <= 0.5){
			$score++;
			map { $genes{ $_->getGene->id }++ } @{ $primer->getTranscripts };
		}
		elsif ( $primer->cnv_score($self) > 1.8 ) {

			$score++;
			map { $genes{ $_->getGene->id }++ } @{ $primer->getTranscripts };
		}
		elsif ( $primer->cnv_score($self) < 0.4 ) {

			$score++;
			map { $genes{ $_->getGene->id }++ } @{ $primer->getTranscripts };
		}

	}
	$self->{is_multiplex_ok}->{$multi} = 1;

	#if ($score  > $nb/4){
	#	$self->{is_multiplex_ok}->{$multi} =  undef;
	#}
	#if (scalar keys %genes > 3){
	#	$self->{is_multiplex_ok}->{$multi} = undef;
	#}
	return $self->{is_multiplex_ok}->{$multi};
}

sub hotspot {
	my ( $self) = @_;

	my $capture = $self->getCapture();
	my $file =$capture->hotspots_filename; 
	my $bam = $self->getBamFile();
	my $sambamba = $self->buffer()->software("sambamba");
	my @lines = `$sambamba depth base $bam -L $file 2>/dev/null`;
	

	chomp(@lines);
	#REF	POS	COV	A	C	G	T	DEL	REFSKIP	SAMPLE
	my @header = split(" ",shift @lines);
	my $res;
	foreach my $l  (@lines){
		my @data = split(" ",$l);
			my $hash ={};
		foreach (my $i=0;$i<@header;$i++){
			
			$hash->{$header[$i]} = $data[$i];
		}
		
	#	my $chromosome = $self->project->getChromosome($hash->{REF});
	#	$hash->{POS} ++;
	#	my $t = $chromosome->genesIntervalTree->fetch($hash->{POS},$hash->{POS}+1);
		my $id = $hash->{REF}.":".$hash->{POS};
		my $h = $capture->hotspots->{$id};
		$hash->{A_REF} = $h->{ref};
		$hash->{A_ALT} = $h->{alt};
		$hash->{NAME} = $h->{name};
		$hash->{ID} = $hash->{REF}.":".($hash->{POS}+1);
		$hash->{GENBO_ID} =$h->{genbo_id};
		$hash->{PROT} =$h->{protid};
		push(@{$res->{$h->{gene}}},$hash);	
	}
	return $res;
}

sub hotspot2 {
	my ( $self, $motif ) = @_;
	confess();
	my $bam = $self->getBamFile();
	my $sam = Bio::DB::Sam->new(
		-fasta => $self->project->genomeFasta,
		-bam   => $bam
	);
	my @alignments = $sam->get_features_by_location(
		-seq_id => $motif->{chromosome},
		-start  => $motif->{start},
		-end    => $motif->{end}
	);

	my $cpt = 0;
	my $res;
	for my $a (@alignments) {
		my $seqid   = $a->seq_id;
		my $starta  = $a->start;
		my $enda    = $a->end;
		my $strand  = $a->strand;
		my $ref_dna = $a->dna;
		next if $enda < $motif->{end};
		next if $starta > $motif->{start};

		my $query_dna = $a->query->dna;

		my ( $ref, $matches, $query ) = $a->padded_alignment();
		my @scores = $a->qscore;
		my @tref   = split( "", $ref );

		my $start          = -1;
		my $motif_sequence = $motif->{sequence};
		my $z              = 0;
		for ( my $i = 0 ; $i < @tref ; $i++ ) {
			next if index( "ATCG", $tref[$i] ) == -1;
			$z++;

			#next if $tref[$i] ne "T";
			my $string = substr( $ref, $i );
			$string =~ s/\-//g;

			if ( $string =~ /^$motif_sequence/ ) {
				$start = $i;
				last;
			}
		}
		next if $start == -1;
		my $end      = -1;
		my $temp_seq = "";
		$z = $start;
		for ( my $i = $start ; $i < @tref ; $i++ ) {
			next if index( "ATCG", $tref[$i] ) == -1;
			$z++;
			$temp_seq .= $tref[$i];

			#	next if $tref[$i] ne "A";
			if ( $temp_seq eq $motif_sequence ) {
				$end = $i;
				last;
			}
		}
		next if $end == -1;
		my @tquery = split( '', $query_dna );

		my $haplo = substr( $query, $start, abs( $start - $end ) + 1 );
		$haplo =~ s/\-//g;

		my $index = index( $query_dna, $haplo );
		my $find  = 0;
		foreach ( my $i = 1 ; $i < length($haplo) - 2 ; $i++ ) {
			my $m =
			  ( $scores[ $i + $index ] +
				  $scores[ $i + $index - 1 ] +
				  $scores[ $i + $index + 1 ] ) / 3;
			$find++ if $m < 15;

			#print $tquery[$i+$index]." ".$scores[$i+$index]."\n";
		}
		if ( $find < 3 ) {
			$cpt++;
			$res->{$haplo}->{p}  += 0;
			$res->{$haplo}->{m}  += 0;
			$res->{$haplo}->{dp} += 1;
			$res->{$haplo}->{p}++ if $a->strand == 1;
			$res->{$haplo}->{m}++ if $a->strand == -1;
			unless ( exists $res->{$haplo}->{align} ) {

				$res->{$haplo}->{align}->[0] =
				  substr( $ref, $start, abs( $start - $end ) + 1 );
				$res->{$haplo}->{align}->[1] =
				  "";    #substr($matches,$start,abs($start-$end)+1);
				$res->{$haplo}->{align}->[2] =
				  substr( $query, $start, abs( $start - $end ) + 1 );
				my @t1 = split( '', $res->{$haplo}->{align}->[0] );
				my @t2 = split( '', $res->{$haplo}->{align}->[2] );
				my @to;
				for ( my $i = 0 ; $i < @t1 ; $i++ ) {
					if ( $t1[$i] eq $t2[$i] ) {
						$res->{$haplo}->{align}->[1] .= "|";
						push( @to, "|" );    #616464
						$t1[$i] =
						  "<td  style='background-color:#48B04B;color:white' >"
						  . $t1[$i] . "</td>";
						$t2[$i] = "<td style='background-color:#E7E7E7'>"
						  . $t2[$i] . "</td>";
					}

					#					elsif ($t1[$i] eq '-' || $t2[$i] eq '-'){
					#						$res->{$haplo}->{align}->[1] .=" ";
					#					}
					else {
						$t1[$i] =
						  "<td style='background-color:#FF4136;color:white' >"
						  . $t1[$i] . "</td>";
						$t2[$i] =
						  "<td style='background-color:#FF4136;color:white' >"
						  . $t2[$i] . "</td>";
						$res->{$haplo}->{align}->[1] .= "X";
					}
				}
				$res->{$haplo}->{table_align}->[0] = \@t1;

				#	$res->{$haplo}->{array_align}->[1] = \@to;
				$res->{$haplo}->{table_align}->[1] = \@t2;
			}
		}

	}

	my $limit = ( $cpt * 10 ) / 100;
	foreach my $k ( keys %$res ) {
		$res->{$k}->{pourcent} = int( ( $res->{$k}->{dp} / $cpt ) * 100 );
		delete( $res->{$k} ) if ( $res->{$k}->{p} + $res->{$k}->{m} ) < $limit;

		#	delete($res->{$k}) if ($res->{$k}->{p} +  $res->{$k}->{m}) < $limit;

	}
	$motif->{results}->{ $self->name() } = $res;
	return $res;

}

sub getFamily {
	my $self = shift;
	return unless ( $self->family() );
	return $self->project->getFamily( $self->family() );
}

sub getSomaticGroup {
	my $self = shift;
	return unless ( $self->somatic_group() );
	return $self->project->getSomaticGroup( $self->somatic_group() );
}

sub getWindowGvcf {
	my ( $self, $window ) = @_;
	my $dir_out      = $self->project->getCallingPipelineDir("gvcf");
	my $dir_out_gvcf = $dir_out . "/" . $self->name;
	unless ( -e $dir_out_gvcf ) {
		system("mkdir -p $dir_out_gvcf && chmod a+rw $dir_out_gvcf");

	}
	my $outfile =
		$dir_out_gvcf . "/"
	  . $self->name() . "."
	  . $window->{chromosome} . "."
	  . $window->{start} . "."
	  . $window->{end}
	  . ".g.vcf";
	return $outfile;
}

sub getCallingPipelineDir {
	my ( $self, $type ) = @_;
	return $self->{_cdir}->{$type} if exists  $self->{_cdir}->{$type};
	my $path = $self->project->getCallingPipelineDir($type);
	return $self->project->makedir($path."/".$self->name);
}
sub getDragenDirName {
	my ( $self, $type ) = @_;
	return $self->{_dir}->{$type} if exists  $self->{_dir}->{$type};
	my $dir_out      = $self->project->project_dragen_pipeline_path_name();
	$self->{_dir}->{$type} = $dir_out. "/$type/" . $self->name;#."/".$type."/";
	return $self->{_dir}->{$type};
}
sub getDragenDir {
	my ( $self, $type ) = @_;
	$self->project->project_dragen_pipeline_path();
	my $dir = $self->getDragenDirName($type);
	unless ( -e $dir ) {
		system("mkdir -p ".$dir." && chmod a+rw ".$dir);
	}
	return $dir;
}
sub get_lmdb_primers {
	my ( $self, $mode ) = @_;
	$mode = "r" unless $mode;
	my $dir_root = $self->project->getCoverageDir();
	my $dir_out  = $dir_root . "/primers";
	unless ( -e $dir_out ) {
		system("mkdir  $dir_out");
		system("chmod a+rwx $dir_out");
	}
	my $no2 = GenBoNoSqlLmdb->new(
		dir     => $dir_out,
		mode    => $mode,
		name    => $self->name,
		vmtouch => $self->buffer->vmtouch
	);
	return $no2;
}

sub get_lmdb_cache {
	my ( $self, $mode ) = @_;
	
	$mode = "r" unless $mode;
	my $dir_out = $self->project->getCacheDir();
	unless (-e $dir_out."/".$self->name.".cache"){
		$mode = "c";
	}
	my $no2     = GenBoNoSqlLmdbCache->new(
		dir     => $dir_out,
		mode    => $mode,
		name    => $self->name.".cache",
		is_compress => 1,
		vmtouch => $self->buffer->vmtouch
	);
	if ( $mode eq "c"){
		#confess();
		$no2->put("cdate",time);
		system("chmod a+w ".$no2->filename);
	}
	return $no2;
}
 
 
sub get_lmdb_patients {
	my ( $self, $mode ) = @_;
	$mode = "r" unless $mode;
	my $dir_out = $self->project->lmdb_cache_patients_dir();
	my $no2     = GenBoNoSqlLmdb->new(
		dir     => $dir_out,
		mode    => $mode,
		name    => $self->name,
		vmtouch => $self->buffer->vmtouch
	);
	return $no2;
}

sub countPublicVariations {
	my ( $self, $chr ) = @_;
	if ($chr) {
		confess();

	}
	else {
		my $nb = 0;
		map { $nb++ if $_->isPublic() } @{ $self->getStructuralVariations };
		return $nb;

	}
}

sub countVariations {
	my ( $self, $chr ) = @_;
	if ($chr) {
		confess();
	}
	else {
		return scalar( @{ $self->getStructuralVariations } );

	}
}

sub countSubstitutions {
	my ( $self, $chr ) = @_;
	if ($chr) {
		confess();
	}
	else {
		return scalar( @{ $self->getVariations } );

	}
}

sub countIndels {
	my ( $self, $chr ) = @_;
	if ($chr) {
		my $v = $self->getVectorIndels($chr);
		return $self->countThisVariants($v);
	}
	else {
		return (
			scalar( @{ $self->getDeletions } ) +
			  scalar( @{ $self->getInsertions } ) );
	}
}

sub countHomozygote {
	my ( $self, $chr ) = @_;
	if ($chr) {
		confess();
	}
	else {
		my $nb = 0;
		map { $nb++ if $_->isHomozygote($self) }
		  @{ $self->getStructuralVariations };
		return $nb;

		#return $nb;
	}
}

sub countHeterozygote {
	my ( $self, $chr ) = @_;
	if ($chr) {
		confess();
	}
	else {
		my $nb = 0;
		map { $nb++ if $_->isHeterozygote($self) }
		  @{ $self->getStructuralVariations };
		return $nb;
	}
}

sub return_icon {
	my ($self) = @_;
	my $icon_sex = "";
	my $type = "s";
		$type = "d" if $self->isIll;
	if ( $self->isChild ) {
		
		$icon_sex =
		  qq{<img src='images/polyicons//iconfinder_baby-boy_$type.png' style='padding-right:1px'>}
		  ;    #qq{<img src='https://img.icons8.com/color/24/000000/boy.png'>};
		$icon_sex =
		  qq{<img src='images/polyicons//iconfinder_baby-girl_$type.png' style='padding-right:1px'>}
		  unless $self->isMale();
	}
	elsif ( $self->isMother ) {
		$icon_sex =
		  qq{<img src='images/polyicons/icons8-person-female-24_$type.png' style='padding-right:1px'>};
	}
	elsif ( $self->isFather ) {
		$icon_sex = qq{<img src='images/polyicons/icons8-person-24_$type.png' style='padding-right:1px'>};
	}
	return $icon_sex;
}

has small_icon => (
	is => 'ro',

	#isa		=> 'ArrayRef[Str]',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $icon;
		if ( $self->isChild ) {
				my $type = "s";
				$type = "d" if $self->isIll;
				my $sex = "girl";
				 $sex = "boy" if $self->isMale();
				$icon = "<img src='/icons/Polyicons/baby-$sex-$type.png' style='filter: drop-shadow(0px 0px 1px grey)'>"
		}
		elsif ( $self->isMother ) {
			$icon = qq{<img src='/icons/Polyicons/female-s.png' style='filter: drop-shadow(1px 1px 1px grey)'>};
			$icon = qq{<img src='/icons/Polyicons/female-d.png'style='filter: drop-shadow(1px 1px 1px grey)'>} if ($self->isIll);
		}
		elsif ( $self->isFather ) {
				$icon = qq{<img src='/icons/Polyicons/male-s.png'style='filter: drop-shadow(1px 1px 1px grey)'>};
				$icon = qq{<img src='/icons/Polyicons/male-d.png'style='filter: drop-shadow(1px 1px 1px grey)'>} if ($self->isIll);
		}
		return $icon;
	}
);

has nb_reads => (
	is => 'ro',

	#isa		=> 'ArrayRef[Str]',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $bam  = $self->getBamFile();
		my $samtools = $self->buffer->software('samtools');
		my $cmd  = qq{$samtools idxstats $bam | cut -f 1,3};
		my @sums = `$cmd`;
		chomp(@sums);
		my $h;
		foreach my $l (@sums) {
			my ( $a, $b ) = split( " ", $l );
			$h->{all} += $b;
			next if $a eq "*";
		
				#$h->{all} += $b;
			eval {
				my $chr = $self->project->getChromosome($a,1);
				if ($chr){
					
				#next unless $chr;
					$h->{ $chr->name } = $b if $chr;
				}
			};

		}
		warn $self->name();
		
		$h->{norm1} = 1/$h->{all};#/1_000_000_000;
#		warn $h->{all};
		$h->{norm} = $h->{all}/1_000_000_000;
	#	warn $h->{all};
		#if (!$self->project->isGenome){
		#	my $max = $project
		#	$h->{norm_exome} = $h->{all}/100_000_000;
		#}
		return $h;
	},
);

has normalized_reads => (
	is => 'ro',

	#isa		=> 'ArrayRef[Str]',
	lazy    => 1,
	default => sub {
	my $self = shift;
	my $value = $self->project->mean_amount_reads();
	return ($self->nb_reads->{all}/$value);
	},
	);
sub sd_value_dude {
	my ( $self, $chr_name, $start, $end ) = @_;
	my $chr_name_control = $chr_name;
 	$chr_name_control.="_".$self->sex if $chr_name eq "X";
 	confess() unless $self->project->isGenome();
	return $self->getCapture->sd_controls_dude->getMean( $chr_name_control, $start, $end );
}
sub sd_value_patient {
	my ( $self, $chr_name, $start, $end ) = @_;
	my $chr_name_control = $chr_name;
	 my $data = $self->getNoSqlDepth->getDepth( $chr_name_control, $start, $end );
	 	my $stat = new Statistics::Descriptive::Full;
	$stat->add_data(@$data);
	warn $stat->trimmed_mean(0.1);
	#warn $data->[0];
	return $stat->standard_deviation();
}
has wisecondor_bins => (
	is => 'ro',

	#isa		=> 'ArrayRef[Str]',
	lazy    => 1,
	default => sub {
		my $self = shift;
	my $dir = $self->project->getVariationsDir("wisecondor");	
	my $filein_pat = $dir."/".$self->name."_bins.bed.gz";
	return  $filein_pat
	},
);
sub ploidy_value {
	my ( $self, $chr ) = @_;
	my $dir = $self->project->getVariationsDir("wisecondor");
	my $filein_pat = $self->wisecondor_bins;
	my $tabix = $self->buffer->software("tabix");
	my $chr_name = $chr->name();
	my $mean = `$tabix  $filein_pat  $chr_name | grep -v 'NaN'| cut -f 5,6 | awk \'\{ total2+= \$2;total += \$1; count++ \} END \{  if(count>0){print total/count,total2/count}else { print 0,0} \}\' `;
	chomp($mean);
	
	return split(" ",$mean);
}

sub ploidy_value2 {
	my ( $self, $chr ) = @_;
	my $chr_name = $chr->name();
	confess() unless $self->project->isGenome();
	my $chr_name_control = $chr_name;
 	$chr_name_control.="_".$self->sex if $chr_name eq "X";
	my $norm1 =  $self->nb_reads->{norm1};
	my $no = $self->getCapture()->infos_controls_dude();
	#my $infos = $lmdb2->get("infos_control");
	my $infos = $no->get("infos_control");
	my $chr_total;
	my $nb_chr;;
	my $nb;
	my $data;
	foreach my $p (keys %{$infos}){
	if ( $chr_name eq "X"){
			next if $infos->{$p}->{sex} ne $self->sex;
	}
	my $total = $infos->{$p}->{reads}->{$chr_name};
	#warn ;
	my $norm = 1/$infos->{$p}->{reads}->{all};#/1_000_000_000;
	
	$chr_total = ( $infos->{$p}->{reads}->{$chr_name} * $norm);
	#warn $chr_total;
	push(@$data,$chr_total);
	#warn int( $infos->{$p}->{reads}->{$chr_name} * $norm);
	}
	my $stat = new Statistics::Descriptive::Full;
	$stat->add_data(@$data);
	my $m = $stat->mean();
	
	my $r2 = ($self->nb_reads->{$chr_name}*$norm1);
	#warn $m." ".$r2;
	#die();
	$chr_total = int($chr_total/scalar(keys %{$infos}));
	
	return  ($r2/$m);
	#warn Dumper $reads;
	die();
}


sub cnv_region_ratio_norm {
	my ( $self, $chr_name, $start, $end ) = @_;
	my $sum  = $self->getNoSqlDepth->getMean( $chr_name, $start, $end );
	my $total = $self->nb_reads->{norm};
	my $ratio1 = int($sum / $total);
	return $ratio1;
}

sub cnv_value_dude {
	my ( $self, $chr_name, $start, $end ) = @_;
	confess() unless $self->project->isGenome();
 	my $chr_name_control = $chr_name;
 	$chr_name_control.="_".$self->sex if $chr_name eq "X";
	# this patient
	#my $data = $self->getCapture->ratio_controls_dude->getDepth( $chr_name_control, $start, $end );
	#warn "1";
	#my $data2 = $self->getCapture->ratio_controls_dude_new->getDepth( $chr_name_control, $start, $end );
	#for (my$i =0;$i<@$data;$i++){
	#	next if int($data->[$i]) eq int($data2->[$i]);
	#	next if $data2->[$i] == -1;
	#	print $i." ".$data->[$i]." ".$data2->[$i]."\n";
	#}
	#warn Dumper @$data2;
#	warn ""
	my $ratio =	int($self->getCapture->ratio_controls_dude->getMean( $chr_name_control, $start, $end ));
	my $ratio1 = $self->cnv_region_ratio_norm($chr_name, $start, $end);
	if ($ratio1 <=5 && $ratio <10 ){
		return 0;
	} 
	if ($ratio <10) {
		my $data = $self->getCapture->ratio_controls_dude->getDepth($chr_name_control, $start, $end);
	#	warn Dumper $data;
	#	die();
		my @neg = grep {$_ == -1 || $_ == 65535} @$data;
		return -1 if scalar(@neg)> scalar(@$data)*0.5;
	}
	return -1 if $ratio  == 0 ;
	return $ratio1/$ratio;
}


sub sr_raw {
	my ( $self, $chr, $start,$debug ) = @_;
	my $samblaster = $self->buffer->software("samblaster");
	my $extractSplitReads_BwaMem =  $self->buffer->software("extractSplitReads_BwaMem");
	my $samtools = $self->buffer->software("samtools");
	my $bamfile = $self->getBamFile();
	my $chr_name = $chr->fasta_name;
	my $end = $start;
	$start -=5;
	$end += 5;
	my @res = `$samtools view  -q 30 -h $bamfile $chr_name:$start-$end |  $samblaster --excludeDups --addMateTags --maxSplitCount 1 --minNonOverlap 20 --ignoreUnmated 2>/dev/null  | grep -v "^\@" | cut -f 6  `;
	
	warn qq{$samtools view  -q 30 -h $bamfile $chr_name:$start-$end |  $samblaster --excludeDups --addMateTags --maxSplitCount 1 --minNonOverlap 20 --ignoreUnmated 2>/dev/null  | grep -v "^\@" | cut -f 6 | grep "[[:digit:]]*S"} if $debug;
	chomp(@res);
	warn Dumper @res if $debug;
	my @tr = grep {$_=~/\d+S$/} @res;
	my @tf = grep {$_=~/^\d+S/} @res;
	my $value = scalar(@res);
	return(0,0,0) if $value ==0;
	warn $value." ".scalar(@tf)." ".scalar(@tr) if $debug;
	return ($value,scalar(@tf),scalar(@tr));
	
}

sub mean_align_quality {
	my ( $self, $chr, $start, $end ) = @_;
	#my $samblaster = $self->buffer->software("samblaster");
	#my $extractSplitReads_BwaMem =  $self->buffer->software("extractSplitReads_BwaMem");

	my $samtools = $self->buffer->software("samtools");
	my $bamfile = $self->getBamFile();
	my $chr_name = $chr->fasta_name;
	$start -=5;
	$end += 5;
	my @res = ` $samtools view $bamfile $chr_name:$start-$end | awk '{sum+=\$5} END { print sum,NR}'`;
	chomp(@res);
	my $value = $res[0];
	return 0 unless $value;
	my ($a,$b) = split(" ",$res[0]);
	return 0 unless $b;
	return 0 if $b ==0;
	return int($a/$b);
}

## Validations Methods

has validations => (
	is => 'ro',

	#isa		=> 'ArrayRef[Str]',
	lazy    => 1,
	default => sub {
		my ($self) = @_;
		my $all_validations;

		#warn $self->project->validations_query()->db;
		#die($self->project->validations_query(1));

		if ( $self->project->isDefidiag ) {
			$all_validations = $self->project->validations_query()
			  ->getAllValidationsForPatient( $self, $self->project->cgi_user );
		}
		else {
			$all_validations =
			  $self->project->validations_query->getAllValidationsForPatient(
				$self);
		}
		return $all_validations;

	},
);

has validations_cnv => (
	is => 'ro',

	#isa		=> 'ArrayRef[Str]',
	lazy    => 1,
	default => sub {
		my ($self) = @_;
		my $all_validations;
		if ( $self->project->isDefidiag ) {
			$all_validations =
			  $self->project->validations_query()
			  ->getAllCnvsValidationsForPatient( $self,
				$self->project->cgi_user );
		}
		else {
			$all_validations = $self->project->validations_query()
			  ->getAllCnvsValidationsForPatient($self);
		}
		return $all_validations;

	},
);

has validation_status => (
	is => 'ro',

	#isa		=> 'ArrayRef[Str]',
	lazy    => 1,
	default => sub {
		my ($self) = @_;
		my $u;
		my $res;
		if ( $self->project->isDefidiag ) {
			$res = $self->project->validations_query()
			  ->getLatestStatusPatients( $self, $self->project->cgi_user );
		}
		else {
			$res = $self->project->validations_query()
			  ->getLatestStatusPatients($self);
		}
		return $res;
	},
);

sub getLatestValidationStatus {
	my ($self) = @_;
	return $self->validation_status->{ $self->id }->[0]
	  if exists $self->validation_status->{ $self->id };
	return undef;
}

sub getDudeFiles {
	my ($self) = shift;
	my @lFiles;
	foreach my $fileName ( @{ $self->getVariationsFiles } ) {
		push(@lFiles, $fileName) if ($fileName =~ /dude/);
	}
	return \@lFiles;
}

sub get_string_identification {
	my ($self) = @_;
	my $stv =  $self->name.":".$self->id."-".$self->identity_vigilance_vcf()."-".$self->identity_vigilance."-".$self->sex."-".$self->status;
	$stv =~ s/_/-/g;
	return $stv;
}
sub get_string_validations {
	my ($self) = @_;
		my $h1 = $self->validations;
		my $h2 = $self->validations_cnv;
		my $h3 = $self->validation_status();
		my $stv = "";
  			foreach my $k (sort keys %$h1) {
			$stv .= $k.":";
			foreach my $v (@{$h1->{$k}}){
				$stv .= $v->{validation_id}.";";
			}
	
			}
			foreach my $k (sort keys %$h2) {
				$stv .= $k.":";
				foreach my $v (@{$h2->{$k}}){
					$stv .= $v->{validation_id}.";";
				}
			}	
			foreach my $k (sort keys %$h3) {
			$stv .= $k.":";
				foreach my $v (@{$h3->{$k}}){
					$stv .= $v->{id}.";".$v->{modification_date}.";".$v->{date}.";".$v->{term};
				}
			}
			$stv =~ s/_/-/g;
			return $self->name.":".$stv if $stv;
			return $stv;
}

sub getJunctionsAnalysePath {
	my ($self) = @_;
	my $path_analisys_root = $self->getProject->get_path_rna_seq_junctions_root();
	confess("\n\nERROR: PATH $path_analisys_root not found. Die.\n\n") unless (-d $path_analisys_root);
	my $path_analisys;
	opendir my ($dir), $path_analisys_root;
	my @found_files = readdir $dir;
	closedir $dir;
	my $pat_name = $self->name();
	foreach my $file (@found_files) {
		next if $file eq '.';
		next if $file eq '..';
		if ($file =~ /$pat_name/) {
			$path_analisys = $path_analisys_root.'/'.$file;
			last;
		}
	}
	confess("\n\nERROR: PATH RNA JUNCTION not found. Die.\n\n") unless ($path_analisys);
	confess("\n\nERROR: PATH RNA JUNCTION not found. Die.\n\n") unless (-d $path_analisys);
	$path_analisys .= '/AllRes/';
	return $path_analisys;
}

sub getPatients_used_control_rna_seq_junctions_analyse {
	my $self = shift;
	my $h_desc;
	$h_desc = $self->getProject->get_hash_patients_description_rna_seq_junction_analyse();
	return unless $h_desc;
	return if not exists $h_desc->{$self->name()}->{used_ctrl};
	my @lPat;
	foreach my $other_pat (@{$self->getProject->getPatients()}) {
		push(@lPat, $other_pat) if exists $h_desc->{$self->name()}->{used_ctrl}->{$other_pat->name()};
	}
	return \@lPat;
} 

has use_not_filtred_junction_files => (
	is => 'rw',
	lazy    => 1,
	default => 1,
);

has junction_RI_file => (
	is => 'ro',
	lazy    => 1,
	default => sub {
		my ($self) = @_;
		return unless ($self->use_not_filtred_junction_files());
		my $path_analyse = $self->getJunctionsAnalysePath();
		my $path_RI_file = $path_analyse.'/AllresAll_RI.txt';
		return $path_RI_file if (-e $path_RI_file);
		return;
	},
);

has junction_RI_file_filtred => (
	is => 'ro',
	lazy    => 1,
	default => sub {
		my ($self) = @_;
		my $path_analyse = $self->getJunctionsAnalysePath();
		my $path_RI_file = $path_analyse.'/AllresRI_f.txt';
		return $path_RI_file if (-e $path_RI_file);
		return;
	},
);

has junction_SE_file => (
	is => 'ro',
	lazy    => 1,
	default => sub {
		my ($self) = @_;
		return unless ($self->use_not_filtred_junction_files());
		my $path_analyse = $self->getJunctionsAnalysePath();
		my $path_SE_file = $path_analyse.'/AllresAll_SE.txt';
		return $path_SE_file if (-e $path_SE_file);
		return;
	},
);

has junction_SE_file_filtred => (
	is => 'ro',
	lazy    => 1,
	default => sub {
		my ($self) = @_;
		my $path_analyse = $self->getJunctionsAnalysePath();
		my $path_SE_file = $path_analyse.'/AllresSE_f.txt';
		return $path_SE_file if (-e $path_SE_file);
		return;
	},
);

sub setJunctions {
	my ($self) = @_;
	my $h_ids;
	foreach my $junction (@{$self->getProject->getJunctions()}) {
		if (exists $junction->{annex}->{$self->name}) {
			$h_ids->{$junction->id()} = undef;
		}
	}
	return $h_ids;
}

sub getFiltredJunctionsRI {
	my ($self) = shift;
	my @lObj;
	foreach my $obj (@{$self->getJunctions()}) {
		next unless $obj->isRI($self);
		next unless $obj->is_filtred_results($self);
		push (@lObj, $obj);
	}
	return \@lObj;
}

sub getFiltredJunctionsSE {
	my ($self) = shift;
	my @lObj;
	foreach my $obj (@{$self->getJunctions()}) {
		next unless $obj->isSE($self);
		next unless $obj->is_filtred_results($self);
		push (@lObj, $obj);
	}
	return \@lObj;
}

sub getJunctionsRI {
	my ($self) = shift;
	my @lObj;
	foreach my $obj (@{$self->getJunctions()}) {
		next unless $obj->isRI($self);
		push (@lObj, $obj);
	}
	return \@lObj;
}

sub getJunctionsSE {
	my ($self) = shift;
	my @lObj;
	foreach my $obj (@{$self->getJunctions()}) {
		next unless $obj->isSE($self);
		push (@lObj, $obj);
	}
	return \@lObj;
}

sub getSampleProfile {
	my ($self) = shift;
	my $query    = $self->getProject()->buffer->getQuery();
	return $query->getProfileSample($self->id);
}

1;
