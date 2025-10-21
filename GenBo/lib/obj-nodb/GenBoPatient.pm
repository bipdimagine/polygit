package GenBoPatient;
use File::Basename;
use strict;
use Moo;
use Data::Dumper;
use Config::Std;
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
use QueryPbsv;
use QuerySniffles;

use List::MoreUtils qw{ natatime };
use QueryDragenSv;
use GenBoCapture;

use List::MoreUtils qw(firstidx );
extends "GenBo";

has patient_id => (
	is => 'ro',

	#isa		=> 'Str',
	reader   => 'getPatientId',
	required => 1,
);

has type => (
	is => 'ro',
	#isa		=> 'Str',
	required => 1,
);

sub isRna {
	my $self = shift;
	return  lc($self->type) eq "rna";
}

has bar_code => (
	is     => 'ro',
	reader => 'barcode',

	#isa		=> 'Str',
);

has bar_code2 => (
	is     => 'ro',
	reader => 'barcode2',

	#isa		=> 'Str',
);

has capture_id => (
	is => 'rw',

	#isa		=> 'Str',
	reader   => 'getCaptureId',
	required => 1,
);

has isGenome => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		foreach my $c ( @{ $self->getCaptures } ) {
			return 1 if lc( $c->analyse ) =~ /genome/;
		}
		return undef;
	},
);

has isExome => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		foreach my $c ( @{ $self->getCaptures } ) {
			return 1 if lc( $c->analyse ) =~ /exome/;
		}
		return undef;
	},
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

has person_id => (
	is     => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $person_id = $self->getProject->buffer->getQuery->getPersonIdFromPatientId($self->id());
		return $person_id;
	}
);

has person_infos => (
	is     => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $person_id = $self->getProject->buffer->getQuery->getPersonInfos($self->person_id());
		return $person_id;
	}
);

has control => (
	is     => 'ro',
	reader => 'is_control',
);

has identity_vigilance_vcf => ( is => 'ro', 
default => sub {
	return "";
}
);

has identity_vigilance => ( is => 'ro',
default => sub {
	return "xx";
}
 );

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
		my $path        = $self->buffer()->config_path("root","sequences");
		return $self->getRun()->fastq_dir();

	},
);


has constructor => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self        = shift;
		my $run         = $self->getRun();
		return  "short-read";
	}
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
	
		my $samtools = $self->buffer->software("samtools");
		my $intspan_capture =
		  $self->project->getChromosome('Y')->getCapturesGenomicSpan();
		my $intSpan = Set::IntSpan::Fast::XS->new(
			"6736370-6736371,6736442-6736443,6737844-6737845,2654896-2655740");
		$intSpan =  Set::IntSpan::Fast::XS->new(
			"6868329-6868330,6868401-6868402,6869803-6869804,2786855-2787699") if $self->project->version =~/HG38/;
		my $zint = $intspan_capture->intersection($intSpan);
		
		return -1 if $zint->is_empty;
		my $iter = $zint->iterate_runs();
		my @tt;
		my $chr_name = $self->project->getChromosome('Y')->fasta_name();
		my $align_file      = $self->getAlignFileName();
		my $max_mean = -50;
		while ( my ( $from, $to ) = $iter->() ) {
			my $cmd = qq{$samtools depth  $align_file -r $chr_name:$from-$to 2>/dev/null   | cut -f 3  | awk \'NR==1 \|\| \$1 > max {max=\$1} END {print max}\'};
			#my $chr = $self->project->getChromosome('Y');
			#warn Dumper $self->depth("Y",$from,$to);
			#warn $from." ".$to;
			#warn $self->meanDepth($self->project->getChromosome('Y'),$from,$to);
			my ($res) = `$cmd`;
			chomp($res);
			#my $res = $self->meanDepth($self->project->getChromosome('Y'),$from,$to);
			$res = 0 unless $res;
			

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
#		warn "start";
		my $covm = $self->coverage_SRY();
		warn  $covm;
#		warn $covm;
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
			project_id   => $self->getProject->id,
			type=> "SNP"
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
		my $lMethods = $query->getCallingMethods(
			patient_name => $self->name() . "",
			project_id   => $self->getProject->id,
			type => "SV"
		);
		return $lMethods;
	},
);


has callingCNVMethods => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my @lMethods;
		my $query    = $self->getProject()->buffer->getQuery();
		my $lMethods = $query->getCallingMethods(
			patient_name => $self->name() . "",
			project_id   => $self->getProject->id,
			type => "CNV"
		);
		return $lMethods;
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
	my $nn =0;
	my $hfiles = $self->callingFiles;
	foreach my $method ( keys %{$hfiles} ) {
		#next unless $method eq 'manta';
		#next unless $method eq 'haplotypecaller4';
		my $vcfFile = $hfiles->{$method};
		next if exists $already_parse->{ $reference->name }->{$vcfFile};
		$already_parse->{ $reference->name }->{$vcfFile}++;
		$self->{queryVcf}->{$vcfFile} = $self->getQueryVcf( $vcfFile, $method ) unless exists $self->{queryVcf}->{$vcfFile};
		my $queryVcf = $self->{queryVcf}->{$vcfFile};
		my $z        = $queryVcf->parseVcfFileForReference($reference);
		foreach my $type ( keys %$z ) {
			foreach my $id ( keys %{ $z->{$type} } ) {
			$nn++;
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
	my $o = [];
	$o = $self->myflushobjects2( $hashRes, $cursor );
	
	return $o;
}

has RnaseqSEA_SE => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self  = shift;
		my $align_method = $self->alignmentMethods->[0];
		my $path = $self->getProject->buffer()->config_path("root","project_data")."/".$self->getProjectType()."/".$self->name()."/".$self->version()."/junctions/$align_method/rnaseqsea/";
		my $file = $path.'/'.$self->name().'_SE.txt.gz';
		return $file if -e $file;
		return;
	}
);

has RnaseqSEA_RI => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self  = shift;
		my $align_method = $self->alignmentMethods->[0];
		my $path = $self->getProject->buffer()->config_path("root","project_data")."/".$self->getProjectType()."/".$self->name()."/".$self->version()."/junctions/$align_method/rnaseqsea/";
		my $file = $path.'/'.$self->name().'_RI.txt.gz';
		return $file if -e $file;
		return;
	}
);

has regtools_file => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self  = shift;
		my $file = $self->project->getJunctionsDir('regtools').'/'.$self->alignmentMethod().'/'.$self->name().'.tsv.gz';
		return $file if -e $file;
		return;
	}
); 

sub setJunctionsForReference {
	my ($self, $reference) = @_;
	my @l_obj;
	my $path = $self->getProject->get_path_rna_seq_junctions_analyse_all_res();
	my $se_file = $self->project->RnaseqSEA_SE;
	my $ri_file = $self->project->RnaseqSEA_RI;
	my ($hash_junc, $hash_junc_sj_ids);
	if ($ri_file and -e $ri_file ) {
		foreach my $hres (@{$self->getProject->getQueryJunction($ri_file, 'RI')->parse_file($reference)}) {
			my $id = $hres->{id};
			$hash_junc->{$id} = $hres;
			$hash_junc_sj_ids->{$hres->{sj_id}}->{$id} = undef;
			#push(@ares,$hres);
		}
	}
	if ($se_file and -e $se_file) {
		foreach my $hres (@{$self->getProject->getQueryJunction($se_file, 'SE')->parse_file($reference)}) {
			my $id = $hres->{id};
			$hash_junc->{$id} = $hres;
			$hash_junc_sj_ids->{$hres->{sj_id}}->{$id} = undef;
			
		}
	}
	
 	my $hSJ= {};
 	my $regtools_file = $self->regtools_file();
	if ($regtools_file and -e $regtools_file) {
		foreach my $hres (@{$self->getProject->getQueryJunction($regtools_file, 'regtools')->parse_file($reference)}) {
			my $id = $hres->{id};
			$hash_junc->{$id} = $hres;
			$hash_junc_sj_ids->{$hres->{sj_id}}->{$id} = undef;
		}
	}
	
	foreach my $id (keys %{$hash_junc} ){
		my @ps = keys %{$hash_junc->{$id}->{annex}};
		foreach my $p (@ps) {
			
			#die() unless exists $hash_junc->{$id}->{annex}->{$p}->{canonic_count};
			
			if (int($hash_junc->{$id}->{annex}->{$p}->{canonic_count}) < 10) {
				unless (exists $hash_junc->{$id}->{annex}->{$p}->{is_sj}){
					delete $hash_junc->{$id}->{annex}->{$p};
					next;
				}
				else {
					my $genes = $reference->getChromosome->getGenesByPosition($hash_junc->{$id}->{start},$hash_junc->{$id}->{end});
					my $max ={};
					$max->{$p} = [];
					foreach my $g (@$genes){
						my $tree = $g->tree_junctions();
						my $res = $tree->fetch($hash_junc->{$id}->{start},$hash_junc->{$id}->{end});
						unless(@$res) {
							$res =[];
							my $r = $tree->fetch_nearest_down($hash_junc->{$id}->{start});
							push(@$res,$r) if $r;
							 $r = $tree->fetch_nearest_up($hash_junc->{$id}->{start});
							push(@$res,$r) if $r;
							 $r = $tree->fetch_nearest_down($hash_junc->{$id}->{end});
							push(@$res,$r) if $r;
							 $r = $tree->fetch_nearest_up($hash_junc->{$id}->{end});
							push(@$res,$r) if $r;
						}
							push( @{$max->{$p}}, map{$_->{count}->{$p} } grep {exists $_->{count}->{$p}} @$res);
					}
	
					if (@{$max->{$p}}){
							$hash_junc->{$id}->{annex}->{$p}->{canonic_count} = max(@{$max->{$p}});
							$hash_junc->{$id}->{annex}->{$p}->{junc_normale_count} = $hash_junc->{$id}->{annex}->{$p}->{canonic_count};
						}
						else {
							if (exists $hash_junc->{$id}->{annex}->{$p}) {
								
							$hash_junc->{$id}->{annex}->{$p}->{canonic_count} = 0;
							$hash_junc->{$id}->{annex}->{$p}->{junc_normale_count} = 0;
							
							delete $hash_junc->{$id}->{annex}->{$p} if ($hash_junc->{$id}->{annex}->{$p}->{alt_count}+0.01) < 5;
							}
						}
					
					
					}
				}
			}
		next unless keys %{$hash_junc->{$id}->{annex}};
		
		#delete $hash_junc->{$id}->{genes};
		my $obj = $self->getProject->flushObject( 'junctions', $hash_junc->{$id});
		warn ref($obj).' -> '.$obj->id();
		push(@l_obj, $obj);
		
	}
	return \@l_obj;	
}

has tempArray => (
	is => 'rw',

	#isa		=> 'ArrayRef[Str]',
	lazy    => 1,
	default => sub {
		[],;
	},
);


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
				$varObj->patients_object()->{ $self->id() } = undef;

				#here add all method for this variant in method calling hash

				$varObj->{annex}->{ $self->id } =  $valHash->{annex}->{ $self->id };
				
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
	push( @$lRes, @{ $self->getLargeInsertions( $chrName, $start, $end ) } );
	push( @$lRes, @{ $self->getInversions( $chrName, $start, $end ) } );
	#push( @$lRes, @{ $self->getBoundaries( $chrName, $start, $end ) } );
	return $lRes;
}

sub getQueryVcf {
	my ( $self, $fileName, $method ) = @_;
	my %args;
	$args{patient} = $self;
	$args{file}    = $fileName;
	$args{method}  = $method;
	
	my $queryVcf;
	
	if ($method eq "pbsv"){
		 $queryVcf = QueryPbsv->new( \%args );
	}
	elsif ($method eq "Sniffles2"){
		 $queryVcf = QuerySniffles->new( \%args );
	}
	elsif ($method eq "Spectre"){
		 $queryVcf = QuerySniffles->new( \%args );
	}
	elsif ($method eq "dragen-sv"){
		 $queryVcf = QueryDragenSv->new( \%args );
	}
	else {
		$queryVcf = QueryVcf->new( \%args );
	}
	 
	
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

#TODO: menage a faire indels / 
has callingFiles => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		my $methods = $self->getCallingMethods();
		my $files   = {};
		foreach my $method_name (@$methods) {
			my $file = $self->_getCallingFileWithMethodName( $method_name, "variations" );
			$files->{$method_name} = $file if $file;
		}
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

has callingSRMethods => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my @lMethods;
		my $query    = $self->getProject()->buffer->getQuery();
		my $lMethods = $query->getCallingMethods(
			patient_name => $self->name() . "",
			project_id   => $self->getProject->id,
			type => "SR"
		);
		return $lMethods;
	},
);

has callingSRFiles => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		my $methods = $self->callingSRMethods();
		my $files   = {};
		foreach my $method_name (@$methods) {
			my $file = $self->_getCallingSVFileWithMethodName( $method_name,"variations" );
			$files->{$method_name} = $file if $file;
		}
		return $files;
	},
);


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
	return $self->_getFile($method,$self->callingSVFiles()->{sv});
}

sub getSRFiles {
	my $self      = shift;
	my @t  = values %{$self->callingSRFiles};
	warn( " warn I was unable to find variation vcf file SR for :  ". $self->name() ) unless scalar @t;
	return \@t;
}
sub getSRFile {
	my ( $self, $method ) = @_;
	return $self->_getFile($method,$self->callingSRFiles());
}

sub _getFile {
	my ( $self, $method, $hash ) = @_;
	if ($method) {
		confess( "can t find vcf for $method " . $self->name ) unless exists( $hash->{$method} );
		return $hash->{$method};
	}
	my @all = values %$hash;
	return "" if scalar(@all) == 0;
	return $all[0] if scalar(@all) == 1;
	confess($self->name
		  . "you have exactly "
		  . scalar(@all)
		  . " methods defined on your project "
		  . $self->getProject->name() );
}


has alignmentUrl => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		my $methods = $self->alignmentMethods();
#		die() if ( scalar( @$methods > 1 ) );
		my $method_name = $methods->[0];
		my $bam_dir     = $self->getProject->getAlignmentUrl($method_name);
		my $bamf        = $self->getAlignFileName($method_name);
		my (@t) = split( "/", $bamf );

		#warn $bam_dir;
		return $bam_dir . "/" . $t[-1];
	},
);

has bamUrl => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		return $self->alignmentUrl();
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

has SJFiles => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		#warn $self->name."----";
		my $methods = $self->alignmentMethods();
		my $files   = {};
		foreach my $method_name (@$methods) {
			my $sjdir = $self->getProject->getSJDir($method_name);
			my $file = $sjdir."/".$self->name.".SJ.tab";
			if (-e $file ){
			my $bgzip = $self->buffer->software("bgzip");
			my $tabix = $self->buffer->software("tabix");
			system("$bgzip -f $file && tabix -f -p bed $file.gz");
			die() unless -e "$file.gz";
			die() unless -e "$file.gz.tbi";
		}
			$file .=".gz";
			next unless $file;
			$files->{$method_name} = $file;
		}
		return $files;
	},

);
sub fastqFiles{
	my ($self) = @_;
	return $self->{fastq} if exists $self->{fastq};
	my $pm =   dirname (__FILE__) . "/packages/fastq/fastq.pm";
	require "$pm";
	my $files = fastq::find_file_pe($self);
	$self->{fastq} = $files;
	return $self->{fastq};
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

sub getAlignFileName {
	my ( $self, $method_name ) = @_;
	my $bam = $self->getBamFileName($method_name);
	return $bam if -e  $bam;
	my $cram = $self->getCramFileName($method_name);
	return $cram;
}

sub getBamFileName {
	my ( $self, $method_name ) = @_;
	my $bam_dir;
	if($method_name){
		 $bam_dir = $self->getProject->getAlignmentDir( $method_name );
		 
	}
	else {
	my $methods = $self->alignmentMethods();
	confess( $self->project->name." ".Dumper($methods)) if scalar(@$methods) > 1;
	 $bam_dir = $self->getProject->getAlignmentDir( $methods->[0] );
	}
	die() unless $bam_dir;
	my $bam     = $bam_dir . "/" . $self->name . ".bam";
	return $bam;
}

sub getPhysicalFilesDir {
	my ( $self, $method_name,$version ) = @_;
	return $self->{files_dir} if exists $self->{cram_dir};
	$self->{files_dir} = $self->buffer->config_path("root","processed-data")."/".$self->id."/";
	return $self->project->makedir($self->{files_dir});
	return $self->{files_dir};
}

sub getPhysicalFileName {
	my ( $self, $method_name,$version,$type ) = @_;
	confess("bam or cram" ) unless $type;
	
	my $bam_dir = $self->getPhysicalFilesDir();
	unless ($method_name){
		confess("miss method name ")  unless $type;
	}
	die() unless $bam_dir;
	my $bam     = $bam_dir . "/" . $self->name .".$version.".$method_name.".$type";
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
			warn "NO BAM FILES " . $self->name." methods :".Dumper  $self->alignmentMethods();
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

sub getAlignmentFile {
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
}

sub getBamFileForPipeline {
	my ($self,$text,$fork) = @_;
	 $fork= 5 unless $fork;
	$text =time unless $text;
	if ($self->isCram){
		my $bam_tmp = $self->project->getCallingPipelineDir("$text-".$self->name)."/".$self->name."-".time.".$text.bam";
		my $samtools = $self->buffer->software("samtools");
		my $ref = $self->project->genomeFasta();
		my $cram_prod = $self->getBamFile();
		system("$samtools view -@ $fork  -T $ref  $cram_prod  -o $bam_tmp --write-index && $samtools index $bam_tmp -@ $fork ");
		return $bam_tmp;
	}
	else {
		return $self->getBamFile();
	}
}
sub isCram {
	my ($self,$method) = @_;
	my $file = $self->getAlignmentFile($method);
	return $file =~ /\.cram/;
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
	my @lVcfFiles = values %{ $self->callingFiles() };
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
			  unless exists( $self->callingFiles()->{$method} );
		}
		confess( "can t find vcf for $method" . Dumper $self->callingFiles() )
		  unless exists( $self->callingFiles()->{$method} );
		return $self->callingFiles()->{$method};
	}
	my @all = values %{ $self->callingFiles() };
	return "" if scalar(@all) eq 0;

	confess($self->name
		  . "you have exactly "
		  . scalar(@all)
		  . " !$method! methods defined on your project "
		  . $self->getProject->name() )
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
	elsif ( $type eq 'junctions' ) {
		$dir = $project->getJunctionsDir($method_name).'/'.$self->alignmentMethod().'/';
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

has SV_extension => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $t = {
			"canavas" => ".vcf.gz",
			"manta"	=> ".vcf.gz",
			"wisecondor" => "_aberrations.bed.gz",
			"dragen-sv" =>  ".sv.vcf.gz",
			"dragen-cnv" => ".cnv.vcf.gz",
		};
		return $t;

	}
);

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
	my $ext = ".vcf.gz";
	$ext = $self->SV_extension->{$method_name} if exists $self->SV_extension->{$method_name};
	
	my $file =  $dir . "/" .$self->name() . $ext;
	confess("no vcf or bed file  $file  "
		  . $self->name()
		  . " method is "
		  . $method_name )
	  unless -e $file;
	return $file;
}

sub _getFileByExtention {
	my ( $self, $dir, $type ) = @_;
	my $name = $self->name();
	confess() unless $dir;
	my $exts = {
		"variation" => [
			"vcf.gz", "gz",  "vcf", "gff3", "sam", "txt",
			"casava", "tab", "casava2", "tsv", "tsv.gz"
		],
		"align" => ["bam","cram"],
		"SJ" => ["gz"]
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
has isKestrel => (
is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return -1 unless $self->kestrel;
		return 0 if $self->kestrel->{data}->{Confidence} eq "Negative";
		return 1;
	}
);
has isadVNTR => (
is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		
		return -1 unless $self->adVNTR;
		return 0 if $self->adVNTR->{data}->[0]->[2] eq "Negative";
		return 1;
	}
);
sub kestrel {
	my ($self) = @_;
	return $self->{kestrel} if exists $self->{kestrel};
	$self->vntyperResults();
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
		my $file = $self->project->getVariationsDir("vntyper")."/".$self->name.".json";
		unless (-e $file){
			$self->{kestrel} = undef;
			$self->{adVNTR} = undef;
			return;
		}
		my $json_text = `cat $file`;
		my $h = decode_json $json_text;
		$self->{kestrel} = $h->{kestrel};
		$self->{adVNTR} = $h->{adVNTR};
		return;
		my $json = JSON::XS->new->utf8->decode($json_text);
		$self->{kestrel} = 
		$self->{kestrel} =[];
		unless( -e $file){
			$self->{kestrel} =[];
			#$self->{adVNTR} =[];
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
			
			push(@t,["kestrel",$date,split(" ",$lines[$i])] ) ;
		}
		push(@t,["kestrel",$date]) unless @t; 
		$self->{header_adVNTR} = [];
		my $type = 
		$self->{header_adVNTR} = ["date",split(" ",$lines[$advntr+1])] if $advntr > 0;
		$self->{header_kestrel} = ["date",split(" ",$lines[$kestrel+1])]; 
		$self->{kestrel} =\@t;
		$limit2 = scalar(@lines) if $advntr >= $kestrel ;
		my @t2;
		
		for (my $i = $advntr+2 ; $i<$limit2;$i++) {
			
			
			my @tt = split(" ",$lines[$i]);
			my @aa = split("&",$tt[1]);
			my @bb = split("_",$aa[0]);
			my $repeat = $bb[1];
			$tt[0] = "Insertion" ;
			$tt[0] = "Deletetion" if $aa[0] =~ /^D/;
			push(@t2,["adVNTR",$date,$repeat,$tt[0],$tt[1],"-",$tt[2],$tt[3],$tt[4]] ) ;
		}
		$date = "-" if $advntr == -1;
		push(@{$self->{run}},"kestrel") if $kestrel;
		#push(@{$self->{run}},"adVNTR") if $advntr;
		#push(@t2,["adVNTR",$date]) unless @t2; 
		#$self->{adVNTR} =\@t2;
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

sub mean_normalize_depth {
	my ( $self, $chr_name, $start, $end ) = @_;
	$end = $start +1 if $start == $end;
	return $self->getNoSqlDepth->getMeanNormalize( $chr_name, $start, $end, $self );
}
	
sub normalize_depth {
	my ( $self, $chr_name, $start, $end ) = @_;
	my $chr   = $self->project->getChromosome($chr_name);
	my $array = $self->getNoSqlDepth->getDepth( $chr->name, $start, $end );
	my $res;
	foreach my $s1 (@$array){
		push(@$res, int($s1/$self->normalized_reads));# *100)/100);
	}
	return $res;
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
				$self->{statistics}->{$key}->{corrected_mean} *= 2 if $self->isMale();
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
#sub hotspot_genome {
#	my ( $self) = @_;
#	my $bam = $self->getBamFile();
#	my $samtools = $self->buffer()->software("samtools");
#	my $chr = $self->project->getChromosome("chrM");
#	my $start = "3243";
#	my $end = "3243";
#	my $ref ="A";
#	my $alt ="G";
#	my $region = $chr->fasta_name.":".$start."-".$end;
#	my $pileup = `$samtools mpileup $bam -r $region | cut -f 5`;
#	$pileup =~ s/[+-]\d+[ACGTNacgtn]+//g; 
#	my $nref = ($pileup =~ tr/$ref//i);
#	my $nalt = ($pileup =~ tr/$alt//i);
#	my $res
#	my $hash;
#	$hash->{A_REF} = $h->{ref};
#		$hash->{A_ALT} = $h->{alt};
#		$hash->{NAME} = $h->{name};
#		$hash->{ID} = $hash->{REF}.":".($hash->{POS}+1);
#		$hash->{GENBO_ID} =$h->{genbo_id};
#		$hash->{PROT} =$h->{protid};
#	push(@{$res->{"MT-TL1"}},$hash);	
#}



sub hotspot {
	my ($self, $fork) = @_;
	$fork = 1 if not $fork;
	my $res;
	my $file = $self->getCapture->hotspots_filename;
	return if not -e $file;
	my $bam = $self->getAlignmentFile();
	my $h = $self->getCapture->hotspots;
	
	my $tmp = "/data-beegfs/tmp/";
	my $t =time;
	
	my $nb   = int( scalar(keys %{$h}) / $fork + 1 );
	my $pm   = new Parallel::ForkManager($fork);
	
	my $iter = natatime( $nb, keys %{$h} );
	my $hrun;
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $h ) = @_;

			unless ( defined($h) or $exit_code > 0 ) {
				print
				  qq|No message received from child process $exit_code $pid!\n|;
				die();
				return;
			}
			foreach my $gid (sort keys %{$h->{res}}) {
				foreach my $hgid (@{$h->{res}->{$gid}}) {
					push(@{$res->{$gid}}, $hgid);
				}
			}
		}
	);
	$self->project->buffer->dbh_deconnect();
	my $bam_obj = Bio::DB::HTS->new(-bam => $bam,-fasta => $self->project->genomeFasta);
	
	
	
			my $t1   = time;
	while ( my @tmp = $iter->() ) {
		$self->project->disconnect();
	#	my $pid = $pm->start and next;
		
		my $hres;

		foreach my $id (@tmp) {
			my @ltmp = split('_', $h->{$id}->{genbo_id});
			my $chr_name = $ltmp[0];
			my $chr = $self->getProject->getChromosome($ltmp[0]);
			my $start = $ltmp[1];
			
			my $hash = {};
			$hash->{'A'} = 0;
			$hash->{'T'} = 0;
			$hash->{'G'} = 0;
			$hash->{'C'} = 0;
			$hash->{'INS'} = 0;
			$hash->{'DEL'} = 0;
			
			
			
			
			
			#ATTENTION: le start presente un decalage d une base (+1) apr rapport au BED. Samtools + bed = demarage a 0 donc decalage d une base
			my $end = $start+1;
			my $region = $chr->fasta_name.":".($start-1)."-".$end; 
			my $t = time;
			$bam_obj->pileup($region, sub {
	    		my ($seqid, $pos, $pileups) = @_;
				return if $pos != $start;
				warn $pos;
	    		foreach my $p (@$pileups) {
	        	next if $p->is_refskip;  # ignorer les sauts de type 'N' (RNA-seq)
	
	       		if ($p->is_del) {
	            	$hash->{'DEL'}++;
	        	} elsif ($p->indel > 0) {
	            	$hash->{'INS'}++;
	        	} else {
	        		my $aln = $p->alignment;
					my $qpos = $p->qpos;
					my $base = substr($aln->qseq, $qpos, 1);
	            	$hash->{$base}++ if $base =~ /^[ATCG]$/;
	        	}
	    		}
			});
			$hash->{'REF'} = $chr_name;
			$hash->{'COV'} = $hash->{'A'} + $hash->{'T'} + $hash->{'G'} + $hash->{'C'} + $hash->{'INS'} + $hash->{'DEL'};
			$hash->{'A_REF'} = $h->{$id}->{'ref'};
			$hash->{'A_ALT'} = $h->{$id}->{'alt'};
			$hash->{'NAME'} = $chr_name.':'.$start;
			$hash->{'ID'} = $h->{$id}->{'genbo_id'};
			$hash->{'GENBO_ID'} = $h->{$id}->{'genbo_id'};
			$hash->{'PROT'} = $h->{$id}->{'protid'};
			my $gid = $h->{$id}->{gene};
			push(@{$res->{res}->{$gid}}, $hash);	
		}
		#$pm->finish( 0, $hres );
	}
	#$pm->wait_all_children();
	return $res->{res};
}
sub _hotspot {
	my $self = shift;
	my @lRegions;
	my $file = $self->getCapture->hotspots_filename;
	return if not -e $file;
	my $h = $self->getCapture->hotspots;
	return undef if not $file;
	my $cram = $self->getBamFile();
	my $samtools = $self->buffer()->software("samtools");
	my $res;
	foreach my $id (keys %{$h}) {
		my @ltmp = split('_', $h->{$id}->{genbo_id});
		my $chr_name = $ltmp[0];
		my $chr = $self->getProject->getChromosome($ltmp[0]);
		my $start = $ltmp[1];
		
		
		#ATTENTION: le start presente un decalage d une base (+1) apr rapport au BED. Samtools + bed = demarage a 0 donc decalage d une base
		my $end = $start;
		my $region = $chr->fasta_name.":".$start."-".$end;
		my $cmd = "$samtools mpileup $cram -r $region | cut -f 5";
		my $pileup = `$cmd`;
		
		my $hash = {};
		$hash->{'A'} = 0;
		$hash->{'T'} = 0;
		$hash->{'G'} = 0;
		$hash->{'C'} = 0;
		$hash->{'INS'} = 0;
		$hash->{'DEL'} = 0;
		foreach my $s (split('', $pileup)) {
			$hash->{'A'}++ if lc($s) eq 'a';
			$hash->{'T'}++ if lc($s) eq 't';
			$hash->{'G'}++ if lc($s) eq 'g';
			$hash->{'C'}++ if lc($s) eq 'c';
			$hash->{'INS'}++ if lc($s) eq '+';
			$hash->{'DEL'}++ if lc($s) eq '-';
		}
		$hash->{'REF'} = $chr_name;
		$hash->{'COV'} = $hash->{'A'} + $hash->{'T'} + $hash->{'G'} + $hash->{'C'} + $hash->{'INS'} + $hash->{'DEL'};
		$hash->{'A_REF'} = $h->{$id}->{'ref'};
		$hash->{'A_ALT'} = $h->{$id}->{'alt'};
		$hash->{'NAME'} = $chr_name.':'.$start;
		$hash->{'ID'} = $h->{$id}->{'genbo_id'};
		$hash->{'GENBO_ID'} = $h->{$id}->{'genbo_id'};
		$hash->{'PROT'} = $h->{$id}->{'protid'};
		push(@{$res->{$h->{$id}->{gene}}},$hash);	
	}
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
sub getEpi2meDirName {
	my ( $self, $type ) = @_;
	return $self->{_dir}->{$type} if exists  $self->{_dir}->{$type};
	my $dir_out      = $self->project->project_epi2me_pipeline_path_name();
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

sub getEpi2meDir {
	my ( $self, $type ) = @_;
	$self->project->project_epi2me_pipeline_path();
	my $dir = $self->getEpi2meDirName($type);
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
sub get_lmdb_cache_beegfs {
	my ( $self, $mode ) = @_;
	confess();
#	$mode = "r" unless $mode;
#	my $dir_out = $self->project->rocks_directory_beegfs();
#	unless (-e $dir_out."/".$self->name.".cache"){
#		$mode = "c";
#	}
#	my $no2     = GenBoNoSqlLmdbCache->new(
#		dir     => $dir_out,
#		mode    => $mode,
#		name    => $self->name.".cache",
#		is_compress => 1,
#		vmtouch => $self->buffer->vmtouch
#	);
#	if ( $mode eq "c"){
#		#confess();
#		$no2->put("cdate",time);
#		system("chmod a+w ".$no2->filename);
#	}
#	return $no2;
}
sub get_lmdb_cache {
	my ( $self, $mode ) = @_;
	
	$mode = "r" unless $mode;
	my $dir_out = $self->project->rocks_directory();
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
 sub get_rocks_cache {
	my ( $self, $mode ) = @_;
	
	$mode = "r" unless $mode;
	my $dir_out = $self->project->rocks_directory();
	unless (-e $dir_out."/".$self->name.".cache"){
		$mode = "c";
	}
	my $no2     = GenBoNoSqlRocks->new(
		dir     => $dir_out,
		mode    => $mode,
		name    => $self->name.".cache",
	);
	if ( $mode eq "c"){
		#confess();
		$no2->put("cdate",time);
		system("chmod a+w ".$no2->filename);
	}
	return $no2;
}
 sub get_lmdb_cache_polydiag {
	my ( $self, $mode ) = @_;
	
	$mode = "r" unless $mode;
	my $dir_out = $self->project->getCacheDir();
	unless (-e $dir_out."/".$self->name.".polydiag.cache"){
		$mode = "c";
	}
	my $no2     = GenBoNoSqlLmdbCache->new(
		dir     => $dir_out,
		mode    => $mode,
		name    => $self->name.".polydiag.cache",
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

has small_icon_url => (
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
				$icon =  "/icons/Polyicons/baby-$sex-$type.png";
		}
		elsif ( $self->isMother ) {
			$icon = qq{/icons/Polyicons/female-s.png};
			$icon = qq{/icons/Polyicons/female-d.png} if ($self->isIll);
		}
		elsif ( $self->isFather ) {
				$icon = qq{/icons/Polyicons/male-s.png};
				$icon = qq{/icons/Polyicons/male-d.png} if ($self->isIll);
		}
		return $icon;
	}
);
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
		my $bam  = $self->getAlignmentFile();
		my $cmd;
		if ($self->isCram){
			my $f = $bam;
			$f =~ s/cram/idxstats/;
			die($f) unless -e $f;
			$cmd = qq{cat $f | cut -f 1,3}
			
		}
		else {
		my $samtools = $self->buffer->software('samtools');
		 $cmd  = qq{$samtools idxstats $bam | cut -f 1,3};
		}
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
		$h->{norm1} = 1/$h->{all};#/1_000_000_000;
#		warn $h->{all};
		$h->{norm} = $h->{all}/1_000_000_000;
		$h->{norm} = $h->{all}/50_000_000;
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
 	return 1 unless ($self->isGenome);
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
	
	my $ratio =	int($self->getCapture->ratio_controls_dude->getMean( $chr_name_control, $start, $end ));
	my $ratio1 = $self->cnv_region_ratio_norm($chr_name, $start, $end);
	if ($ratio1 <=5 && $ratio <10 ){
		return 0;
	} 
	if ($ratio <10) {
		my $data = $self->getCapture->ratio_controls_dude->getDepth($chr_name_control, $start, $end);
		my @neg = grep {$_ == -1 || $_ == 65535} @$data;
		return -1 if scalar(@neg)> scalar(@$data)*0.5;
	}
	return -1 if $ratio  == 0 ;
	return $ratio1/$ratio;
}


sub sr_raw {
	my ( $self, $chr, $start,$debug ) = @_;
	confess();
	my $samblaster = $self->buffer->software("samblaster");
	my $extractSplitReads_BwaMem =  $self->buffer->software("extractSplitReads_BwaMem");
	my $samtools = $self->buffer->software("samtools");
	my $bamfile = $self->getBamFile();
	my $chr_name = $chr->fasta_name;
	my $end = $start;
	$start -=5;
	$end += 5;
	my @res = `$samtools view  -q 30 -h $bamfile $chr_name:$start-$end |  $samblaster --excludeDups --addMateTags --maxSplitCount 1 --minNonOverlap 20 --ignoreUnmated 2>/dev/null  | grep -v "^\@" | cut -f 6  `;
#	warn qq{$samtools view  -q 30 -h $bamfile $chr_name:$start-$end |  $samblaster --excludeDups --addMateTags --maxSplitCount 1 --minNonOverlap 20 --ignoreUnmated 2>/dev/null  | grep -v "^\@" | cut -f 6};
	warn qq{$samtools view  -q 30 -h $bamfile $chr_name:$start-$end |  $samblaster --excludeDups --addMateTags --maxSplitCount 1 --minNonOverlap 20 --ignoreUnmated 2>/dev/null  | grep -v "^\@" | cut -f 6 | grep "[[:digit:]]*S"} if $debug;
	chomp(@res);
	my @tr = grep {$_=~/\d+S$/} @res;
	my @tf = grep {$_=~/^\d+S/} @res;
	my $value = scalar(@res);
	return(0,0,0) if $value ==0;
	warn $value." ".scalar(@tf)." ".scalar(@tr) if $debug;
	return ($value,scalar(@tf),scalar(@tr));
	
}

sub mean_align_quality {
	my ( $self, $chr, $start, $end ) = @_;
	confess();
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
	my $idv = $self->identity_vigilance;
	$idv = "" unless $idv;
	my $idvcf = $self->identity_vigilance_vcf();
	 $idvcf ="" unless  $idvcf;
	my $stv =  $self->name.":".$self->id."-".$idvcf."-".$idv."-".$self->sex."-".$self->status;
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

sub getSJFile {
	my ( $self, $method, $nodie ) = @_;
	if ($method) {
		if ($nodie) {
			return ""
			  unless exists( $self->SJFiles()->{$method} );
		}
		confess( "can t find vcf for $method" . Dumper $self->SJFiles() )
		 unless exists( $self->callingFiles()->{$method} );
		return $self->SJFiles->{$method};
	}
	my @all = values %{ $self->SJFiles};
	return "" if scalar(@all) eq 0;

	confess($self->name
		  . "you have exactly "
		  . scalar(@all)
		  . " !$method! methods defined on your project "
		  . $self->getProject->name() )
	  if scalar(@all) ne 1;
	return $all[0];
	
	
	
	
	return $self->{file}->{SJ} if exists $self->{file}->{SJ};
	my $path = $self->project->getSJPath();
	my $file = $path."/".$self->name.".SJ.tab";
	if (-e $file ){
		my $bgzip = $self->buffer->software("bgzip");
		my $tabix = $self->buffer->software("tabix");
		system("$bgzip $file && tabix -p bed $file.gz");
		die()unless -e "$file.gz";
	}
	$self->{file}->{SJ} = "$file.gz";
	return $self->{file}->{SJ};
}
sub getJunctionsAnalysePath {
	my ($self) = @_;
	my $path_analisys_root = $self->getProject->get_path_rna_seq_junctions_root();
	confess("\n\nERROR: PATH $path_analisys_root not found. Die.\n\n") unless (-d $path_analisys_root);
	if (-d $path_analisys_root.'/AllRes/') {
		$path_analisys_root .= '/AllRes/';
		return $path_analisys_root;
	}
	my $path_analisys;
	opendir my ($dir), $path_analisys_root;
	my @found_files = readdir $dir;
	closedir $dir;
	die;
	my $pat_name = $self->name();
	foreach my $file (@found_files) {
		next if $file eq '.';
		next if $file eq '..';
		warn $path_analisys_root.'/'.$file;
		if ($file =~ /$pat_name/) {
			$path_analisys = $path_analisys_root.'/'.$file;
			last;
		}
	}
	confess("\n\nERROR: PATH RNA ($path_analisys) JUNCTION not found. Die.\n\n") unless ($path_analisys);
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

sub update_software_version {
	my ($self,$name,$cmd,$real_version) = @_;
	my $query    = $self->getProject()->buffer->getQuery();
	my ($vid,$v) = $query->getLatestSoftwareVersion($name);
	if (lc($name) eq "dragen"){
		die("problem version $name $real_version vs $v")if ($v ne $real_version);
	}
	confess() unless $vid;
	#my ($svid,$sv) = $query->getLatestSoftwareVersionByPatient($name,$self->id);
	#return 1 if $svid eq $vid;
	$query->update_software_version($vid,$self->id,$cmd);
}
sub upd_file {
	my ($self) = @_;
	my $dir = $self->project->getSVDir("UPD");
	return $dir."/".$self->name.".json";
}

has fastq_screen_path => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir = $self->getProject->fastq_screen_path();
		$dir .= '/fastq_screen_'.$self->name().'/';
		unless (-d $dir) {
			$self->getProject->makedir($dir);
		}
		return $dir;
	},
);

has fastq_screen_file_html => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $h;
		my $fastq_screen_pat_dir = $self->fastq_screen_path();
		if (-d $fastq_screen_pat_dir) {
			opendir my ($dir), $fastq_screen_pat_dir;
			my @found_files = readdir $dir;
			closedir $dir;
			my (@lFiles);
			foreach my $file (@found_files) {
				next if $file eq '.';
				next if $file eq '..';
				next if not $file =~ /\.html/;
				return $fastq_screen_pat_dir.'/'.$file;
			}
		}
		return;
	},
);

has fastq_screen_file_html_url => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $h;
		my $fastq_screen_html = $self->fastq_screen_file_html();
		return if not $fastq_screen_html or not -e $fastq_screen_html;
		$fastq_screen_html =~ s/\/\//\//g;
		$fastq_screen_html =~ s/\/data-isilon\/sequencing\/ngs\///g;
		my $polyweb_url = $self->getProject->buffer->config->{polyweb_url}->{polyweb_NGS};
		return $polyweb_url.'/'.$fastq_screen_html;
	},
);

has fastq_screen_file_specie => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $file = $self->fastq_screen_path().'/'.$self->name().'_screen_nom_espece.txt';
		return $file if -e $file;
		return;
	},
); 

has fastq_screen_found_specie => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $file = $self->fastq_screen_file_specie();
		return if not $file or not -e $file;
		open (F, $file);
		my $specie = <F>;
		chomp($specie);
		close (F);
		return $specie;
	},
);

has fastq_screen_perc_contaminants => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $html_file = $self->fastq_screen_file_html();
		return if not $html_file;
		my $txt_file = $html_file;
		$txt_file =~s /\.html/\.txt/;
		return if not -e $txt_file;
		my $value;
		open (F, $txt_file);
		while (<F>) {
			my $line = $_;
			chomp($line);
			my @lTmp = split(' ', $line);
			next if scalar(@lTmp) < 2;
			next if $lTmp[0] ne 'Contaminants';
			my $total = $lTmp[1];
			my $reads = $lTmp[4] + $lTmp[6];
			$value = $reads / $total;
		}
		close (F);
		return $value;
	},
);

has fastq_screen_has_contaminants => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $perc_contaminants = $self->fastq_screen_perc_contaminants();
		return if not $perc_contaminants;
		return 1 if $perc_contaminants >= 1;
		return;
	},
);

1;
