package GenBoChromosome;

use strict;
use Moo;

use Carp;
use Data::Dumper;
use Config::Std;
use GenBoGenomic;
use GenBoNoSqlLmdb;
use GenBoNoSqlLmdbPipeline;  
use GenBoNoSqlRocksAnnotation;
use GenBoNoSqlRocksVariation;
use GenBoNoSqlRocksVector;
use GenBoNoSqlRocksGenome;
use GenBoNoSqlRocksPolyviewerVariant;
use Carp;

#use GenBoNoSqlRocksVariation;
#use GenBoNoSqlRocksPolyviewerVariant;
extends "GenBoGenomic";
use Storable qw/thaw retrieve/;
use Set::IntervalTree;
use List::MoreUtils qw(natatime);
use List::Util qw[min max];
use Storable qw(store retrieve freeze dclone thaw);
has name => (
	is       => 'ro',
	reader   => 'name',
	required => 1,
);
has fasta_name => (
	is       => 'ro',
	reader   => 'fasta_name',
	required => 1,
);

sub retreiveDejaVuIdSV {
	my ( $self, $start, $end, $type ) = @_;
	#return [] if exists $self->{dejavuIdTree}->{$type};
	return $self->{dejavuIdTree}->{$type}->fetch( $start, $end )
	  if exists $self->{dejavuIdTree}->{$type};
#		warn "$type ".$self->name;
		
	my $tree = $self->project->dejavuSVIntervalTree->get( "sv_dv_interval_tree", $self->name . "_" . $type );
	$self->{dejavuIdTree}->{$type} = $tree;
	
	my $tab_id = $tree->fetch( $start, $end );
#	warn " END $type ".$self->name;
	return $tab_id;
}

sub setStructuralVariants {
	my ( $self, $typeVar ) = @_;

	my @objs;
	foreach my $method ( @{ $self->project->callingSVMethods } ) {
		foreach my $p ( sort { $a->name cmp $b->name }
			@{ $self->getProject()->getPatients } )
		{
			foreach my $o (
				@{ $p->setStructuralVariantsForReference( $self, $method ) } )
			{
				$o->{references_object}->{ $self->id } = undef;
				$self->{ $o->type_object }->{ $o->id } = undef;
			}
		}
	}
}

sub setJunctions {
	my ($self) = @_;
	my $h_ids;
	my $path = $self->getProject->get_path_rna_seq_junctions_analyse_all_res();
	my $se_file = $self->project->RnaseqSEA_SE;
	my $ri_file = $self->project->RnaseqSEA_RI;
	my ($hash_junc, $hash_junc_sj_ids, $check_sj_junc);
	if ($ri_file and -e $ri_file ) {
		foreach my $hres (@{$self->getProject->getQueryJunction($ri_file, 'RI')->parse_file($self)}) {
			my $id = $hres->{id};
			$hash_junc->{$id} = $hres;
			$hash_junc_sj_ids->{$hres->{sj_id}}->{$id} = undef;
			#push(@ares,$hres);
		}
		$check_sj_junc = 1;
	}
	if ($se_file and -e $se_file) {
		foreach my $hres (@{$self->getProject->getQueryJunction($se_file, 'SE')->parse_file($self)}) {
			my $id = $hres->{id};
			$hash_junc->{$id} = $hres;
			$hash_junc_sj_ids->{$hres->{sj_id}}->{$id} = undef;
		}
		$check_sj_junc = 1;
	}
	my $regtools_file = $self->getProject->get_path_rna_seq_junctions_analyse_all_res().'/regtools/'.$self->getProject->name().'_regtools.tsv.gz';
#	my $regtools_file = '/data-isilon/bipd-src/mbras/test_GOSR2.tsv.gz';
	if (-e $regtools_file) {
		foreach my $hres (@{$self->getProject->getQueryJunction($regtools_file, 'regtools')->parse_file($self)}) {
			my $id = $hres->{id};
			$hash_junc->{$id} = $hres;
			$hash_junc_sj_ids->{$hres->{sj_id}}->{$id} = undef;
			#push(@ares,$hres);
		}
	}
	
#	warn "step 1";
 my $hSJ= {};
if ($check_sj_junc) {
	foreach my $patient (@{$self->getProject->getPatients()}) {
		my $file = $patient->getSJFile();
		next unless $file;
		my $queryJunction = QueryJunctionFile->new( file=>$file );
		$queryJunction->parse_SJ_file($patient,$self,$hash_junc,$hash_junc_sj_ids);
	}
}
	
#warn "step 2";
foreach my $id (keys %{$hash_junc} ){
	
	if ($check_sj_junc) {
		my @ps = keys %{$hash_junc->{$id}->{annex}};
		foreach my $p (@ps) {
		#die() unless exists $hash_junc->{$id}->{annex}->{$p}->{canonic_count};
		if ($hash_junc->{$id}->{annex}->{$p}->{canonic_count} < 10) {
			unless (exists $hash_junc->{$id}->{annex}->{$p}->{is_sj}){
				delete $hash_junc->{$id}->{annex}->{$p};
				next;
			}
			else {
				my $genes = $self->getGenesByPosition($hash_junc->{$id}->{start},$hash_junc->{$id}->{end});
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
	}
	
	next unless keys %{$hash_junc->{$id}->{annex}};
	
	#delete $hash_junc->{$id}->{genes};
	my $obj = $self->getProject->flushObject( 'junctions', $hash_junc->{$id});
	$h_ids->{$obj->id()} = undef;
}
return $h_ids;	

# 
#foreach my $id (keys %{} ){
#	my $id = $hres->{sj_id};
#	if ($hres->{canonic_count} == 0 ){
#		next unless  exists $hSJ->{$id};
#	}
#	
#	my $obj = $self->getProject->flushObject( 'junctions', $hres );
#	$obj->{isCanonique} = 1 if (exists $hres->{isCanonique} and $hres->{isCanonique} == 1);
#	$h_ids->{$obj->id()} = undef;
#
# unless  (exists $hSJ->{$id}){
# 	warn Dumper $hres;
# 	#die();
# }
# 
#next unless  (exists $hSJ->{$id});
#	
#}
#	my $path3 = $self->getProject->getJunctionsDir('star');
#	foreach my $patient (@{$self->getProject->getPatients()}) {
#		my $star_file = $path3.'/'.$patient->name().'.SJ.tab';
#		my $star_file_bz = $star_file.'.gz';
#		if (-e $star_file) {
#			my $cmd1 = "bgzip $star_file";
#			`$cmd1`;
#			my $cmd2 = "tabix -p bed $star_file_bz";
#			`$cmd2`;
#		}
#		next if not -e $star_file_bz;
#		foreach my $hres (@{$self->getProject->getQueryJunction($star_file_bz,'STAR')->parse_dragen_file($patient, $self)}) {
#			my $obj = $self->getProject->flushObject( 'junctions', $hres );
#			$h_ids->{$obj->id()} = undef;
#		}
#	}
#	return $h_ids;

}

has ucsc_name => (
	is => 'ro',

	#reader	=> 'getUcscName',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $name = $self->name();
		confess() if $name =~ /chr/;

		$name = 'M' if $name eq 'MT';
		$name = "chr" . $name;
		return $name;
	},
);

has karyotypeId => (
	is => 'ro',

	#reader	=> 'getUcscName',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $name = $self->name();
		return 23                if $name eq "X";
		return 24                if $name eq "Y";
		return 25                if $name eq "MT";
		return int( rand(1000) ) if ( $name =~ /-/ );
		return $name;

	},
);

has chromosome => (
	is      => 'ro',
	reader  => 'chrom',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->id();
	}
);

has isChromosomeCache => (
	is      => 'rw',
	default => undef,
);

has tree_primers => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $tree = Set::IntervalTree->new;
		foreach my $a ( @{ $self->getPrimers } ) {
			$tree->insert( $a->id, $a->start, $a->end + 1 );
		}
		return $tree;

	}
);

has tree_annotations_genes => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $tree = Set::IntervalTree->new;
		foreach my $a ( @{ $self->tree->{transcripts} } ) {
			$tree->insert(@$a);
		}
		delete $self->tree->{transcript};
		return $tree;

	}
);

has start => (
	is      => 'rw',
	default => 1,
);

has end => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->length();
	}
);

has fai => (
	is       => 'ro',
	required => 1,

);
has length => (
	is       => 'ro',
	required => 1,
);

has intspan_captures => (
	is      => 'ro',
	reader  => 'getCapturesGenomicSpan',
	lazy    => 1,
	default => sub {
		my $self       = shift;
		
		
		my $intSpanRes = Set::IntSpan::Fast::XS->new();
		
		if ($self->project->isGenome){
			$intSpanRes = Set::IntSpan::Fast::XS->new("1-".$self->length);
		}
		else {
			foreach my $capture ( @{ $self->getProject()->getCaptures() } ) {
				$intSpanRes = $intSpanRes->union( $capture->getIntSpanForChromosome($self) );
			}
		}
		return $intSpanRes;
	},
);

sub getExtentedGenomicSpan {
	my ( $self, $ext ) = @_;
	my $span = Set::IntSpan::Fast::XS->new();
	my $iter = $self->getCapturesGenomicSpan()->iterate_runs();
	while ( my ( $start, $end ) = $iter->() ) {
		$start -= $ext;
		$start = 1 if $start < 1;
		$end += $ext;
		$end = $self->length() if $end > $self->length;
		$span->add_range( $start, $end );
	}
	return $span;
}

has intspan_captures_extended => (
	is      => 'ro',
	reader  => 'getCapturesGenomicSpan_extended',
	lazy    => 1,
	default => sub {
		my $self       = shift;
		my $intSpanRes = Set::IntSpan::Fast::XS->new();
		foreach my $capture ( @{ $self->getProject()->getCaptures() } ) {
			$intSpanRes = $intSpanRes->union(
				$capture->getIntSpanForChromosome_extended($self) );
		}
		return $intSpanRes;
	},
);

has intspan_captures_referenceExtended => (
	is      => 'ro',
	reader  => 'getCapturesGenomicSpan_referenceExtended',
	lazy    => 1,
	default => sub {
		my $self       = shift;
		my $intSpanRes = Set::IntSpan::Fast::XS->new();
		foreach my $capture ( @{ $self->getProject()->getCaptures() } ) {
			$intSpanRes = $intSpanRes->union(
				$capture->getIntSpanForChromosome_referenceExtended($self) );
		}
		return $intSpanRes;
	},
);

has arrayIntSpanReference => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $foo  = Array::IntSpan->new();
		foreach my $ref ( @{ $self->getReferences } ) {
			$foo->set_range( $ref->start, $ref->end, $ref->id );
		}
		return $foo;
	},
);

has genomeFasta => (
	is      => 'rw',
	lazy    => 1,
	reader  => 'getGenomeFasta',
	default => sub {
		my $self = shift;
		return $self->getProject()->getGenomeFasta();
	},
);

sub findReferences {
	my ( $self, $obj ) = @_;

	#die($self->name." ".$self->getType) unless $obj->isGene();
	my $foo =
	  $self->arrayIntSpanReference()->get_range( $obj->start, $obj->end );

	confess() if scalar(@$foo) == 0;
	return $foo;
}

###### SET OBJECTS #####

sub getReferences {
	my ( $self, $start, $end, $nb ) = @_;
	my $hThisRefId;
	if ($nb) {
		my $ech = int( $self->length / $nb ) + 1;
		$start = 1;
		my $end;
		for ( my $i = 0 ; $i < $nb ; $i++ ) {
			$end = $start + $ech;

			foreach
			  my $refId ( keys( %{ $self->setReferences( $start, $end ) } ) )
			{
				$hThisRefId->{$refId} = undef;
				$self->references_object()->{$refId} = undef;
			}
			$start = $end + 1;
		}
		if ( $start < $self->length ) {

			foreach my $refId (
				keys( %{ $self->setReferences( $start, $self->length ) } ) )
			{
				$hThisRefId->{$refId} = undef;
				$self->references_object()->{$refId} = undef;
			}
		}

		my $lRef =
		  $self->getProject()->myflushobjects( $hThisRefId, "references" );

		my (@find) = sort { $a->start() <=> $b->start } @{$lRef};
		return \@find if scalar(@find) > 0;
		confess('\n\nNo reference found... Exit\n\n');

		#return $self->getProject()->myflushobjects($hThisRefId, "references");

	}
	if ( $start and $end ) {
		foreach my $refId ( keys( %{ $self->setReferences( $start, $end ) } ) )
		{
			$hThisRefId->{$refId} = undef;
			$self->references_object()->{$refId} = undef;
		}
		my $lRef =
		  $self->getProject()->myflushobjects( $hThisRefId, "references" );
		my (@find) =
		  grep { ( $_->start() eq $start ) and ( $_->end() eq $end ) } @{$lRef};
		return \@find if scalar(@find) > 0;
		confess('\n\nNo reference found... Exit\n\n');
	}
	foreach my $refId ( keys( %{ $self->setReferences() } ) ) {
		$hThisRefId->{$refId} = undef;
		$self->references_object()->{$refId} = undef;
	}
	return $self->getProject()->myflushobjects( $hThisRefId, "references" );
}

sub setReferences {
	my ( $self, $start, $end ) = @_;
	my $lRefIds = {};
	my $args;
	$args->{id}    = 'reference_' . $self->ucsc_name() . '_1_' . $self->length;
	$args->{name}  = $args->{id};
	$args->{start} = 0;
	$args->{end}   = $self->length;
	if ( ($start) and ($end) ) {
		$args->{id} =
		  'reference_' . $self->ucsc_name() . '_' . $start . '_' . $end;
		$args->{name}  = $args->{id};
		$args->{start} = $start;
		$args->{end}   = $end;
	}
	$args->{chromosomes_object} = { $self->id() => undef };
	$args->{project}            = $self->getProject();
	$self->getProject()->flushObject( 'references', $args );
	$lRefIds->{ $args->{id} } = undef;
	return $lRefIds;
}

sub __setReferences {
	my $self           = shift;
	my $lRefIds        = {};
	my $intSpan_refExt = $self->getCapturesGenomicSpan_referenceExtended();
	my $iter           = $intSpan_refExt->iterate_runs();
	while ( my ( $start, $end ) = $iter->() ) {
		my $args;
		$args->{id} =
		  'reference_' . $self->ucsc_name() . '_' . $start . '_' . $end;
		$args->{name}               = $args->{id};
		$args->{start}              = $start;
		$args->{end}                = $end;
		$args->{chromosomes_object} = { $self->id() => undef };
		$args->{project}            = $self->getProject();
		$self->getProject()->flushObject( 'references', $args );
		$lRefIds->{ $args->{id} } = undef;
	}
	return $lRefIds;
}

##### METHODS #####
has cytobandTree => (
	is      => 'ro',
	lazy    => 1,
	default => sub {

		my $self  = shift;
		my $tree  = Set::IntervalTree->new;
		my $dir   = $self->project->get_cytology_directory();
		
		my $tabix = Bio::DB::HTS::Tabix->new( filename => $dir . "/cytoband.bed.gz" );
		my $res = $tabix->query( $self->ucsc_name );    # if $start;
		return $tree if not $res;
		
		while ( my $line = $res->next ) {
			chomp($line);
			my ( $chr, $start, $end, $name, $color ) = split( " ", $line );
			my @v = ( "$name;$color", $start + 1, $end + 2 );
			$tree->insert(@v);
		}
		return $tree;
	}
);

sub getCytoband {
	my ( $self, $start, $end ) = @_;
	my $res = {};
	my $aa  = $self->cytobandTree->fetch( $start, $end );
	foreach my $a (@$aa) {
		my ( $k, $c ) = split( ";", $a );
		$res->{name}->{$k} = $c;
		$res->{type}->{$c} = $k;
	}
	return $res;

}

sub sequence {
	my ( $self, $start, $end, $debug ) = @_;
	my $ucscName = $self->ucsc_name();
	my $fai      = $self->fai();
	$start = 1 unless $start;

	#$end = $fai->{length} unless $end;
	$end = $self->length() unless $end;
	warn $self unless $end;
	my $l = $fai->{length_line};

	#warn "problem genomic position $start $end ". $fai->{length};

	#$end = $fai->{length} if $ucscName eq "chrM" && $end>$fai->{length};
	#$end = $fai->{length} if $ucscName eq "chrM" && $end>$fai->{length};
	$end = $self->length() if $ucscName eq "chrM" && $end > $self->length();

	#	if ($self->name() eq 'MT') { warn 'Im MT !!!'; }
	#	unless ($end) {
	#		warn 'chr name: '.$self->name();
	#		warn 'fai: '.$fai;
	#		warn Dumper $fai;die;
	#	}

#confess("problem genomic position $start $end ". $fai->{length}) if $end > $fai->{length};
	my $start_f     = $start - 1;
	my $file_start  = $fai->{start} + ($start_f) + int( $start_f / $l );
	my $dec         = abs( int( ($start_f) / $l ) - int( ($end) / $l ) );
	my $file_length = abs( $end - $start_f ) + $dec;
	my $fastaFile   = $self->getGenomeFasta();
	open( GENOME, $fastaFile ) or die("$!");
	binmode GENOME;
	seek( GENOME, $file_start, 0 );
	my $seq2;
	read( GENOME, $seq2, $file_length, 0 );
	$seq2 =~ s/\n//g;
	close GENOME;

	#	warn 'seq: '.$seq2;
	return $seq2;
}

sub getSequence {
	my ( $self, $start, $end, $debug ) = @_;
	$self->sequence( $start, $end );
}

has intergenic_intspan => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $z    = $self->project->liteIntervalTree->get_data( "genes_padding",
			$self->name );
		my $intspan = Set::IntSpan::Fast::XS->new();
		foreach my $p (@$z) {
			$intspan->add_range( $p->[1], $p->[2] - 1 );
		}
		my $intspan2 = Set::IntSpan::Fast::XS->new( "1-" . $self->length );
		my $intspan3 = $intspan2->diff($intspan);
		return $intspan3;
	}
);
has intronic_intspan => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;

		return $self->intergenic_intspan;
	}
);

sub getFastGenesIds {
	my ( $self, $obj, $debug ) = @_;

	my $start = $obj->start();
	my $end   = $obj->end();

	#	warn Dumper $self->fastGenes() if $debug;
	my $t;
	$self->intergenic_intspan();
	if ( $start == $end ) {

		$t = $self->return1position( $start, $self->fastGenes );
	}
	else {
		$t = $self->returnIntervalPosition( $start, $end, $self->fastGenes );
	}
	return $t if $self->name() eq "MT";

	if (   ( $obj->isVariant && scalar(@$t) == 0 )
		|| ( $obj->isVariant && $obj->project->isDiagnostic() ) )
	{
		$start -= 200;
		$end += 200;
		$t = $self->returnIntervalPosition( $start, $end, $self->fastGenes );
	}
	return $t;

}

sub return1position {
	my ( $self, $start, $hash ) = @_;

	my $v = $hash->lookup($start);
	my @value;

	foreach my $name (@$v) {
		my ( $a, $b ) = split( "_", $name );
		push( @value, $a . "_" . $self->name() );
	}
	return \@value;
}

sub returnIntervalPosition {
	my ( $self, $start, $end, $hash ) = @_;
	my $v = $hash->get_range( $start, $end );
	my %seen;
	foreach my $l (@$v) {
		foreach my $name ( @{ $l->[2] } ) {
			my ( $a, $b ) = split( "_", $name );
			$seen{ $a . "_" . $self->name() }++;
		}
	}
	return [ keys %seen ];
}

sub getTranscriptsByPosition {
	my ( $self, $start, $end ) = @_;
	my $ids         = $self->transcriptsIntervalTree->fetch( $start, $end + 1 );
	my $transcripts = [];
	foreach my $id (@$ids) {
		my $tr = $self->project->newTranscript($id);
		push( @$transcripts, $tr );
	}
	return $transcripts;
}

sub getGenesByPosition {
	my ( $self, $start, $end ) = @_;
	my $ids   = $self->genesIntervalTree->fetch( $start, $end + 1 );
	my $genes = [];
	foreach my $id (@$ids) {
		my $gene = $self->project->newGene($id);
		push( @$genes, $gene );
	}
	return $genes;
}

sub getFastGenesByPosition {
	my ( $self, $start, $end ) = @_;

	return $self->return1position( $start, $self->fastGenes ) if $start == $end;
	return $self->returnIntervalPosition( $start, $end, $self->fastGenes );
}

sub getIntSpanCapture {
	my ( $self, $span_limit ) = @_;

	#confess() if $patient;
	$span_limit = 0 unless defined $span_limit;
	return $self->{intspan_cpature}->{$span_limit}
	  if exists $self->{intspan_cpature}->{$span_limit};
	my $project = $self->getProject();
	my $span    = Set::IntSpan::Fast::XS->new();
	foreach my $capture ( @{ $project->getCaptures } ) {
		$span = $span->union(
			$capture->getIntSpanForChromosome( $self, $span_limit ) );
	}
	$self->{intspan_cpature}->{$span_limit} = $span;
	return $self->{intspan_cpature}->{$span_limit};

}

sub getWindowCaptureForCalling {
	my ( $self, $span_limit, $window ) = @_;
	my $res;
	my $from = 1;
	my $to   = $self->length;
	if ( $self->name eq 'MT' ) {
		my $intspan = Set::IntSpan::Fast->new( $from . "-" . $to );
		push(
			@$res,
			{
				start      => 1,
				end        => $self->length,
				intspan    => $intspan,
				ext_gvcf   => $self->ucsc_name . ".$from.$to.g.vcf",
				chromosome => $self->ucsc_name
			}
		);
		return $res;
	}
	return $self->getWindowCaptureForCallingGenome( $span_limit, $window )
	  if $self->project->isGenome();
	$window = 300_000 unless $window;

	#die();
	return $self->getWindowCaptureForCallingCapture( $span_limit, $window )
	  ;    # if $self->project->isGenome();

}

sub getWindowCaptureForCallingCapture {
	my ( $self, $span_limit, $window ) = @_;
	$span_limit = 250       unless $span_limit;
	$window     = 1_000_000 unless $window;

	#warn $span_limit;
	#warn $window;
	#die();
	my $res;

	my $intspan = $self->getIntSpanCaptureForCalling( $span_limit, $window );

	#my $nb_jobs = 50;
	#$window = int($nb/$nb_jobs);

	my $iter      = $intspan->iterate_runs();
	my $len       = 0;
	my $temp_span = Set::IntSpan::Fast->new();

	while ( my ( $from, $to ) = $iter->() ) {
		$temp_span->add_range( $from, $to );
		my $l = $self->buffer->Intspan_length($temp_span);

		if ( $l >= $window ) {
			my ( $start, $end ) = $self->buffer->Intspan_start_end($temp_span);
			push(
				@$res,
				{
					start      => $start,
					end        => $end,
					intspan    => $temp_span,
					ext_gvcf   => $self->ucsc_name . ".$start.$end.g.vcf",
					chromosome => $self->ucsc_name
				}
			);
			$temp_span = Set::IntSpan::Fast->new();

		}
	}
	unless ( $temp_span->is_empty ) {
		my ( $start, $end ) = $self->buffer->Intspan_start_end($temp_span);
		push(
			@$res,
			{
				start      => $start,
				end        => $end,
				intspan    => $temp_span,
				ext_gvcf   => $self->ucsc_name . ".$start.$end.g.vcf",
				chromosome => $self->ucsc_name
			}
		);
		$temp_span = Set::IntSpan::Fast->new();
	}
	return $res;
}

sub getWindowCaptureForCallingCapture1 {
	my ( $self, $span_limit, $window ) = @_;
	$span_limit = 250       unless $span_limit;
	$window     = 1_000_000 unless $window;

	#warn $span_limit;
	#warn $window;
	#die();
	my $res;

	warn "couou";
	my $intspan = $self->getIntSpanCaptureForCalling( $span_limit, $window );
	warn "end";
	my @array = $intspan->as_array();
	my $nb    = scalar(@array);

	#my $nb_jobs = 50;
	#$window = int($nb/$nb_jobs);
	my $iter = natatime $window, @array;

	my $final_intspan = Set::IntSpan::Fast->new();
	while ( my @tmp = $iter->() ) {

		my $start          = $tmp[0];
		my $end            = $tmp[-1];
		my $intspan_region = Set::IntSpan::Fast->new( $start . "-" . $end );
		$final_intspan->add_range( $start - 5, $end + 5 );
		my $inter = $intspan->intersection($intspan_region);
		push(
			@$res,
			{
				start      => $start,
				end        => $end,
				intspan    => $inter,
				ext_gvcf   => $self->ucsc_name . ".$start.$end.g.vcf",
				chromosome => $self->ucsc_name
			}
		);

	}

	return $res;
}

sub chunk {
	my ( $self, $size ) = @_;
	return $self->{chunk}->{$size} if exists $self->{chunk}->{$size};
	my $from = 1;
	my $to   = $self->length;
	my $regions;
	my $start;
	my $end;
	while ( $from < $to ) {
		$start = $from;
		$end   = $from + $size;
		if ( $end > $to ) {
			$end = $to;
		}
		my $intspan_region = Set::IntSpan::Fast->new( $from . "-" . $end );

		#my $inter = $intspan->intersection($intspan_region);
		#warn $inter->as_string;
		#warn $inter->as_string;
		#warn "coucou" if $inter->is_empty();
		push( @$regions,
			{ start => $from, end => $end - 1, chr => $self->name } );

#  push(@$res,{start=>$from,end=>$end,intspan=>$intspan_region,ext_gvcf=>$self->ucsc_name.".$from.$end.g.vcf",chromosome=>$self->ucsc_name}) unless $intspan_region->is_empty;

		#print chrom_name + ":" + str(region_start) + "-" + str(end)
		$from = $end;
	}

	$regions->[-1]->{end} = $self->length;
	$self->{chunk}->{$size} = $regions;
	return $regions;
}

sub getWindowCaptureForCallingGenome {
	my ( $self, $span_limit, $window ) = @_;

	#$span_limit = 250 unless $span_limit;
	$window = 1_000_000 unless $window;
	my $res;

	#my $intspan = $self->getIntSpanCaptureForCalling($span_limit,$window);
	my $regions = $self->chunk($window);
	foreach my $r (@$regions) {
		my $from           = $r->{start};
		my $end            = $r->{end};
		my $intspan_region = Set::IntSpan::Fast->new( $from . "-" . $end );
		push(
			@$res,
			{
				start      => $from,
				end        => $end,
				intspan    => $intspan_region,
				ext_gvcf   => $self->ucsc_name . ".$from.$end.g.vcf",
				chromosome => $self->fasta_name
			}
		) unless $intspan_region->is_empty;
	}
	return $res;
}

sub getWindow {
	my ( $self, $from, $end, $span_limit ) = @_;
	my $intspan        = $self->getIntSpanCaptureForCalling($span_limit);
	my $intspan_region = Set::IntSpan::Fast->new( $from . "-" . $end );
	my $inter          = $intspan->intersection($intspan_region);
	return {
		start      => $from,
		end        => $end,
		intspan    => $inter,
		ext_gvcf   => $self->ucsc_name . ".$from.$end.g.vcf",
		chromosome => $self->ucsc_name
	};

}

sub getIntSpanCaptureForCalling {
	my ( $self, $span_limit ) = @_;
	#confess() if $patient;
	my $span_chr = Set::IntSpan::Fast::XS->new( "1-" . $self->length );
	my $project  = $self->getProject();
	my $span     = Set::IntSpan::Fast::XS->new();
	foreach my $capture ( @{ $project->getCaptures } ) {
		$span = $span->union(
			$capture->getIntSpanForChromosome( $self, $span_limit ) );
	}
	if ( $project->isDiagnostic ) {
		#my $span3 = Set::IntSpan::Fast::XS->new() ;
		my $trs;
		foreach my $capture ( @{ $project->getCaptures } ) {
			push( @$trs, @{ $capture->transcripts_name() } );
		}
		my $span_transcript = Set::IntSpan::Fast::XS->new();
		
		my $nb_ok = 0;
		my $nb_ok_1 = 0;
		my $nb_ok_2 = 0;
		my $nb_error = 0;
		
		foreach my $tr (@$trs) {
			my $t = $project->newTranscript($tr);
			next if $t->getChromosome()->name ne $self->getChromosome()->name;
			my $g = $t->getGene();

#$span_transcript->{$t->getChromosome()->name} =  Set::IntSpan::Fast::XS->new() unless exists $span_transcript->{$tr->getChromosome()->name};
			$span_transcript->add_range( $g->start() - 1000, $g->end + 1000 );
		}
		$span = $span->union($span_transcript);

		my $primers      = $self->getPrimers();
		my $span_primers = Set::IntSpan::Fast::XS->new();
		foreach my $primer (@$primers) {

		#next if $primer->getChromosome()->name ne $self->getChromosome()->name;

			my $chr   = $primer->getChromosome();
			my $start = $primer->start - 1000;
			my $end   = $primer->end + 1000;

			$span_primers->add_range( $start, $end );

		}

		#	warn $span_primers->as_string();
		$span = $span->union($span_primers);
	}
	else {
		#my $span3 = $self->project->liteIntervalTree->get_intspan( "transcripts_padding",$self->name );
		#$span = $span->union($span3);
	}
	return $span_chr->intersection($span);

}

sub setPrimersForPatient {
	my ( $self, $patients ) = @_;
	confess();
	my $captures = [];
	foreach my $patient (@$patients) {
		push( @$captures, $patient->getCapture() );
	}

	#$captures = [$patient->getCapture()] if ($patient) ;
	my %hchrs;
	my @objs;
	my $cv;
	foreach my $c (@$captures) {
		next if exists $cv->{ $c->id };
		$cv->{ $c->id }++;
		foreach my $p ( @{ $c->parsePrimersForChromosome($self) } ) {
			$self->{primers_object}->{ $p->id } = 0;
		}

		#	$self->{primers_object}->{}
		#	last;
	}

	#	warn scalar (keys %{$self->{primers_object}});
	
	return $self->{primers_object};

	#return \@objs;

}

sub _constructPrimersFromIntspan {
	my ($self,$intspan,$captures) = @_;
	my $iter = $intspan->iterate_runs();
	my $primers;
	while (my ( $from, $to ) = $iter->()) {
		my $size = 	abs($from-$to);
		foreach (my $i=$from;$i<$to;$i+=1000){
			my $start = $i;
			my $end = $i+1000;
			$end=$to if $end > $to;
			next if abs($start-$to) <50;
			my $hpos;
	 
			
			$hpos->{chromosomes_object}->{$self->id} =undef;
			$hpos->{gstart} = $start - 15;
			$hpos->{gend} = $end+15;
			$hpos->{start} = $start;
			$hpos->{end}   = $end;
			confess() if ($hpos->{start} > $hpos->{end});
			($hpos->{end}, $hpos->{start}) = ($hpos->{start},$hpos->{end})  if ($hpos->{start} > $hpos->{end});
			$hpos->{start_forward}   = $start -15;
			$hpos->{end_forward}   = $start;
			$hpos->{id}="primer".$self->name."_$start-$end";
			$hpos->{id2}="primer".$self->name."_$start-id";
			$hpos->{name}= $self->name."_$start-$end";
			$hpos->{multiplex}= 1;
			$hpos->{length}   = ($end-$start)+1;
			$hpos->{intspan} =  Set::IntSpan::Fast::XS->new($hpos->{start}."-".$hpos->{end} );
		
			$hpos->{cnv} ={};
			foreach my $c (@$captures){
				$hpos->{$c->type_object}->{$c->id} =undef;
			}
			#$self->{primer_size}->{$hpos->{id}} = abs($startf-$startr)-1; 
			push(@$primers,$hpos);
		}
		}
	my $objs = $self->getProject()->flushObjects("primers",$primers);
	foreach my $o (@$objs){
		$self->{primers_object}->{ $o->id } = undef;
		$o->{chromosomes_object}->{$self->id} = undef;	
	}
	return $objs;
		
	
}


sub setPrimers {
	my ( $self, $patient ) = @_;
	my $captures = $self->project->getCaptures();
	#$captures = [ $patient->getCapture() ] if ($patient);
	if (scalar(@$captures) > 1) {
		my $intspan;# = Set::IntSpan::Fast::XS->new();
		foreach my $c (@$captures) {
			my $span = $c->getIntSpanForChromosome($self,50);
			$intspan = $span unless $intspan;
			$intspan = $intspan->intersection($span);
		}
		$self->_constructPrimersFromIntspan($intspan,$captures);
		return $self->{primers_object};
	}
	my %hchrs;
	my @objs;
	#if ($self->project->isExome){
	#	warn "EXOME";
#		warn $self->project->getCaptures()->[-1]->name;
#		$captures = [$self->project->selectCapture ] ;
#	}
	my $tt =0;
	foreach my $c (@$captures) {
		foreach my $p ( @{ $c->parsePrimersForChromosome($self) } ) {
			$tt ++;
			$self->{primers_object}->{ $p->id } = 0;
		}

		#	$self->{primers_object}->{}
		#	last;
	}
	#warn scalar(keys %{$self->{primers_object}});
	return $self->{primers_object};

}

sub set_ordered_primers {
	my ($self) = @_;
	return $self->{ordered_primers} if exists $self->{ordered_primers};
	my $objs = $self->getProject()->myflushobjects( $self->primers_object(), "primers" );
	my @primers = map { $_->id } sort { $a->start <=> $b->start } @{$objs};

	$self->{ordered_primers} = \@primers;
	return $self->{ordered_primers};
}

sub getPrimers {
	my ($self) = @_;
	return $self->getProject()->myflushobjects( $self->set_ordered_primers(), "primers" );
}

sub lmdb_image_transcripts_uri {
	my ( $self, $mode ) = @_;
	my $fname = $self->name;
	my $hindex = "image_transcripts_uri".$fname;
	return $self->{lmdb}->{$hindex} if exists $self->{lmdb}->{$hindex};

	my $dir_out = $self->getProject->getCoverageDir() . "/lmdb_images_uri/";
	system("mkdir $dir_out && chmod a+rwx $dir_out") unless -e $dir_out;

	$mode = "r" unless $mode;
	if ( $mode eq "r" && -z "$dir_out/$fname" ) {
		my $no2 = GenBoNoSqlLmdb->new(
			dir         => $dir_out,
			mode        => "c",
			name        => $fname,
			is_compress => 1,
			vmtouch     => $self->buffer->vmtouch
		);
		$no2->create();
		$no2->close();

	}

	my $no2 = GenBoNoSqlLmdb->new(
		dir         => $dir_out,
		mode        => $mode,
		name        => $fname,
		is_compress => 1,
		vmtouch     => $self->buffer->vmtouch
	);

	$no2->clean_files() if ( $mode eq "c" );
	$no2->create()      if ( $mode eq "c" );
	 $self->{lmdb}->{$hindex} = $no2;
	return  $self->{lmdb}->{$hindex} ;

}

## lmdb cache

has lmdb_cache_dir => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self     = shift;
		my $dir_root = $self->project->lmdb_cache_dir();
		my $dir_out  = $dir_root . "/" . $self->name();
		unless ( -e $dir_out ) {
			system("mkdir  $dir_out");
			system("chmod a+rwx $dir_out");
		}
		return $dir_out;
	}
);

sub get_lmdb_ncboost_chromosomes {
	my ( $self, $mode ) = @_;
	$mode = "r" unless $mode;
	
	my $hindex = "ncboost";
	return $self->{lmdb}->{$hindex} if exists $self->{lmdb}->{$hindex};
	my $no = GenBoNoSqlLmdb->new(
		dir         => $self->project->lmdb_ncboost_path() . '/vector/',
		mode        => $mode,
		name        => $self->id() . '.vector',
		is_compress => 1,
		vmtouch     => $self->buffer->vmtouch
	);
	$no->clean_files() if ( $mode eq "c" );
	$self->{lmdb}->{$hindex} = $no;
	return $no;
}

sub get_lmdb_ncboost_chromosomes_scores {
	my ( $self, $mode ) = @_;
	$mode = "r" unless $mode;
	my $hindex = "ncboost_score";
	return $self->{lmdb}->{$hindex} if exists $self->{lmdb}->{$hindex};
	my $no = GenBoNoSqlLmdbScore->new(
		dir         => $self->project->lmdb_ncboost_path(),
		mode        => $mode,
		name        => $self->id() . '.uc',
		is_compress => 1,
		vmtouch     => $self->buffer->vmtouch
	);
	$no->clean_files() if ( $mode eq "c" );
	$self->{lmdb}->{$hindex} = $no;
	return $no;
}

sub getGenBoBitVectorNcboostCategory {
	my ( $self, $category ) = @_;
	confess("\n\nERROR: need category argument in method "
		  . ref($self)
		  . "::get_lmdb_ncboost_vector_category(). Die...\n\n" )
	  unless ($category);
	confess("\n\nERROR: $category not defined in genbo config in method "
		  . ref($self)
		  . "::get_lmdb_ncboost_vector_category(). Die...\n\n" )
	  unless (
		exists $self->project->buffer->config->{scaled_score_ncboost}
		->{$category} );
	my $no = $self->get_lmdb_ncboost_chromosomes( $self, 'r' );
	my $v  = $no->get($category);
	return $v;
}

sub get_lmdb_hotspot {
	my ( $self, $type ) = @_;
	return $self->buffer->get_lmdb_database( "hotspot", $self->name, $type );
}

has intspan_pseudo_autosomal => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;

		return Set::IntSpan::Fast::XS->new()
		  if $self->name ne "X" && $self->name ne "Y";
		my $pos;

		if ( $self->project->getVersion =~ /HG19/ ) {
			if ( $self->name eq "X" ) {
				$pos = "60001-2699520,154931044-155260560";
			}
			elsif ( $self->name eq "Y" ) {
				$pos = "10001-2649520,59034050-59363566";
			}
		}
		elsif ( $self->project->getVersion =~ /HG38/ ) {
			if ( $self->name eq "X" ) {
				$pos = "10001-2781479,155701383-156030895";
			}
			elsif ( $self->name eq "Y" ) {
				$pos = "10001-2781479,56887903-57217415";
			}
		}
		else {
			confess();
		}

		return Set::IntSpan::Fast::XS->new($pos);
	}
);

sub isAutosomal {
	my ( $self, $start, $end ) = @_;

	#return 1 unless $self->name ne "X" or $self->name ne "Y";
	my $intspan = $self->intspan_pseudo_autosomal;
	return 1 if $intspan->is_empty;
	my $intspan2 = Set::IntSpan::Fast::XS->new( $start, $end );
	my $intspan3 = $intspan->intersection($intspan2);
	return $intspan3->is_empty;
}

sub ploidy {
	my ( $self, $sample, $start, $end ) = @_;
	if ( ( $self->name eq "Y" or $self->name eq "X" ) && $sample->sex == 1 ) {
		return 1 if $self->isAutosomal( $start, $end );
	}
	return 2;

	#return 1 unless $self->name ne "X" or $self->name ne "Y";
}

sub isPseudoAutosomal {
	my ( $self, $start, $end ) = @_;
	$end = $start + 1 unless $end;
	return if $self->name ne "X" && $self->name ne "Y";
	my $intspan = $self->intspan_pseudo_autosomal;
	return 1 if $intspan->is_empty;
	my $intspan2 = Set::IntSpan::Fast::XS->new( $start, $end );
	my $intspan3 = $intspan->intersection($intspan2);
	return !( $intspan3->is_empty );

}

sub _get_lmdb {
	my ( $self, $modefull, $fname, $no_index ) = @_;
	$modefull = "r" unless $modefull;
	my $hindex = $modefull."_".$fname;
	return $self->{lmdb}->{$hindex} if exists $self->{lmdb}->{$hindex};
	my ( $mode, $pipeline );
	( $mode, $pipeline ) = split( '', $modefull ) if ($modefull);
	die() unless $fname;
	$mode = "r" unless $mode;
	my $dir_out = $self->lmdb_cache_dir();

	if ( $mode eq "r" && -z "$dir_out/$fname" ) {
		return undef;
	}
	my $no2;
	if ( $modefull and $pipeline eq 'p' ) {

		#$dir_out = $self->project->lmdb_cache_variations_dir();
		confess() if $mode ne 'c';
		confess();
		$dir_out = $self->project->lmdb_pipeline_dir() . "/$fname/";
		$no2     = GenBoNoSqlLmdbPipeline->new(
			dir_prod    => $self->lmdb_cache_dir(),
			dir         => $dir_out,
			mode        => $mode,
			is_index    => 1,
			name        => $fname,
			is_compress => 1,
			vmtouch     => $self->buffer->vmtouch
		) unless $no_index;
		$no2 = GenBoNoSqlLmdbPipeline->new(
			dir_prod    => $self->lmdb_cache_dir(),
			dir         => $dir_out,
			mode        => $mode,
			name        => $fname,
			is_compress => 1,
			vmtouch     => $self->buffer->vmtouch
		) if $no_index;
		$self->{lmdb}->{$hindex} = $no2;
		return $no2;
	}
	else {
		$no2 = GenBoNoSqlLmdb->new(
			dir         => $dir_out,
			mode        => $mode,
			name        => $fname,
			is_compress => 1,
			is_index    => 1,
			vmtouch     => $self->buffer->vmtouch
		) unless $no_index;
		$no2 = GenBoNoSqlLmdb->new(
			dir         => $dir_out,
			mode        => $mode,
			name        => $fname,
			is_compress => 1,
			vmtouch     => $self->buffer->vmtouch
		) if $no_index;
		my $sname = $no2->filename;

	}
	$no2->clean_files() if ( $mode eq "c" );
	$no2->create()      if ( $mode eq "c" );
	$self->{lmdb}->{$hindex} = $no2;
	return $no2;
}

sub lmdb_hash_variants {
	my ( $self, $mode ) = @_;
	confess();
	if ( $mode eq "close" ) {
		$self->{lmdb}->{hash_variant}->close
		  if exists $self->{lmdb}->{hash_variant};
		delete $self->{lmdb}->{hash_variant};
		return;
	}
	return $self->{lmdb}->{hash_variant}
	  if exists $self->{lmdb}->{hash_variant};
	$self->{lmdb}->{hash_variant} =
	  $self->_get_lmdb( $mode, "hashes_variants", 1 );
	return $self->{lmdb}->{hash_variant};
}

sub lmdb_polyviewer_mask {
	my ( $self, $mode ) = @_;
	my $hname = "mask";
	confess();
	if ( $mode && $mode eq "close" ) {
		$self->{lmdb}->{$hname}->close if exists $self->{lmdb}->{$hname};
		delete $self->{lmdb}->{$hname};
		return;
	}
	return $self->{lmdb}->{$hname} if exists $self->{lmdb}->{$hname};
	die() unless $mode;
	if ( $mode eq "close" ) {
		$self->{lmdb}->{hash_variant}->close if exists $self->{lmdb}->{$hname};
		delete $self->{lmdb}->{$hname};
		return;
	}
	return $self->{lmdb}->{$hname} if exists $self->{lmdb}->{$hname};
	$self->{lmdb}->{$hname} = $self->_get_lmdb( $mode, "$hname.polyviewer", 1 );
	return $self->{lmdb}->{$hname};
}

sub _lmdb_polyviewer {
	my ( $self, $type, $patient, $mode ) = @_;
	my $hname;
	unless ( ref($patient) =~ /Patient/ ) {
		$hname = $patient . ".$type";
	}
	else {
		$hname = $patient->name . ".$type";
	}
	if ( $mode && $mode eq "close" ) {
		my @t = grep {$_=~ /$hname/ } keys %{$self->{lmdb}};
		foreach my $k (@t){
			$self->{lmdb}->{$k}->close if exists $self->{lmdb}->{$k};
			delete $self->{lmdb}->{$k};
		}
	
		return;
	}
	return $self->{lmdb}->{$hname} if exists $self->{lmdb}->{$hname};
	$self->{lmdb}->{$hname} = $self->_get_lmdb( $mode, "$hname.polyviewer", 1 );
	return $self->{lmdb}->{$hname};

}

sub lmdb_polyviewer_variants {
	my ( $self, $patient, $mode ) = @_;
	my $no = $self->_lmdb_polyviewer( "variants", $patient, $mode );
	#$no->test(1);
	return $no;

}
###################
## Vector
###################




# !!!!!!!! rocks 
sub rocks_vector {
	my ( $self, $mode ) = @_;
	$mode = "r" unless $mode;
	my $name = $self->name."-vector-".$mode.$$;
	
	return $self->project->{rocks}->{$name} if exists $self->project->{rocks}->{$name};
	 $self->project->{rocks}->{$name} = $self->flush_rocks_vector($mode);
	 return $self->project->{rocks}->{$name};

}

sub flush_rocks_vector {
	my ( $self, $mode ) = @_;
	$mode = "r" unless $mode;
	return GenBoNoSqlRocksVector->new(chromosome=>$self->name,dir=>$self->project->rocks_directory("vector"),mode=>$mode,name=>$self->name);
}

#------ lmdb 
sub get_lmdb_patients {
	my ( $self, $modefull ) = @_;
	my $hindex = "patient_";
	$hindex = "patient_".$modefull if ($modefull);
	return $self->{lmdb}->{$hindex} if exists $self->{lmdb}->{$hindex};
	my $dir_out = $self->project->lmdb_cache_patients_dir();
	my ( $mode, $pipeline );
	( $mode, $pipeline ) = split( '', $modefull ) if ($modefull);
	$mode = "r" unless $mode;
	my $no2;
	if ( $pipeline and $pipeline eq 'p' ) {

		confess() if $mode ne 'c';
		my $dir_out2 = $self->project->lmdb_pipeline_dir() . "/lmdbd_patients/";
		$no2 = GenBoNoSqlLmdbPipeline->new(
			dir_prod    => $dir_out,
			dir         => $dir_out2,
			mode        => $mode,
			name        => $self->name,
			is_compress => 1
		);

	}
	else {
		$no2 = GenBoNoSqlLmdb->new(
			dir         => $dir_out,
			mode        => $mode,
			name        => $self->name,
			is_compress => 1,
			vmtouch     => $self->buffer->vmtouch
		);
	}
	$self->{lmdb}->{$hindex} = $no2;
	return $no2;
}

############
# ROCKS 
##############





sub rocks_polyviewer_directory {
	my ($self) = @_;
	return $self->project->rocks_directory("polyviewer");
}

sub rocks_polyviewer_variants {
	my ( $self, $mode ) = @_;
	$mode ="r" unless $mode;
	my $hmode = $mode.$$;
	return $self->project->{rocks}->{"polyviewer".$hmode} if exists $self->project->{rocks}->{"polyviewer".$hmode};
	 $self->project->{rocks}->{"polyviewer".$hmode} = GenBoNoSqlRocksPolyviewerVariant->new(dir=>$self->rocks_polyviewer_directory,mode=>$mode,name=>$self->name.".polyviewer_variant");
	
	  #$self->{rockssss}->{"polyviewer".$mode}->rocks;
	 return $self->project->{rocks}->{"polyviewer".$hmode};

}

sub rocks_polyviewer_variants_genes {
	my ( $self,$mode ) = @_;
	return $self->project->{rocks}->{"variants_genes".$mode} if exists $self->project->{rocks}->{"variants_genes".$mode};
	 $self->project->{rocks}->{"variants_genes".$mode} = GenBoNoSqlRocks->new(dir=>$self->rocks_polyviewer_directory,mode=>$mode,name=>$self->name.".polyviewer_variant_genes");
	 return $self->project->{rocks}->{"variants_genes".$mode};

}

sub rocks_polyviewer_patient {
	my ( $self,$patient, $mode ) = @_;
	$mode ="r" unless $mode;
	confess();
	my $hmode = $patient->id."-".$mode ."$$";
	return $self->project->{rocks}->{"polyviewer".$hmode} if exists $self->project->{rocks}->{"polyviewer".$hmode};
	my $dir = $self->rocks_polyviewer_directory."patients/".$patient->name."/";
	$self->project->makedir($dir);
	 $self->project->{rocks}->{"polyviewer".$hmode} = GenBoNoSqlRocks->new(dir=>$dir,mode=>$mode,name=>$self->name.".polyviewer_patient");
	  #$self->{rockssss}->{"polyviewer".$mode}->rocks;
	 return $self->project->{rocks}->{"polyviewer".$hmode};

}


sub get_polyviewer_variants_genes {
		my ( $self,$patient,$id ) = @_;
		if ($self->project->isRocks){
			return  $self->rocks_polyviewer_patient($patient)->get($id);#$self->rocks_polyviewer_patient($patient)->get_patient_variant_genes($id);
		}
		else {
			warn $self->lmdb_polyviewer_variants_genes($patient,"r")->filename;
			warn $id;
			return $self->lmdb_polyviewer_variants_genes($patient,"r")->get($id);
		}
}

sub lmdb_polyviewer_variants_genes {
	my ( $self, $patient, $mode ) = @_;
	$mode ="r" unless $mode;
	
	return $self->_lmdb_polyviewer( "variants-genes", $patient, $mode );
}


sub get_polyviewer_genes {
	my ( $self, $patient, $gene_id) = @_;
	
	if ($self->project->isRocks){
		return $self->rocks_polyviewer_genes->get_patient_gene($patient,$gene_id);
	}
	confess();
	return $self->lmdb_polyviewer_genes($patient)->get($gene_id);
	
}

sub rocks_polyviewer_genes {
	my ($self) = @_;
	return $self->rocks_polyviewer_variants();
}

sub lmdb_polyviewer_genes {
	my ( $self, $patient, $mode ) = @_;
	confess() if $self->project->isRocks();
	return $self->_lmdb_polyviewer( "hgenes", $patient, $mode );
}

sub getPolyviewer_score {
	my ( $self, $patient, $key, $mode ) = @_;
	return $self->_lmdb_polyviewer( "variants-genes", $patient, $mode )
	  ->get($key);

	#return $no->
	my $name    = $patient->name;
	my $dir_out = $self->lmdb_cache_dir();
	if ( -e $dir_out . "/$name.genes.polyviewer" ) {
		return $self->lmdb_polyviewer_genes($patient)->get($key);
	}
	else {
		$key = $key . "@" . $name;
		return $self->lmdb_hash_variants($mode)->get($key);
	}
}

sub get_lmdb_variations_object {
	my ( $self, $mode ) = @_;
	return $self->_get_lmdb( $mode, "variations_objects" );
}



sub lmdb_variations {
	my ( $self, $mode ) = @_;
	$mode = "r" unless $mode;
	if ( $mode eq "close" ) {
		$$self->{lmdb}->{variations}->close();
		delete $self->{lmdb}->{variations};
	}
	if ( exists$self->{lmdb}->{variations} ) {
		return $self->{lmdb_variations};
	}
	$self->{lmdb}->{variations} = $self->get_lmdb_variations($mode);
	
}
#sub get_lmdb_annex {
#	my ( $self, $mode ) = @_;
#	return $self->{lmdb}->{annex} if exists $self->{lmdb}->{annex};
#	my $dir_out = $self->project->lmdb_cache_variations_dir();
#	my $no2 = GenBoNoSqlLmdb->new(
#			dir         => $dir_out,
#			mode        => $mode,
#			name        => $self->name.".annex",
#			is_compress => 3,
#			is_integer=>1,
#			vmtouch     => $self->buffer->vmtouch
#			);
#	$self->{lmdb}->{annex} = $no2;
#	return $no2;
#}

sub get_old_lmdb_variations {
	my ($self,$mode,$dir_out) = @_;
		return  GenBoNoSqlLmdb->new(
			dir         => $self->project->lmdb_cache_variations_dir(),
			mode        => $mode,
			is_index    => 1,
			name        => $self->name,
			is_compress => 1,
			vmtouch     => $self->buffer->vmtouch
			);
}


sub get_lmdb_junctions_canoniques {
	my ( $self, $modefull) = @_;
	my $hindex = "dv_junctions_canoniques_";
	$hindex = "dv_junctions_canoniques_".$modefull if ($modefull);
	return $self->{lmdb}->{$hindex} if exists $self->{lmdb}->{$hindex};
	$modefull = "r" unless $modefull;
	my ( $mode, $pipeline ) = split( '', $modefull );
	my $dir_out = $self->project->get_dejavu_junctions_path('canoniques');
	return  GenBoNoSqlLmdb->new(
		dir         => $dir_out,
		mode        => $mode,
		is_index    => 1,
		name        => $self->name,
		is_compress => 1,
		vmtouch     => $self->buffer->vmtouch
	);
	return $self->{lmdb}->{$hindex};
}



sub get_lmdb_variations {
	my ( $self, $modefull,$rocks) = @_;
	my $hindex = "variations_";
	$hindex = "variations_".$modefull if ($modefull);
	return $self->{lmdb}->{$hindex} if exists $self->{lmdb}->{$hindex};
	$modefull = "r" unless $modefull;
	my ( $mode, $pipeline ) = split( '', $modefull );
	my $dir_out = $self->project->lmdb_cache_variations_dir();

	
	
	if ($mode eq "c"){
		if ($rocks) {
			my $dir_out_rocks = $self->project->rocks_cache_dir;
			system ("mkdir $dir_out && chmod a+rwx $dir_out" ) unless -e  $dir_out_rocks;
			$self->{lmdb}->{$hindex} =  $self->get_rocks_variations($mode);

		}
		else {
			$self->{lmdb}->{$hindex} =  $self->get_old_lmdb_variations($mode,$dir_out);
		}
	}
	else {
		#if ( -e  $dir_out_rocks){
		#		$self->{lmdb}->{$hindex} =  $self->get_rocks_variations($mode);
		#}
		#else {
			$self->{lmdb}->{$hindex} =  $self->get_old_lmdb_variations($mode,$dir_out);
		#}
	}
	
	return $self->{lmdb}->{$hindex};
}




sub get_rocks_variations {
	my ( $self, $modefull,$rocks) = @_;
	my $hindex = "variations_";
	$hindex = "variations_".$modefull.'_' if ($modefull);
	$hindex .= $self->name.'_';
	$hindex .= $$;
	return $self->project->{rocks}->{$hindex} if exists $self->project->{rocks}->{$hindex};
	my $rg;
	$self->project->{rocks}->{$hindex} = GenBoNoSqlRocksVariation->new(dir=>$self->project->rocks_directory("genbo"),mode=>"r",name=>$self->name.".genbo.rocks");
	return $self->project->{rocks}->{$hindex};

}
sub get_rocks_variations_index {
	my ( $self, $modefull,$rocks) = @_;
	my $hindex = "variations_";
	$hindex = "variations_".$modefull if ($modefull);
	$hindex .= $$;
	return $self->project->{rocks}->{$hindex} if exists $self->project->{rocks}->{$hindex};
	my $rg;
	$self->project->{rocks}->{$hindex} = GenBoNoSqlRocksVariation->new(dir=>$self->project->rocks_directory("index"),mode=>"r",name=>$self->name);
	return $self->project->{rocks}->{$hindex};

}

sub get_lmdb_cnvs {
	my ( $self, $mode ) = @_;
	return $self->{lmdb}->{cnvs} if exists $self->{lmdb}->{cnvs};
	$mode = "r" unless $mode;
	my $dir_out = $self->project->noSqlCnvsDir;
	$self->{lmdb}->{cnvs} = GenBoNoSqlLmdb->new(
		dir         => $dir_out,
		mode        => "$mode",
		name        => $self->name,
		is_compress => 1,
		vmtouch     => $self->buffer->vmtouch
	);
	return $self->{lmdb}->{cnvs};
}







#sub get_lmdb_patients_variations {
#	my ( $self, $mode, $patient ) = @_;
#	confess(
#"\n\nERROR: need GenBoPatient object in argument for GenBoChromosome->get_lmdb_patients_variations method. Die...\n\n"
#	) unless ($patient);
#	my $hindex = "patients_variations".$patient->id;
#	return $self->{lmdb}->{$hindex} if exists $self->{lmdb}->{$hindex};
#	$mode = "r" unless $mode;
#	my $dir_out = $self->project->lmdb_cache_patients_dir();
#	my $no2     = GenBoNoSqlLmdb->new(
#		dir     => $dir_out,
#		mode    => $mode,
#		name    => $self->name . "_" . $patient->name,
#		vmtouch => $self->buffer->vmtouch
#	);
#	$self->{lmdb}->{$hindex} = $no2;
#	return $no2;
#}


sub get_lmdb_categories {
	my ( $self, $mode ) = @_;
	return $self->_get_lmdb( $mode, "categories_annotations" );
}

sub get_lmdb_genes {
	my ( $self, $mode ) = @_;
	return $self->_get_lmdb( $mode, "genes" );
}

#sub get_lmdb_calling_methods {
#	my ( $self, $mode ) = @_;
#	return $self->_get_lmdb( $mode, "calling_methods" );
#}
#

#sub lmdb_score_impact {
#	my ( $self, $mode ) = @_;
#	my $buffer = $self->buffer();
#	$mode = "r" unless $mode;
#	#confess() unless $mode;
#	#$mode ="" unless $mode;
#
#	my $hindex = "lmdb_score_impact_".$mode;
#	return $self->{lmdb}->{$hindex} if exists $self->{lmdb}->{$hindex};
#		
#	my $dir_out = $self->project->dir_lmdb_score_impact();
#	$self->{lmdb}->{$hindex} =
#	  GenBoNoSqlLmdb->new(
#		dir         => $dir_out,
#		mode        => $mode,
#		name        => $self->name,
#		is_compress => 1,
#		vmtouch     => $self->buffer->vmtouch
#	  );    #GenBoNoSql->new(dir=>$output,mode=>$param);
#	 # $self->{lmdb}->{$hindex}->{vmtouch} = 1;
#	$self->{lmdb}->{$hindex}->clean_files() if ( $mode eq "c" );
#	return $self->{lmdb}->{$hindex};
#}

sub isVectorScore {
	my ( $self, $mode ) = @_;
	return -e $self->project->dir_lmdb_score_impact() . "/" . $self->name;
}


sub prepare_rocks_vector {
		my ( $self, $ids, $nodie ) = @_;
		$self->rocks_vector->prepare_vector($ids);
} 
sub getVectorScore {
	my ( $self, $id, $nodie ) = @_;
	confess() unless $id;
	if ( $id eq "all" ) {
		my $v = $self->getNewVector();
		$v->Fill();
		return $v;
	}
	if ($self->project->isRocks){
		return $self->getNewVector() if $self->size_vector == 0;
		return $self->rocks_vector->get_vector_chromosome($id);
		return $self->rocks_vector->get($id);
	}
	
	
	confess();
	my $v = $self->lmdb_score_impact("r")->get($id);
	if ($nodie) {
		return undef if ref($v) ne "Bit::Vector";
	}
	confess("miss vector => $id "
		  . $self->name . " "
		  . $self->project->dir_lmdb_score_impact() )
	  if ref($v) ne "Bit::Vector";

	#$v = $self->getNewVector() unless defined $v;
	return $v;
}

sub lmdb_vector_transcripts {
	my ( $self, $mode ) = @_;
	my $buffer = $self->buffer();
	my $hindex = "lmdb_vector_transcripts".$mode;
	return $self->{lmdb}->{$hindex} if exists $self->{lmdb}->{$hindex};
	$mode = "r" unless $mode;
	my $dir_out = $self->project->dir_lmdb_bitvector_transcripts();
	$self->{lmdb}->{$hindex} =
	  GenBoNoSqlLmdb->new(
		dir         => $dir_out,
		mode        => $mode,
		name        => $self->name,
		is_compress => 1,
		vmtouch     => $self->buffer->vmtouch
	  );    #GenBoNoSql->new(dir=>$output,mode=>$param);
	$self->{lmdb}->{$hindex}->clean_files()
	  if ( $mode eq "c" );
	return $self->{lmdb}->{$hindex};
}

sub get_lmdb_hgmd_class_ids {
	my ( $self, $mode, $class ) = @_;
	confess(
"\n\nERROR: need class in argument for GenBoChromosome->get_lmdb_hgmd_class_ids_intspan method. Die...\n\n"
	) unless ($class);
	$mode = "r" unless $mode;
	my $hindex = "hgmd_class_ids".$mode;
	return $self->{lmdb}->{$hindex} if exists $self->{lmdb}->{$hindex};


	my $dir_out =
		$self->project->get_public_data_directory("hgmd")
	  . '/class_'
	  . $class . '/';

	#my $dir_out = '/data-xfs/dev/mbras/TMP/hgmd/'.'/class_'.$class.'/';

	unless ( -e $dir_out ) {
		system("mkdir  $dir_out");
		system("chmod a+rwx $dir_out");
	}
	my $no2 = GenBoNoSqlLmdb->new(
		dir         => $dir_out,
		mode        => $mode,
		name        => $self->name,
		is_compress => 1,
		vmtouch     => $self->buffer->vmtouch
	);
	$no2->clean_files() if ( $mode eq "c" );
	$self->{lmdb}->{$hindex} = $no2;
	return $no2;
}

has list_hgmd_class_DM_var_ids => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $no   = $self->get_lmdb_hgmd_class_ids( 'r', 'DM' );
		my @lIds = @{ $no->get_keys() };
		$no->close();
		return \@lIds;
	},
);

has hash_hgmd_class_DM_var_ids => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my %h    = map { $_ => 1 } ( @{ $self->list_hgmd_class_DM_var_ids() } );
		return \%h;
	},
);

has hash_hgmd_class_DM_genes_name => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return
		  $self->project->buffer->queryHgmd->getDataHGMDPro_genes_for_class(
			$self->id(), 'DM' );
	},
);

has list_hgmd_DM_postions => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->project->buffer->queryHgmd->getDataHGMDPro_positions_for_class( $self->id(), 'DM' );
	},
);

#table avec les var ids stocke PENDANT le cache
has cache_hash_get_var_ids => (
	is      => 'rw',
	lazy    => 1,
	default => sub { {} },
);

#table avec les genes ids stocke PENDANT le cache
has cache_hash_get_gene_ids => (
	is      => 'rw',
	lazy    => 1,
	default => sub { {} },
);

has intspan_duplicate_region => (
	is => 'ro',

	#isa => 'ArrayRef[Str]',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $hash_dup;
		my $span = Set::IntSpan::Fast::XS->new;
		foreach my $patient ( @{ $self->getProject->getPatients } ) {
			my $filebed =
				$self->project->getVariationsDir("duplicate_region_calling")
			  . "/regions/"
			  . $patient->name()
			  . ".dup.bed";
			if ( -e $filebed ) {
				open( BED, $filebed );
				while (<BED>) {
					chomp();
					my ( $chr, $start, $end ) = split(" ");
					next if $chr ne $self->ucsc_name();
					$span->add_range( $start, $end );
				}
			}
		}
		return $span;
	},
);
has genesIntervalTree => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;

		if ( exists $self->{set_genes} ) {
			my $z = $self->project->liteIntervalTree->get_data( "genes_padding",$self->name );
			my $tree = Set::IntervalTree->new;

			foreach my $a (@$z) {
				next unless exists $self->{hset_genes}->{ $a->[0] };
				$tree->insert(@$a);
			}
			return $tree;
		}
		return $self->project->liteIntervalTree->get( "genes_padding",$self->name );

		#return  $self->project->liteIntspan->get("intspan_genes",$self->name);

	}
);

has regulatoryRegionsIntervalTree => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $tree = $self->project->litregionsIntervalTree->get("regulatory_region",$self->name);
		confess() unless $tree;
		return $tree;
	}
);
has annotatedRegionsIntervalTree => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;

		if ( exists $self->{set_genes} ) {
			my $z = $self->project->liteIntervalTree->get_data( "genes_padding",
				$self->name );
			my $tree = Set::IntervalTree->new;

			foreach my $a (@$z) {
				next unless exists $self->{hset_genes}->{ $a->[0] };
				$tree->insert(@$a);
			}
			return $tree;
		}
		my $tree = $self->project->liteIntervalTree->get( "genes_padding", $self->name );
		return $self->project->litregionsIntervalTree->append("regulatory_region",$self->name,$tree);
	}
);


has transcriptsIntervalTree => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->project->liteIntervalTree->get( "transcripts_padding",
			$self->name );

		#return  $self->project->liteIntspan->get("intspan_genes",$self->name);
	}
);

########
# lmdb SCORE
########

sub get_lmdb_score_db {
	my ( $self, $database ) = @_;
	warn $database;
	#my $database = "cadd";
	my $buffer = $self->buffer;
	my $chr    = $self->name;
	
	my $hindex = "lmdb_score_db_".$database;
	return $self->{lmdb}->{$hindex} if exists $self->{lmdb}->{$hindex};
	

	$self->{lmdb}->{$hindex} = GenBoNoSqlLmdbScore->new(
		dir        => $self->project->get_public_data_directory($database),
		mode       => "r",
		name       => $chr,
		is_integer => 1,
		vmtouch    => $self->buffer->vmtouch
	);
	return $self->{lmdb}->{$hindex};
}

sub close_lmdb_score {
	my ($self) = @_;
	$self->close_lmdb();
	my $buffer = $self->buffer;
	foreach my $database ( keys %{ $self->buffer->{lmdb_score} } ) {

		if ( exists $buffer->{lmdb_score}->{$database}->{ $self->name } ) {
			$buffer->{lmdb_score}->{$database}->{ $self->name }->close();
			delete $buffer->{lmdb_score}->{$database}->{ $self->name };
		}

	}

}

sub close_lmdb {
	my ($self) = @_;
	foreach my $lmdb (keys %{$self->{lmdb}}){
				$self->{lmdb}->{$lmdb}->close();
				delete $self->{lmdb}->{$lmdb};
	}
}
sub get_lmdb_spliceAI {
	my ($self)   = @_;
	confess();
	my $buffer   = $self->buffer;
	my $chr      = $self->name;
	my $database = "spliceAI";
	my $hindex = "lmdb_score_db_".$database;
	return $self->{lmdb}->{$hindex} if exists $self->{lmdb}->{$hindex};
	$self->{lmdb}->{$hindex} = GenBoNoSqlLmdb->new(
		dir         => $self->project->get_public_data_directory($database),
		mode        => "r",
		name        => $chr,
		is_integer  => 1,
		is_compress => 1,
		vmtouch     => $self->buffer->vmtouch
	);
	return $self->{lmdb}->{$hindex};
}

has intspan_spliceAI => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->get_lmdb_spliceAI()->get(0);
	}
);



sub score_gene_spliceAI {
	my ( $self, $v, $gene_name ) = @_;
	return undef unless defined $v;
	return { "AG" => 0, "AL" => 0, "DG" => 0, "DL" => 0 } if $v == 0;
	return undef unless exists $v->{$gene_name};
	my @data = unpack( "W4 C4", $v->{$gene_name} );

	#warn Dumper(@data);
	my %score;
	my @type = ( "AG", "AL", "DG", "DL" );
	for ( my $i = 0 ; $i < 4 ; $i++ ) {
		$score{ $type[$i] } = $data[$i] / 100;

		#push(@score,$data[$i]/100);
	}
	return \%score;
}

sub get_lmdb_score {
	my ( $self, $database, $variation ) = @_;
	
	my $position = $variation->start();
	my $allele   = $variation->sequence();
	my $db       = $self->get_lmdb_score_db($database);
	my $value    = $db->get($position);
	my $value2   = $db->get_score($position);
	return unless $value2;
	return $value2 if ( $value2 and $database eq 'ncboost' );
	return $value2->{$allele};
}

sub get_lmdb_database {
	my ( $self, $database, $type ) = @_;
	return $self->buffer->get_lmdb_database( $database, $self->name, $type );
}

#TODO: hgms DM for gene to do for HG38
sub is_hgmd_DM_for_gene {
	my ($self, $hgmd_id, $gene) = @_;
	
	return;
	
	my $db =  $self->getChromosome->get_lmdb_database("hgmd",'relation_variant_gene');
	my $pub = $db->get($hgmd_id);
	return 1 if ($pub and exists $pub->{$gene->external_name()});
	return if ($pub);
	confess("\n\nERROR GenBoChromosome::is_hgmd_DM_for_gene -> $hgmd_id not found. Die.\n\n");
}

sub is_clinvar_pathogenic_for_gene {
	my ($self, $clinvar_id, $gene) = @_;
	my $db =  $self->getChromosome->get_lmdb_database("clinvar",'relation_variant_gene');
	my $pub = $db->get($clinvar_id);
	return 1 if ($pub and exists $pub->{$gene->external_name()});
	return if ($pub);
	confess("\n\nERROR GenBoChromosome::is_clinvar_pathogenic_for_gene -> $clinvar_id not found. Die.\n\n");
}

sub purge_lmdb_score {
	my ($self) = @_;
	my $buffer = $self->buffer;
	my $chr    = $self->name;
	foreach my $db ( keys %{ $buffer->{lmdb_score} } ) {
		$buffer->{lmdb_score}->{$db}->{$chr}->close
		  if exists $buffer->{lmdb_score}->{$db}->{$chr};
		delete $buffer->{lmdb_score}->{$db}->{$chr}
		  if exists $buffer->{lmdb_score}->{$db}->{$chr};
	}
}

sub rocksdb {
	my($self,$db) = @_;
	confess unless $db;
	return $self->project->{rocks}->{$db} if exists $self->project->{rocks}->{$db};
	my $dir = $self->buffer->get_index_database_directory($db);
	$self->project->{rocks}->{$db} =  GenBoNoSqlRocksAnnotation->new(dir=>$dir,mode=>"r",name=>$self->name);
	return $self->project->{rocks}->{$db};
	
}

sub rocks_predictions {
	my ( $self, $mode ) = @_;
	$mode = "r" unless $mode;
	my $name = "prediction-".$mode.$$;
	return $self->project->{rocks}->{$name} if exists $self->project->{rocks}->{$name};
	 $self->project->{rocks}->{$name} =  GenBoNoSqlRocksAnnotation->new(dir=>$self->get_public_data_directory("prediction_matrix"),mode=>$mode,name=>$self->name);
	 return $self->project->{rocks}->{$name};
}

sub rocks_dejavu {
	my ( $self, $mode ) = @_;
	$mode = "r" unless $mode;
	my $name = "dejavu-".$mode.$$;
	return $self->project->{rocks}->{$name} if exists $self->project->{rocks}->{$name};
	 $self->project->{rocks}->{$name} = GenBoNoSqlRocksGenome->new(dir=>$self->project->deja_vu_rocks_dir,mode=>$mode,genome=>$self->project->genome_version_generic,index=>"genomic",chromosome=>$self->name);
	 return $self->project->{rocks}->{$name};
}

sub getDejaVuInfos {
	my ( $self, $id ) = @_;
	my $hres;
	confess();
	my $no = $self->rocks_dejavu();
	
	my $no1 = $self->project->lite_deja_vu2();
	my $string_infos = $no->dejavu($id);
	#warn $id;
	#return;
	return $hres unless ($string_infos);
	foreach my $string_infos_project (split('!', $string_infos)) {
		my @lTmp = split(':', $string_infos_project);
		my $project_name = 'NGS20'.$lTmp[0];
		my $nb_ho = $lTmp[1];
		my $nb_all = $lTmp[2];
		my $nb_he = $nb_all - $nb_ho;
		my $patients_info = $lTmp[3];
		$hres->{$project_name}->{nb} = $nb_all;
		$hres->{$project_name}->{ho} = $nb_ho;
		$hres->{$project_name}->{he} = $nb_he;
		$hres->{$project_name}->{patients} = "";
		$hres->{$project_name}->{string}   = "";
		my $hpat = $no1->get( "projects", $project_name );
		next;
		my $nb_pat = 0;
		foreach my $pat_id (split(',', $patients_info)) {
			$nb_pat++;
			my $pat_name = $hpat->{$pat_id};
			next unless ($pat_name);
			if ($nb_pat <= $nb_ho) { $hres->{$project_name}->{string} .= $pat_name.":1;"; }
			else { $hres->{$project_name}->{string} .= $pat_name.":2;"; }
			$hres->{$project_name}->{patients} .= $pat_name.";";
		}
		chop( $hres->{$project_name}->{patients} );
		chop( $hres->{$project_name}->{string} );
	}
	return $hres;
}

sub getShortResumeDejaVuInfosForDiagforRocksId {
	my ($self, $rocks_id) = @_;
	
}

sub getDejaVuInfosForDiagforVariant{
	my ($self, $v) = @_;
	my $chr = $self;
#	my $in_this_run_patients =  $self->project->{in_this_run_patients};
#	$in_this_run_patients->{total} =0 unless $in_this_run_patients->{total};
	#my $no1 = $self->project->lite_deja_vu2();
	my $no = $self->rocks_dejavu();
	my $h = $no->dejavu($v->rocksdb_id);
	#my $h = $no1->get($v->getChromosome->name,$v->id);
	#warn Dumper $h;
	#die();
	my $similar = $self->project->similarProjects();
	my $exomes = $self->project->exomeProjects();
	my $pe =  $self->project->countExomePatients();
	my $ps =  $self->project->countSimilarPatients();
	my $res;
	$res->{similar_projects} = 0;
	$res->{other_projects} = 0;
	$res->{exome_projects} = 0;
	$res->{other_patients} = 0;
	$res->{exome_patients} = 0;
	$res->{similar_patients} = 0;
	$res->{other_patients_ho} = 0;
	$res->{exome_patients_ho} = 0;
	$res->{similar_patients_ho} = 0;
	$res->{total_in_this_run_patients} = 0;
	#$res->{total_in_this_run_patients} = $in_this_run_patients->{total} + 0;
	if ($res->{total_in_this_run_patients} == 0 ){
		$res->{total_in_this_run_patients} = scalar(@{$self->project->getPatients});
	}
	$res->{in_this_run_patients} = 0;
	$res->{in_this_run_patients} += scalar(@{$v->getPatients});
	return $res unless ($h);
	
	foreach my $l (split("!",$h)) {
		my($p,$nho,$nhe,$info) = split(":",$l);
		$p = "NGS20".$p;
		next if $p eq $self->name();
		
		#TODO: here ! a faire
#		if (exists $in_this_run_patients->{$p}){
#			
#			my (@samples) = split(",",$info);
#			foreach my $s (@samples){
#				if (exists $in_this_run_patients->{$p}->{$s}){
#					$res->{in_this_run_patients} ++;
#				}
#			}
#		}
		#IN EXOME 	
		if (exists $exomes->{$p}){
			$res->{exome_projects}  ++; 
			$res->{exome_patients}   += $nhe;
			$res->{exome_patients_ho}   += $nho;
		}
		#in similar;
		if (exists $similar->{$p}){
			$res->{similar_projects}  ++;
			$res->{similar_patients} += $nhe;
			$res->{similar_patients_ho} += $nho;
		}
		else {
			$res->{other_projects} ++;
			$res->{other_patients}+= $nhe;
			$res->{other_patients_ho}+= $nho;
			$res->{other_projects_ho} ++ if $nho > 0;
		}
	}	
	$res->{total_exome_projects} =  scalar(keys %{$self->project->exomeProjects()});
	$res->{total_exome_patients} =  $self->project->countExomePatients();
	$res->{total_similar_projects} =  scalar(keys %{$self->project->similarProjects()});
	$res->{total_similar_patients} =  $self->project->countSimilarPatients();
	return $res;
	
}

1;
