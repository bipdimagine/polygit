package GenBoProject;
use FindBin qw($Bin);
use strict;
use POSIX qw(strftime);
use Moo;
use Carp;
use Data::Dumper;
extends 'GenBo';
use GenBoFamily;
use GenBoPatient;
use GenBoSomaticGroup;
use GenBoChromosome;
use GenBoVariation;
use GenBoMnp;
use decode_prediction_matrix;
use GenBoDeletion;
use GenBoInsertion;
use GenBoLargeDeletion;
use GenBoInversion;
use GenBoBoundary;
use GenBoLargeDuplication;
use GenBoLargeInsertion;
use GenBoCnv;
use GenBoMei;

#use GenBoComplex;
use GenBoRegulatoryRegion;
use GenBoReference;
use GenBoGene;
use GenBoGeneCache;
use GenBoExon;
use GenBoIntron;
use GenBoTranscript;
use GenBo;
use GenBoProtein;
use GenBoRun;
use GenBoPrimer;
use GenBoPanel;
use GenBoPhenotype;
use GenBoBundle;
use GenBoSVDel;
use GenBoSVDup;
use BioTools;
use File::Path qw(make_path);
use base 'Exporter';
use Storable qw(freeze thaw );
use JSON::XS;
use GenBoNoSql;
use GenBoNoSqlText;
use GenBoNoSqlDejaVu;
use GenBoNoSqlDejaVuCNV;
use GenBoNoSqlDejaVuSV;
#use GenBoNoSqlDejaVuSV_agregate;
use GenBoNoSqlDejaVuJunctions;
use GenBoNoSqlDejaVuJunctionsResume;
use GenBoNoSqlDejaVuJunctionsCanoniques;
use GenBoNoSqlDejaVuJunctionsPhenotype;
use GenBoNoSqlAnnotation;
use GenBoNoSqlLmdbInteger;
use GenBoJunction;
use GenBoJunctionCache;
use GenBoNoSqlRocks;
use Storable qw(store retrieve freeze dclone thaw);

#use LMDB_File qw(:flags :cursor_op);
use GenBoNoSqlLmdb;
use GenBoNoSqlIntervalTree;
use QueryValidationAcmg;
use GenBoNoSqlLmdbScore;
use Time::HiRes qw (time);
use Sys::Hostname;
use File::Temp;
use Bio::DB::HTS::Faidx;


sub hasHgmdAccess {
	my ( $self, $user ) = @_;
	return 1 if ( $self->buffer->queryHgmd->getHGMD($user) == 1 );
	return;
}

has extAlignFile => (
is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		return "cram" if $self->getVersion() =~ /HG38/;
		return "bam";
		return 1;
	},
);

# just stupid mehtod to print on stdout a dot each  modulo $max  it's for waiting for CGI
# method to write dot on stdout for cgi waiting
sub print_dot {
	my ( $self, $max ) = @_;
	return unless ( $self->cgi_object() );
	$self->{count_dot}++;
	if ( $self->count_dot() >= $max ) {
		$self->count_dot(0);
		print '.';
	}
}

has count_dot => (
	is      => 'rw',
	lazy    => 1,
	default => 0,
);
has start_timer1 => (
	is      => 'rw',
	lazy    => 1,
	default => sub { return time; }
);
has timer => (
	is      => 'rw',
	lazy    => 1,
	default => 0,
);

sub start_timer {
	my ($self) = @_;
	$self->start_timer1(time);
}

sub add_timer {
	my ($self) = @_;
	my $a = $self->timer() + abs( time - $self->start_timer1 );
	$self->timer($a);
}

has cgi_object => (
	is      => 'rw',
	lazy    => 1,
	default => undef,
);

has is_lmdb_variation_cache => (
	is      => 'rw',
	lazy    => 1,
	default =>sub {
		my $self = shift;
		my $chr = $self->getChromosomes->[0];
		my $dir = $self->lmdb_cache_variations_dir();#$self->project->lmdb_pipeline_dir() . "/lmdb_variations/";
		return -e $dir."/".$chr->name.".lmdb_var";
		#return 0;
	}
);
####################
# end method waiting
####################

has no_cache => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		return 0;
	}
);

has get_list_emails => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my @lUsers;
		foreach my $h (
			@{ $self->buffer->getQuery()->getOwnerProject( $self->id() ) } )
		{
			push( @lUsers, $h->{email} );
		}
		
		foreach my $g (@{ $self->buffer->getQuery()->getGroups( $self->id() ) } ){
				push( @lUsers, $g );
		}
		
		return \@lUsers;
	},
);

has test => (
	is      => 'ro',
	default => 0,
);
has print_waiting => (
	is      => 'rw',
	default => undef,
);
has id => (
	is      => 'ro',
	lazy    => 1,
	reader  => 'id',
	default => sub {
		my $self = shift;
		return $self->infosProject->{id};
	},
);

has isFamilial => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my @t    = grep {
			exists $self->pedigree_details()->{$_}->{mother}
			  or exists $self->pedigree_details()->{$_}->{father}
		} keys %{ $self->pedigree_details() };
		return 1 if scalar(@t) >= 1;
		foreach my $fam_name (keys %{$self->pedigree_details()}) {
			my $nb_childs = 0;
			$nb_childs = scalar(keys %{$self->pedigree_details()->{$fam_name}->{children}}) if (exists $self->pedigree_details()->{$fam_name}->{children});
			return 1 if ($nb_childs > 1);
		}
		return;
#	return 1 if (scalar( grep {exists $self->pedigree_details()->{$_}->{mother} }keys %{$self->pedigree_details()}) ne scalar(@{$self->getPatients}));
#	return;
	},
);

has isFamilialStudy => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;

		#confess();
		return 1 if ( $self->isFamilial() );
		return;
	},
);

has isSomatic => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return 1 if ( $self->infosProject->{is_somatic} == 1 );
		return;
	},
);

has isSomaticStudy => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->isSomatic();
	},
);

has biotools => (
	is   => 'ro',
	lazy => 1,

	default => sub {
		my $self = shift;
		return BioTools->new();
	},
);

has project => (
	is      => 'ro',
	lazy    => 1,
	reader  => 'getProject',
	default => sub {
		my $self = shift;
		return $self;
	},
);

has year => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my ( $a, $b ) = split( "_", $self->name );
		my $y = substr $a, -4;

		return $y;
	},
);

has isDiagnostic => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return undef if $self->isExome() || $self->isGenome();
		foreach my $c ( @{ $self->getCaptures } ) {
			return $c->analyse
			  if lc( $c->analyse ) !~ /exome/
			  && lc( $c->analyse ) !~ /genome/
			  && lc( $c->analyse ) !~ /rna/;

			#			return "diagnostic" if $c->analyse eq "diagnostic";
			#			return "ciliome" if $c->analyse eq "ciliome";
			#			return "callosome" if $c->analyse eq "callosome";
		}
		return 1;
	},
);

has validations => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $zen  = $self->validations_query->getAllValidations(1);
		return $zen;
	},
);


sub getValidationVariation {
	my ( $self, $val_id, $patient ) = @_;
	return undef unless exists $self->validations->{$val_id};
	
	my $tvid = $self->validations->{$val_id};
	if ( $self->isDefidiag && $patient ) {
		my (@test)
		  ; # = grep {$_->{sample_id} ne $patient->id && $_->{sample_name} ne $patient->name} @{$tvid};
		foreach my $t ( @{$tvid} ) {
			if (   $t->{sample_id} ne $patient->id
				&& $t->{sample_name} eq $patient->name )
			{
				next;
			}
			push( @test, $t );
		}
		return undef unless @test;
		$tvid = \@test;
	}
	return $tvid->[0];
}

#has validations_query => (
#	is      => 'ro',
#	lazy    => 1,
#	default => sub {
#		my $self = shift;
#		my $vquery;
#		if ( $self->isExome() or $self->isGenome() ) {
#			$vquery = QueryValidationAcmg->new(
#				dbh      => $self->buffer->dbh,
#				database => "test"
#			);
#		}
#		else {
#			$vquery = validationQuery->new(
#				dbh          => $self->buffer->dbh,
#				capture_name => $self->getCapture->validation_db()
#			);
#		}
#		return $vquery;
#	}
#
#);

sub validations_query {
	my ($self,$type) = @_;
	return $self->buffer->validations_query();
	 return $self->buffer->{validations_query1} if exists   $self->{validations_query1};
	  $self->{validations_query1} = QueryValidationAcmg->new(
				dbh      => $self->buffer->dbh,
				database => "ACMG"
			);
#			
#		if ( $type or $self->isExome() or $self->isGenome() ) {
#			 $self->{validations_query1} = QueryValidationAcmg->new(
#				dbh      => $self->buffer->dbh,
#				database => "test"
#			);
#		}
#		else {
#			 $self->{validations_query1} = validationQuery->new(
#				dbh          => $self->buffer->dbh,
#				capture_name => $self->getCapture->validation_db()
#			);
#		}
		return   $self->{validations_query1};
}

has validation_db => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		
		my $self = shift;
		return $self->infosProject->{validation_db}
		  if $self->infosProject->{validation_db};
		
		if ( $self->isExome() or $self->isGenome() ) {
			return  $self->infosProject->{validation_db};
		}
		
		my %db;
		foreach my $c ( @{ $self->getCaptures } ) {
			$db{ $c->validation_db }++ if $c->validation_db;
		}
		my @res = keys %db;
		confess if scalar(@res) > 1;
		return $res[0];
	},
);

has isPolyQuest => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return 1 if $self->isExome() or $self->isGenome() or $self->isCiliome();

	}

);
has isRnaseq => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		foreach my $c ( @{ $self->getCaptures } ) {
			return 1 if lc( $c->analyse ) =~ /rna/;
		}
		return undef;
	},
);
has isRna => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $find ;
		foreach my $p ( @{ $self->getPatients } ) {
			$find ++ if $p->isRna;
		}
		return $find;
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
has isCiliome => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		foreach my $c ( @{ $self->getCaptures } ) {
			return 1 if lc( $c->analyse ) =~ /ciliome/;
		}
		return undef;
	},
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
has cgi_user => (
	is      => 'rw',
);

has isDefidiag => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		if ( $self->validation_db() and $self->validation_db() eq 'defidiag' ) {
			return 1;
		}
		return undef;
	},
);
has isDefidiagSolo => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return undef unless $self->isDefidiag();
		my $find = 1;
		foreach my $f ( @{ $self->getFamilies } ) {
			$find = undef if scalar( @{ $f->getMembers } ) > 1;
			last unless $find;
		}

		return $find;
	},
);

has isRnaSeq => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		foreach my $c ( @{ $self->getCaptures } ) {

			#warn $c->analyse;
			return 1 if lc( $c->analyse ) =~ /rna/;
		}
		return undef;
	},
);
has name => (
	is       => 'ro',
	required => 1,
);

has release => (
	is       => 'rw',
	required => 1,
);

has buffer => (
	is       => 'rw',
	required => 1,
	weak_ref => 1,
);

has db => (
	is      => 'ro',
	lazy    => 1,
	default => 'genbo',
);

has infosProject => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self  = shift;
		my $query = $self->buffer->getQuery();
		my $res   = $query->getProjectByName( $self->name() );
		return $res;
	},
);

has local_config_file => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return
			$self->getCacheBitVectorDir()
		  . '/genbo_'
		  . $self->name() . '.cfg';
	}
);

has is_pipeline_dir => (
	is      => 'rw',
	default => sub {
		return undef;
	}
);

has lmdb_pipeline_dir => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		my $name    = $self->name;
		my $dir_out = File::Temp->newdir(
			TEMPLATE => "$name.temporary.XXXXXXX",
			DIR      => '/data-beegfs/tmp/',
		);
		$self->is_pipeline_dir(1);
		return $dir_out;
	}
);

sub delete_lmdb_pipeline_dir {
	my $self = shift;
	if ( $self->is_pipeline_dir ) {
		warn "coucou " . $self->lmdb_pipeline_dir;
		if ( -e $self->lmdb_pipeline_dir
			&& $self->lmdb_pipeline_dir =~ /temporary/ )
		{
			system( "rm -r " . $self->lmdb_pipeline_dir . '/' );
		}

	}
}



sub get_lmdb_cache_summary {
	my ( $self, $mode ) = @_;
	
	$mode = "r" unless $mode;
	my $dir_out = $self->getCacheDir();
	unless (-e $dir_out."/".$self->name.".summary.cache"){
		$mode = "c";
	}
	my $no2     = GenBoNoSqlLmdbCache->new(
		dir     => $dir_out,
		mode    => $mode,
		name    => $self->name.".summary.cache",
		is_compress => 1,
		vmtouch => $self->buffer->vmtouch
	);
	if ( $mode eq "c"){
		$no2->put("cdate",time);
		system("chmod a+w ".$no2->filename);
	}
	return $no2;
}

sub get_lmdb_cache_cnv {
	my ( $self, $mode ) = @_;
	
	$mode = "r" unless $mode;
	my $dir_out = $self->getCacheDir();
	unless (-e $dir_out."/".$self->name.".table_cnv.cache"){
		$mode = "c";
	}
	my $no2     = GenBoNoSqlLmdbCache->new(
		dir     => $dir_out,
		mode    => $mode,
		name    => $self->name.".table_cnv.cache",
		is_compress => 1,
		vmtouch => $self->buffer->vmtouch
	);
	if ( $mode eq "c"){
		$no2->put("cdate",time);
		system("chmod a+w ".$no2->filename);
	}
	return $no2;
}

has lmdb_cache_dir => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->getCacheBitVectorDir() . "/lmdb_cache";
		return $self->makedir($dir);
	}
);
has rocks_cache_dir => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
	my $genome_version = $self->genome_version();
	my $annot_version = $self->annotation_version();

	my $name = $self->buffer->config_path("root","cache").'/'.$genome_version;#$self->buffer()->getDataDirectory("cache")."/rocks/".$genome_version;
	$name .= '.'.$annot_version if ($annot_version and $annot_version ne '.');
	$name .= "/".$self->name();
	$self->makedir($name);
	return $name;
	
	}
);


sub rocks_pipeline_directory {
	my ($self,$type) = @_;
	
	my $path = $self->project_pipeline_path . "/rocks_tmp/";
	$path .="$type/" if $type;
	return $self->makedir($path);
}
sub rocks_directory_path {
	my ($self,$type,$create) = @_;
	my $path =  $self->rocks_cache_dir();
	$path .="/$type/" if $type;
	return $path;
	
}

sub rocks_directory {
	my ($self,$type,$create) = @_;
	my $path =  $self->rocks_directory_path($type);
	#unless ($create){
	#confess() unless -e $path;  
	#return $path;
	#}
	return $self->makedir($path);
}

sub rocksdb {
	my($self,$db) = @_;
	confess unless $db;
	return $self->{rocks}->{$db} if exists $self->{rocks}->{$db};
	my $dir = $self->buffer->get_index_database_directory($db);
	$self->{rocks}->{$db} =  GenBoNoSqlRocksAnnotation->new(dir=>$dir,mode=>"r",name=>$db);
	return $self->{rocks}->{$db};
}

has lmdb_cache_variations_dir => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;

		#test before debug
		my $dir = $self->lmdb_cache_dir() . "/variations";
		return $self->makedir($dir);
	}
);

has lmdb_cache_patients_dir => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->lmdb_cache_dir() . "/patients";
		return $self->makedir($dir);
	}
);

has lmdb_cache_dejavu_dir => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->lmdb_cache_dir() . "/dejavu";
		return $self->makedir($dir);

	}
);

has lmdb_cache_dejavuho_dir => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->lmdb_cache_dir() . "/dejavu_ho";
		return $self->makedir($dir);

	}
);
has quality_dir => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->project_path() . "/quality";
		return $self->makedir($dir);

	}
);

has lmdb_ncboost_path => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->annotation_public_path() . '/ncboost/';
		return $self->makedir($dir);

	}
);

has lmdb_gnomad_path => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->annotation_public_path() . '/gnomad/';
		return $self->makedir($dir);

	}
);

sub get_lmdb_ncboost_chromosomes_vectors {
	my ( $self, $chr, $mode ) = @_;
	confess(
"\n\nERROR: GenBoProject::get_lmdb_ncboost_chromosomes_vectors need chr object. Die.\n\n"
	) unless ($chr);
	confess(
"\n\nERROR: GenBoProject::get_lmdb_ncboost_chromosomes_vectors need chr object. Die.\n\n"
	) unless ( ref($chr) =~ /GenBoChromosome/ );
	$mode = "r" unless $mode;
	my $no = GenBoNoSqlLmdb->new(
		dir         => $self->lmdb_ncboost_path() . '/vector/',
		mode        => $mode,
		name        => $chr->id() . '.vector',
		is_compress => 1
	);
	$no->clean_files() if ( $mode eq "c" );
	return $no;
}

sub get_lmdb_gnomad_chromosomes_vectors {
	my ( $self, $chr, $mode ) = @_;
	confess(
"\n\nERROR: GenBoProject::get_lmdb_ncboost_chromosomes_vectors need chr object. Die.\n\n"
	) unless ($chr);
	confess(
"\n\nERROR: GenBoProject::get_lmdb_ncboost_chromosomes_vectors need chr object. Die.\n\n"
	) unless ( ref($chr) =~ /GenBoChromosome/ );
	$mode = "r" unless $mode;
	my $no = GenBoNoSqlLmdb->new(
		dir         => '/data-xfs/dev/mbras/gnomad_test_vector/',
		mode        => $mode,
		name        => $chr->id() . '.vector',
		is_compress => 1
	);

	$no->clean_files() if ( $mode eq "c" );
	return $no;
}

has bundle_infos => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self  = shift;
		my $query = $self->buffer->getQuery();
		if ( $self->isExome ) {
			my $res = $query->getBundleTranscripts( $self->id );
			return $res;
		}
		else {
			my $captures = $self->getCaptures();

			my $res;
			my %temp_bundle;
			foreach my $capture (@$captures) {
				my $hquery = $query->getCaptureTranscripts( $capture->id );
				%temp_bundle = ( %temp_bundle, %{ $hquery->{bundle} } );
				foreach my $tr_id ( @{ $capture->transcripts_name() } ) {
					die($tr_id) unless $self->getProject->rocksGenBo->synonym($tr_id);
					push(
						@{ $res->{transcripts_name} },
						$self->getProject->rocksGenBo->synonym($tr_id)
					);

				}
			}
			$res->{bundle} = \%temp_bundle;
			foreach my $b ( keys %temp_bundle ) {
				foreach my $bc ( @{ $temp_bundle{$b} } ) {
					my $tr = $bc->{ENSEMBL_ID};
					if (not $tr =~ /_.+/) {
						foreach my $chr (@{$self->getChromosomes()}) {
							my $id2 = $tr.'_'.$chr->id();
							if ($self->rocksGenBo->synonym($id2)) {
								$tr = $id2;
								last;
							}
						}
					}
					push(
						@{
							$res->{transcripts}->{$self->getProject->rocksGenBo->synonym($tr)}
						},
						$b
					);
				}
			}
			return $res;
		}
	},
);

sub return_bundle {
	my ( $self, $tr ) = @_;
	return $self->bundle_infos()->{transcripts}->{$tr};
}

has bundles => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		$self->bundle_infos->{bundle};
	},
);

has bundle_transcripts => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;

		#my @t = keys %{$self->bundle_infos->{transcripts}};

		return $self->bundle_infos->{transcripts_name};
	},
);

has _getListTranscripts => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		if ( $self->isDiagnostic ) {

			return $self->bundle_transcripts();
		}
		else {
			my @a =
			  grep { $_ !~ /GL/ }
			  grep { $_ =~ /ENST/ }
			  @{ $self->lmdbMainTranscripts->get('transcripts') };
			return \@a if @a;
			my @z = grep { $_ !~ /GL/ }
			  grep { $_ =~ /ENST/ } @{ $self->lmdbMainTranscripts->get_keys };
			die();
			return \@z;
		}
	},
);

sub getListTranscripts {
	my ( $self, $args ) = @_;
	my $t = $self->_getListTranscripts;
	return $t unless ($args);
	my @res = @$t;
	if ( $args->{transcripts} && $args->{transcripts} ne "all" ) {
		my @ids = split( ",", $args->{transcripts} );
		my $trs = $self->newTranscripts( \@ids );
		my @res = map { $_->id } @$trs;
		return \@res;
	}
	return $t;

	confess();

}

sub getListGenBoTranscripts {
	my ( $self, $args ) = @_;
	my $t   = $self->_getListTranscripts;
	my $trs = $self->newTranscripts($t);
	return $trs;

}

has mean_amount_reads => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $nb;
		foreach my $patient (@{$self->getPatients()}){
			$nb+= $patient->nb_reads->{all};
		}
		return int($nb/(scalar @{$self->getPatients()}));
	},
);
has phenotypes => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $h;
		foreach my $phenotype ( @{ $self->getPhenotypes() } ) {
			$h->{ $phenotype->name() } = undef;
		}
		my @lPhenotypes = keys %$h;
		return \@lPhenotypes;
	},
);


has similarProjectsId => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $results;
		my $query = $self->buffer->getQuery();
	#	
			my $phenotypes = $self->getPhenotypes();
			return {} unless $phenotypes;
			return {} unless @$phenotypes;
			
			foreach my $ph (@$phenotypes){
				map { $results->{$_}++ } @{ $query->getSimilarProjectsIdByPhenotype($ph)};
			}
			 
		if ( $self->isExome or $self->isGenome() ) { 
			return {} unless $results;
			return $results;
		}
		
		foreach my $c ( @{ $self->getCaptures() } ) {
			my $analyse = $c->infos->{analyse};
			my $vdb     = $c->validation_db;
			if ( $self->isExome or $self->isGenome() ) {
				my $phenotypes = $self->phenotypes();
				foreach my $ph (@$phenotypes){
					map { $results->{$_}++ }  @{ $query->getSimilarProjectsIdByPhenotype( $ph )};
				}
			}
			elsif ($analyse) {
				map { $results->{$_}++ }  @{ $query->getSimilarProjectsIdByAnalyse($analyse) };
			}
			else {
				confess();
				map { $results->{$_}++ } @{ $query->getSimilarProjects( $c->id ) };
			}
			if ($vdb) {
				map { $results->{$_}++ }
				  @{ $query->getSimilarProjectsIdByValidation_db($vdb) };
			}

		}
		delete $results->{ $self->id };
		return $results;
	},
);

has similarProjects => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $results;
		my $query = $self->buffer->getQuery();
	#	
			my $phenotypes = $self->getPhenotypes();
			return {} unless $phenotypes;
			return {} unless @$phenotypes;
			map { $results->{$_}++ } @{ $phenotypes->[0]->projects_name } if $phenotypes;
			 # warn Dumper $results;
			 
		if ( $self->isExome or $self->isGenome() ) { 
			return {} unless $results;
			return $results;
		}
		foreach my $c ( @{ $self->getCaptures() } ) {
			my $analyse = $c->infos->{analyse};
			my $vdb     = $c->validation_db;
			if ( $self->isExome or $self->isGenome() ) {
				my $phenotypes = $self->phenotypes();
				map { $results->{$_}++ }
				  @{ $query->getSimilarProjectsByPhenotype( $phenotypes->[0] )
				  };
			}
			elsif ($analyse) {
				map { $results->{$_}++ }
				  @{ $query->getSimilarProjectsByAnalyse($analyse) };

				#next;
			}
			else {
				map { $results->{$_}++ } @{ $query->getSimilarProjects( $c->id ) };
			}
			if ($vdb) {
				map { $results->{$_}++ }
				  @{ $query->getSimilarProjectsByValidation_db($vdb) };
			}

		}
		delete $results->{ $self->name };
		return $results;
	},
);




has countSimilarPatients => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		my $similar = $self->similarProjects();

		my $query = $self->buffer->getQuery();
		my $nb    = $query->countPatients( [ keys %$similar ] );
		return $nb;

	},
);

has exomeProjectsId => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self  = shift;
		my @lists = $self->buffer()->listProjectsExomeForDejaVu();
		my %hash;
		map { $hash{$_}++ } @{ $self->buffer()->getQuery->listProjectsIdExomeForDejaVu() };
		return \%hash;
	},
);
has exomeProjects => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self  = shift;
		my @lists = $self->buffer()->listProjectsExomeForDejaVu();
		my %hash;
		map { $hash{$_}++ } @{ $self->buffer()->getQuery->listProjectsExomeForDejaVu() };
		return \%hash;
	},
);

has countExomePatients => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		my $similar = $self->exomeProjects();

		my $query = $self->buffer->getQuery();
		my $nb    = $query->countPatients( [ keys %$similar ] );
		return $nb;

	},
);
has description => (
	is => 'ro',

	lazy    => 1,
	default => sub {
		my $self = shift;

		return $self->infosProject->{description};
	},
);
has creation_date => (
	is => 'ro',

	lazy    => 1,
	default => sub {
		my $self = shift;

		return $self->infosProject->{creation_date};
	},
);

has projectTypeId => (
	is      => 'ro',
	lazy    => 1,
	reader  => 'getProjectTypeId',
	default => sub {
		my $self = shift;

		return $self->infosProject->{projectTypeId};
	},
);

has projectType => (
	is      => 'ro',
	lazy    => 1,
	reader  => 'getProjectType',
	default => sub {
		my $self = shift;
		return $self->infosProject->{projectType};

	},
);

has is_ngs => (
	is      => 'ro',
	lazy    => 1,
	reader  => 'is_ngs',
	default => sub {
		my $self = shift;
		return $self->getProjectType eq "ngs";
	},
);

has is_classic => (
	is      => 'ro',
	lazy    => 1,
	reader  => 'is_classic',
	default => sub {
		my $self = shift;
		return $self->getProjectType eq "classic";
	},
);

has is_cnv => (
	is      => 'ro',
	lazy    => 1,
	reader  => 'is_cnv',
	default => sub {
		my $self = shift;
		return $self->getProjectType eq "cnv";
	},
);

has is_array => (
	is      => 'ro',
	lazy    => 1,
	reader  => 'is_array',
	default => sub {
		my $self = shift;
		return $self->getProjectType eq "array";
	},
);

has geneWindow => (
	is      => 'rw',
	lazy    => 1,
	reader  => 'getGeneWindow',
	default => sub {
		my $self = shift;
		return ( 5000 * 2 ) if $self->is_cnv;
		return 500;
	},
);

has is_human_genome => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return 1 if $self->getVersion() =~ /HG/;
		return;
	},
);

has version => (
	is      => 'rw',
	lazy    => 1,
	reader  => 'getVersion',
	default => sub {
		my $self = shift;
		return $self->infosProject->{version};
	},
);

has fastq_screen_path => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir = $self->getProjectRootPath().'/fastq_screen/';
		unless (-d $dir) {
			$self->makedir($dir);
		}
		return $dir;
	},
);

has project_root_path => (
	is      => 'rw',
	lazy    => 1,
	reader  => 'getProjectRootPath',
	default => sub {
		my $self     = shift;
		my $dirNgs   = $self->buffer->config->{project_data}->{ngs};
		my $pathRoot = $self->buffer->config_path("root","project_data");
		my $path1    = $pathRoot . "/" . $dirNgs . "/" . $self->name() . '/';
		$self->makedir($path1);

		# LN -S vers XFS
		if ( exists $self->buffer->config->{project_data}->{alias}
			and $self->buffer->config->{project_data}->{alias} ne '' )
		{
			my $pathAliasRoot = $self->buffer->config->{project_data}->{alias};
			my $path2 = $pathAliasRoot . "/" . $dirNgs . "/" . $self->name();
			unless ( -l $path2 ) {
				if ( -d $path2 ) {
					confess( $self->name
						  . " \n\nERROR: [GenBoProject -> getProjectRootPath] Directory found but symbolic link expected... Die...\n\n  $path2 $pathAliasRoot\n\n"
					);
					die;
				}
				#system("ln -s $path1 $path2");
			}
		}
		return $path1;
	},
);

has project_path => (
	is      => 'rw',
	lazy    => 1,
	reader  => 'getProjectPath',
	default => sub {
		my $self = shift;
		my $path = $self->getProjectRootPath() . "/" . $self->genome_version . "/";
		return $self->makedir($path);
	},
);

has project_log => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $path = $self->getProjectRootPath() . "/log/";
		return $self->makedir($path);
	},
);
has project_pipeline_path => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self     = shift;
		my $pathRoot = $self->buffer->config_path("root","project_pipeline");
		my $path     = $pathRoot . "/tmp." . $self->name() . "/";
		$self->makedir($path);
		$path .= $self->getVersion . "/";
		return $self->makedir($path);
	},
);
has project_dragen_pipeline_path_name => (
is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self     = shift;
		my $pathRoot = $self->buffer->config_path("root","dragen_pipeline");
		my $path     = $pathRoot . "/" . $self->name() . "/";
		$path .= $self->getVersion . "/";
		return $path;
	},
);
has project_epi2me_pipeline_path_name => (
is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self     = shift;
		my $pathRoot = $self->buffer->config_path("root","pipeline");
		my $path     = $pathRoot . "/" . $self->name() . "/";
		$path .= $self->getVersion . "/";
		return $path;
	},
);
has dragen_fastq => (
is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self     = shift;
		my $pathRoot = $self->buffer->config_path("root","dragen_pipeline");
		my $path     = $pathRoot . "/" . $self->getRun->name()."/";
		return ($self->makedir($path));
	},
);
has project_dragen_pipeline_path => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self     = shift;
		return $self->makedir($self->project_dragen_pipeline_path_name);
	},
);

has project_epi2me_pipeline_path => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self     = shift;
		return $self->makedir($self->project_epi2me_pipeline_path_name);
	},
); 

has project_dragen_demultiplex_path => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self     = shift;
		my $pathRoot = "/data-dragen/fastq/";
		my $run = $self->getRuns->[0];
		die() if scalar(@{$self->getRuns}) > 1;
		
		my $path     = $pathRoot . "/" . $run->name().".".$self->name.time."/";
		$self->makedir($path);
		return $path;
	},
);
has project_metrics_path => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $path = $self->getAlignmentRootDir();
		$path = $path . "/stats";
		return $self->makedir($path);
	},
);

sub projectType {
	my $self = shift;
	return "ngs";
}

has pipeline_tmp_dir => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $tmp1 = "/tmp/pipeline/";
		system("mkdir /tmp/pipeline;chmod a+rwx  /tmp/pipeline") unless -e $tmp1;
		my $tmp = File::Temp::tempdir(
			TEMPLATE => $self->name . ".XXXXXXXX",
			DIR      => "$tmp1"
		);
		system("date > $tmp/test.out");
		$self->{tmp_dir} = 1;
		warn $tmp;
		return $tmp;
	},
);

sub destroy_tmp_dir {
	my $self = shift;
	if ( exists $self->{pipeline_tmp_dir} ) {
		my $dir = $self->{pipeline_tmp_dir};
		return unless -e $dir;
		confess() unless $dir =~ /pipeline/;
		system("find $dir -type f  -exec rm -f {} \ ; rmdir $dir ");
		delete $self->{pipeline_tmp_dir};
	}
}

sub getRootDir {
	my $self = shift;
	return $self->getProjectPath();
}

has sequenceRootDir => (
	is      => 'rw',
	lazy    => 1,
	reader  => 'getSequencesRootDir',
	default => sub {
		my $self = shift;
		confess();
		return;
		my $path = $self->buffer()->config_path("root","project_data");
		$path =
			$path . "/"
		  . $self->getProjectType() . "/"
		  . $self->name()
		  . "/sequences/";
		return $self->makedir($path);
	},
);

has sequencesDir => (
	is      => 'ro',
	lazy    => 1,
	reader  => 'getSequencesDirectories',
	default => sub {
		my $self = shift;
		confess();
		my $sequence_dir       = "sequences";
		my $root_dir           = $self->getSequencesRootDir();
		my $methods_sequencing = $self->getSequencingMachines();
		my @dirs;
		foreach my $m (@$methods_sequencing) {
			my $dir = $root_dir . "/" . lc($m) . "/";
			push( @dirs, $dir );
			$self->makedir($dir);
		}
		return \@dirs;
	},
);

has genome_version => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		my $version = $self->getVersion();
		my ( $g, $v ) = split( /\./, $version );
		return $g;
	},
);

sub update_public_database_version {
	my ($self) = @_;

	#read gencode from capture ?
	my $query   = $self->buffer->getQuery();
	my $version = $query->getMaxPublicDatabaseVersion( $self->id );
	confess() unless $version;
	$query->insertPublicDatabaseVersion( $self->id, $version );

	return $query->getPublicDatabaseVersion( $self->id );

}
has public_database_version => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;

		#my $self = shift;
		my $query = $self->buffer->getQuery();
		my $id    = $query->getPublicDatabaseVersion( $self->id );
		return $id if defined $id;
		$id = $self->update_public_database_version();
		die() unless $id;
		return $id;

	},
);

has public_data_root => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->buffer->public_data_root;
	},
);

sub get_gencode_from_capture {
	my ($self) = @_;

	my $captures = $self->getCaptures();
	return -1 unless scalar( @{ $self->getCaptures() } );
	my $vs;
	map { $vs->{ $_->gencode_version }++ } @$captures;
	my @t = keys %$vs;
	confess() if ( scalar(@t) > 1 );
	my $gencode = 0;
	if ( scalar(@t) == 1 and $t[0] > 0 ) {
		$gencode = $t[0];
	}
	return $gencode;
}

sub is_gencode_upgradable {
	my ($self) = @_;
	return $self->get_gencode_from_capture() == 0;
}

sub update_gencode_version {
	my ($self) = @_;

	#read gencode from capture ?
	my $query   = $self->buffer->getQuery();
	my $gencode = $self->get_gencode_from_capture();
	if ( $gencode == -1 ) {

   #c'est le cas ou il n'y a pas encore de cpature ni de patient dans le projet;
		return -1;

		#$query->insertGencodeVersion($self->id,$gencode);
		#return $query->getGencodeVersion($self->id);
	}
	if ( $gencode == 0 ) {
		$gencode = $query->getMaxGencodeVersion( $self->id );
	}
	confess() unless $gencode;
	$query->insertGencodeVersion( $self->id, $gencode );
	return $query->getGencodeVersion( $self->id );
}

has gencode_version => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self  = shift;
		my $query = $self->buffer->getQuery();
		my $id    = $query->getGencodeVersion( $self->id );
		
		return $id if $id =~ /M/;
		
		return $id if defined $id && $id > -1;

		return $self->update_gencode_version();
	},
);

sub get_gencode_directory {
	my ( $self, $version ) = @_;
	if ($self->annotation_genome_version() =~ /MM38/) {
		return $self->buffer->config_path("root","public_data")."/repository/MM38/annotations/gencode.vM25/lmdb/";
	}
	if ($self->annotation_genome_version() =~ /MM39/) {
		return $self->buffer->config_path("root","public_data")."/repository/MM39/annotations/gencode.vM32/lmdb/";
	}
	
	my $database = "gencode";
	$version = $self->gencode_version unless $version;
	return $self->{directory}->{$version}->{$database} if exists $self->{directory}->{$version}->{$database};
	confess() unless exists $self->buffer->gencode->{$version}->{directory};
	$self->{directory}->{$version}->{$database} = $self->public_data_root . "/". $self->annotation_genome_version . "/". $self->buffer->gencode->{$version}->{directory};
	confess( "\n\nERROR: score:$database " . $self->{directory}->{$version}->{$database}."\n/n" ) unless -e $self->{directory}->{$version}->{$database};

	#$self->{directory}->{$version}->{$database} = "/tmp/pipeline/annot/lmdb/";
	#warn $self->{directory}->{$version}->{$database};
	return $self->{directory}->{$version}->{$database};
}


sub get_wisecondor_reference {
	my ( $self, $version ) = @_;
	$version = "5kb" unless ($version);
	return $self->{directory}->{wisecondor}->{$version} if exists $self->{directory}->{wisecondor}->{$version};
	$self->{directory}->{wisecondor}->{$version} = $self->public_data_root . "/". $self->annotation_genome_version . "/wisecondor/$version/ref.npz";
	return $self->{directory}->{wisecondor}->{$version};
}

sub get_cytology_directory {
	my ( $self ) = @_;
	return $self->{directory}->{cyto} if exists $self->{directory}->{cyto};
	$self->{directory}->{cyto} = $self->public_data_root . "/". $self->annotation_genome_version . "/cytology/";
	return $self->{directory}->{cyto};
}


sub get_public_data_directory {
	my ( $self, $database, $version ) = @_;
	return $self->{directory}->{$database} if exists $self->{directory}->{$database} ;
	$version = $self->public_database_version unless $version;
	#warn $database." ".$version;
	
	if (exists $self->buffer->public_data->{$version}->{$database}->{config}->{semantic}){
		
		$self->{directory}->{$database} = $self->buffer->config_path("root","public_data")."/repository/semantic/".$self->buffer->public_data->{$version}->{$database}->{config}->{directory};
	}
	else {
		$self->{directory}->{$database} = $self->public_data_root . "/". $self->annotation_genome_version . "/". $self->buffer->public_data->{$version}->{$database}->{config}->{directory};
	}
	confess( "\n\nERROR: Public data :\nDatabase: $database\nDir: " . $self->{directory}->{$database}."\n\n".Dumper ($self->buffer->public_data->{$version}->{$database}->{config}) ) unless -e $self->{directory}->{$database};
	return $self->{directory}->{$database};
}

sub get_gencode_description {
	my ( $self, $version ) = @_;
	return $self->get_public_data_description( 'gencode', $version );
}

sub get_public_data_description {
	my ( $self, $database, $version ) = @_;
	if ( $database eq 'gencode' or $database eq 'annotations' ) {
		$version = $self->gencode_version unless $version;
	}
	else {
		$version = $self->public_database_version unless $version;
	}

	return $self->{json}->{$version}->{$database}
	  if exists $self->{json}->{$version}->{$database};
	my $dir;
	if ( $database eq 'gencode' or $database eq 'annotations' ) {
		$dir = $self->get_gencode_directory($version);
	}
	else {
		$dir = $self->get_public_data_directory( $database, $version );
	}
	my $f = "$dir/description.json";
	$f = "$dir/version.json" unless -e $f;
	confess($f) unless -e $f;
	open( JSON, $f );
	my $desc = decode_json <JSON>;
	$self->{json}->{$version}->{$database} = $desc;
	return $self->{json}->{$version}->{$database};
}

sub get_public_data_version {
	my ( $self, $database, $version ) = @_;
	unless($version){
		$version = $self->public_database_version unless $version;
	}
	return $self->get_public_data_description( $database, $version )->{version};
}

#LMDB public gene based

has lmdbHPO_genes_to_phenotypes => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return GenBoNoSqlLmdb->new(
			name        => "hpo_genes_to_phenotypes",
			dir         => $self->get_public_data_directory("hpo"),
			mode        => "r",
			is_compress => 1,
			vmtouch=>$self->buffer->vmtouch
		);

#       return GenBoNoSqlLmdb->new(name=>"pLI",dir=>$sqliteDir,mode=>"r",is_compress=>1);
	}
);



has lmdbpolyScore => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return GenBoNoSqlLmdb->new(
			name        => "polyScore",
			dir         => $self->get_public_data_directory("polyScore"),
			mode        => "r",
			is_compress => 1,
			vmtouch=>$self->buffer->vmtouch
		);

#	return GenBoNoSqlLmdb->new(name=>"pLI",dir=>$sqliteDir,mode=>"r",is_compress=>1);
	}
);



has lmdbOmim => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		#warn  $self->get_public_data_directory("omim");
		
		return GenBoNoSqlLmdb->new(
			name        => "omim",
			dir         => $self->get_public_data_directory("omim"),
			mode        => "r",
			is_compress => 1,
			vmtouch=>$self->buffer->vmtouch
		);
	}
);

has getRockPartialTranscriptDir  => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->get_gencode_directory();
	}
);
has getFastaGencodeDir  => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		#confess() unless -e $self->get_gencode_directory()."/fasta";
		return $self->get_gencode_directory()."/../fasta";
	}
);
sub fastaProteinIndex {
	my($self,$version) = @_;
	die() if $version ne "HG38" && $version ne "HG19";
	return $self->{fasta_index}->{$version} if exists $self->{fasta_index}->{$version};
	my $file = $self->get_gencode_directory()."/../fasta/$version/proteins.fa.gz";
	return undef unless -e $file;
	$self->{fasta_index}->{$version} = Bio::DB::HTS::Faidx->new($file);
	return $self->{fasta_index}->{$version};
}
has getFastaProteinIndexHG38  => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return unless -e $self->get_gencode_directory()."/fasta";
	#	confess() unless -e $self->get_gencode_directory()."/fasta";
		return $self->get_gencode_directory()."/fasta";
	}
);

has lmdbPartialTranscripts => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return unless -e $self->getRockPartialTranscriptDir().'/partial_transcripts_lmdb';
		my $no = GenBoNoSqlLmdb->new(
			dir         => $self->getRockPartialTranscriptDir(),
			mode        => 'r',
			is_index    => 1,
			name        => 'partial_transcripts_lmdb',
			is_compress => 1,
		);	
		return $no;
	}
);


sub rocksPartialTranscripts{
	my $self = @_;
	return $self->{rocks}->{partial} if exists $self->{rocks}->{partial};
	$self->{rocks}->{partial} = GenBoNoSqlRocks->new(
			dir         => $self->getRockPartialTranscriptDir(),
			mode        => 'r',
			is_index    => 1,
			name        => 'partial_transcripts',
			is_compress => 1,
		);	
		return $self->{rocks}->{partial};
}





has annotation_version_history => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self  = shift;
		my $query = $self->buffer->getQuery();
		my (@list_annot, $h_annot);
		foreach my $l ( @{ $query->ListCacheHistoryVersion( $self->id() ) } ) {
			next unless ( $self->isCacheDone( $l->[0] ) );
			my ($gencode, $db_annot) = split('\.', $l->[0]);
			$h_annot->{$gencode}->{$db_annot} = undef;
		}
		foreach my $gencode (sort {$a <=> $b} keys %$h_annot) {
			foreach my $db_annot (sort {$a <=> $b} keys %{$h_annot->{$gencode}}) {
				push(@list_annot, $gencode.'.'.$db_annot);
			}
		}
		return \@list_annot;
	},
);

has annotation_version_current => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $annot_version;
		return $self->gencode_version . "."
		  . $self->buffer->getQuery()->getMaxPublicDatabaseVersion( $self->id );
	},
);

has annotation_version => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $annot_version;
		return $self->gencode_version . "." . $self->public_database_version;
	},
);

has annotation_public_path => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;

		#		confess();
		confess(
"\n\nERROR: no annotation version found for this project in DataBase. Die\n\n"
		) unless ( $self->annotation_version );
		my $dir =
			$self->buffer()->config_path("root","public_data")
		  . '/annotations/'
		  . $self->annotation_version . "/";
		confess("public_data annotation $dir") unless -e $dir;

		return $dir;

	},
);


sub get_dejavu_junctions_path {
	my ($self, $phenotype_name) = @_;
	my $dir = $self->buffer->config_path("root","dejavu") . $self->annotation_genome_version . "/" . $self->buffer()->config->{'deja_vu_JUNCTION'}->{junctions};
	$dir .= '/'.$phenotype_name.'/' if ($phenotype_name);
	
	confess("junction dejavu $dir") unless -e $dir;
	return $dir;
}

has DejaVuJunction_path => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->get_dejavu_junctions_path();
	},
);

has DejaVuCNV_path => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;

#		confess();
#confess("\n\nERROR: no annotation version found for this project in DataBase. Die\n\n") unless ($self->annotation_version);
return $self->buffer()->deja_vu_public_dir($self->annotation_genome_version,"CNV");

		#return $dir;

	},
);

has DejaVuProjectsCNV_path => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;

#		confess();
#confess("\n\nERROR: no annotation version found for this project in DataBase. Die\n\n") unless ($self->annotation_version);

		my $dir =  $self->buffer()->deja_vu_public_dir($self->annotation_genome_version,"CNV")."/projects/";
		system("mkdir $dir && chmod a+rwx $dir ") unless -e $dir;
		confess("cnv dejavu $dir") unless -e $dir;

		return $dir;

	},
);

has DejaVuCNVFile => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $file =
		  $self->DejaVuProjectsCNV_path . "/" . $self->name . ".dejavu";
		return $file;
	},
);


has DejaVuSVeq_path => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;

#		confess();
#confess("\n\nERROR: no annotation version found for this project in DataBase. Die\n\n") unless ($self->annotation_version);
		my $dir =  $self->buffer()->deja_vu_public_dir($self->annotation_genome_version,"SVeq");
		confess("sveq dejavu $dir") unless -e $dir;

		return $dir;

	},
);

has DejaVuProjectsSVeq_path => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;

#		confess();
#confess("\n\nERROR: no annotation version found for this project in DataBase. Die\n\n") unless ($self->annotation_version);
 	my $dir =$self->buffer()->deja_vu_public_dir($self->annotation_genome_version,"SVeq")."/projects";
 
confess("sveq dejavu $dir") unless -e $dir;

		return $dir;
	},
);

has DejaVuSVeqFile => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $file =
		  $self->DejaVuProjectsSVeq_path . "/" . $self->name . ".SVeqDejavu";
		return $file;
	},
);

sub changeAnnotationVersion {
	my ( $self, $annot_version, $force ) = @_;

	if ($force) {
		my @lTmp = split( '\.', $annot_version );
		$self->gencode_version( $lTmp[0] );
		$self->public_database_version( $lTmp[1] );
		$self->annotation_version($annot_version);
		delete $self->{'directory'};
		delete $self->{'annotation_public_path'};
		delete $self->{'annotationsDirectory'};
		delete $self->{'lmdb_hgmd'};
		delete $self->{'lmdbPli'};
		delete $self->{'omimDirectory'};
		delete $self->{'litePredictionMatrix'};
		return;
	}
	foreach my $history_version ( @{ $self->annotation_version_history() } ) {
		warn $history_version." ".$annot_version;
		if ( $history_version eq $annot_version ) {
			my @lTmp = split( '\.', $annot_version );
			$self->gencode_version( $lTmp[0] );
			$self->public_database_version( $lTmp[1] );
			$self->annotation_version($annot_version);
			delete $self->{'directory'};
			delete $self->{'annotation_public_path'};
			delete $self->{'annotationsDirectory'};
			delete $self->{'lmdb_hgmd'};
			delete $self->{'lmdbPli'};
			delete $self->{'omimDirectory'};
			delete $self->{'litePredictionMatrix'};
			return;
		}
	}
	my $l_annotation_version_history =
	  join( ', ', sort @{ $self->annotation_version_history() } );
	  warn Dumper $self->annotation_version_history();
	confess("\n\nERROR: Project "
		  . $self->name().join(";",@{$self->annotation_version_history() })
		  . " never analised with annotation version $annot_version. Die.\n\nAnnotation Project History: $l_annotation_version_history\n\n"
	);
}

has dirGenome => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->buffer()->config_path("root","public_data") . "/genome/"
		  . $self->genome_version . "/";
		return $dir;

	},
);

has dirCellranger => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->buffer()->config_path("root","public_data")
		  . "/cellranger";
		return $dir;

	},

);

sub getCellRangerIndex {
	my ( $self, $method ) = @_;
	my $dir;
	$dir = $self->dirCellranger . "/" . $self->genome_version . "/";
	confess() unless -e $dir;

	return $dir;
}
has annotation_genome_version => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		my $version = $self->getVersion();
		$version = "HG19" if $version =~ /HG19/;
		$version = "HG38" if $version =~ /HG38/;
		$version = "HG38" if $version =~ /MT/;
		return $version;
	}
);
has current_genome_version => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		return $self->annotation_genome_version();
	}
);
has lift_genome_version => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		return "HG38" if $self->current_genome_version() eq "HG19";
		return "HG19" if $self->current_genome_version() eq "HG38";
		confess($self->current_genome_version);
	}
);

has genome_version_generic => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		return $self->annotation_genome_version();
	#	my $version = $self->getVersion();
	#	$version = "HG19" if $version =~ /HG19/;
	#	$version = "HG38" if $version =~ /HG38/;
	#	return $version;
	}
);
has capture_dir => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->buffer()->config_path("root","public_data") . "/capture/"
		  . $self->genome_version_generic . "/";
		return $dir;

	},

);
has chain_dir => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->buffer()->config_path("root","public_data") . "/chain/";
		return $dir;

	},

);
sub liftover_chain_file {
	my ($self,$vto) = @_;
	my $vfrom = $self->genome_version_generic;
	my $file = $self->chain_dir.$self->buffer()->config->{'public_data'}->{"liftover_chain_".$vfrom."_".$vto};
	confess($file." : chain file not found ") unless -e $file;
	return $file;
	
}
has dirGenomeGeneric => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->buffer()->config_path("root","public_data") . "/genome/"
		  . $self->genome_version_generic . "/";
		return $dir;

	},

);

has gnomad_rsname_dir => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self   = shift;
		my $dir = $self->buffer->get_lmdb_database_directory('gnomad-genome').'/rsname/';
		return $dir;
	},
);






has rds_gencode_file => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $path = my $version = $self->getVersion();
		my $file = $self->buffer()->config_path("root","public_data").'repository/'.$self->annotation_genome_version.'/annotations/'.'/gencode.v'.$self->gencode_version."/rds/".$self->annotation_genome_version."_gencode".$self->gencode_version.".rds";
		die($file) unless -e $file;
		return $file;
	}
);

has rds_junctions_canoniques_gencode_file => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $path = my $version = $self->getVersion();
		my $file = $self->buffer()->config_path("root","public_data").'repository/'.$self->annotation_genome_version.'/annotations/'.'/gencode.v'.$self->gencode_version."/rds/Junc_".$self->annotation_genome_version."_gencode".$self->gencode_version.".rds";
		die($file) unless -e $file;
		return $file;
	}
);

has gtf_file => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $path = my $version = $self->getVersion();
		my $file = $self->buffer()->config_path("root","public_data").'repository/'.$self->annotation_genome_version.'/annotations/'.'/gencode.v'.$self->gencode_version."/gtf/annotation.gtf";
		$file = $self->buffer()->config_path("root","public_data").'/repository/'.$version.'/'.$self->buffer()->config->{'public_data'}->{gtf} unless -e $file;
		die($file) unless -e $file;	
		return $file;
	},
);


has gtf_file_dragen => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $path = my $version = $self->getVersion();
		my $file =
			$self->buffer()->config_path("root","public_data") . '/repository/'
		  .  $self->annotation_genome_version  . '/'
		  . $self->buffer()->config->{'public_data'}->{gtf_dragen};
		return $file;
	},
);


has gtf_file_star => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $path = my $version = $self->getVersion();
		my $file =
			$self->buffer()->config_path("root","public_data") . '/repository/'
		  . $version . '/'
		  . $self->buffer()->config->{'public_data'}->{gtf_star};
		return $file;
	},
);


has refFlat_file => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $path = my $version = $self->getVersion();
		my $file = $self->buffer()->config_path("root","public_data") 
			. 'repository/'.$self->annotation_genome_version  .'/annotations/'
		 	.   '/gencode.v'.$self->gencode_version."/refFlat/refFlat.txt";
		warn $file;
		$file = $self->buffer()->config_path("root","public_data") . '/repository/'
		 	. $version
		  	. '/refFlat/refFlat.txt' unless -e $file;
		
		return $file;
	},
); 

has refFlat_file_star => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $path = my $version = $self->getVersion();
		my $file =
			$self->buffer()->config_path("root","public_data") . '/'
		  . $version
		  . '/refFlat/refFlat_no_chr.txt';
		return $file;
	},
);

has rRNA_file => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $path = my $version = $self->getVersion();
		my $file =
			$self->buffer()->config_path("root","public_data") . '/'
		  . $version
		  . '/refFlat/rRNA.interval_list';
		return $file;
	},
);

#has ensembl_version => (
#	is      => 'rw',
#	lazy    => 1,
#	reader  => 'getEnsemblDirectory',
#	default => sub {
#		my $self    = shift;
#		my $version = $self->getVersion();
#		my $ensembl = $self->buffer()->config->{'ensembl_versions'}
#		  ->{ "ensembl_" . $version };
#		return $ensembl;
#	},
#);

has metricsHeaderFile => (
	is      => 'rw',
	lazy    => 1,
	reader  => 'getMetricsHeader',
	default => sub {
		my $self              = shift;
		my @patients          = @{ $self->getPatients() };
		my $onePatientBamFile = $patients[0]->getBamFile();
		my $metricsHeaderFile =
		  $self->getAlignmentPipelineRootDir() . "/header.txt";
		my $cmd    = "samtools view -H $onePatientBamFile > $metricsHeaderFile";
		my $header = `$cmd`;

		#my $header = "/data-xfs/public-data/HG19/capture/headerMetrics.txt";
		return $metricsHeaderFile;
	},

);

has sequencing_machines => (
	is      => 'rw',
	lazy    => 1,
	reader  => 'getSequencingMachines',
	default => sub {
		my $self      = shift;
		my $query     = $self->buffer->getQuery();
		my $projectId = $self->id();
		my $res       = $query->getSequencingMachines($projectId);
		return $res;
	},
);

has _genomeFasta => (
	is      => 'rw',
	lazy    => 1,
	reader  => 'getGenomeFasta',
	default => sub {
		my $self = shift;

		my $f = $self->dirGenome() . "/fasta/all.fa";
		$f =~ s/\/+\//\//;
		confess($f) unless -e $f;
		return $f;
	},
);

sub genomeFasta {
	my ( $self, $bam ) = @_;
	my $ref = $self->_genomeFasta();
	return $ref unless $bam;
	my $samtools = $self->buffer->software("samtools");
	my @header = `$samtools view -H  $bam`;
	chomp(@header);
	my ($pangenome) =  grep{$_ =~ /chr6_cox_hap2/} @header;
	if ($pangenome ){
		# Remplacer "HG19_*" par "HG19_DRAGEN" en utilisant des délimiteurs "|"
		$ref =~ s|/HG19_[^/]+|/HG19_DRAGEN|;
		$ref =~ s|/HG38_[^/]+|/HG38_DRAGEN|;
	}
	return $ref;
}
sub getGenomeIndex {
	my ( $self, $method ) = @_;
	my $dir;
	confess("\n\nERROR: path "
		  . $self->dirGenome()
		  . "$method doesnt exist. Die.\n\n" )
	  unless -e $self->dirGenome() . "/$method";
	$dir = $self->dirGenome() . "/$method";
	$dir .= "/hg19"   if ( $method eq "bowtie2" );
	$dir .= "/genome" if ( $method eq "hisat2" );
	$dir .= "/latest" if ( $method eq "cellranger-star");
	return $dir;
}

has genomeFai => (
	is      => 'rw',
	lazy    => 1,
	reader  => 'getGenomeFai',
	default => sub {
		my $self      = shift;
		my $fastaFile = $self->getGenomeFasta();
		
		my $file_fai = $fastaFile . ".fai";
		confess("can't find fai file for genome") unless -e $file_fai;
		#open(FAI1, $file_fai ) or die ("Can't open, $!");
		my @data = `cat $file_fai`;
		unless (@data) {
			die("problem: ".$file_fai) unless @data;
		}
		#chomp(@data);
		my $arrayChrs = [];
		#while ( my $line = <FAI> ) {
		foreach my $line (@data){
			chomp($line);
			my @data = split( " ", $line );
			my $chrfai;
			my $chr  = $data[0];
			my $ochr = $chr;
			$chr = $self->buffer()->ucsc2ensembl($chr);
			next if $chr =~ /^GL/;
			next if $chr =~ /^Un_/;
			next if $chr =~ /^hs37d5/;
			next if $chr =~ /^NC_007605/;
			next if $chr =~ /KI/;
			next if $chr =~ /GL/;
			next if $chr =~ /EBV/;
			next if $chr =~ /JH/;
			next if $chr =~ /MU0/;
			next if $chr =~ /_random/;
			next if $chr =~ /_hap/i;
			next if $chr =~ /HLA/i;
			next if $chr =~ /KB/i;
			$chrfai->{id}                 = $chr;
			$chrfai->{name}               = $chr;
			$chrfai->{fasta_name}         = $ochr;
			$chrfai->{chromosomes_object} = { $chr => undef };
			$chrfai->{length}             = $data[1];
			$chrfai->{fai}->{start}       = $data[2];
			$chrfai->{fai}->{length_line} = $data[3];
			$chrfai->{fai}->{length_char} = $data[4];
			push( @$arrayChrs, $chrfai );
		}
		
		return $arrayChrs;
	},
);

sub get_only_list_chromosomes {
	my ( $self, $chr_name ) = @_;
	my $chr = $self->getChromosomes();
	my %keep;
	map { $keep{$_}++ } split( ",", $chr_name );
	my @t;
	foreach my $c (@$chr) {

		unless ( exists $keep{ $c->name } ) {
			die() unless exists $self->{chromosomes_object}->{ $c->id };
			delete $self->{chromosomes_object}->{ $c->id };
		}
	}
	return $self->getChromosomes();

}

sub isChromosomeName {
	my ( $self, $name ) = @_;
	confess() unless $name;
	my $names;
	foreach my $chr ( @{ $self->getChromosomes } ) {
		$names->{ $chr->name }++;
		$names->{ $chr->ucsc_name }++;
	}
	return exists $names->{$name};

}

sub setChromosomes {
	my $self = shift;
	my %chrIds;
	foreach my $hchr ( @{ $self->getGenomeFai() } ) {
		my $obj = $self->flushObject( 'chromosomes', $hchr );
		$chrIds{ $obj->id } = undef;
	}
	$self->{genomeFai} = {};
	return \%chrIds;
}

sub getChromosomes {
	my $self       = shift;
	my $chr        = $self->SUPER::getChromosomes();

	my @chr_sorted = sort { $a->karyotypeId <=> $b->karyotypeId } @$chr;
	return \@chr_sorted;
}

has callingMethods => (
	is      => 'rw',
	lazy    => 1,
	reader  => 'getCallingMethods',
	default => sub {
		my $self = shift;
		my %dejavu;
		my $methods;
		foreach my $p ( @{ $self->getPatients() } ) {
			foreach my $m ( @{ $p->getCallingMethods } ) {
				$methods->{$m} = undef;
			}
		}

		#my $query = $self->buffer->getQuery();
		return [ keys %$methods ];

		#return $query->getOriginMethods( $self->id(), "SNP" );
	},
);
has callingSVMethods => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my %dejavu;
		my $methods;
		foreach my $p ( @{ $self->getPatients() } ) {
			foreach my $m ( @{ $p->callingSVMethods } ) {
				$methods->{$m} = undef;
			}
		}

		#my $query = $self->buffer->getQuery();
		return [ keys %$methods ];

		#return $query->getOriginMethods( $self->id(), "SNP" );
	},
);
has callingMethod => (
	is      => 'rw',
	lazy    => 1,
	reader  => 'getCallingMethod',
	default => 'None',
);

has variationsMethods => (
	is      => 'rw',
	lazy    => 1,
	reader  => 'getVariationMethods',
	default => sub {
		my $self  = shift;
		my $query = $self->buffer->getQuery();
		return $query->getOriginMethods( $self->id(), "SNP" );
	},
);

has alignmentMethods => (
	is      => 'rw',
	lazy    => 1,
	reader  => 'getAlignmentMethods',
	default => sub {
		my $self = shift;
		my $hAlignMeth;
		foreach my $pat ( @{ $self->getPatients() } ) {
			foreach my $methName ( @{ $pat->alignmentMethods() } ) {
				$hAlignMeth->{$methName} = undef;
			}
		}
		my @lMethods = keys(%$hAlignMeth);
		return \@lMethods;
	},
);

has patientsNameVcfParsed => (
	is      => 'rw',
	lazy    => 1,
	reader  => 'getPatientsNameVcfParsed',
	default => sub { },
);

has allPath => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self  = shift;
		my $types = $self->buffer()->config->{project_datatypes}
		  ->{ $self->getProjectType() };

		my $rootDir        = $self->getProjectPath();
		my $methodsCalling = $self->getVariationMethods();
		foreach my $type (@$types) {
			foreach my $method (@$methodsCalling) {
				$self->{dir}->{$type}->{$method} =
				  $rootDir . "" . lc($type) . "/" . lc($method) . "/";
			}
		}

		my $methodsAlign = $self->getAlignmentMethods();
		my $align        = "align";
		foreach my $m (@$methodsAlign) {
			$self->{dir}->{align}->{$m} =
			  $rootDir . "" . $align . "/" . lc($m) . "/";
		}

		my $sequenceDir       = $self->{buffer}->config_path("root","project_data");
		my $methodsSequencing = $self->getSequencingMachines();
		my @dirs;
		foreach my $m (@$methodsSequencing) {
			my $dir =
				$sequenceDir . "/"
			  . $self->getProjectType() . "/"
			  . $self->name()
			  . "/sequences/"
			  . lc($m) . "/";
			$dir =~ s/\/\//\//;
			$self->{dir}->{sequence}->{$m} = $dir;
		}

#		$self->{dir}->{pipeline}->{root} = $self->project_pipeline_path;
#		foreach my $m (@$methodsAlign){ $self->{dir}->{pipeline_align}->{$m}= $self->project_pipeline_path . "align/" . lc($m) . "/"; }
#		foreach my $m (@$methodsCalling){ $self->{dir}->{pipeline_calling}->{$m} = $self->project_pipeline_path . "calling/" . lc($m) . "/"; }

		$self->{dir}->{data}->{deja_vu} = $rootDir . "deja_vu/";

		return $self->{dir};
	},
);

has pipelineDir => (
	is      => 'ro',
	reader  => 'getPipelineDir',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->project_pipeline_path;
	},
);

has pipelineDragen => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->project_dragen_pipeline_path;
	},
);

has pipelineEpi2me => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->project_epi2me_pipeline_path;
	},
);
has metricsDir => (
	is      => 'ro',
	reader  => 'getMetricsDir',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->project_metrics_path;
	},
);

has dejavuDirTest => (
	is      => 'rw',
	lazy    => 1,
	default => undef,
);

has dejavuDir => (
	is => 'ro',

	reader  => 'getDejaVuDir',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->allPath()->{data}->{deja_vu};
		return $self->makedir($dir);

	},
);
has pedigree => (
	is      => 'rw',
	lazy    => 1,
	default => -1,
);

has pedigree_details => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $listHash =
		  $self->buffer->getQuery->getPatientProjectInfo( $self->id() );
		my ( $hTmp, $hPed );
		foreach my $hFam (@$listHash) {
			unless ( $hFam->{family} ) {
				$hFam->{family} = 'ByDefault_' . $hFam->{name};
			}
			if ( not $hFam->{status} or $hFam->{status} eq '0' ) {
				$hFam->{status} = 2;
			}
			if ( not $hFam->{sex} or $hFam->{sex} eq '0' ) { $hFam->{sex} = 1; }
			my $fam_name = $hFam->{family};
			$hTmp->{$fam_name}->{father} = $hFam->{father}
			  if ( $hFam->{father} );
			$hTmp->{$fam_name}->{mother} = $hFam->{mother}
			  if ( $hFam->{mother} );
			$hTmp->{$fam_name}->{'all'}->{ $hFam->{name} }->{'sex'} =
			  $hFam->{sex};
			$hTmp->{$fam_name}->{'all'}->{ $hFam->{name} }->{'status'} =
			  $hFam->{status};
		}

		foreach my $fam_name ( keys %$hTmp ) {
			if ( exists $hTmp->{$fam_name}->{father} ) {
				$hPed->{$fam_name}->{father} = $hTmp->{$fam_name}->{father};
			}
			if ( exists $hTmp->{$fam_name}->{mother} ) {
				$hPed->{$fam_name}->{mother} = $hTmp->{$fam_name}->{mother};
			}
			foreach my $pat_name ( keys %{ $hTmp->{$fam_name}->{all} } ) {
				if ( exists $hTmp->{$fam_name}->{father}
					and $hTmp->{$fam_name}->{father} eq $pat_name )
				{
					$hPed->{$fam_name}->{'all'}->{$pat_name}->{'sex'} =
					  $hTmp->{$fam_name}->{all}->{$pat_name}->{sex};
					$hPed->{$fam_name}->{'all'}->{$pat_name}->{'status'} =
					  $hTmp->{$fam_name}->{all}->{$pat_name}->{status};
				}
				elsif ( exists $hTmp->{$fam_name}->{mother}
					and $hTmp->{$fam_name}->{mother} eq $pat_name )
				{
					$hPed->{$fam_name}->{'all'}->{$pat_name}->{'sex'} =
					  $hTmp->{$fam_name}->{all}->{$pat_name}->{sex};
					$hPed->{$fam_name}->{'all'}->{$pat_name}->{'status'} =
					  $hTmp->{$fam_name}->{all}->{$pat_name}->{status};
				}
				else {
					$hPed->{$fam_name}->{'children'}->{$pat_name} = undef;
					$hPed->{$fam_name}->{'all'}->{$pat_name}->{'sex'} =
					  $hTmp->{$fam_name}->{all}->{$pat_name}->{sex};
					$hPed->{$fam_name}->{'all'}->{$pat_name}->{'status'} =
					  $hTmp->{$fam_name}->{all}->{$pat_name}->{status};
				}
			}
		}
		$self->pedigree(1);
		return $hPed;
	},
);

has stats => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return {} unless -e $self->getStatsFile();
		open( my $fh, '<', $self->getStatsFile );
		my $json_text = <$fh>;
		return decode_json($json_text);

	},
);

has all_families => (
	is      => 'rw',
	lazy    => 1,
	default => sub { [] },
);

#has isFamilialStudy => (
#	is		=> 'ro',
#	lazy	=> 1,
#	default => sub {
#		my $self = shift;
#		return -e $self->getPedigreeFile;
#	},
#);

has families => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self         = shift;
		my $listPatients = $self->getPatients();
		my $hashFamily;
		foreach my $pat (@$listPatients) {
			my $familyName = $pat->family;
			push( @{ $hashFamily->{$familyName}->{members} }, $pat->name() );
			if ( $pat->isHealthy ) {
				push(
					@{ $hashFamily->{$familyName}->{members_healthy} },
					$pat->name()
				);
			}
			else {
				push(
					@{ $hashFamily->{$familyName}->{members_ill} },
					$pat->name()
				);
			}
			if ( $pat->isChild ) {
				push( @{ $hashFamily->{$familyName}->{child} }, $pat->name() );
				if ( $pat->isHealthy ) {
					push(
						@{ $hashFamily->{$familyName}->{child_healthy} },
						$pat->name()
					);
				}
				else {
					push(
						@{ $hashFamily->{$familyName}->{child_ill} },
						$pat->name()
					);
				}

			}
			else {
				push( @{ $hashFamily->{$familyName}->{parents} },
					$pat->name() );
				$hashFamily->{$familyName}->{father} = $pat->name()
				  if ( $pat->isFather() );
				$hashFamily->{$familyName}->{mother} = $pat->name()
				  if ( $pat->isMother() );
			}

		}
		return $hashFamily;
	},
);

has softwares_path => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->buffer->config->{software};
	}
);

has gatk_indels_gold_file => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;

		my $dir  = $self->dirGenome();
		my $file = $dir . $self->buffer->config->{gatk_files}->{indels_gold};
		return $file if ( -e $file );
		confess();
	}
);

has gatk_genomes1k_file => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->dirGenome();
		my $file = $dir . $self->buffer->config->{gatk_files}->{genomes1k};
		return $file if ( -e $file );
		confess();
		return;
	}
);

has gatk_genomes1k_indels_file => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->dirGenome();
		my $file = $dir . $self->buffer->config->{gatk_files}->{genomes1k_indels};
		return $file if ( -e $file );
		confess();
		return;
	}
);

has gatk_genomes1k_snps_phase1_file => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->dirGenome();
		my $file =
		  $dir . $self->buffer->config->{gatk_files}->{genomes1k_snps_phase1};
		return $file if ( -e $file );
		confess();
		return;
	}
);

has gatk_hapmap_file => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->dirGenome();
		my $file = $dir . $self->buffer->config->{gatk_files}->{hapmap};
		return $file if ( -e $file );
		confess();
		return;
	}
);

has gatk_dbsnp_file => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->dirGenome();
		my $file = $dir . $self->buffer->config->{gatk_files}->{dbsnp};
		return $file if ( -e $file );
		confess();
		return;
	}
);

has gatk_omni_file => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->dirGenome();
		my $file = $dir . $self->buffer->config->{gatk_files}->{omni};
		return $file if ( -e $file );
		return;
	}
);

has gatk_mills_indels_file => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->dirGenome();
		my $file = $dir . $self->buffer->config->{gatk_files}->{mills_indels};
		return $file if ( -e $file );
		confess();
		return;
	}
);

has ensembl_annotations => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->buffer->config->{ensembl_annotations};
	}
);

has cache_verbose => (
	is      => 'rw',
	lazy    => 1,
	default => 1,
);

has impact_factors => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $h;
		foreach my $impact ( keys %{ $self->buffer->config->{impact_factors} } )
		{
			my @lTmp =
			  split( ',', $self->buffer->config->{impact_factors}->{$impact} );
			foreach my $annot (@lTmp) { $h->{$impact}->{$annot} = undef; }
		}
		return $h;
	},
);

has impacts_ensembl_annotations => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $h;
		foreach my $impact (
			keys %{ $self->buffer->config->{impacts_ensembl_annotations} } )
		{
			my @lTmp = split( ';',
				$self->buffer->config->{impacts_ensembl_annotations}->{$impact}
			);
			foreach my $annot (@lTmp) { $h->{$impact}->{$annot} = undef; }
		}
		return $h;
	},
);


has getTsoAdaptors => (
	is 		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $dir = $self->buffer->config_path("root","public_data");
		my $file = $dir . "/".$self->buffer->config->{adaptor_flexbar}->{tso};
		return $file if (-e $file);
		confess();
	}
);

has getIlluminaAdaptors => (
	is 		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $dir = $self->buffer->config_path("root","public_data");
		my $file = $dir."/".$self->buffer->config->{adaptor_flexbar}->{illumina};
		warn $file;
		return $file if (-e $file);
		confess();
	}
);



###### SET OBJECTS #####

sub getSoftware {
	my ( $self, $name ) = @_;
	return $self->softwares_path->{$name};
}

sub getRunFromId {
	my ( $self, $id ) = @_;
	$self->setRuns();
	return $self->{objects}->{runs}->{$id}
	  if exists $self->{objects}->{runs}->{$id};
	confess($id);
}

sub setPhenotypes {
	my $self = shift;
	my %hids;
	my $query = $self->buffer->queryPhenotype();
	my %hIds;
	map { $hIds{$_}++ } @{ $query->getPhenotypesIdFromProjectId( $self->id ) };
	return \%hIds;
}

sub setPanels {
	my $self = shift;
	my %hids;
	my $query = $self->buffer->queryPanel();
	my %hIds;
	map { $hIds{$_}++ } @{ $query->getAllPanelsIds() };
	return \%hIds;
}

sub setBundles {
	my $self = shift;
	my %hIds;
	foreach my $panel ( @{ $self->getPanels() } ) {
		foreach my $bundle ( @{ $panel->getBundles() } ) {
			$hIds{ $bundle->id() } = undef;
		}
	}
	return \%hIds;
}

sub setRuns {
	my $self  = shift;
	my $query = $self->buffer->getQuery();
	my $res   = $query->getRuns( $self->id );
	my %hids;
	foreach my $h (@$res) {
		$h->{name} = "run_" . $h->{id};
		$hids{ $h->{id} } = undef;
		$h->{project} = $self;
		next if exists $self->{objects}->{runs}->{ $h->{id} };
		
		$self->{objects}->{runs}->{ $h->{id} } = new GenBoRun($h);
	}
	return \%hids;
}

sub getFamiliesFromListNames {
	my ( $self, $lNames ) = @_;
	return \[] unless ($lNames);
	my @lFam;
	foreach my $fam_name (@$lNames) {
		push( @lFam, $self->getFamily($fam_name) );
	}
	return \@lFam;
}

sub getPatientsFromListNames {
	my ( $self, $lNames ) = @_;
	return \[] unless ($lNames);
	my @lPat;
	foreach my $pat_name (@$lNames) {
		push( @lPat, $self->getPatient($pat_name) );
	}
	return \@lPat;
}

sub getPatientsControl {
	my $self = shift;
	my @z    = grep { $_->is_control } @{ $self->getPatientsAndControl };
	return \@z;
}

sub getPatientsAndControl {
	my $self = shift;
	return $self->myflushobjects( $self->patients_object(), "patients" );
}

sub getPatientOrControl {
	my ( $self, $patient_name ) = @_;
	my $patient;
	foreach my $this_pat ( @{ $self->getPatientsAndControl() } ) {
		if ( $this_pat->name() eq $patient_name ) {
			$patient = $this_pat;
			last;
		}
	}
	confess("\n\nERROR: $patient_name not found... Die...\n\n")
	  unless ($patient);
	return $patient;
}

sub species_id {
	my ($self) = @_;
	$self->getPatients;
	return $self->{species_id} ;
}

sub setPatients {
	my $self  = shift;
	my $query = $self->buffer->getQuery();
	my $res   = $query->getPatients( $self->id );
	my %names;
	my $spec;
	foreach my $h (@$res) {
		$h->{id} = $h->{patient_id};
		$names{ $h->{id} } = undef;
		unless ( $h->{family} ) { $h->{family} = 'ByDefault_' . $h->{name}; }
		$h->{vstatus} = $h->{status} ;
		if ( not $h->{status} or $h->{status} eq '0' ) { $h->{status} = 2; }
		if ( not $h->{sex}    or $h->{sex} eq '0' )    { $h->{sex}    = 1; }
		$h->{project} = $self;
		$spec->{$h->{species_id}} ++;
		next if exists $self->{objects}->{patients}->{ $h->{id} };
		
		$self->{objects}->{patients}->{ $h->{id} } =
		  $self->flushObject( 'patients', $h );
		  $self->{species_id} = $h->{species_id};
		 
	}
	return \%names unless (%names);
	confess("problem species for project ".$self->name()) if scalar(keys %$spec) ne 1;
	return \%names;
}

has somatic_details => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return unless ( $self->isSomatic() );
		my $listHash =
		  $self->buffer->getQuery->getPatientSomaticProjectInfo( $self->id() );
		my $hSom;
		if ( scalar(@$listHash) > 0 ) {
			foreach my $hGroup (@$listHash) {
				$hSom->{ $hGroup->{somatic_group} }->{ $hGroup->{name} }
				  ->{group} = $hGroup->{somatic_group};
				$hSom->{ $hGroup->{somatic_group} }->{ $hGroup->{name} }
				  ->{status} = $hGroup->{status};
				if ( $hGroup->{status} == 1 ) {
					$hSom->{ $hGroup->{somatic_group} }->{ $hGroup->{name} }
					  ->{tissue} = 'C';
				}
				elsif ( $hGroup->{status} == 2 ) {
					$hSom->{ $hGroup->{somatic_group} }->{ $hGroup->{name} }
					  ->{tissue} = 'T';
				}
			}
		}
		else {
			foreach my $patient ( @{ $self->getPatients() } ) {
				my $pat_name = $patient->name();
				$hSom->{$pat_name}->{$pat_name}->{group}  = $pat_name;
				$hSom->{$pat_name}->{$pat_name}->{status} = $patient->status();
				if ( $patient->status() == 1 ) {
					$hSom->{$pat_name}->{$pat_name}->{tissue} = 'C';
				}
				elsif ( $patient->status() == 2 ) {
					$hSom->{$pat_name}->{$pat_name}->{tissue} = 'T';
				}
			}
		}
		return $hSom;
	},
);

# methode pour creer les objest GenBoGroupCache
sub setSomaticGroups {
	my $self = shift;
	my $hIds;
	my $h = $self->somatic_details();
	foreach my $group_name ( keys %$h ) {
		my $hArgs;
		$hArgs->{id} = $group_name;
		foreach my $patient_name ( keys %{ $h->{$group_name} } ) {
			my $patient = $self->getPatient($patient_name);
			$hArgs->{patients}->{$patient_name} = undef;
		}
		$hArgs->{project}         = $self->project();
		$hArgs->{chromosome}      = $self;
		$hArgs->{name}            = $group_name;
		$hArgs->{somatic_details} = $h->{$group_name};
		$hArgs->{hash_samples}    = $h->{$group_name};
		$self->{objects}->{somatic_groups}->{ $hArgs->{'id'} } =
		  $self->flushObject( 'somatic_groups', $hArgs );
		$hIds->{ $hArgs->{'id'} } = undef;
	}
	return $hIds;
}

# methode pour creer les objest GenBoFamily
sub setFamilies {
	my $self = shift;
	my $hPed = $self->pedigree_details();
	my %names;
	return \%names unless ($hPed);
	foreach my $famName ( keys %{$hPed} ) {
		my $hashArgs;
		my $hPed_fam = $hPed->{$famName};
		next if ( scalar( keys %{ $hPed_fam->{all} } ) == 0 );
		$hashArgs->{project}          = $self;
		$hashArgs->{id}               = $famName;
		$hashArgs->{name}             = $famName;
		$hashArgs->{hash_samples}     = $hPed_fam->{all};
		$hashArgs->{pedigree_details} = $hPed_fam;
		$names{ $hashArgs->{id} }     = undef;
		$self->{objects}->{families}->{ $hashArgs->{'id'} } =
		  $self->flushObject( 'families', $hashArgs );
	}
	return \%names;
}

sub getGroups {
	my $self  = shift;
	my $query = $self->buffer->getQuery();
	my $res   = $query->getGroups( $self->id );
	my $groups;
	foreach my $h (@$res) {
		my $p = $self->getPatient( $h->{pname} );
		push( @{ $groups->{ $h->{name} } }, $p );
	}
	return $groups;
}

sub setVariations {
	my $self = shift;
	my $varIds;
	foreach my $patient ( @{ $self->getPatients() } ) {
		foreach my $var ( @{ $patient->getVariations() } ) {
			$varIds->{ $var->id() } = undef;
		}
	}
	return $varIds;
}

sub setInsertions {
	my $self = shift;
	my $insIds;
	foreach my $patient ( @{ $self->getPatients() } ) {
		foreach my $ins ( @{ $patient->getInsertions() } ) {
			$insIds->{ $ins->id() } = undef;
		}
	}
	return $insIds;
}

sub setLargeInsertions {
	my $self = shift;
	my $insIds;
	foreach my $patient ( @{ $self->getPatients() } ) {
		foreach my $ins ( @{ $patient->getLargeInsertions() } ) {
			$insIds->{ $ins->id() } = undef;
		}
	}
	return $insIds;
}

sub setDeletions {
	my $self = shift;
	my $delIds;
	foreach my $patient ( @{ $self->getPatients() } ) {
		foreach my $del ( @{ $patient->getDeletions() } ) {
			$delIds->{ $del->id() } = undef;
		}
	}
	return $delIds;
}

sub setLargeDeletions {
	my $self = shift;
	my $delIds;
	foreach my $patient ( @{ $self->getPatients() } ) {
		foreach my $del ( @{ $patient->getLargeDeletions() } ) {
			$delIds->{ $del->id() } = undef;
		}
	}
	return $delIds;
}
sub setInversions {
	my $self = shift;
	my $delIds;
	foreach my $patient ( @{ $self->getPatients() } ) {
		foreach my $del ( @{ $patient->getInversions() } ) {
			$delIds->{ $del->id() } = undef;
		}
	}
	return $delIds;
}
sub setBoundaries {
	my $self = shift;
	my $delIds;
	foreach my $patient ( @{ $self->getPatients() } ) {
		foreach my $del ( @{ $patient->getBoundaries() } ) {
			$delIds->{ $del->id() } = undef;
		}
	}
	return $delIds;
}
sub setMnps {
	my $self = shift;
	my $delIds;
	foreach my $patient ( @{ $self->getPatients() } ) {
		foreach my $del ( @{ $patient->getMnps() } ) {
			$delIds->{ $del->id() } = undef;
		}
	}
	return $delIds;
}

sub setCaptures {
	my $self = shift;
	my %captIds;
	foreach my $patient ( @{ $self->getPatients() } ) {
		foreach my $thisCapture ( @{ $patient->getCaptures() } ) {
			$captIds{ $thisCapture->id() } = undef;
		}
	}
	return \%captIds;
}


sub setGenes {
	my $self    = shift;
	my $lAlLVar = $self->getStructuralVariations();
	my $hGenesId;
	foreach my $var (@$lAlLVar) {
		my $lGenes = $var->getGenes();
		foreach my $gene (@$lGenes) {
			$hGenesId->{ $gene->{genbo_id} } = undef;
		}
	}
	return $hGenesId;
}

sub setTranscripts {
	my $self = shift;
	my $hRes;
	my $lGenesId = $self->getGenes();
	foreach my $gene (@$lGenesId) {
		my $hTranscriptsId = $gene->setTranscripts();
		foreach my $id ( keys(%$hTranscriptsId) ) { $hRes->{$id} = undef; }
	}
	return $hRes;
}

###### CREATE OBJECTS #####

sub flushObjects {
	my ( $self, $type, $lHashArgs ) = @_;
	my @lObjects;
	my $print = $self->print_waiting();
	my $nb    = 1;
	foreach my $thisHashArgs (@$lHashArgs) {
		$nb++;
		my $object = $self->flushObject( $type, $thisHashArgs );
		push( @lObjects, $object );
		print "+" if $print && $nb % 2000 == 0;
	}

	return \@lObjects;
}

###### CREATE OBJECTS StructuralVariants #####

sub flushObjectSV {
	my ( $self, $hashArgs ) = @_;
	my $type = $hashArgs->{type};
	my $object;
	die() unless exists $hashArgs->{chromosome};
	my $chromosome = $self->getChromosome( $hashArgs->{chromosome} );
	my $intspan    = Set::IntSpan::Fast::XS->new(
		$hashArgs->{start} . "-" . $hashArgs->{end} );
	my $l1 = abs( $hashArgs->{start} - $hashArgs->{end} ) + 1;
	foreach my $obj ( values %{ $self->{objects}->{$type} } ) {
		next if $chromosome->name ne $obj->getChromosome->name();
		my $inter1 = $intspan->intersection( $obj->getGenomicSpan );
		next if $inter1->is_empty;
		my $li = $self->buffer->Intspan_length($inter1);
		next if ( $li < ( 0.75 * $l1 ) ) and ( $li < 0.75 * $obj->length );
		return $obj;

	}
	$self->createObject( $type, $hashArgs );

	return $self->{objects}->{$type}->{ $hashArgs->{id} };
}

sub flushObject {
	my ( $self, $type, $hashArgs ) = @_;
	my $object;
	unless ( exists $self->{objects}->{$type}->{ $hashArgs->{id} } ) {
		$self->createObject( $type, $hashArgs );
	}
	
	$object = $self->{objects}->{$type}->{ $hashArgs->{id} };
	return $object;
}

sub define_filter {
	my ( $self, $filters ) = @_;

	if ( exists $filters->{panels} ) {
		foreach my $panel ( @{ $filters->{panels} } ) {
			$self->addPanel($panel);
		}
	}
	if ( exists $filters->{transcripts} ) {
		$self->addTranscripts( $filters->{transcripts} );
	}
	if ( exists $filters->{patients} ) {
		$self->get_only_list_patients( @{ $filters->{patients} } );
	}

}

sub genes_panels {
	my ( $self, $chr ) = @_;
	die() unless exists $self->{genes_panels};
	my $os = [];
	foreach my $p ( keys %{ $self->{genes_panels} } ) {
		next unless exists $self->{genes_panels}->{$p}->{ $chr->name };
		push(
			@$os,
			@{
				$self->myflushobjects(
					$self->{genes_panels}->{$p}->{ $chr->name }, "genes" )
			}
		);
	}
	return $os;
}

sub addTranscripts {
	my ( $self, $ids ) = @_;
	$self->getChromosomes();
	my $rids = [];
	my $name = "panel";
	foreach my $id (@$ids) {
		my $tid = $self->liteAnnotations->get( "synonyms", $id );
		push( @$rids, $tid );
	}
	my $objs = $self->myflushobjects( $rids, "transcript" );
	foreach my $o (@$objs) {
		my $chr = $o->getChromosome();
		push(
			@{ $self->{genes_panels}->{$name}->{ $chr->id } },
			$o->getGene->id
		);
		push( @{ $chr->{set_genes} },       $o->getGene->id );
		push( @{ $chr->{set_transcripts} }, $o->id );
	}

}

sub setPanel {
	my ( $self, @panel ) = @_;
	my $ids;

	foreach my $name (@panel) {
		my $id = $self->buffer->getQuery()->getCaptureId($name);
		push( @$ids, $id );
	}
	foreach my $p ( @{ $self->getPatients } ) {
		$p->captures_object($ids);
	}
	$self->captures_object($ids);
	$self->isExome(undef);
}

sub addPanel {
	my ( $self, $name ) = @_;
	my $query = $self->buffer->getQuery();
	my $ids;
	foreach my $t (
		keys %{ $query->getCaptureTranscriptsbyName($name)->{transcripts} } )
	{
		my $id = $self->liteAnnotations->get( "synonyms", $t );
		push( @$ids, $id );
	}

	#die();
	$self->getChromosomes();
	my $objs = $self->myflushobjects( $ids, "transcripts" );
	foreach my $o (@$objs) {
		my $chr = $o->getChromosome();
		push(
			@{ $self->{genes_panels}->{$name}->{ $chr->id } },
			$o->getGene->id
		);
		push( @{ $chr->{set_genes} }, $o->getGene->id );
		$chr->{hset_genes}->{ $o->getGene->id }++;
		push( @{ $chr->{set_transcripts} }, $o->id );
	}
	return $ids;

}

sub newTranscripts {
	my ( $self, $ids ) = @_;
	$self->getChromosomes();
	return $self->myflushobjects( $ids, "transcripts" );
}

sub newRegulatoryRegions {
	my ( $self, $ids ) = @_;
	$self->getChromosomes();
	return $self->myflushobjects( $ids, "regulatory_regions" );
}

sub newGenes {
	my ( $self, $ids ) = @_;
	$self->getChromosomes();

	return $self->myflushobjects( $ids, "genes" );
}

sub newPhenotype {
	my ( $self, $id ) = @_;
	return $self->myflushobjects( [$id], "phenotypes" )->[0];
}
sub newRegulatoryRegion {
	my ( $self, $id ) = @_;
	$self->getChromosomes();
	return $self->myflushobjects( [$id], "regulatory_regions" )->[0];
}
sub newTranscript {
	my ( $self, $id ) = @_;
	$self->getChromosomes();
	if (not $id =~ /_.+/) {
		foreach my $chr (@{$self->getChromosomes()}) {
			my $id2 = $id.'_'.$chr->id();
			if ($self->rocksGenBo->synonym($id2)) {
				$id = $id2;
				last;
			}
		}
	}
	#	my $id1 =  $self->getGenBoId($id);
	#	return undef unless $id1;
	my $s = $self->myflushobjects( [$id], "transcripts" );
	return $s->[0] if $s ;
	return;

}

sub newGene {
	my ( $self, $id ) = @_;
	confess if not $id;
	return $self->{objects}->{genes}->{$id} if (exists $self->{objects}->{genes}->{$id});
	$self->getChromosomes();
	my $id1 = $self->rocksGenBo->synonym($id);
	return undef unless $id1;
	return $self->myflushobjects( { $id1 => undef }, "genes" )->[0];
}

sub newIntergenic {
	my ( $self, $chr_name, $start, $end ) = @_;
	my $name = 'intergenic_' . $chr_name . '_' . $start . '_' . $end;
	my $hashArgs;
	$hashArgs->{id}            = $name;
	$hashArgs->{name}          = $name;
	$hashArgs->{strand}        = '1';
	$hashArgs->{external_name} = $name;
	$hashArgs->{genbo_id}      = $name;
	$hashArgs->{kyotoId}       = $name;
	$hashArgs->{project}       = $self->project();
	$hashArgs->{chromosome}    = $self;
	$hashArgs->{description} =
	  'intergenic from ' . $start . ' to ' . $end . ' in chr' . $chr_name;
	$hashArgs->{start}       = $start;
	$hashArgs->{end}         = $end;
	$hashArgs->{transcripts} = [];
	return $self->flushObject( 'genes', $hashArgs );
}

sub newVariant {
	my ( $self, $id ) = @_;
	my @lFieldsId = split( '_', $id );
	my $chrName   = $lFieldsId[0];
	my $refAll    = $lFieldsId[2];
	my $varAll    = $lFieldsId[3];
	my $toto;
	$toto->{$id} = undef;
	$self->getChromosomes();
	my $typeObject = 'variations';
	confess($id) unless $refAll;
	if    ( length($refAll) > length($varAll) ) { $typeObject = 'deletions'; }
	elsif ( length($refAll) < length($varAll) ) { $typeObject = 'insertions'; }
	return $self->myflushobjects( $toto, $typeObject )->[0];
}

sub _newVariantFromRockdbId {
	my ($self, $chr, $rocksid, $patient) = @_;
	my ($pos, $alt) = split('!', $rocksid);
	my $ref = $chr->sequence($pos, $pos);
	my $var_id = $chr->id.'_'.int($pos);
	if ($alt =~ /\+/) {
		$alt = $ref.$alt;
		$alt =~ s/\+//;
	}
	if ($alt =~ /[0-9]+/) {
		my $alt_del = $chr->sequence($pos, $pos+int($alt));
		$var_id .= '_'.$alt_del.'_'.$ref;
	}
	else {
		$var_id .= '_'.$ref.'_'.$alt;
	}
	return $self->_newVariant($var_id, $patient);
}

sub _newVariant {
	my ( $self, $id, $patient ) = @_;
	my $hash;
	my $type;
	my $strucType;
	my ( $chr_name, $start, $ref, $alt ) = split( '_', $id );
	$alt ="" unless $alt;
	$hash->{id}         = $id;
	$hash->{annex}      = undef;
	#$hash->{line_infos} = "";

	$hash->{start}  = $start;
	$hash->{end}    = $hash->{start};
	$hash->{strand} = 1;
	$hash->{vcf_id} = $id;

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
		$hash->{end}        = $start;
	}
	$hash->{structuralTypeObject} = $type;
	$hash->{structuralType}       = $strucType;
	my $chr = $self->getChromosome($chr_name);

	$hash->{chromosomes_object} = { $chr->id() => undef };
	if ($patient) {
		my $y;
		$y->{del}->{$id} = $hash;
		$y->{del}->{$id}->{annex}->{ $patient->id } = undef;

		delete $patient->{deletions_object}->{none} if $type eq "deletions";

		delete $patient->{insertions_object}->{none} if $type eq "insertions";
		delete $patient->{variations_object}->{none} if $type eq "variations";
		my $o = $patient->myflushobjects($y)->[0];
		$patient->{variations_object}->{ $o->id } = undef if $o->isVariation();
		$patient->{deletions_object}->{ $o->id }  = undef if $o->isDeletion();
		$patient->{insertions_object}->{ $o->id } = undef if $o->isInsertion();
		$patient->{variants_object}->{ $o->id }   = undef;

		#$o->annex()->{$patient->id} = 1;
		$o->annex()->{ $patient->id }->{ho}         = "-";
		$o->annex()->{ $patient->id }->{nb_all_ref} = "?";
		$o->annex()->{ $patient->id }->{nb_all_mut} = "?";
		$o->{references_object}->{ $self->id }      = undef;

		$chr->getReference->{variations_object}->{ $o->id } = undef
		  if $o->isVariation();
		$chr->getReference->{deletions_object}->{ $o->id } = undef
		  if $o->isDeletion();
		$chr->getReference->{insertions_object}->{ $o->id } = undef
		  if $o->isInsertion();
		return $o;
	}
	my $obj = $self->flushObject( $type, $hash );
}

sub getVariantFromId {
	my ( $self, $id ) = @_;
	my @lFieldsId = split( '_', $id );
	my $chrName   = $lFieldsId[0];
	my $refAll    = $lFieldsId[2];
	my $varAll    = $lFieldsId[3];
	my @find;
#	confess($refAll." ".$varAll);

	if ( length($refAll) == length($varAll) ) {
		my $refObj =
		  $self->getChromosome($chrName)
		  ->getReferences( int( $lFieldsId[1] ), int( $lFieldsId[1] ) + 1 )
		  ->[0];
		(@find) = grep { $_->id() eq $id } @{ $refObj->getVariations() };
	}
	elsif ( length($refAll) < length($varAll) ) {
		my $refObj = $self->getChromosome($chrName)->getReferences(
			int( $lFieldsId[1] ) - length($varAll),
			int( $lFieldsId[1] ) + length($varAll)
		)->[0];
		my @lIns;
		push( @lIns, @{ $refObj->getInsertions() } );
		push( @lIns, @{ $refObj->getLargeDuplications() } );
		(@find) = grep { $_->id() eq $id } @lIns;
	}
	elsif ( length($refAll) > length($varAll) ) {
		my $refObj = $self->getChromosome($chrName)->getReferences(
			int( $lFieldsId[1] ) - length($refAll),
			int( $lFieldsId[1] ) + length($refAll)
		)->[0];
		my @lDel;
		push( @lDel, @{ $refObj->getDeletions() } );
		push( @lDel, @{ $refObj->getLargeDeletions() } );
		(@find) = grep { $_->id() eq $id } @lDel;
	}
	return $find[0] if scalar(@find) == 1;
	warn "\n\nERROR: no variant found with ID $id. Exit...\n\n";
	#return;
	confess("\n\nERROR: no variant found with ID $id. Exit...\n\n");
}

has hash_patients_name => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		my $h;
		my $id = 1;
		foreach my $patient (@{$self->getPatients()}) {
			$h->{$patient->name()}->{id} = $id;
			$h->{$patient->name()}->{fam} = $patient->family();
			$h->{by_id}->{$id} = $patient->name();
			$id++;
		}
		return $h;
	},
);

sub myflushobjects {
	my ( $self, $ids, $type ) = @_;
	my $array_ids;
	if ( ref($ids) eq 'HASH' ) {
		if ( exists $ids->{none} ) {
			return [];
		}
		$array_ids = [ keys %$ids ];
	}
	elsif ( ref($ids) eq 'ARRAY' ) {
		$array_ids = $ids;
	}
	else {
		confess($ids);
	}
	my @objs = ();

	#confess if $ids =~ /ENSG0/;
	#unless (exists $ids->{none}) {
	foreach my $id (@$array_ids) {
		#warn Dumper  $self->{objects}->{$type};
		
		unless ( exists $self->{objects}->{$type}->{$id} ) {
			if (   $type eq 'genes'
				or $type eq 'transcripts'
				or $type eq 'proteins'
				or $type eq 'exons'
				or $type eq 'introns' )
			{
				
				my $obj;
#				if ($self->getVersion() =~ /HG38/) {
					$obj = $self->rocksGenBo->genbo($id);
#				}
#				else {
#					$obj = $self->lmdbGenBo->get($id);
#					unless ($obj) {
#						if    ( $id =~ /_X/ ) { $id =~ s/_X/_Y/; }
#						elsif ( $id =~ /_Y/ ) { $id =~ s/_Y/_X/; }
#						$obj = $self->lmdbGenBo->get(
#							$self->liteAnnotations->get( "synonyms", $id ) );
#					}
#					unless ($obj) {
#						my $sid = $self->getGenBoId($id);
#						confess( $type . " " . $sid ) unless ($sid);
#						$obj = $self->lmdbGenBo->get($id);
#						$id  = $sid;
#					}
#					unless ($obj) {
#						my @t = split( "_", $id );
#						$id  = $t[0] . "_" . $t[-1];
#						$obj = $self->lmdbGenBo->get($id);
#					}
#				}
				
				confess( $id . " " . $type ) unless $obj;
				$obj->{project}                  = $self;
				$obj->{buffer}                   = $self->buffer;
				$self->{objects}->{$type}->{$id} = $obj;
			}
			elsif ( $type eq 'variations' ) { $self->getVariantFromId($id); }
			elsif ( $type eq 'deletions' )  { $self->getVariantFromId($id); }
			elsif ( $type eq 'indels' )     { die(); }
			elsif ( $type eq 'insertions' ) { $self->getVariantFromId($id); }
			elsif ( $type eq 'boundaries' ) { $self->getVariantFromId($id); }
			elsif ( $type eq 'large_duplications' ) {
				$self->getVariantFromId($id);
			}
			elsif ( $type eq 'large_deletions' ) {
				$self->getVariantFromId($id);}
			elsif ( $type eq 'large_insertions' ) {
				warn Dumper keys %{$self->{objects}->{$type}};
				confess($id);

				$self->getVariantFromId($id);
			}
			elsif ( $type eq 'mnps' ) {
				warn $id . "-";
				$self->getVariantFromId($id);
			}
			elsif ( $type eq 'runs' ) { $self->getRunFromId($id); }
			elsif ( $type eq 'patients' ) {
				$self->setPatients();
				confess() unless exists $self->{objects}->{$type}->{$id};
			}
			elsif ( $type eq 'captures' ) {

				#$self->setCaptures();
				$self->createObject( $type, { id => $id } );
				confess() unless exists $self->{objects}->{$type}->{$id};
			}
			elsif ( $type eq 'panels' ) {
				$self->createObject( $type, { id => $id } );
				confess() unless exists $self->{objects}->{$type}->{$id};
			}
			elsif ( $type eq 'phenotypes' ) {
				$self->createObject( $type, { id => $id } );
				confess() unless exists $self->{objects}->{$type}->{$id};
			}
			elsif ( $type eq 'bundles' ) {
				$self->createObject( $type, { id => $id } );
				confess() unless exists $self->{objects}->{$type}->{$id};
			}
			elsif ( $type eq 'regulatory_regions' ) {
				my $obj = $self->createObject( $type, { id => $id } );
				$id = $obj->id;
				confess($type." ".$id) unless exists $self->{objects}->{$type}->{$obj->id};
			}
			elsif ( $type eq 'junctions' ) {
				my $obj = $self->createObject( $type, { id => $id } );
				$id = $obj->id;
				confess($type." ".$id) unless exists $self->{objects}->{$type}->{$obj->id};
			}
			else {
				warn "\n\nproblem... no $type";
				confess("je fais quoi ici $type");
			}

		}

		#confess($id) unless (exists  $self->{objects}->{$type}->{$id});
		push( @objs, $self->{objects}->{$type}->{$id} );

	}

	#}

	return \@objs;
}

#sub myflushobjects_old {
#	my ( $self, $ids, $type ) = @_;
#	my $array_ids;
#	if ( ref($ids) eq 'HASH' ) {
#		if ( exists $ids->{none} ) {
#			return [];
#		}
#		$array_ids = [ keys %$ids ];
#	}
#	elsif ( ref($ids) eq 'ARRAY' ) {
#		$array_ids = $ids;
#	}
#	else {
#		confess($ids);
#	}
#	my @objs = ();
#
#	#confess if $ids =~ /ENSG0/;
#	#unless (exists $ids->{none}) {
#	foreach my $id (@$array_ids) {
#
#		unless ( exists $self->{objects}->{$type}->{$id} ) {
#			if    ( $type eq 'genes' )       { $self->setKyotoGene($id); }
#			elsif ( $type eq 'transcripts' ) { $self->setKyotoTranscript($id); }
#			elsif ( $type eq 'proteins' )    { $self->setKyotoProtein($id); }
#			elsif ( $type eq 'variations' )  { $self->getVariantFromId($id); }
#			elsif ( $type eq 'deletions' )   { $self->getVariantFromId($id); }
#			elsif ( $type eq 'indels' )      { die(); }
#			elsif ( $type eq 'insertions' )  { $self->getVariantFromId($id); }
#			elsif ( $type eq 'large_duplications' ) {
#				$self->getVariantFromId($id);
#			}
#			elsif ( $type eq 'large_deletions' ) {
#				$self->getVariantFromId($id);
#			}
#			elsif ( $type eq 'mnps' ) {
#				warn $id . "-";
#				$self->getVariantFromId($id);
#			}
#			elsif ( $type eq 'runs' ) { $self->getRunFromId($id); }
#			elsif ( $type eq 'patients' ) {
#				$self->setPatients();
#				confess() unless exists $self->{objects}->{$type}->{$id};
#			}
#			elsif ( $type eq 'captures' ) {
#
#				#$self->setCaptures();
#				$self->createObject( $type, { id => $id } );
#				confess() unless exists $self->{objects}->{$type}->{$id};
#			}
#			else {
#				confess("je fais quoi ici $type");
#			}
#
#		}
#
#		#confess($id) unless (exists  $self->{objects}->{$type}->{$id});
#		push( @objs, $self->{objects}->{$type}->{$id} );
#
#	}
#
#	#}
#
#	return \@objs;
#}

sub flushMemoryForObject {
	my ( $self, $type, $id ) = @_;
	delete $self->{objects}->{$type}->{$id};
}

has hashTypeObject => (
	is      => 'ro',
	default => sub {
		my $self           = shift;
		my $hashTypeObject = {
			'phenotypes'         => 'GenBoPhenotype',
			'panels'             => 'GenBoPanel',
			'bundles'            => 'GenBoBundle',
			'captures'           => 'GenBoCapture',
			'variations'         => 'GenBoVariation',
			'deletions'          => 'GenBoDeletion',
			'insertions'         => 'GenBoInsertion',
			'large_deletions'    => 'GenBoLargeDeletion',
			'large_duplications' => 'GenBoLargeDuplication',
			'large_insertions'	 => 'GenBoLargeInsertion',
			'inversions' 	 	 => 'GenBoInversion',
			'boundaries' 	 	 => 'GenBoBoundary',
			'references'         => 'GenBoReference',
			'transcripts'        => 'GenBoTranscript',
			'proteins'           => 'GenBoProtein',
			'positions'          => 'Position',
			'exons'              => 'GenBoExon',
			'introns'            => 'GenBoIntron',
			'primers'            => 'GenBoPrimer',
			'chromosomes'        => 'GenBoChromosome',
			'genes'              => 'GenBoGene',
			'mnps'               => 'GenBoMnp',
			'multiplexes'        => 'GenBoPcrMultiplex',
			'part_chromosomes'   => 'GenBoPartChromosome',
			'complex'            => 'GenBoComplex',
			'families'           => 'GenBoFamily',
			'patients'           => 'GenBoPatient',
			'somatic_groups'     => 'GenBoSomaticGroup',
			'svduplications'     => 'GenBoSVDup',
			'svdeletions'        => 'GenBoSVDel',
			'regulatory_regions' => 'GenBoRegulatoryRegion',
			'junctions'			 => 'GenBoJunction',
			'meis'				 =>'GenBoMei',
		};
		return $hashTypeObject;
	}
);

has time_test => (
	is      => 'rw',
	default => sub {
		return 0;
	}
);

has check_ped_db_only => (
	is      => 'rw',
	lazy    => 1,
	default => 0,
);

sub createObject {
	my ( $self, $type, $hash ) = @_;
	my $hashTypeObject = $self->hashTypeObject();
	confess("\n\nERROR: No type defined to create object -$type- !! die...\n\n".Dumper $hash)  unless exists $hashTypeObject->{$type};
	$hash->{project} = $self;

	my $typeObj = $hashTypeObject->{$type};
	if ( $type eq 'captures' ) {
		my $query = $self->buffer->getQuery();
		my $hash1 = $query->getCaptureInfos( $hash->{id} );

		#die();
		$hash1->{method} = "capture" unless exists $hash1->{method};
		$hash1->{method} = "capture" unless defined $hash1->{method};

		if ( lc( $hash1->{method} ) eq 'pcr' ) {
			$typeObj = $hashTypeObject->{multiplexes};
		}
		my $hash2 = {};
		if ( exists $hash1->{analyse} ) {
			if ( $hash1->{analyse} ne "exome" && $hash1->{analyse} ne "genome" )
			{
				$hash2 = $query->getCaptureTranscripts( $hash->{id} );
			}
		}

		#my $hash2 = $query->getCaptureTranscripts($hash->{id});
		my %newHash = ( %$hash1, %$hash2 );
		$hash->{infos} = \%newHash;

	}
	elsif ( $type eq 'regulatory_regions' ) {
	 	$hash = $self->liteRegulations->get( "annotations", $hash->{id} );
	 	$hash->{project} = $self;
	 	$hash->{buffer} = $self->buffer;
		
	}
	
	
	my $z = 0;
	$z = time;
	my $object;
	if ($type eq "insertions" && exists $hash->{isMei}) {
		$object = $hashTypeObject->{meis}->new($hash);
	#	$object = $self->get_void_object($type);
	#	die();
	}
	elsif ($type eq "variations"
		or $type eq "deletions"
		or $type eq "insertions"
		or $type eq "large_duplication"
		or $type eq "large_deletions" 
		or $type eq "large_insertions"
		or $type eq "inversions" )
	{
		$object = $self->get_void_object($type);
		foreach my $k ( keys %$hash ) {
			$object->{$k} = $hash->{$k};
		}
	}
	else {
		$object = $typeObj->new($hash);
	}

	#	#my $zz = dclone($object);
	#if ( $type eq "variations" ) {
	#	my $a = abs( time - $z );
	#	$self->time_test( $self->time_test + $a );
	#}

	#$object->{project} = undef;
	
	$self->{objects}->{$type}->{ $hash->{id} } = $object;
	
	delete $hash->{project};
	delete $hash->{buffer};
	return $object;
	#$object->{project} = $self->encoder->encode($object);
}

sub get_void_object {
	my ( $self, $type ) = @_;
	
	return dclone( $self->{void}->{$type} ) if exists $self->{void}->{$type};
	my $typeObj = $self->hashTypeObject()->{$type};
	$self->{void}->{$type} = $typeObj->new( id => "titi");
	return dclone( $self->{void}->{$type} );
}

###### METHODS #####

my $maskDB = {
	'dbsnp'       => 1,
	'1000genomes' => 2,
	'pheno_snp'   => 4,
	'evs'         => 8,
	'none'        => 16,
	'1percent'    => 32,
	'new'         => 1024
};

sub setNoPedigree {
	my $self    = shift;
	my $lPatObj = $self->getPatients();
	my $nb_fam  = 1;
	foreach my $p (@$lPatObj) {

		#if ($p->compute_sex() ne -1)
		#$p->sex($p->compute_sex());
		$p->status(1);
		$p->family( "F" . $nb_fam );
		$nb_fam++;
		$p->children(1);
	}
}

sub testPedigree {
	my $self = shift;
	my $hPatNames;
	my @errors;
	foreach my $pat ( @{ $self->getPatients() } ) {
		$hPatNames->{ $pat->name() }->{indb} = 1;
	}
	my $ped_file;
	unless ($ped_file) { $ped_file = $self->getPedigreeFile(); }
	die() unless -e $ped_file;
	open( PED, $ped_file );
	my $line = 1;
	while (<PED>) {
		chomp($_);
		my $line    = $_;
		my @lFields = split( " ", $line );
		my $patName = $lFields[1];
		$hPatNames->{$patName}->{nb}++;
		$hPatNames->{$patName}->{father} = $lFields[2] if $lFields[2] ne '0';
		$hPatNames->{$patName}->{mother} = $lFields[3] if $lFields[3] ne '0';
		$hPatNames->{$patName}->{sex}    = $lFields[4];
		$hPatNames->{$patName}->{fam}    = $lFields[0];
		$hPatNames->{$patName}->{indb}   = undef;
		$hPatNames->{$patName}->{line} .= $line . " ";

		#$line++;
	}
	close(PED);
	my $errors;
	my $type;
	foreach my $patName ( keys(%$hPatNames) ) {

		if ( $hPatNames->{$patName}->{indb} ) {
			$type = "Not in pedigree";
			push( @{ $errors->{$type} }, $patName );

		}

		$type = "Duplicate name ";

		if ( $hPatNames->{$patName}->{nb} > 1 ) {
			push( @{ $errors->{$type} }, $patName );
		}
	}
	$type = "Pedigree Structure Error ";

	foreach my $patName ( keys(%$hPatNames) ) {
		if (   !$hPatNames->{$patName}->{father}
			&& !$hPatNames->{$patName}->{mother} )
		{
			next;
		}
		my $father = $hPatNames->{$patName}->{father};
		my $mother = $hPatNames->{$patName}->{mother};
		my $fam    = $hPatNames->{$patName}->{fam};
		my $text1  = "problem with Fam: - $fam -";
		if ( !( exists $hPatNames->{$father} )
			|| $hPatNames->{$father}->{fam} ne $hPatNames->{$patName}->{fam} )
		{
			push(
				@{ $errors->{$type} },
" $fam : \t  -$patName-  \t FATHER : '$father' not in the pedigree "
			);
		}
		if ( !( exists $hPatNames->{$mother} )
			|| $hPatNames->{$mother}->{fam} ne $hPatNames->{$patName}->{fam} )
		{
			push(
				@{ $errors->{$type} },
				"$fam \t  -$patName- \t  MOTHER : '$mother' not in the pedigree"
			);
		}

		if ( $hPatNames->{$father}->{sex} ne 1 ) {
			push(
				@{ $errors->{$type} },
				"$fam -> '$father' declared as FATHER but sex is 2"
			);
		}
		if ( $hPatNames->{$mother}->{sex} ne 2 ) {
			push(
				@{ $errors->{$type} },
				"$fam -> '$mother' declared as MOTHER but sex is 1 "
			);
		}

	}
	return $errors;
}

sub getSomaticFile {
	my $self = shift;
	my $file = $self->getProjectRootPath() . '/' . $self->name() . '.somatic';
	return $file;
}

sub getPedigreeFile {
	my $self = shift;
	my $file = $self->getProjectRootPath() . '/' . $self->name() . '.ped';
	return $file;
}

sub getStatsFile {
	my $self = shift;
	my $file = $self->getProjectRootPath() . '/' . $self->name() . '.stat.json';
	$file =~ s/\/\//\//;
	return $file;
}





sub makePath {
	my $self = shift;
	$self->allPath unless exists $self->{dir};

	my $dir = $self->buffer()->config_path("root","project_data");
	my $dd =
		$dir . ""
	  .  $self->getProjectType . "/"
	  . $self->name();

	my @dirs;
	foreach my $k ( keys %{ $self->{dir} } ) {
		push( @dirs, values %{ $self->{dir}->{$k} } );
	}
	make_path( @dirs, { mode => 0777 } );
	foreach my $dd ( sort { length($a) <=> length($b) } @dirs ) {
		chmod( 0777, $dd );
		system("chmod -R a+rwx $dd");
	}
	system("chmod -R a+rwx $dd");
	foreach my $p ( @{ $self->getPatients() } ) {
		my $methods  = $p->getCallingMethods();
		my   $methods2  = $p->callingSVMethods();
		foreach my $method_name (@$methods) {

			$self->getVariationsDir($method_name);
			$self->getIndelsDir($method_name);
		}
		
		foreach my $method_name (@$methods2) {
			$self->getVariationsDir($method_name);
		}
	}
}

sub getAlignmentUrl {
	my ( $self, $method_name ) = @_;
	my $path = "/NGS/"
	  . $self->name . "/"
	  . $self->genome_version
	  . "/align/"
	  . $method_name . '/';
	return $path;
}

sub getCellRangerDir {
	my ($self) = @_;
	my $path = $self->project_path . "/cellranger/";
	return $self->makedir($path);
}

sub getPipelineTrackingDir {
	my ($self) = @_;
	my $path = $self->getProjectRootPath() . "/tracking/";
	return $self->makedir($path);
}

sub getSVDir {
	my ( $self, $method_name ) = @_;
	my $path = $self->project_path . "/SV/";
	$self->makedir($path);
	if ($method_name){
		$path .= $method_name . '/';
		return $self->makedir($path);
	}
	return $path;
}

sub getCNVDir {
	my ( $self, $method_name ) = @_;
	my $path = $self->getSVDir . "/CNV/";
	return $self->makedir($path);
}

sub getSVeqDir {
	my ( $self, $method_name ) = @_;
	my $path = $self->getSVDir . "/SVeq/";
	return $self->makedir($path);
}

sub project_path_version {
	my ($self,$version) = @_;
	my $path = $self->getProjectRootPath() . "/" .$version . "/";
	return $self->makedir($path);
}

sub getAlignmentRootDirWithVersion {
	my ($self,$version) = @_;
	my $path = $self->getProjectRootPath() . "/".$version. "/align/";
	return $self->makedir($path);
}
##
sub getAlignmentStatsDir {
	my ( $self, $method_name,$version ) = @_;
	confess() unless $method_name;
	my $path = $self->getAlignmentDir($method_name) . "/stats/";
	return $self->makedir($path);
}



sub getAlignmentDir {
	my ( $self, $method_name,$version ) = @_;
	confess() unless $method_name;
	my $path = $self->getAlignmentRootDir($version) . "/" . $method_name . '/';
	return $self->makedir($path);
}


sub getAlignmentDirName {
	my ($self) = @_;
	my $path = $self->project_path . "/align/";
	return $path;

	#return $self->makedir($path);
}


sub existsAlignmentDir {
	my ($self) = @_;
	my $path = $self->project_path . "/align/";
	return -e $path;

	#return $self->makedir($path);
}



has CellRangerDir => (
	is => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $path = $self->project_path . "/cellranger/";
		return $path;

	},
);

sub getCellRangerRootDir {
	my ($self) = @_;
	my $path = $self->CellRangerDir;
	return $self->makedir($path);
}
sub getAlignmentRootDir {
	my ($self,$version) = @_;
	my $path;
	if ($version){
	 $path = $self->getProjectRootPath() . "/".$version. "/align/";
	}
	else {
	 $path = $self->project_path . "/align/";
	}
	return $self->makedir($path);
}

sub getRecalDir {
	my ( $self, $method_name ) = @_;
	my $path = $self->getAlignmentRootDir . "/" . $method_name . '/recal';
	return $self->makedir($path);
}

sub getAlignmentPipelineDir {
	my ( $self, $method_name ) = @_;
	my $path = $self->getAlignmentPipelineRootDir . "/" . $method_name . '/';
	return $self->makedir($path);
}

sub getAlignmentPipelineRootDir {
	my ( $self, $method_name ) = @_;
	my $path = $self->project_pipeline_path . "/align/";
	return $self->makedir($path);
}

#TODO: etrange ici...
sub getCountingDir {
	my ( $self, $method_name ) = @_;
	my $path = $self->project_path . "/count/";
	$self->makedir($path);
	$path .= $method_name . '/';
	return $self->makedir($path);
}

sub getVariationsDir {
	my ( $self, $method_name ) = @_;
	my $path = $self->project_path . "/variations/";
	$self->makedir($path);
	$path .= $method_name . '/';
	return $self->makedir($path);
}
sub getJunctionsDir {
	my ( $self, $method_name ) = @_;
	confess() unless $method_name;
	my $path = $self->project_path . "/junctions/";
	$self->makedir($path);
	$path .= $method_name . '/';
	return $self->makedir($path);
}
sub getAnnotSVDir {
	my ( $self, $method_name ) = @_;
	my $path = $self->project_path . "/variations/";
	$self->makedir($path);
	$path .= $method_name . '/annotSV/';
	return $self->makedir($path);
}

sub getGvcfDir {
	my ( $self, $method_name ) = @_;
	my $path = $self->project_path . "/gvcf/";
	$self->makedir($path);
	$path .= $method_name . '/';
	return $self->makedir($path);
}

sub getDejaVuProjectDir {
	my $self = shift;
	my $path = $self->project_path . "/deja_vu/";
	confess("\n\nERROR: $path (dejavu dir) doesn't exist... Die...\n\n")
	  unless ( -d $path );
	return $path;
}

sub getRawCallingDir {
	my ($self) = @_;
	my $path = $self->project_path . "/raw_calling/";
	return $self->makedir($path);
}

sub getIndelsDir {
	my ( $self, $method_name ) = @_;
	my $path = $self->project_path . "/indels/";
	$self->makedir($path);
	$path .= $method_name . '/';
	return $self->makedir($path);
}

sub getLargeIndelsDir {
	my ( $self, $method_name ) = @_;
	my $path = $self->project_path . "/large_indels/";
	$self->makedir($path);
	$path .= $method_name . '/';
	return $self->makedir($path);
}

sub getBedPolyQueryDir {
	my ($self) = @_;
	my $path = $self->rocks_cache_dir . "/bed_polyquery/";
	return $self->makedir($path);
}

sub makedir {
	my ( $self, $dir ) = @_;
	$dir =~ s/\/\//\//g;
	return $dir if -e $dir;
	system("mkdir -p $dir");
	system("chmod a+rwx $dir");
	if ( -d $dir ) { return $dir; }
	return undef;
	confess("\n\nERROR: $dir directory doesn't exist !!\n\n");
}

sub getCallingPipelineDir {
	my ( $self, $method_name ) = @_;
	my $path = $self->project_pipeline_path . "/calling/" . $method_name . '/';
	$self->makedir($path);
	$path =~ s/\/\//\//;
	return $self->makedir($path);
}

sub getNbVarHoRegionsFromType {
	my ($self, $type) = @_;
	my $methodName = 'ho_regions_'.$type.'_value';
	return $self->$methodName();
}

has ho_regions_short_value => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->buffer->config->{ho_regions_genomes}->{short}
		  if ( $self->isGenome() );
		return $self->buffer->config->{ho_regions_exomes}->{short};
	},
);

has ho_regions_medium_value => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->buffer->config->{ho_regions_genomes}->{medium}
		  if ( $self->isGenome() );
		return $self->buffer->config->{ho_regions_exomes}->{medium};
	},
);

has ho_regions_large_value => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->buffer->config->{ho_regions_genomes}->{large}
		  if ( $self->isGenome() );
		return $self->buffer->config->{ho_regions_exomes}->{large};
	},
);

has ho_regions_extra_large_value => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->buffer->config->{ho_regions_genomes}->{extra_large}
		  if ( $self->isGenome() );
		return $self->buffer->config->{ho_regions_exomes}->{extra_large};
	},
);

has getEnsemblVersion => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $v    = $self->buffer->config->{ensembl}->{ $self->getVersion() };
		unless ($v) {
			my $hg = "";
			$hg = "HG19" if $self->getVersion() =~ /HG19/;
			$hg = "HG38" if $self->getVersion() =~ /HG38/;
			die() if $hg eq "";
			$v = $self->buffer->config->{ensembl}->{$hg};
		}
		die unless $v;
		return $v;
	},
);

sub getEnsemblDir {
	my $self = shift;
	my $dir  = $self->getKyotoDir() . "/" . $self->getEnsemblVersion();
	return $self->makedir($dir);

}

sub getEnsemblFreezeDir {
	my $self = shift;
	my $dir  = $self->getEnsemblDir() . '/freeze';
	return $self->makedir($dir);
}

sub statsFromVcfFilesParsed {
	my ( $self, $verbose ) = @_;
	my %arrayRefPatient;
	my $listPatientObject = $self->getPatients();
	foreach my $patObj (@$listPatientObject) {
		my $patName = $patObj->name();
		if ( exists $self->{'patientsNameVcfParsed'}->{$patName} ) {
			$arrayRefPatient{'PATIENTS'}->{$patName}->{'SNPs'} = 0;
			$arrayRefPatient{'PATIENTS'}->{$patName}->{'INS'}  = 0;
			$arrayRefPatient{'PATIENTS'}->{$patName}->{'DEL'}  = 0;
		}
	}
	my $lVar = $self->getObjects('variations');
	my $lIns = $self->getObjects('insertions');
	my $lDel = $self->getObjects('deletions');
	$arrayRefPatient{'TOTAL'}->{'SNPs'} = 0;
	$arrayRefPatient{'TOTAL'}->{'INS'}  = 0;
	$arrayRefPatient{'TOTAL'}->{'DEL'}  = 0;
	foreach my $obj (@$lVar) {
		my @lPatForThisObj = keys( %{ $obj->patients_object() } );
		$arrayRefPatient{'TOTAL'}->{'diffSNPs'} = scalar(@$lVar);
		foreach my $patName (@lPatForThisObj) {
			$arrayRefPatient{'PATIENTS'}->{$patName}->{'SNPs'}++;
		}
		my $lenPat = scalar(@lPatForThisObj);
		unless ( exists( $arrayRefPatient{'DETAILS'}->{'SNPs'}->{$lenPat} ) ) {
			$arrayRefPatient{'DETAILS'}->{'SNPs'}->{$lenPat} = 0;
		}
		$arrayRefPatient{'DETAILS'}->{'SNPs'}->{$lenPat}++;
		$arrayRefPatient{'TOTAL'}->{'SNPs'} =
		  int( $arrayRefPatient{'TOTAL'}->{'SNPs'} ) + $lenPat;
		$arrayRefPatient{'TOTAL'}->{'SNPs'} += $lenPat;
	}
	foreach my $obj (@$lIns) {
		$arrayRefPatient{'TOTAL'}->{'diffINS'} = scalar(@$lIns);
		my @lPatForThisObj = keys( %{ $obj->patients_object() } );
		foreach my $patName (@lPatForThisObj) {
			$arrayRefPatient{'PATIENTS'}->{$patName}->{'INS'}++;
		}
		my $lenPat = scalar(@lPatForThisObj);
		unless ( exists( $arrayRefPatient{'DETAILS'}->{'INS'}->{$lenPat} ) ) {
			$arrayRefPatient{'DETAILS'}->{'INS'}->{$lenPat} = 0;
		}
		$arrayRefPatient{'DETAILS'}->{'INS'}->{$lenPat}++;
		$arrayRefPatient{'TOTAL'}->{'INS'} += $lenPat;
	}
	foreach my $obj (@$lDel) {
		$arrayRefPatient{'TOTAL'}->{'diffDEL'} = scalar(@$lDel);
		my @lPatForThisObj = keys( %{ $obj->patients_object() } );
		foreach my $patName (@lPatForThisObj) {
			$arrayRefPatient{'PATIENTS'}->{$patName}->{'DEL'}++;
		}
		my $lenPat = scalar(@lPatForThisObj);
		unless ( exists( $arrayRefPatient{'DETAILS'}->{'DEL'}->{$lenPat} ) ) {
			$arrayRefPatient{'DETAILS'}->{'DEL'}->{$lenPat} = 0;
		}
		$arrayRefPatient{'DETAILS'}->{'DEL'}->{$lenPat}++;
		$arrayRefPatient{'TOTAL'}->{'DEL'} += $lenPat;
	}
	if ($verbose) { $self->printStatsFromVcfFilesParsed( \%arrayRefPatient ); }
	return \%arrayRefPatient;
}

sub printStatsFromVcfFilesParsed {
	my ( $self, $hashStats ) = @_;
	my $projectName = $self->name();
	print
"\n\n##### STATS SNPs / Insertions / Deletions already found in project $projectName #####\n\n";
	print "### TOTAL\n";
	print "  -> " . $$hashStats{"TOTAL"}->{"SNPs"} . " SNPs\n";
	print "  -> " . $$hashStats{"TOTAL"}->{"diffSNPs"} . " unique SNPs\n";
	print "  -> " . $$hashStats{"TOTAL"}->{"INS"} . " insertions\n";
	print "  -> " . $$hashStats{"TOTAL"}->{"diffINS"} . " unique insertions\n";
	print "  -> " . $$hashStats{"TOTAL"}->{"DEL"} . " deletions\n";
	print "  -> " . $$hashStats{"TOTAL"}->{"diffDEL"} . " unique deletions\n\n";

#	foreach my $nbPat (sort{$a <=> $b}(keys($$hashStats{"DETAILS"}->{"SNPs"}))) {
#		print "  -> ". $$hashStats{"DETAILS"}->{"SNPs"}->{$nbPat} ." SNPs found in ". $nbPat ." patient(s)\n";
#	}
#	print "\n";
#	foreach my $nbPat (sort{$a <=> $b}(keys($$hashStats{"DETAILS"}->{"INS"}))) {
#		print "  -> ". $$hashStats{"DETAILS"}->{"INS"}->{$nbPat} ." insertions found in ". $nbPat ." patient(s)\n";
#	}
#	print "\n";
#	foreach my $nbPat (sort{$a <=> $b}(keys($$hashStats{"DETAILS"}->{"DEL"}))) {
#		print "  -> ". $$hashStats{"DETAILS"}->{"DEL"}->{$nbPat} ." deletions found in ". $nbPat ." patient(s)\n";
#	}
#	foreach my $patName (sort(keys($$hashStats{'PATIENTS'}))) {
#		print "\n### Patient $patName\n";
#		print "  -> ". $$hashStats{'PATIENTS'}->{$patName}->{"SNPs"} ." SNPs\n";
#		print "  -> ". $$hashStats{'PATIENTS'}->{$patName}->{"INS"} ." insertions\n";
#		print "  -> ". $$hashStats{'PATIENTS'}->{$patName}->{"DEL"} ." deletions\n";
#	}
	print "\n\n";
}

sub getCaptureFiles {
	my $self = shift;
	my @lFilesFound;
	foreach my $patObj ( @{ $self->getPatients() } ) {
		push( @lFilesFound, @{ $patObj->getCaptureFiles() } );
	}
	my $hCapture;
	foreach my $file (@lFilesFound) { $hCapture->{$file} = undef; }
	my @res = keys(%$hCapture);
	if ( scalar(@res) > 1 ) {
		confess(
"\n\nERROR: Different capture files declared in this project !! Exit.\n"
		);
	}
	return \@res;
}

sub getCaptureFile {
	my $self = shift;
	return $self->getCaptureFiles();
}

sub getCaptureBedFiles {
	my $self = shift;
	my @lFilesFound;
	foreach my $patObj ( @{ $self->getPatients() } ) {
		push( @lFilesFound, @{ $patObj->getCaptureBedFiles() } );
	}
	my $hCapture;
	foreach my $file (@lFilesFound) { $hCapture->{$file} = undef; }
	my @res = keys(%$hCapture);
	if ( scalar(@res) > 1 ) {
		confess(
"\n\nERROR: Different capture files declared in this project !! Exit.\n"
		);
	}
	return \@res;
}

sub getMaskDatabase {
	my ( $self, $value ) = @_;
	confess("$value") unless exists $maskDB->{$value};
	return $maskDB->{$value};
}

#sub getEnsemblAnnotDir {
#	my $self = shift;
#	mkdir $self->getEnsemblDir()."/annotations" unless -e  $self->getEnsemblDir()."/annotations/";
#	#die() unless -e
#	return $self->getEnsemblDir()."/annotations";
#}

sub existsObjects {
	my ( $self, $type, $id ) = @_;
	return exists $self->{objects}->{$type}->{$id};
}


sub liteAnnotations {
	my ($self) = @_;
	my $hid = "liteAnnotations".$$;
	return $self->{$hid} if exists $self->{$hid};
	$self->{$hid} = GenBoNoSqlAnnotation->new(
			dir  => $self->get_gencode_directory,
			mode => "r"
		);
	return $self->{$hid};
}
#has liteAnnotations => (
#	is      => 'rw',
#	lazy    => 1,
#	default => sub {
#		my $self = shift;
#		return GenBoNoSqlAnnotation->new(
#			dir  => $self->get_gencode_directory,
#			mode => "r"
#		);
#	}
#);
has liteRegulations => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return GenBoNoSqlAnnotation->new(
			dir  => "/data-xfs/dev/pnitschk/svn-genbo/GenBo/script/ngs_exome/index_public_annotations/genes/0_add_local_annotation/toto",
			mode => "r"
		);
	}
);
sub transcriptsCoverageLite {
	my ( $self, $mode, $reopen ) = @_;
	$mode = 'r' unless $mode;
	return $self->{transcriptsCoverageSqlite}
	  if exists $self->{transcriptsCoverageSqlite};
	my $sqliteDir = $self->getProject()->getCoverageDir() . "/lite/";
	$self->{transcriptsCoverageSqlite} =
	  GenBoNoSql->new( dir => $sqliteDir, mode => $mode );
	return $self->{transcriptsCoverageSqlite};
}

has lmdbGenBo => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self      = shift;
		my $sqliteDir = $self->get_gencode_directory;
		die( "you don t have the directory : " . $sqliteDir )
		  unless -e $sqliteDir;
		return GenBoNoSqlLmdb->new(
			name        => "genbo",
			dir         => $sqliteDir,
			mode        => "r",
			is_compress => 1,
			vmtouch=>$self->buffer->vmtouch
		);
	}
);

sub rocksGenBo {
	my ( $self, $mode ) = @_;
	$mode = "r" unless $mode;
	my $name = "genbo-".$mode.$$;
	return $self->{rocks}->{$name} if exists $self->{rocks}->{$name};
	
	$self->{rocks}->{$name}  = GenBoNoSqlRocksAnnotation->new(
			name        => "genbo",
			dir         => $self->get_gencode_directory,
			mode        => "r",
		);
		return $self->{rocks}->{$name};
	
}
sub random {
	my ( $self ) = @_;
	return $self->{rocks}->{random} if exists $self->{rocks}->{random};
	$self->{rocks}->{random} = rand(time);
	return $self->{rocks}->{random};
}
sub rockspLI {
	my ( $self, $mode ) = @_;
	$mode = "r" unless $mode;
	my $name = "pli-".$mode.$self->random();
	return $self->{rocks}->{$name} if exists $self->{rocks}->{$name};
	$self->{rocks}->{$name} = GenBoNoSqlRocksAnnotation->new(
			name        => "pLI",
			dir         => $self->get_public_data_directory("pLI"),
			mode        => "r",
			debug		=> $name,
		);
		
		return $self->{rocks}->{$name};
}


has lmdbMainTranscripts => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self      = shift;
		my $sqliteDir = $self->get_gencode_directory;
		die( "you don t have the directory : " . $sqliteDir )
		  unless -e $sqliteDir;
		confess( $self->get_gencode_directory . "/main_transcripts" )
		  unless -e $self->get_gencode_directory . "/main_transcripts";
		return GenBoNoSqlLmdb->new(
			name        => "main_transcripts",
			dir         => $sqliteDir,
			mode        => "r",
			is_compress => 1,
			vmtouch=>$self->buffer->vmtouch
		);
	}
);  
# TEST 
has liteIntervalTree => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self      = shift;
		my $sqliteDir = $self->get_gencode_directory;
		die( "you don t have the directory : " . $sqliteDir ) unless -e $sqliteDir;
		return GenBoNoSqlIntervalTree->new( dir => $sqliteDir, mode => "r" );    #->new(dir=>$sqliteDir,mode=>"r");
	}
);

has litregionsIntervalTree => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self      = shift;
		my $sqliteDir = "/data-xfs/dev/pnitschk/svn-genbo/GenBo/script/ngs_exome/index_public_annotations/genes/0_add_local_annotation/toto/";
		die( "you don t have the directory : " . $sqliteDir )
		  unless -e $sqliteDir;
		return GenBoNoSqlIntervalTree->new( dir => $sqliteDir, mode => "r" )
		  ;    #->new(dir=>$sqliteDir,mode=>"r");
	}
);

#has dejavuSVDirectory => (
#	is      => 'ro',
#	lazy    => 1,
#	default => sub {
#		my $self = shift;
#		my $dir  = "/data-xfs/Manue/Test_SV/DejaVu/newVersion/";
#		return $dir;
#	}
#);
has dejavuCNV => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $sqliteDir =   $self->DejaVuCNV_path;
#		warn $sqliteDir;
		die("you don t have the directory : ".$sqliteDir) unless -e $sqliteDir;
		return  GenBoNoSqlDejaVuCNV->new( dir => $sqliteDir, mode => "r" );
		#return GenBoNoSqlIntervalTree->new(dir=>$sqliteDir,mode=>"r");#->new(dir=>$sqliteDir,mode=>"r");
	}
);
sub localdejavuCNV {
	my ($self,$mode) = @_;
	return $self->{lite}->{"ldv".$mode} if exists  $self->{lite}->{"ldv".$mode};
	my $sqliteDir =  $self->getCNVDir()."/dejavu/";
	$self->{lite}->{"ldv".$mode} = GenBoNoSqlDejaVuCNV->new( dir => $sqliteDir, mode => $mode );
}

has dejavuSV => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $sqliteDir =   $self->DejaVuSVeq_path;
		warn $sqliteDir;
		die("you don t have the directory : ".$sqliteDir) unless -e $sqliteDir;
		return  GenBoNoSqlDejaVuSV->new( dir => $sqliteDir, mode => "r" );
		#return GenBoNoSqlIntervalTree->new(dir=>$sqliteDir,mode=>"r");#->new(dir=>$sqliteDir,mode=>"r");
	}
);
has dejavuSVAgregate => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $sqliteDir =   $self->DejaVuSVeq_path;
		die("you don t have the directory : ".$sqliteDir) unless -e $sqliteDir;
		warn $sqliteDir;
		return  GenBoNoSqlDejaVuSV_agregate->new( dir => $sqliteDir, mode => "r" );
		#return GenBoNoSqlIntervalTree->new(dir=>$sqliteDir,mode=>"r");#->new(dir=>$sqliteDir,mode=>"r");
	}
);


has dejavuJunctions => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $release = $self->annotation_genome_version();
		$release = 'HG19' if ($release =~ /HG19/);
		my $sqliteDir = $self->DejaVuJunction_path();
		die("you don t have the directory : ".$sqliteDir) unless -e $sqliteDir;
		return  GenBoNoSqlDejaVuJunctions->new( dir => $sqliteDir, mode => "r" );
	}
);

has dejavuJunctionsResume => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $release = $self->annotation_genome_version();
		$release = 'HG19' if ($release =~ /HG19/);
		$release = 'HG38' if ($release =~ /HG38/);
		my $sqliteDir = $self->DejaVuJunction_path();
		die("you don t have the directory : ".$sqliteDir) unless -e $sqliteDir;
		
		confess();
		
		return  GenBoNoSqlDejaVuJunctionsResume->new( dir => $sqliteDir, mode => "r" );
	}
);

has dejavuJunctionsCanoniques => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $release = $self->annotation_genome_version();
		$release = 'HG19' if ($release =~ /HG19/);
		$release = 'HG38' if ($release =~ /HG38/);
		my $sqliteDir = $self->get_dejavu_junctions_path('canoniques');
		die("you don t have the directory : ".$sqliteDir) unless -e $sqliteDir;
		
		confess();
		
		return  GenBoNoSqlDejaVuJunctionsCanoniques->new( dir => $sqliteDir, mode => "r" );
	}
);

sub dejavuJunctionsPhenotype {
	my ($self, $phenotype_name) = @_;
	$phenotype_name =~ s/ /_/g;
	$phenotype_name = lc($phenotype_name);
	confess("\n\nERROR: DV Junctions Phenotype mandatory\n\n") unless ($phenotype_name);
	return $self->{'dv_junctions_'.$phenotype_name} if (exists $self->{'dv_junctions_'.$phenotype_name});
	my $release = $self->annotation_genome_version();
	$release = 'HG19' if ($release =~ /HG19/);
	$release = 'HG38' if ($release =~ /HG38/);
	my $sqliteDir = $self->get_dejavu_junctions_path($phenotype_name);
	die("you don t have the directory : ".$sqliteDir) unless -e $sqliteDir;
	
	confess();
		
	$self->{'dv_junctions_'.$phenotype_name} = GenBoNoSqlDejaVuJunctionsPhenotype->new( phenotype_name => $phenotype_name, dir => $sqliteDir, mode => "r" );
	return $self->{'dv_junctions_'.$phenotype_name};
}


has dejavuSVIntervalTree => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $sqliteDir =   $self->DejaVuCNV_path;
		die("you don t have the directory : ".$sqliteDir) unless -e $sqliteDir;
		return GenBoNoSqlIntervalTree->new(dir=>$sqliteDir,mode=>"r");#->new(dir=>$sqliteDir,mode=>"r");
	}
);


#WARNING

sub noSqlCoverage {
	my ( $self, $mode ) = @_;

	#confess();
	return $self->{noSqlCoverage} if exists $self->{noSqlCoverage};
	$mode = "w" unless $mode;

	#	my $output   =$self->getCacheDir() . "/coverage_lite_test";
	my $output = $self->getCacheDir() . "/coverage_lite";
	$self->{noSqlCoverage} = GenBoNoSql->new( dir => $output, mode => "$mode" );
	
	return $self->{noSqlCoverage};

}

has transcriptsCacheDir => (
	is => 'ro',

	#isa		=> 'ArrayRef[Str]',
	lazy    => 1,
	default => sub {
		my $self   = shift;
		my $output = $self->rocks_cache_dir() . "/transcripts";
		return $self->makedir($output);

	},
);

has transcriptsCoverageDir => (
	is => 'ro',

	#isa		=> 'ArrayRef[Str]',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->getProject->transcriptsCacheDir() . "/coverage/";
		return $self->makedir($dir);
	},
);

has transcriptsDudeDir => (
	is => 'ro',

	#isa		=> 'ArrayRef[Str]',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->getProject->transcriptsCacheDir() . "/dude/";
		return $self->makedir($dir);

	},
);

has noSqlCnvsDir => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->getVariationsDir("dude") . "/cnv_lite";
		return $self->makedir($dir);

	}
);

sub noSqlCnvs {
	my ( $self, $mode ) = @_;
	return $self->{noSqlCnvs} if exists $self->{noSqlCnvs};
	$mode = "w" unless $mode;

	#my $output   =$self->getCacheDir() . "/coverage_lite";
	my $output = $self->noSqlCnvsDir();

	$self->{noSqlCnvs} = GenBoNoSql->new( dir => $output, mode => "$mode" );
	return $self->{noSqlCnvs};

}

#has noSqlCoverage =>(
#	is		=> 'ro',
#	lazy	=> 1,
#	default => sub {
#		my $self = shift;
#
#		my $output   =$self->getCacheDir() . "/coverage_lite_test";
#		warn $output;
#		my $no = GenBoNoSql->new(dir=>$output,mode=>"w");
#	}
#);

has noSqlCoverage_readOnly => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;

		my $output = $self->getCacheDir() . "/coverage_lite";
		return GenBoNoSql->new( dir => $output, mode => "r" );

	}
);

sub existsnoSqlQuality {
	my ($self) = @_;
	return -e $self->quality_dir() . "/" . $self->name . ".lite";
}

sub noSqlQuality {
	my ( $self, $param ) = @_;
	return $self->{noSqlQuality} if exists $self->{noSqlQuality};
	$param = "r" unless $param;

	my $output = $self->quality_dir;   #$self->getCacheDir() . "/check_quality";
	if ( $param eq "r" ) {
		if ( -z "$output/" . $self->name . ".lite"
			or !( -e "$output/" . $self->name . ".lite" ) )
		{
			my $no = GenBoNoSql->new( dir => $output, mode => "c" );
			$no->put( $self->name, "_void_", undef );
			$no->close();
		}

	}
	$self->{noSqlQuality} = GenBoNoSql->new( dir => $output, mode => $param );
	return $self->{noSqlQuality};
}

sub getEnsgIDs {
	my ( $self, $name ) = @_;
	return [$name] if $name =~ /^ENSG/;
	confess();
	my $syno = $self->liteAnnotations->get( "synonyms", $name );
	return $syno;

}

sub if_exists_liteObject {
	my ( $self, $id, $type ) = @_;
	confess();
#	foreach
#	  my $text ( @{ $self->liteAnnotations->get_text( "annotations", $id ) } )
#	{
#		if ($type) {
#			my @lFields = split( ' ', $text );
#			foreach my $field (@lFields) {
#				if ( $field eq $type ) { return 1; }
#			}
#		}
#		else { return 1; }
#	}
#	return;
}

sub liteObject {
	my ( $self, $id, $type ) = @_;
	confess();
	return $self->liteAnnotations->get( "annotations", $id );
}



sub getGenBoId {
	my ( $self, $geneId ) = @_;
	#delete $self->{liteAnnotations};
	my $syno = $self->liteAnnotations->get( "synonyms", $geneId );
	return $syno;

	#my ($toto,$titi) = split("_",$geneId);
	return $geneId unless $syno;
	return $syno->[0] if scalar(@$syno) == 1;
	confess("more than one syno $geneId");
}

sub getGenBoIds {
	my ( $self, $geneId ) = @_;
	#delete  $self->{liteAnnotations};
	
	my $syno = $self->liteAnnotations->get_key_values( "synonyms", $geneId );
	return $syno;
}

sub setPrimers {
	my ($self) = @_;
	my %hash;
	foreach my $c ( @{ $self->getCaptures } ) {
		map { $hash{ $_->id }++ } @{ $c->getPrimers };
	}
	return \%hash;
}

sub setKyotoGene {
	my ( $self, $geneId ) = @_;

	#		confess();
	#my $hgene = thaw $self->kyotoGenes->{$geneId};
	my $hgene = $self->liteObject( $geneId, "genes" );

	confess("can't find $geneId") unless ($hgene);

	#	unless ($hgene){
	#		my ($tid,$tchr) = split("_",$geneId);
	#
	#		 $hgene = thaw $self->kyotoGenes->{$tid};
	#
	#
	#	}

  #confess("this gene ($geneId) doesn't exist in kyoto file !") unless ($hgene);
	$hgene->{id}      = $geneId;
	$hgene->{kyotoId} = $geneId;
	my ( $n, $c ) = split( "_", $geneId );
	$hgene->{name} = $n;
	$hgene->{chromosomes_object}->{ $hgene->{chromosome} } = undef;
	my $obj = $self->flushObject( 'genes', $hgene );

	return 1;
}

sub getKyotoTranscripts {
	my ( $self, $tr_id ) = @_;
	return $self->liteObject($tr_id);

	#return  thaw $self->kyotoTranscriptsdb->get($tr_id);
}
my %htest;

sub setKyotoTranscript {
	my ( $self, $transcriptId ) = @_;
	confess();

	#my $hgene = thaw $self->kyotoGenes->{$geneId};
	$htest{$transcriptId}++;
	my $hTranscript = $self->liteObject( $transcriptId, "transcripts" );
	confess("can't find transcript $transcriptId") unless ($hTranscript);

	confess(
		"WARN: this transcript ($transcriptId) doesn't exist in kyoto file !")
	  unless ($hTranscript);
	my ( $id, $c ) = split( "_", $transcriptId );
	my $id2 .= $id . "_" . $hTranscript->{chromosome};
	$hTranscript->{id} = $id2;
	$hTranscript->{chromosomes_object} =
	  { $hTranscript->{chromosome} => undef };
	$hTranscript->{kyotoId} = $transcriptId;
	$hTranscript->{name}    = $id;
	my $obj = $self->flushObject( 'transcripts', $hTranscript );

	return 1;
}

sub setKyotoProtein {
	my ( $self, $id ) = @_;
	confess();
	my $hprotein = $self->liteObject( $id, "proteins" );
	confess("can't find $id") unless ($hprotein);

#my $hprotein = thaw $self->kyotoProteins->{$id};
#confess ("WARN: this protein ($id) doesn't exist in kyoto file !") unless ($hprotein);
	$hprotein->{id} = $id;

	unless ( exists $hprotein->{chromosome} ) {
		my ( $i, $chr ) = split( "_", $hprotein->{genbo_id} );
		$hprotein->{chromosome} = $chr;
	}

	unless ( exists $hprotein->{chromosome} ) {
		my ( $i, $chr ) = split( "_", $hprotein->{genbo_id} );
		$hprotein->{chromosome} = $chr;
	}

	$hprotein->{chromosomes_object} = { $hprotein->{chromosome} => undef };
	my $tname = $hprotein->{transcript};
	$tname .= "_" . $hprotein->{chromosome} unless ( $tname =~ /_/ );
	$hprotein->{transcripts_object} = { $tname => undef };
	$hprotein->{genes_object} =
	  { $hprotein->{gene_kyoto_id} => undef };    #."_".$hprotein->{chromosome};
	my $obj = $self->flushObject( 'proteins', $hprotein );
	return 1;
}

### mask coding

my $maskCoding = {
	intergenic                             => 1,
	intronic                               => 2,
	"non_coding_transcript_intron_variant" => 2,
	splice_site                            => 4,
	essential_splicing                     => 8,
	pseudogene                             => 16,
	"non_coding_transcript_exon_variant"   => 32,

	ncrna                 => 32,
	utr                   => 64,
	stop                  => 128,
	phase                 => 256,
	silent                => 512,
	nonsynonymous         => 1024,
	frameshift            => 2048,
	nonframeshift         => 4096,
	maturemirna           => 8192,
	upstream              => 16384,
	downstream            => 32768,
	predicted_splice_site => 65536
};


has maskImpact => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		return {
			intergenic            => 1,
			intronic              => 2,
			splice_site           => 4,
			essential_splicing    => 8,
			pseudogene            => 16,
			ncrna                 => 32,
			utr                   => 64,
			stop                  => 128,
			phase                 => 256,
			silent                => 512,
			nonsynonymous         => 1024,
			frameshift            => 2048,
			nonframeshift         => 4096,
			maturemirna           => 8192,
			upstream              => 16384,
			downstream            => 32768,
			predicted_splice_site => 65536,
	#		duplication           => 131072,
	#		deletion           => 262144,
		};
	}
);

has maskImpactTextForLegend => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		return {
			intergenic            => "intergenic",
			intronic              => "intronic",
			splice_site           => "splice region",
			essential_splicing    => "splice donor/acceptor",
			pseudogene            => "pseudogene",
			ncrna                 => "non coding transcript",
			utr                   => "utr",
			stop                  => "stop gained",
			phase                 => "stop/start lost",
			silent                => "synonymous",
			nonsynonymous         => "missense",
			frameshift            => "frameshift",
			nonframeshift         => "inframe insertion/deletion",
			maturemirna           => "mature miRNA",
			upstream              => "upstream",
			downstream            => "downstream",
			predicted_splice_site => "predicted splice region",
		#	duplication           => "predicted splice region",
		#	deletion           	  => "predicted splice region",
		};
	}
);

has maskImpactText => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $maskCodingText;

#$maskCodingText->{exonic} = ['utr', 'stop', 'phase', 'silent', 'nonsynonymous', 'frameshift', 'nonframeshift'];
#$maskCodingText->{coding} = ['stop', 'phase', 'silent', 'nonsynonymous', 'frameshift', 'nonframeshift'];
		$maskCodingText->{high} =
		  [ 'maturemirna', 'stop', 'phase', 'frameshift',
			'essential_splicing' ];
		$maskCodingText->{medium} = [
			'maturemirna',        'stop',
			'phase',              'frameshift',
			'essential_splicing', 'nonsynonymous',
			'nonframeshift',      'splice_site',
			'predicted_splice_site'
		];

		#$maskCodingText->{moderate} = $maskCodingText->{medium};
		$maskCodingText->{low} = [
			'maturemirna',         'stop',
			'phase',               'frameshift',
			'essential_splicing',  'ncrna',
			'nonsynonymous',       'nonframeshift',
			'splice_site',         'predicted_splice_site',
			'pseudogene',          'intronic',
			'silent',              'utr',
			'upstream_downstream', 'upstream',
			'downstream'
		];
		return $maskCodingText;
	}
);

$maskCoding->{exonic} =
  $maskCoding->{utr} | $maskCoding->{stop} | $maskCoding->{phase} |
  $maskCoding->{silent} | $maskCoding->{nonsynonymous} |
  $maskCoding->{frameshift} | $maskCoding->{nonframeshift};
$maskCoding->{coding} =
  $maskCoding->{stop} | $maskCoding->{phase} | $maskCoding->{silent} |
  $maskCoding->{nonsynonymous} | $maskCoding->{frameshift} |
  $maskCoding->{nonframeshift};

$maskCoding->{high} =
  $maskCoding->{maturemirna} | $maskCoding->{stop} | $maskCoding->{phase} |
  $maskCoding->{frameshift} | $maskCoding->{essential_splicing};
$maskCoding->{moderate} =
  $maskCoding->{high} | $maskCoding->{nonsynonymous} |
  $maskCoding->{nonframeshift} | $maskCoding->{splice_site} |
  $maskCoding->{predicted_splice_site};
$maskCoding->{medium} = $maskCoding->{moderate};
$maskCoding->{low} =
  $maskCoding->{moderate} | $maskCoding->{ncrna} | $maskCoding->{pseudogene} |
  $maskCoding->{intronic} | $maskCoding->{silent} | $maskCoding->{utr} |
  $maskCoding->{upstream} | $maskCoding->{downstream};

#$maskCoding->{all} = $maskCoding->{low} | $maskCoding->{intronic} ;

my $maskCodingDiag;
$maskCodingDiag->{exonic} =
  $maskCoding->{utr} | $maskCoding->{stop} | $maskCoding->{phase} |
  $maskCoding->{silent} | $maskCoding->{nonsynonymous} |
  $maskCoding->{frameshift} | $maskCoding->{nonframeshift};
$maskCodingDiag->{coding} =
  $maskCoding->{stop} | $maskCoding->{phase} | $maskCoding->{silent} |
  $maskCoding->{nonsynonymous} | $maskCoding->{frameshift} |
  $maskCoding->{nonframeshift};

$maskCodingDiag->{high} =
  $maskCoding->{maturemirna} | $maskCoding->{stop} | $maskCoding->{phase} |
  $maskCoding->{frameshift} | $maskCoding->{essential_splicing};
$maskCodingDiag->{medium} =
  $maskCodingDiag->{high} | $maskCoding->{ncrna} |
  $maskCoding->{nonsynonymous} | $maskCoding->{nonframeshift} |
  $maskCoding->{splice_site};
$maskCodingDiag->{moderate} = $maskCodingDiag->{medium} | $maskCoding->{silent};
$maskCodingDiag->{low} =
  $maskCodingDiag->{moderate} | $maskCoding->{pseudogene} |
  $maskCoding->{intronic} | $maskCoding->{utr} | $maskCoding->{upstream} |
  $maskCoding->{downstream};

sub getMaskCoding {
	my ( $self, $value ) = @_;
	confess("$value") unless exists $maskCoding->{$value};
	return $maskCoding->{$value};
}

#####################
# Sift polyphen prediction
#####################""""

our $maskPrediction = {
	'polyphen_probably damaging' => 64,
	'polyphen_possibly damaging' => 32,
	'sift_deleterious'           => 16,
	'polyphen_benign'            => 8,
	'sift_tolerated'             => 4,
	'polyphen_unknown'           => 2,
	'sift_unknown'               => 1
};

$maskPrediction->{completely_benign} =
  $maskPrediction->{polyphen_benign} | $maskPrediction->{sift_tolerated};
$maskPrediction->{completely_damaging} =
  $maskPrediction->{'polyphen_probably damaging'} |
  $maskPrediction->{sift_deleterious};

our $PREDICTION_VALUES = {
	polyphen => {
		'probably damaging' => 3,
		'possibly damaging' => 2,
		'benign'            => 1,
		'unknown'           => 0,
		3                   => 'probably damaging',
		2                   => 'possibly damaging',
		1                   => 'benign',
		0                   => 'unknown',

	},

	sift => {
		'unknown'     => 0,
		'tolerated'   => 1,
		'deleterious' => 2,

		#'deleterious'   => 2,
		0 => 'unknown',
		1 => 'tolerated',
		2 => 'deleterious',
		3 => 'tolerated',
		4 => 'deleterious',
	},
};

sub returnPredictionValue {
	my ( $self, $method, $value ) = @_;
	confess("$value") unless exists $PREDICTION_VALUES->{$method}->{$value};
	return $PREDICTION_VALUES->{$method}->{$value};
}

sub getMaskPrediction {
	my ( $self, $value ) = @_;
	confess("$value") unless exists $maskPrediction->{$value};
	return $maskPrediction->{$value};
}



sub getPredictions {
	my ( $self, $chr, $prot_id, $pos, $aa ) = @_;
	confess() unless $chr;
	confess() unless $prot_id;
	my $no = $self->litePredictionMatrix();

	#	die():
	my $id     = $prot_id . "_" . $pos;
	my $matrix = $no->get( $chr->name(), $prot_id, 1 );
	my $types  = {
		"sift" => "sift",

		#"polyphen_humdiv" => "polyphen",
		"polyphen_humvar" => "polyphen"
	};
	my $res;
	foreach my $type ( keys %$types ) {
		my $stype = $types->{$type};
		$res->{$type}->{score} = 0;

		( $res->{$type}->{pred_all}, $res->{$type}->{score} ) =
		  decode_prediction_matrix::prediction_from_matrix( $matrix->{$stype},
			$pos, $aa );
		$res->{$type}->{pred_all} = 0 unless $res->{$type}->{pred_all};
		$res->{$type}->{pred}     = $res->{$type}->{pred_all};
		if ( $type eq "polyphen_humvar" ) {

			# j'inverse la prediction
			# big warning
			$res->{$type}->{pred} = abs( $res->{$type}->{pred} - 3 );
		}
		elsif ( $type eq "sift" ) {
			$res->{$type}->{pred} = 0 if $res->{$type}->{score} > 1;
			$res->{$type}->{pred}++ if $res->{$type}->{score} <= 1;
		}
		else {
			die();
		}

		my $text_value = $self->returnPredictionValue( $types->{$type},
			$res->{$type}->{pred} );
		$text_value = $types->{$type} . "_" . $text_value;

		#	warn $text_value;
		$res->{mask} |= $self->getMaskPrediction("$text_value");

#(,,$res->{$type}->{pred_all}) = decode_prediction_matrix::decode($val,$types->{$type});

	}

	#my $res2 = $self->getPredictionsOld($chr,$prot_id,$pos,$aa);

	return $res;

}

sub getDataBase {
	my $self = shift;
	return "polydb";
	confess() unless exists $self->buffer->config->{polyprod}->{name};

	return $self->buffer->config->{polyprod}->{name};
}

sub isCacheDone {
	my ( $self, $annot_version ) = @_;
	$annot_version = $self->annotation_version() unless ($annot_version);
	my $cache_dir =
		$self->buffer()->getDataDirectory("cache") . "/"
	  . $self->genome_version() . '.'
	  . $annot_version . "/"
	  . $self->name() . '/';
	my $global_infos = $cache_dir . '/vector/global_infos.freeze';
	return 1 if ( -d $cache_dir and -e $global_infos );
	return;
}

has cache_date => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self           = shift;
		my $genome_version = $self->genome_version();
		my $annot_version  = $self->annotation_version();
		my $dir =
			$self->buffer()->getDataDirectory("cache") . "/"
		  . $genome_version . '.'
		  . $annot_version . "/"
		  . $self->name()
		  . '/vector/';
		if ( -d $dir ) {
			my @lStats     = stat $dir;
			my $now_string = strftime "%e %b %Y", localtime( $lStats[9] );
			return $now_string;
		}
		return '?';
	}
);

sub getCacheDir {
	my $self           = shift;
	return $self->{cache_dir} if (exists $self->{cache_dir} and $self->{cache_dir});
	my $genome_version = $self->genome_version();
	my $annot_version = $self->annotation_version();
	$genome_version .= '.'.$annot_version if ($annot_version and $annot_version ne '.');
	$self->{cache_dir} = $self->buffer()->config_path("root","cache").$genome_version;
	unless (-d $self->{cache_dir}){
		system( "mkdir  " . $self->{cache_dir} );
		system( "chmod a+rwx " . $self->{cache_dir} );
	}
	$self->{cache_dir} .= "/".$self->name();
	$self->{cache_dir} =~ s/ \/\//\//g;
	return $self->{cache_dir} if ( -d $self->{cache_dir} );
	system( "mkdir  " . $self->{cache_dir} );
	system( "chmod a+rwx " . $self->{cache_dir} );
	return $self->{cache_dir};
}

has getCacheCNV => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir_root = $self->getCacheDir();
		my $dir = $dir_root."/CNV/";
		unless (-d $dir){
			system( "mkdir  " . $dir );
			system( "chmod a+rwx " . $dir );
		}
		return $dir;
	},
);

has getCacheSV => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir_root = $self->getCacheDir();
		my $dir = $dir_root."/SV/";
		unless (-d $dir){
			system( "mkdir  " . $dir );
			system( "chmod a+rwx " . $dir );
		}
		return $dir;
	},
);

sub isCacheVectorDone {
	my $self = shift;
	my $cache_dir =
		$self->buffer->getDataDirectory("cache") . "/"
	  . $self->getVersion() . "/"
	  . $self->name()
	  . "/vector/";
	return unless ( -d $cache_dir );
	return unless ( -d $cache_dir . '/dejavu/' );
	return unless ( -e $cache_dir . '/global_infos.freeze' );
	my $hChrDone;
	opendir( DIR, $cache_dir . '/dejavu/' ) or die $!;
	while ( my $file = readdir(DIR) ) {
		next if ( $file eq '.' );
		next if ( $file eq '..' );
		my ( $chr, $ext ) = split( '\.', $file );
		$hChrDone->{$chr} = undef;
	}
	closedir(DIR);
	my @lDirToCheck = (
		'all_global',            'all_variants_ids',
		'dejavu',                'dejavu_ho',
		'exome-intronic_global', 'exome-intronic_variants_ids',
		'genes_ids',             'genes_infos',
		'variants_infos'
	);
	foreach my $cat_name (@lDirToCheck) {
		return unless ( -d $cache_dir . '/' . $cat_name . '/' );
		my $hChrCheck = $hChrDone;
		opendir( DIR, $cache_dir . '/' . $cat_name . '/' ) or die $!;
		while ( my $file = readdir(DIR) ) {
			next if ( $file eq '.' );
			next if ( $file eq '..' );
			my ( $chr, $ext ) = split( '\.', $file );
			delete $hChrCheck->{$chr};
		}
		closedir(DIR);
		if ( scalar keys %$hChrCheck > 0 ) {
			warn "ERROR: "
			  . $self->name()
			  . " missing $cat_name -> "
			  . join( ',', sort keys %$hChrCheck ) . " !\n";
			return;
		}
	}
	return 1;
}

sub getCacheBitVectorDir {
	my $self      = shift;
	my $cache_dir = $self->getCacheDir();
	$self->makedir($cache_dir);
	$cache_dir .= '/vector/';
	return $self->makedir($cache_dir);
}

sub getCacheVariationsKyotoFile {
	my ( $self, $chr_name ) = @_;
	my $dir = $self->getCacheDir();
	confess("give a chromosome name for variations cache") unless $chr_name;
	my $file_name = $dir . "/variations.$chr_name." . $self->name() . ".kct";
}

sub getCachePatientKyotoFile {
	my ( $self, $name ) = @_;
	my $dir = $self->getCacheDir();
	confess("give a chromosome name for variations cache") unless $name;
	my $file_name = $dir . "/$name.kct";
}

sub getCacheGenesKyotoFile {
	my ( $self, $chr_name ) = @_;
	my $dir = $self->getCacheDir();
	confess("give a chromosome name for variations cache") unless $chr_name;
	my $file_name = $dir . "/genes.$chr_name." . $self->name() . ".kct";
	return $file_name;
}

sub getCacheFilterGenesKyotoFile {
	my ( $self, $chr_name ) = @_;
	my $dir = $self->getCacheDir();
	confess("give a chromosome name for variations cache") unless $chr_name;
	my $file_name = $dir . "/genes.filter.$chr_name." . $self->name() . ".kct";
	return $file_name;
}

sub getCacheVariationsFile {
	my ( $self, $chr_name ) = @_;
	my $dir = $self->getCacheDir();
	confess("give a chromosome name for variations cache") unless $chr_name;
	my $file_name =
	  $dir . "/variations.$chr_name." . $self->name() . ".store.gz";
}

sub getCacheGenesFile {
	my $self = shift;
	my $dir  = $self->getCacheDir();

	return $dir . "/genes." . $self->name . ".store.gz";
}

sub get_list_chromosomes {
	my ( $self, $chromosome_name, $separator ) = @_;
	$separator       = ","   unless $separator;
	$chromosome_name = "all" unless $chromosome_name;
	if ( $chromosome_name eq 'all' ) {
		return $self->getChromosomes();
		my @z = sort { $a->name cmp $b->name } @{ $self->getPatients() };
	}
	else {
		my @names = split( $separator, $chromosome_name );
		my $chromosomes;
		foreach my $name (@names) {
			unless ( $self->getChromosome($name) ) {
				warn "$name =>"
				  . join( " ", map { $_->name } @{ $self->getChromosomes() } );
				die();
			}
			push( @$chromosomes, $self->getChromosome($name) );
		}
		return $chromosomes;
	}
}

sub get_list__controls_patients {
	my ( $self, $patients_name, $separator ) = @_;
	$separator = "," unless $separator;
	my $patients;
	if ( $patients_name eq 'all' ) {
		my @z = sort { $a->name cmp $b->name } @{ $self->getPatientsAndControl() };
		$patients = \@z;
	}
	else {
		my @names = split( $separator, $patients_name );

		foreach my $name (@names) {
			unless ( $self->getPatientOrControl($name) ) {
				warn "$name =>"
				  . join( " ", map { $_->name } @{ $self->getPatientsAndControl() } );
				die();
			}
		}
		map { push( @$patients, $self->getPatientOrControl($_) ) }
		  split( ",", $patients_name );
	}
	return $patients;
}


sub get_list_patients {
	my ( $self, $patients_name, $separator ) = @_;
	$separator = "," unless $separator;
	my $patients;
	if ( $patients_name eq 'all' ) {
		my @z = sort { $a->name cmp $b->name } @{ $self->getPatients() };
		$patients = \@z;
	}
	else {
		my @names = split( $separator, $patients_name );

		foreach my $name (@names) {
			unless ( $self->getPatient($name) ) {
				warn "$name =>"
				  . join( " ", map { $_->name } @{ $self->getPatients() } );
				die();
			}
		}
		map { push( @$patients, $self->getPatient($_) ) }
		  split( ",", $patients_name );
	}
	return $patients;
}



sub get_only_list_patients {
	my ( $self, $patients_name, $separator ) = @_;
	$separator = "," unless $separator;
	my $patients;
	$patients_name = "all" unless $patients_name;
	if ( $patients_name eq 'all' ) {
		$patients = $self->getPatients();
	}
	else {
		my %names;
		map { $names{$_}++ } split( $separator, $patients_name );

		foreach my $patient ( @{ $self->getPatients() } ) {
			if ( exists $names{ $patient->name } ) {
				push( @$patients, $patient );
				delete $names{ $patient->name };
			}
			else {
				delete $self->{objects}->{patients}->{ $patient->id };
				delete $self->patients_object->{ $patient->id };
			}

		}

	}

	return $patients;
}

sub getCoverageDir {
	my $self = shift;
	my $dir  = $self->getAlignmentRootDir . "/coverage/";
	return $self->makedir($dir);
}

sub getTargetCountDir {
	my $self = shift;
	my $dir  = $self->getCoverageDir."/target_count/";
	return $self->makedir($dir);
}

sub getCoverageCallable {
	my $self = shift;
	my $dir  = $self->getAlignmentRootDir . "/callable/";
	return $self->makedir($dir);
}

sub purge_memory {
	my ( $self, $limit ) = @_;
	my $delete_types = [ "variations", "insertions", "deletions" ];
	foreach my $id ( keys %{ $self->{objects}->{genes} } ) {
		my $gene = $self->{objects}->{genes}->{$id};
		if ( $gene->end() + 501 < $limit ) {
			foreach my $tr ( @{ $gene->getTranscripts() } ) {
				foreach my $p ( @{ $tr->getProteins } ) {
					$self->{dejavu}->{ $p->id }++;
					delete $self->{objects}->{proteins}->{ $p->id };
					$p = undef;
					undef($p);

				}
				foreach my $e ( @{ $tr->getExons } ) {

					delete $self->{objects}->{exons}->{ $e->id };
					$e = undef;
					undef($e);
				}

				return 1 if exists $self->{dejavu}->{ $tr->id };
				delete $self->{objects}->{transcripts}->{ $tr->id };
				undef($tr);
				$tr = undef;
			}

			delete $self->{objects}->{genes}->{ $gene->id };
			$gene = undef;
			undef($gene);
		}
	}
	#delete $self->{objects};
	warn "pruge";
	
	#delete $self->{objects}->{chromosomes};
	
	#delete $self->{objects}->{variants};
	#warn Dumper keys %{$self->{objects}};
	delete $self->{objects};
}

sub purge_memory_reference {
	my ( $self, $ref ) = @_;
	my $delete_types = [ "variations", "insertions", "deletions", ];

	my $nb;

	my $o = $self->{objects}->{references}->{$ref};
	my $objs;

	#push(@$objs,$o->getChromosome());
	push( @$objs, @{ $self->getPatients() } );

	foreach my $type (@$delete_types) {
		my $vtype = $type . "_object";

		foreach my $id ( keys %{ $o->{$vtype} } ) {

			my $obj = $self->{objects}->{$type}->{$id};
			delete $obj->{project};
			delete $obj->{buffer};
			undef($obj);
			$obj = undef;
		}
		delete $o->{$vtype};
	}

	foreach my $type (@$delete_types) {
		foreach my $id ( keys %{ $self->{objects}->{$type} } ) {
			my $obj = $self->{objects}->{$type}->{$id};
			$nb++;
			delete $obj->{project};
			delete $obj->{buffer};
			undef($obj);
			$obj = undef;
		}
		delete $self->{objects}->{$type};
	}

	#warn "delete $nb";

}

##### DEJA  VU KYOTO

#has deja_vu_positions_file => (
#	is      => 'ro',
#	lazy    => 1,
#	default => sub {
#		my $self = shift;
#		return $self->deja_vu_public_dir . "/positions.tab.gz";
#	}
#);



has deja_vu_public_projects_parquet  => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self   = shift;
		return $self->buffer->config_path("root","dejavu").'/projects_parquet/';
	},
);

has lite_deja_vu_projects => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->deja_vu_lite_dir_projects;
		
		my $no = GenBoNoSqlDejaVu->new( dir => $dir, mode => "r" );
		return $no;
	}
);

has lite_deja_vu_rsname => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $no =
		  GenBoNoSql->new( dir => $self->lite_dir . "/rsname/", mode => "r" );
		return $no;
	}
);


sub get_deja_vu_from_position {
	my ( $self, $chr, $start, $end ) = @_;
	my $no = $chr->rocks_dejavu();
	my $h = $no->dejavu_interval( ($start-1), ($end+1) );
	return $h;
}

has countInThisRunPatients => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		my $similar = $self->similarProjects();

		my $query = $self->buffer->getQuery();
		my $nb    = $query->countPatients( [ keys %$similar ] );
		return $nb;

	},
);

has in_this_run_patients => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $h;
	my $res;
	my $total = 0;
	my $inthisrunp;
	foreach my $run (@{$self->getRuns}){
		my $patients = $run->getAllPatientsInfos();
			foreach my $p (@$patients){
				$inthisrunp->{$p->{project_id} }->{$p->{id}} ++;
				$total ++;
			}
		}
		$inthisrunp->{nb_patients} = $total;
		return $inthisrunp;
	},
);
sub getDejaVuInfosForDiagforVariant{
	my ($self, $v) = @_;
	my $chr = $v->getChromosome()->name();
	return  $v->getChromosome()->getDejaVuInfosForDiagforVariant($v);
	
}
sub getDejaVuInfosForDiag {
	my ($self, $vid) = @_;
	my $v = $self->_newVariant($vid);
	return $self->getDejaVuInfosForDiagforVariant($v);
}



sub getDejaVuThisProject {
	my ( $self, $id ) = @_;
	confess();
	my ( $chr, @t ) = split( "_", $id );
	unless ( exists $self->{dejavu_this}->{$chr} ) {

		#		my $no = $self->lite_deja_vu_projects();
		#		 $self->{dejavu_this}->{$chr} =$no->get($self->name,$chr);
		$self->{dejavu_this}->{$chr} =
		  $self->lite_deja_vu_projects()->get( $self->name, $chr );
	}
	my $v = $self->{dejavu_this}->{$chr}->{$id};

	return 1 unless exists $self->{dejavu_this}->{$chr}->{$id};

	my ( $p1, $p2 ) = split( " ", $self->{dejavu_this}->{$chr}->{$id} );

	return 1 unless $p1;
	my $n   = () = $p1 =~ /,/g;
	my @hho = ();

	@hho = split( ",", $p2 ) if $p2;
	return ( $n + 1, scalar(@hho) );
}

sub getDejaVuInfos {
	my ( $self, $id ) = @_;
	my ( $chr, @t ) = split( "_", $id );
	confess();
	my $hres;
	my $no = $self->lite_deja_vu2();
	my $string_infos = $no->get( $chr, $id );
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
		my $hpat = $no->get( "projects", $project_name );
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

has deja_vu_lite_dir => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		confess();
		my $dir = $self->deja_vu_public_dir();
		#TODO: ajouter makedir si createion automatique
		#return $self->makedir($dir);
		return $dir;
	},
);



has deja_vu_rocks_dir => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;

		#my $dir =$self->buffer->deja_vu_public_dir($self->genome_version_generic)."/rocks/";
		#return $self->makedir($dir);
		my $dir =$self->buffer->deja_vu_public_dir($self->genome_version_generic);
		return $dir;
		

	},
);


sub deja_vu_lite_public_dir {
	my ($self,$version,$type)= @_;
	$type = "variations" unless $type;
	$version = $self->genome_version_generic unless $version;
	return $self->buffer->deja_vu_public_dir($version,$type);
}

sub deja_vu_rocks_project_dir {
	my ($self,$version)= @_;
	$version = $self->genome_version_generic unless $version;
	my $type = "projects";
	my $root =  $self->buffer->deja_vu_project_dir($version,$type);
	return $root;
}

sub deja_vu_rocks_public_dir {
	my ($self,$version,$type)= @_;
	$type = "variations" unless $type;
	$version = $self->genome_version_generic unless $version;
	my $root =  $self->buffer->deja_vu_public_dir($version,$type);
	return $root;
}

has lite_deja_vu => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		warn "\n\n\n\n";
		confess();
		my $dir  = $self->deja_vu_lite_dir;

		#	wanr $dir;
		my $no = GenBoNoSqlDejaVu->new( dir => $dir, mode => "r" );
		return $no;
	}
);
#has lite_deja_vu2 => (
#	is      => 'ro',
#	lazy    => 1,
#	default => sub {
#		my $self = shift;
#		my $dir  = $self->deja_vu_lite_dir;
##		warn $dir;
#		my $no = GenBoNoSql->new( dir => $dir, mode => "r" );
#		return $no;
#	}
#
#);

sub lite_deja_vu_local_validation {
	my ( $self, $rsName ) = @_;
	return $self->{lmdb}->{dv_local_validation} if exists $self->{lmdb}->{dv_local_validation};
	
	my $dir  = $self->deja_vu_lite_dir;
		 $self->{lmdb}->{dv_local_validation} = GenBoNoSqlLmdb->new(
			dir         => $dir,
			mode        => "r",
			name        => "local",
			is_compress => 1,
			vmtouch=>$self->buffer->vmtouch
		);
		return $self->{lmdb}->{dv_local_validation};
}

has lite_deja_vu2 => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->deja_vu_lite_dir;
		
		my $no = GenBoNoSqlDejaVu->new( dir => $dir, mode => "r" );
		return $no;
	}
);

sub convertRsNameToIntervalPositions {
	my ( $self, $rsName ) = @_;
	my $no_rsname = $self->gnomad_rsname_dir();
	$rsName =~ s/rs//;
	foreach my $chr_id ( 1 .. 22, 'X', 'Y' ) {
		my $no = GenBoNoSqlLmdb->new(
			dir         => $no_rsname,
			mode        => "r",
			name        => $chr_id,
			is_compress => 1,
			is_integer  => 1,
			vmtouch=>$self->buffer->vmtouch
		);
		my $res = $no->get($rsName);
		$no->close();
		if ($res) {
			my $region = $chr_id . ':' . int( $res - 1 ) . '-' . int( $res + 1 );
			return $region;
		}
	}
	return;
}

sub getDejaVuHo {
	my ( $self, $chr, $ids, $max ) = @_;
	return [] unless $ids;
	my $no   = $self->lite_deja_vu();
	my $rids = [];

	foreach my $id (@$ids) {
		my $nb = $no->get_lmdb( $id . "_ho" );
		push( @$rids, $id ) if $nb && $nb > $max;
	}
	return $rids;
}

sub getCapturedChromosomes {
	my ($self) = @_;    #methode qui prend en argument un projet
	my @chr_list;
	my %chr_list_uniq;
	my $captures = $self->getCaptures();
	foreach my $capture (@$captures) {
		print "\ncapture file name : \n";
		print my $capture_name = $capture->file_name()
		  . "\n";       # $capture is GenBoCapture=HASH(0x23f15d08)
		my $captureFiles  = $capture->files();
		my $captureGzFile = $captureFiles->{gz};
		die() unless -e $captureGzFile;
		my $cmd =
		  $self->buffer->config->{software}->{tabix} . " -l $captureGzFile";
		my @res = `$cmd`;
		chomp(@res);

		foreach my $name (@res) {
			$chr_list_uniq{$name}++;
		}
	}
	confess() if scalar( keys %chr_list_uniq ) == 0;
	return keys %chr_list_uniq;

}

sub getDiagCacheDir {
	my $self = shift;
	my $dir =
		$self->buffer()->getDataDirectory("diag_cache") . "/"
	  . $self->getVersion() . "/"
	  . $self->name()
	  . '/cache_polydiag/';
	return $self->makedir($dir);
}

sub json_ensembl_annotations {
	my $self = shift;
	my @lHashes;
	my $hLen;
	$hLen->{5}  = '80px';
	$hLen->{10} = '90px';
	$hLen->{15} = '100px';
	$hLen->{20} = '110px';
	$hLen->{25} = '120px';
	$hLen->{30} = '130px';
	my $hAnnot = $self->ensembl_annotations();

	foreach my $key ( keys %$hAnnot ) {
		my $value = $hAnnot->{$key};
		my ( $ensembl, $synonym ) = split( ';', $value );
		my $this_hash;
		$this_hash->{id}          = $key;
		$this_hash->{ensembl}     = $ensembl;
		$this_hash->{synonym}     = $synonym;
		$this_hash->{filter_name} = "filter_" . $key;
		if ( length( $this_hash->{ensembl} ) <= 5 ) {
			$this_hash->{length} = $hLen->{5};
		}
		elsif ( length( $this_hash->{ensembl} ) <= 10 ) {
			$this_hash->{length} = $hLen->{10};
		}
		elsif ( length( $this_hash->{ensembl} ) <= 15 ) {
			$this_hash->{length} = $hLen->{15};
		}
		elsif ( length( $this_hash->{ensembl} ) <= 20 ) {
			$this_hash->{length} = $hLen->{20};
		}
		elsif ( length( $this_hash->{ensembl} ) <= 25 ) {
			$this_hash->{length} = $hLen->{25};
		}
		elsif ( length( $this_hash->{ensembl} ) <= 30 ) {
			$this_hash->{length} = $hLen->{30};
		}
		push( @lHashes, $this_hash );
	}
	return \@lHashes;
	#
}

has sqlite_var_infos => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->getCacheBitVectorDir() . '/variants_infos';
		return undef unless ( -d $dir );
		return GenBoNoSql->new( dir => $dir, mode => 'r' );
	}
);

sub getVariantVcfInfos {
	my ( $self, $varId ) = @_;
	my $nosql = $self->sqlite_var_infos();
	return unless ($nosql);
	my @lTmp = split( '_', $varId );
	my $hash = $nosql->get( $lTmp[0], $varId );
	return unless ($hash);
	return $hash;
}

sub noSqlPolydiag {
	my ( $self, $param ) = @_;
	if ($param && $param eq "close") {
		return unless exists $self->{nosqlpolydiag};
		$self->{nosqlpolydiag}->close();
		delete $self->{nosqlpolydiag};
		return;
	}
	return $self->{nosqlpolydiag} if exists $self->{nosqlpolydiag};
	$param = "r" unless $param;
	my $output = $self->getCacheDir() . "/polydiag_lite";
	$self->{nosqlpolydiag} = GenBoNoSql->new( dir => $output, mode => $param );
	return $self->{nosqlpolydiag};
}

sub setVariants {
	my ( $self, $type ) = @_;

	my $method;
	if    ( $type eq 'variations' )      { $method = 'getVariations'; }
	elsif ( $type eq 'insertions' )      { $method = 'getInsertions'; }
	elsif ( $type eq 'deletions' )       { $method = 'getDeletions'; }
	elsif ( $type eq 'large_deletions' ) { $method = 'getLargeDeletions'; }
	elsif ( $type eq 'large_duplications' ) {$method = 'getLargeDuplications';}
	elsif ( $type eq 'large_insertions' ) {$method = 'getLargeInsertions';}
	elsif ( $type eq 'inversions' ) {$method = 'getInversions';}
	my $h;

	foreach my $chr ( @{ $self->getChromosomes() } ) {
		foreach my $var ( @{ $chr->$method() } ) {
			$h->{ $var->id() } = undef;
		}
	}

	return $h;
}

#new database for transcripts vector

has dir_lmdb_score_impact => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->getCacheBitVectorDir() . '/lmdb_score_impact_intspan';
		return $self->makedir($dir);
	}
);

has dir_hash_variants => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->getCacheBitVectorDir() . '/hash_variants';
		return $self->makedir($dir);
	}
);

has dir_lite_cache_tree_objects => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir  = $self->getCacheBitVectorDir() . '/tree_objects';
		return $self->makedir($dir);

	}
);
has dir_lmdb_bitvector_transcripts => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir =
		  $self->getCacheBitVectorDir() . '/lmdb_cache_bitvector_transcripts/';
		return $self->makedir($dir);
	}
);

has isNoSqlDepth => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		foreach my $p ( @{ $self->getPatients() } ) {
			return undef unless ( $p->isNoSqlDepth );
		}
		return 1;
	}
);

has isUpdate => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		foreach my $chr ( @{ $self->getChromosomes() } ) {
			next
			  unless (
				$chr->countThisVariants( $chr->getVariantsVector() ) > 0 );
			my $var = $chr->getVarObject(1);
			return 1
			  if $var->sum_mask_coding() >=
			  $self->maskImpact->{'predicted_splice_site'};
			return;
		}
		return;
	}
);

# do you have run dude

has isDude => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;

		#	warn $self->tabix_primers_file;
		#	return 1 if $self->tabix_primers_file;
		return -e $self->tabix_primers_file;
	},
);

# fichier généré par dude
has tabix_primers_file => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;

		my $dir_out = $self->noSqlCnvsDir;
		return $dir_out . "/primers.bed.gz";
	},
);

sub tabix_primers {
	my ($self) = @_;
	confess( $self->tabix_primers_file ) unless $self->isDude;
	return Bio::DB::HTS::Tabix->new( filename => $self->tabix_primers_file );
}

has dir_varids_to_vector_ids => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->getCacheBitVectorDir() . '/varids_to_vector_ids/';
	}
);

sub need_update {
	my ($self) = @_;
	unless ( -d $self->dir_varids_to_vector_ids() ) { return 1; }
	return 0;
}

sub disconnect {
	my $self = shift;
	$self->buffer->dbh_deconnect();
	delete $self->{liteAnnotations};
	delete $self->{cosmic_db};
	delete $self->{transcriptsCoverageSqlite};
	delete $self->{lmdbGenBo};
	delete $self->{lmdbMainTranscripts};
	$self->close_rocks;
	
	#foreach my $c (values %{$self->{rocks}}){
	#	$c->close() if $c;
	#}
}
sub close_rocks {
	my $self = shift;
	foreach my $c (values %{$self->{rocks}}){
		#warn $c;
		#$c->close() if $c;
		delete $self->{rocks}->{$c}
	}
	delete $self->{rocks};
}
sub preload_patients {
	my $self = shift;
	$self->getFamilies();
	$self->pedigree_details;
	foreach my $p (@{$self->getPatients()}) {
		$p->callingSVMethods();
		$p->callingMethods();
		$p->getBamFile(undef,1);

}
$self->getRuns();
$self->getPhenotypes();
$self->similarProjects();
$self->similarProjectsId();
$self->in_this_run_patients();
$self->exomeProjects();
$self->countSimilarPatients();
$self->countExomePatients();


}
sub preload{
	my $self = shift;
	$self->preload_patients;
	
}

has get_polybtf_project_path => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir = $self->buffer->get_polybtf_path() . '/' . $self->name() . '/';
		return $dir;
	},
);

has get_polybtf_resume => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->buffer->get_polybtf_project_resume( $self->name() );
	},
);

has dir_controls_dude => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $dir =
"/data-xfs/dev/pnitschk/svn-genbo/polypipeline/scripts/scripts_pipeline/dude//genome//novaseq/";
		return $dir;
	},
);

has get_path_rna_seq_polyrna_root  => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $path = $self->buffer()->config_path("root","project_data")."/".$self->getProjectType()."/".$self->name()."/".$self->version()."/polyRNA/";
		return $path if -d $path;
		my @lPotentialRelease = ('HG19', 'HG19_CNG', 'HG19_MT', 'HG38', 'HG38-ERCC', 'MM38', 'MM39');
		foreach my $rel2 (@lPotentialRelease) {
			my $alt_path = $self->buffer()->config_path("root","project_data")."/".$self->getProjectType()."/".$self->name()."/".$rel2."/polyRNA/";
			return $alt_path if -d $alt_path;
			
		}
		return $path;
	},
);

has get_path_rna_seq_junctions_root  => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $path = $self->buffer()->config_path("root","project_data")."/".$self->getProjectType()."/".$self->name()."/".$self->version()."/analysis/";
		return $path;
	},
);

has get_path_rna_seq_junctions_analyse_all_res  => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $path = $self->get_path_rna_seq_junctions_root()."/AllRes/";
		return $path;
	},
);
has RnaseqSEA_SE  => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $path = $self->get_path_rna_seq_junctions_analyse_all_res;
		my $file = $path.'/allResSE.txt';
		my $filegz = $file.".gz";
		if(-e $filegz) {
			return $filegz;
		}
		return if not -e $file;
		$self->tabix_gzip_rnaseqsea($filegz,$file);
		return $filegz;
	},
);
has RnaseqSEA_RI  => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $path = $self->get_path_rna_seq_junctions_analyse_all_res;
		my $file = $path.'/allResRI.txt';
		my $filegz = $file.".gz";
		if(-e $filegz) {
			return $filegz;
		}
		return if not -e $file;
		$self->tabix_gzip_rnaseqsea($filegz,$file);
		return $filegz;
	},
);

sub tabix_gzip_rnaseqsea {
	my ($self,$filegz,$file) = @_;
	my $chr = 4;
	my $start = 5;
	my $end = 6;
	if ($file =~ /SE/){
		#$chr = 2;
		#$start = 6;
		#$end = 7;
		
	}
	if (-e $file){
		 	my $bgzip = $self->buffer->software("bgzip");
		 	my $tabix = $self->buffer->software("tabix");
		 	my $firstline = `head -n 1 $file`;
		 	$firstline ="#".$firstline;
		 	chomp($firstline);
		 	warn  "(echo \"$firstline\"  && tail -n +2 $file | sort -k$chr,$chr -k$start,$start"."n) | $bgzip -s /dev/stdin  > $filegz && $tabix -f -S 1 -s $chr -b $start  -e $end $file.gz";
		 	system("(echo \"$firstline\"  && tail -n +2 $file | sort -k$chr,$chr -k$start,$start"."n) | $bgzip -s /dev/stdin  > $filegz && $tabix -f  -s $chr -b $start  -e $end $file.gz");
		 }
		 else{
		 	confess("no $file");
		 }
		confess("no $filegz") unless (-e $filegz);
}
has get_path_rna_seq_junctions_analyse_description_root  => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $path = $self->buffer()->config_path("root","project_data")."/".$self->getProjectType()."/".$self->name()."/".$self->version()."/RNAseqSEA/";
		return $path;
	},
);

has get_hash_patients_description_rna_seq_junction_analyse => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $project_name = $self->name;
		my (@lFilesDescription, $hType_patients);
		my $path_rna_description = $self->get_path_rna_seq_junctions_analyse_description_root();
		return $hType_patients if (not -d $path_rna_description);
		
		opendir my ($dir), $path_rna_description;
		my @found_files = readdir $dir;
		closedir $dir;
		foreach my $file (@found_files) {
			next if (not $file =~ /SamplesTypes_/);
			push(@lFilesDescription, $path_rna_description.'/'.$file);
		}
		foreach my $file (@lFilesDescription) {
			next if (not -e $file);
			open (FILE, $file);
			my (@lPat, @lCtrl);
			while (<FILE>) {
				chomp($_);
				my $line = $_;
				next if not $line =~ /$project_name/;
				my ($pat_name, $type, $proj_name) = split(' ', $line);
				$hType_patients->{$pat_name}->{lc($type)} = undef;
				push(@lPat, $pat_name) if (lc($type) eq 'pat');
				push(@lCtrl, $pat_name) if (lc($type) eq 'ctrl');
			}
			close(FILE);
			foreach my $pat_name (@lPat) {
				foreach my $ctrl (@lCtrl) {
					$hType_patients->{$pat_name}->{used_ctrl}->{$ctrl} = undef;
				}
			}
		}
		return $hType_patients;
	},
);


sub get_gtf_genes_annotations_igv {
	my ($self) = @_;
	if ($self->getVersion() =~ /HG19/) {
		my $igv_dir = $self->buffer->config_path("root","public_data").'/igv/';
		if (defined $self->gencode_version() && $self->gencode_version() ne '-1') {
			my $file = $igv_dir.'/gencode.'.$self->gencode_version().'.gtf.gz';
			return $file if (-e $file);
		}
		my $file = $igv_dir.'/gencode.gtf.gz';
		return $file;
	}
	return $self->buffer->config_path("root","public_data")."/".$self->getVersion()."/igv/gencode.gtf.gz";
}

sub getQueryJunction {
	my ($self, $fileName, $method) = @_;
	my %args;
	$args{project} = $self;
	$args{file}    = $fileName;
	if (lc($method) eq 'ri') { $args{isRI} = 1; }
	elsif (lc($method) eq 'se') { $args{isSE} = 1; }
	elsif (lc($method) eq 'dragen') { $args{isDRAGEN} = 1; }
	elsif (lc($method) eq 'star') { $args{isSTAR} = 1; }
	elsif (lc($method) eq 'regtools') { $args{isREGTOOLS} = 1; }
	else { confess(); }
	my $queryJunction = QueryJunctionFile->new( \%args );
	return $queryJunction;
}
sub getSJDir {
	my ( $self, $method_name ) = @_;
	my $path = $self->project_path . "/junctions/";
	$self->makedir($path);
	$path .= $method_name . '/';
	return $self->makedir($path);
}





sub setJunctions {
	my ($self) = @_;
	my @lObj;
	foreach my $chr (@{$self->getChromosomes()}) {
		push(@lObj, @{$chr->getJunctions()});
	}
	return \@lObj;
}
		
sub writeCaptureBedFile {
	my ($self,$span,$file) = @_;
	$span = 0 unless $span;
	open ("BED",">$file");
	foreach my $chr (@{$self->getChromosomes}){
		my $intspan = $chr->getIntSpanCapture($span);
		next if $intspan->is_empty;
		my @line = $self->buffer->intspanToBed($chr,$intspan);
		print BED join("\n",@line)."\n";
	}	
	close(BED);
}

sub return_calling_methods_short_name {
	my ($self,$name) = @_;
	unless (exists $self->{short}){
		foreach my $m (@{$self->calling_methods}){
			$self->{short}->{$m->{name}} = $m->{short_name};
		}
	}
	
	return unless exists $self->{short}->{$name};
	return $self->{short}->{$name};
}

sub is_sv_calling_methods{
	my ($self,$name) = @_;
	unless (exists $self->{isSV}){
		$self->{isSV} = {};
		foreach my $mn (@{$self->callingSVMethods}){
			$self->{isSV}->{$mn} = 1;
			my $sn = $self->return_calling_methods_short_name($mn);
			$self->{isSV}->{$sn} = 1;
			
		}
	}
	return exists $self->{isSV}->{$name};
}
sub calling_methods {
	my ($self) = @_;
	return $self->{calling_methods1} if exists $self->{calling_methods1};
	my $methods = [];
			foreach my $method_name (@{$self->callingSVMethods()}) {
				my $method;
				$method->{short_name} = substr $method_name,0,3;
				$method->{name} = $method_name;
				$method->{sv} = 1;
				push(@$methods,$method);
			}
			foreach my $method_name (@{$self->getCallingMethods()}) {
					my $method;
				$method->{short_name} = substr $method_name,0,3;
				$method->{name} = $method_name;
				push(@$methods,$method);
				
			}
	$self->{calling_methods1} = $methods;
	return $self->{calling_methods1};
}


1;

