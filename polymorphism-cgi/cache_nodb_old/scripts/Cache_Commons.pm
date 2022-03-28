package Cache_Commons;
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Moose;
use Data::Dumper;
use Parallel::ForkManager;
use Bio::DB::Sam;
use Storable qw(store retrieve freeze thaw);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Cwd 'abs_path';
use Digest::MD5 qw(md5_hex);
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Carp;
use Compress::Snappy;



has categories => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $categories = {
			global => {
				frequency => {
					freq_none    => 1,
					freq_1       => 1,
					freq_05      => 1,
					freq_01      => 1,
					freq_001     => 1,
					freq_0001    => 1,
					freq_ho_none => 1,
					freq_ho_1    => 1,
					freq_ho_05   => 1,
					freq_ho_01   => 1,
					freq_ho_001  => 1,
					freq_ho_0001 => 1,
					freq_he_none => 1,
					freq_he_1    => 1,
					freq_he_05   => 1,
					freq_he_01   => 1,
					freq_he_001  => 1,
					freq_he_0001 => 1,
					pheno_snp    => 1,
					ngs_score2   => 1,
					ngs_score1   => 1,
					ngs_score0   => 1,
					cadd_not		   => 1,
					cadd_0			   => 1,
					cadd_5			   => 1,
					cadd_10			   => 1,
					cadd_15			   => 1,
					cadd_20			   => 1,
					cadd_25			   => 1,
					cadd_30			   => 1,
					cadd_35			   => 1,
					cadd_40			   => 1,
		
				},
				variation_type => {
					substitution   => 1,
					insertion      => 1,
					deletion       => 1,
					large_deletion => 1,
					large_duplication => 1,
				},
			},
			genes => {
				annotations => {
					utr                => 1,
					splicing           => 1,
					pseudogene         => 1,
					coding             => 1,
					maturemirna        => 1,
					essential_splicing => 1,
					phase              => 1,
					silent             => 1,
					intergenic         => 1,
					stop               => 1,
					ncrna              => 1,
					frameshift         => 1,
					upstream	       => 1,
					downstream         => 1,
					intronic           => 1,
					"non-frameshift"   => 1,
					prediction_0       => 1,
					prediction_1       => 1,
					prediction_2       => 1,
					prediction_3       => 1,
				},
			},
			patients => {
				all => 1,
				he  => 1,
				ho  => 1,
			  }
		};
		return $categories;
	}
);

has output_files => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $outputs = {
			'cache_store_ids' => {
				'/lmdb_cache/' => {
					'CHR_NAME.dv.freeze' => 1,
				},
				'/lmdb_cache/CHR_NAME/' => {
					'categories_annotations' => 1,
					'categories_annotations_index' => 1,
				},
				'/lmdb_cache/patients/' => {
					'CHR_NAME' => 1,
				},
				'/lmdb_cache/variations/' => {
					'CHR_NAME' => 1,
					'CHR_NAME_index' => 1,
				},
			},
			'cache_store_annotations' => {
				'/lmdb_cache/CHR_NAME/' => {
					'genes' => 1,
					'genes_index' => 1,
				},
			},
			'cache_strict_denovo' => {
				'/strict-denovo/' => {
					'CHR_NAME.lite' => 1,
				},
			},
			'cache_loh' => {
				'/somatic_loh/' => {
					'CHR_NAME.lite' => 1,
				},
			},
			'cache_global_infos' => {
				'/' => {
					'global_infos.freeze' => 1,
				},
			},
		};
		return $outputs;
	}
);

has specific_check_steps => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $outputs = {
			'cache_store_ids' => {
				'estim_nb_var' => '/lmdb_cache/CHR_NAME.dv.freeze',
			},
			'cache_store_annotations' => {
				'estim_nb_var' => '/lmdb_cache/CHR_NAME.dv.freeze',
			},
		};
		return $outputs;
	}
);

sub get_regions {
	my ( $chr, $fork ) = @_;
	if ( $fork == 1 ) {
		my $regions;
		my $hregions;
		$hregions->{chromosome} = $chr->name;
		$hregions->{start}      = 1;
		$hregions->{end}        = $chr->length;
		push( @$regions, $hregions );
		return \@$regions;
	}
	my $tabix = $chr->buffer->software("tabix");
	my %hv;
	foreach my $patient ( @{ $chr->project->getPatients() } ) {
		my $calling_files = $patient->callingFiles();
		my @files;
		foreach my $m ( keys %{$calling_files} ) {
			foreach my $v ( values %{ $calling_files->{$m} } ) {
				push( @files, $v );
			}
		}
		foreach my $file (@files) {
			my $file = $patient->getVariationsFiles->[0];
			next unless -e $file;
			open( TOTO, "$tabix $file " . $chr->ucsc_name() . " | cut -f 2 |" );
			while ( my $pos = <TOTO> ) {
				chomp($pos);
				$hv{$pos}++;
			}
			close TOTO;
		}
	}
	my @snps = sort { $a <=> $b } keys %hv;
	if ( scalar(@snps) == 0 ) {
		my $hregions;
		$hregions->{chromosome} = $chr->name;
		$hregions->{start}      = 1;
		$hregions->{end}        = $chr->length;
		return [$hregions];
	}
	my $nb;
	if ( $fork == 1 ) {
		$nb = scalar(@snps);
	}
	else {
		$nb = int( scalar(@snps) / ( $fork - 1 ) );
	}
	$nb = 7_000 if $nb > 7_000;
	my $regions;
	my $iter = natatime $nb, @snps;
	my $old_end;
	while ( my @tmp = $iter->() ) {
		my $hregions;
		$hregions->{chromosome} = $chr->name;
		if ($old_end) {
			$hregions->{start} = $old_end + 1;
		}
		else {
			$hregions->{start} = 1;
		}
		$hregions->{end} = $tmp[-1] + 100;
		$old_end = $hregions->{end};
		push( @$regions, $hregions );
	}
	$regions->[0]->{start} = 1;
	$regions->[-1]->{end}  = $chr->length + 1000;
	return $regions;
}

1;