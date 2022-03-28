#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use File::Basename;
use  File::Temp;

my $project_name;
my $chr_name;
my $fork;

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


GetOptions(
	"project=s" =>\$project_name,
	"fork=s" =>\$fork,
	"chr=s" =>\$chr_name,
);


my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $chr = $project->getChromosome($chr_name);

store_ids($chr,$fork);

sub store_ids {
	my ( $chr, $fork ) = @_;
	my $patients = $chr->project->getPatients();
	my $regions;
	my @lVarRegion;
	my $new_pos       = 1;
	my $old_pos       = 1;
	my $nb_variants   = 0;
	my $nb_var_region = 0;
	$regions = get_regions( $chr, $fork );
	my $pm = new Parallel::ForkManager($fork);
	my $all;
	$pm->run_on_finish(
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;
			unless (defined($hRes) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
			#append in the txt file all the variations found in this region and return by getIds
			my ($region) = grep { $_->{start} eq $hRes->{region}->{start} and $_->{end} eq $hRes->{region}->{end} } @$regions;
			die() unless $region;
			my $freeze_file = $chr->lmdb_cache_dir."/". $region->{chromosome}.".".$region->{start}."-".$region->{end}.'.freeze';
			$region->{freeze} = $freeze_file;
			unlink $freeze_file if -e $freeze_file;
			store($hRes->{'variants'}, $freeze_file);
			die() unless -e $freeze_file;
			$nb_variants += scalar(@{$hRes->{'variants'}});
			foreach my $hv (@{$hRes->{'variants'}}) {
				$chr->{cache_hash_get_var_ids}->{$hv->{id}}++;
			}
		}
	);
	my $cpt = 1;
	foreach my $region (@$regions) {
		warn "$cpt/".scalar(@$regions) if ( $chr->project->cache_verbose() );
		$cpt++;
		my $pid = $pm->start and next;
		unless ( exists $region->{chromosome} ) {
			$region->{chromosome} = $chr->id();
		}
		my $time_start = time;
		my $all        = get_ids( $chr->getProject->name, $region );
		my $time_end   = time;
		my $hres;
		$hres->{variants} = $all;
		$hres->{region}   = $region;
		$hres->{ttime}    = abs( $time_start - $time_end );
		$pm->finish( 0, $hres );
	}
	$pm->wait_all_children();
	warn "end process" if ( $chr->project->cache_verbose() );

	#end fork now I have a file with all variation and  json  it's time to sort this file and store it in lmdb database
	#construct intspan 	for patient I will store lmdb_id in intspan
	my $project = $chr->project;
	my @patient_names = sort { $a cmp $b } map { $_->name } @{ $project->getPatients };
	my $index_patients = 0;
	my $hpatients;
	my @categorie_patient = ( "all", "he", "ho" );
	foreach my $pname (@patient_names) {
		$hpatients->{$pname}->{index} = $index_patients++;
		$hpatients->{$pname}->{name}  = $pname;
		foreach my $c (@categorie_patient) {
			$hpatients->{$pname}->{intspan}->{$c} = Set::IntSpan::Fast::XS->new();
		}
	}

	#initialisation global categorie
	my @global = ( "substitution", "insertion", "deletion", "large_deletion" );
	my $intspan_global_type;
	foreach my $g ( keys %{ $categories->{global}->{variation_type} } ) {
		$intspan_global_type->{$g} = Set::IntSpan::Fast::XS->new();
	}

	my $no2 = $chr->get_lmdb_variations("c");    #open lmdb database
	#ok sort and read the file
	my $hh;
	my $uniq;
	my $size_variants = 0;
	foreach my $region (@$regions) {
		my $freeze = $region->{freeze};
		warn Dumper $region  unless -e $freeze;
		die($freeze) unless -e $freeze;
		my $hall = retrieve $freeze;
		foreach my $hv ( @{$hall} ) {
			my $var_id = $hv->{id};
			next if exists $uniq->{$var_id};
			$uniq->{$var_id}++;
			my $variation = thaw( decompress( $hv->{obj} ) );
			my $index_lmdb = $no2->put( $var_id, $variation );
			$size_variants++;
			$hh->{$var_id} = $variation->{heho_string};
			foreach my $pn ( keys %{ $variation->{patients_details} } ) {
				my $type = $variation->{patients_details}->{$pn}->{type};
				$hpatients->{$pn}->{intspan}->{all}->add($index_lmdb);
				$hpatients->{$pn}->{intspan}->{$type}->add($index_lmdb);
			}
			unless ( exists $intspan_global_type->{ $variation->type() } ) {
				warn "\n\nERROR: doesn't exists \$intspan_global_type->{".$variation->type() . "}\n\n";
				warn Dumper keys %$intspan_global_type;
				die;
			}
			$intspan_global_type->{ $variation->type }->add($index_lmdb);
			unlink $freeze;
		}
	}
	$no2->close();

	my $chr_name = $chr->name();
	#store htable only fort dejavu by project purpose
	store( $hh, $project->lmdb_cache_dir . "/$chr_name.dv.freeze" ) if $hh;    
	my $no3 = $chr->get_lmdb_categories("c");
	foreach my $k ( keys %{$intspan_global_type} ) {
		my $h;
		$h->{name}    = $k;
		$h->{intspan} = $intspan_global_type->{$k};
		my $bitv = Bit::Vector->new_Enum( $size_variants, join( ',', $intspan_global_type->{$k}->as_array ) );
		$h->{bitvector} = $bitv;
		$no3->put( $k, $h );
	}
	$no3->close();
	my $no4 = $chr->get_lmdb_patients("c");
	foreach my $pname (@patient_names) {
		my $h;
		$h->{name} = $pname;
		foreach my $c (@categorie_patient) {
			my $intspan = $hpatients->{$pname}->{intspan}->{$c};
			$h->{intspan}->{$c} = $intspan;
			my $bitv = Bit::Vector->new_Enum( $size_variants, join( ',', $intspan->as_array ) );
			$h->{bitvector}->{$c} = $bitv;
		}
		$no4->put( $pname, $h );
	}
	$no4->close;
	return $nb_variants;
}