#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use Parallel::ForkManager;
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
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use Cache_Commons;



my $fork = 1;
my ($project_name, $chr_name, $no_verbose, $skip_pseudo_autosomal);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'chr=s'        => \$chr_name,
	'no_verbose=s' => \$no_verbose,
	'skip_pseudo_autosomal=s' => \$skip_pseudo_autosomal,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }

my $nbErrors = 0;
my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name );
if ($no_verbose) {
	$project->cache_verbose(0);
}
warn "\n### CACHE: store ids step\n" if ( $project->cache_verbose() );
my $chr = $project->getChromosome($chr_name);
my $patients = $chr->project->getPatients();
my $regions;
my @lVarRegion;
my $new_pos       = 1;
my $old_pos       = 1;
my $nb_variants   = 0;
my $nb_var_region = 0;
$regions = Cache_Commons::get_regions( $chr, $fork );
my $pm = new Parallel::ForkManager($fork);
my $all;
$pm->run_on_finish(
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;
		unless (defined($hRes) or $exit_code > 0) {
			$nbErrors++;
			print qq|No message received from child process $exit_code $pid!\n|;
			return;
		}
		#append in the txt file all the variations found in this region and return by getIds
		my ($region) = grep { $_->{start} eq $hRes->{region}->{start} and $_->{end} eq $hRes->{region}->{end} } @$regions;
		confess() unless $region;
		my $freeze_file = $chr->lmdb_cache_dir."/". $region->{chromosome}.".".$region->{start}."-".$region->{end}.'.freeze';
		$region->{freeze} = $freeze_file;
		unlink $freeze_file if -e $freeze_file;
		store($hRes->{'variants'}, $freeze_file);
		confess() unless -e $freeze_file;
		$nb_variants += scalar(@{$hRes->{'variants'}});
		foreach my $hv (@{$hRes->{'variants'}}) {
			$chr->{cache_hash_get_var_ids}->{$hv->{id}}++;
		}
	}
);
my $cpt = 1;
foreach my $region (@$regions) {
	warn "$cpt/".scalar(@$regions) if ( $project->cache_verbose() );
	$cpt++;
	my $pid = $pm->start and next;
	unless ( exists $region->{chromosome} ) {
		$region->{chromosome} = $chr->id();
	}
	my $time_start = time;
	my $all        = get_ids( $project_name, $region );
	my $time_end   = time;
	my $hres;
	$hres->{variants} = $all;
	$hres->{region}   = $region;
	$hres->{ttime}    = abs( $time_start - $time_end );
	$pm->finish( 0, $hres );
}
$pm->wait_all_children();
warn "end process" if ( $chr->project->cache_verbose() );

if ($nbErrors > 0) {
	confess("\n\nERRORS: $nbErrors errors found... confess...\n\n");
}
if (scalar keys %{$chr->{cache_hash_get_var_ids}} == 0) {
	my $cmd = "touch ".$project->getCacheBitVectorDir()."/lmdb_cache/".$chr->id().".empty";
	`touch $cmd`;
}

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
my $intspan_global_type;
my $categories = Cache_Commons::categories();
foreach my $g ( keys %{ $categories->{global}->{variation_type} } ) {
	$intspan_global_type->{$g} = Set::IntSpan::Fast::XS->new();
}

warn 'store 1/3: lmdb variations' if ( $project->cache_verbose() );
my $no2 = $chr->get_lmdb_variations("c");    #open lmdb database
#ok sort and read the file
my $hh;
my $uniq;
my $size_variants = 0;
foreach my $region (@$regions) {
	my $freeze = $region->{freeze};
	warn Dumper $region  unless -e $freeze;
	confess($freeze) unless -e $freeze;
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
			confess;
		}
		$intspan_global_type->{ $variation->type }->add($index_lmdb);
		unlink $freeze;
	}
}
$no2->close();

my $nb_from_vcf = scalar(keys %{$chr->{cache_hash_get_var_ids}});
my $nb_from_final = scalar(keys %{$hh});
if ($nb_from_vcf == $nb_from_final) {
	warn "check nb var: nb var VCF parsed = $nb_from_vcf and nb var final = $nb_from_final -> OK\n" if ( $project->cache_verbose() );
}
else {
	warn "\n\nERROR:\n";
	warn "   check nb var: nb var VCF parsed = $nb_from_vcf and nb var final = $nb_from_final -> ERROR\n";
	warn "   DIE...\n\n";
	die();
}

warn 'store 2/3: lmdb chr_name freeze' if ( $project->cache_verbose() );
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

warn 'store 3/3: lmdb patients' if ( $project->cache_verbose() );
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

sub get_ids {
	my ( $project_name, $region ) = @_;
	my $buffer = new GBuffer;
	my $project = $buffer->newProject( -name => $project_name );
	my $chr = $project->getChromosome( $region->{chromosome} );
	my $reference = $chr->getReferences( $region->{start}, $region->{end} )->[0];
	my @all;
	my @patient_names = sort { $a cmp $b } map { $_->name } @{ $project->getPatients };
	my $hpatients;
	for ( my $i = 0 ; $i < @patient_names ; $i++ ) {
		$hpatients->{ $patient_names[$i] } = $i;
	}
	foreach my $variation ( @{ $reference->getStructuralVariations } ) {
		
		#$variation->gnomad;
		#$variation->dejaVuInfosForDiag();
		$variation->name();
		$variation->cosmic();
		$variation->getGenes();
		#$variation->score_clinical_local();
		$variation->score_clinvar();
		$variation->text_clinvar();
		$variation->comment_clinical_local();
		$variation->cadd_score();
		$variation->hgmd_id();
		$variation->get_codon_text(1);
		$variation->get_codon_text(-1);
		
		my $hv;
		my $array_patients;
		my $aho = [];
		my $ap  = [];
		foreach my $pat ( @{ $variation->getPatients() } ) {
			my $pn = $pat->name();
			my $hp;
			$hp->{name} = $pn;
			my $patient_id = $hpatients->{$pn};
			push( @$ap, $patient_id );
			push( @$aho, $patient_id ) if ( $variation->isHomozygote($pat) );
			$variation->{patients_details}->{$pn}->{vcf_infos} = $variation->{check_id};
			$variation->{patients_details}->{$pn}->{he} = $variation->{annex}->{ $pat->id() }->{he};
			$variation->{patients_details}->{$pn}->{ho} = $variation->{annex}->{ $pat->id() }->{ho};
			$variation->{patients_details}->{$pn}->{he_ho_details} = "he";
			$variation->{patients_details}->{$pn}->{he_ho_details} = "ho" if ( $variation->{annex}->{ $pat->id() }->{ho} eq '1' );
			$variation->{patients_details}->{$pn}->{he_ho_details} .= ':'.$variation->{annex}->{ $pat->id() }->{nb_all_ref}.':'.$variation->{annex}->{ $pat->id() }->{nb_all_mut};
			$variation->{patients_details}->{$pn}->{type} = "he";
			$variation->{patients_details}->{$pn}->{type} = "ho" if ( $variation->isHomozygote($pat) );
		}
		delete $variation->{buffer};
		delete $variation->{project};
		$variation->{heho_string} = join( ",", sort { $a <=> $b } @$ap );
		if ( scalar(@$aho) ) {
			$variation->{heho_string} = $variation->{heho_string}." ".join( ",", sort { $a <=> $b } @$aho ) . " HO";
		}
		
		# line to prepare dejavu global;
		$hv->{obj}   = compress( freeze($variation) );
		$hv->{start} = $variation->start;
		$hv->{end}   = $variation->end;
		$hv->{id}    = $variation->id;
		push( @all, $hv );
	}
	my @sort = sort { $a->{start} <=> $b->{start} or $a->{end} <=> $b->{end} } @all;
	return \@sort;
}
