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
use List::Util qw(shuffle);
require "$RealBin/Cache_Commons.pm";
use Sys::Hostname;

 my $host = hostname();


warn "*_*_*_*_*_".$host."*_*_*_*_*_";
#use Cache_Commons;



my $fork = 1;
my ($project_name, $chr_name, $no_verbose, $skip_pseudo_autosomal,$version,$annot_version,$is_rna_junction_analyse);

GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'annot_version=s'    => \$annot_version,
	'chr=s'        => \$chr_name,
	'no_verbose=s' => \$no_verbose,
	'skip_pseudo_autosomal=s' => \$skip_pseudo_autosomal,
	'version=s' => \$version,
	'rna_junction' => \$is_rna_junction_analyse
);

warn "*_*_*_*_*_ fork :".$fork."*_*_*_*_*_";
`ulimit -Su unlimited && echo toto`;
system("ulimit -Su unlimited");
system("ulimit -a >/tmp/test");
my (@z) = `ulimit -a `;

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }

my $nbErrors = 0;
my $buffer = new GBuffer;
$buffer->vmtouch(1);
#warn $buffer->config->{'public_data_annotation'}->{root};
#my $color = $colors[ rand @colors ];
my $project = $buffer->newProject( -name => $project_name );
if ($annot_version) {
	$project->changeAnnotationVersion($annot_version);
}
if ($no_verbose) {
	$project->cache_verbose(0);
}
my $chr = $project->getChromosome($chr_name);
if (-e $project->getCacheBitVectorDir()."/lmdb_cache/".$chr->id().".empty") {
	unlink $project->getCacheBitVectorDir()."/lmdb_cache/".$chr->id().".empty";
}
warn "\n### CACHE: store ids step\n" if ( $project->cache_verbose() );
#warn "------------------------";
#system("/software/bin/vmtouch -t /data-isilon/public-data/repository/HG19/cadd/1.6/lmdb/".$chr->name.".uc "."/data-isilon/public-data/repository/HG19/gnomad-genome/2.1/lmdb//snps/".$chr->name." "."/data-isilon/public-data/repository/HG19/gnomad-exome/2.1/lmdb//snps/".$chr->name);#/software/bin/vmtouch -t ".$no1->filename." "
#system("/software/bin/vmtouch -t /data-isilon/public-data/repository/HG19/cadd/1.6/lmdb/".$chr->name.".uc");#/software/bin/vmtouch -t ".$no1->filename." "
#system("/software/bin/vmtouch -t /data-isilon/public-data/repository/HG19/gnomad-genome/2.1/lmdb//snps/".$chr->name);

my $file_ok = $project->lmdb_cache_variations_dir().'/chr'.$chr_name.'_store_ids.ok';

$project->preload_patients();
$project->buffer->disconnect();
$project->buffer->{dbh} ="-";
#system("/software/bin/vmtouch -t /data-isilon/public-data/repository/HG19/gnomad-exome/2.1/lmdb//snps/".$chr->name);
#warn "/data-isilon/public-data/repository/HG19/gnomad-exome/2.1/lmdb//snps/".$chr->name;
#warn "------------------------";
my $patients = $chr->project->getPatients();
my @lVarRegion;
my $new_pos       = 1;
my $old_pos       = 1;
my $nb_variants   = 0;
my $nb_var_region = 0;

my $process;
my $cpt = 1;
#warn Dumper $regions;
my $hregion;

my $time_start = time;
my ($all,$ids) = get_junctions_ids($project, $chr);
my $time_end   = time;

#$hres->{variants} = $all;
$nb_variants += scalar(@{$ids});
foreach my $hv (@{$ids}) {
	$chr->{cache_hash_get_var_ids}->{$hv}++;
}


if ($nbErrors > 0) {
	confess("\n\nERRORS: $nbErrors errors found... confess...\n\n");
}

my $size_variants = scalar @{$all};

############


my $t = time;

warn 'store 1/4: rocks variations' if ( $project->cache_verbose() );
my $hh;
my $h_junct_index;
my $finalrg = GenBoNoSqlRocksVariation->new(dir=>$project->rocks_directory("genbo"),mode=>"c",name=>$chr_name.".genbo.rocks");
foreach my $hv ( @{$all} ) {
	my $var_id = $hv->{id};
	my $junction = thaw( decompress( $hv->{obj} ) );
	my $index = $finalrg->put_batch_variation($var_id, $junction);
	$h_junct_index->{$junction->id()} = $index;
	$hh->{$var_id} = $junction->{heho_string};
}
$finalrg->write_batch();
$finalrg->close();



my $rocks3 = $chr->rocks_vector("c");
$rocks3->size($size_variants);

if($size_variants == 0 ){
	my $vnull =  Bit::Vector->new(0);
	foreach my $patient (@{$project->getPatients}){	
		$rocks3->put_vector_patient_batch($patient,"all",$vnull);
		$rocks3->put_vector_patient_batch($patient,"he",$vnull);
		$rocks3->put_vector_patient_batch($patient,"ho",$vnull);
		$rocks3->put_vector_patient_batch($patient,"ri",$vnull);
		$rocks3->put_vector_patient_batch($patient,"se",$vnull);
	}
	$rocks3->write_batch();
	$rocks3->close();
	warn "no variants !!!!";
	warn "\n\nEND!\n";
	
	open (FILE, ">$file_ok");
	print FILE "OK";
	close(FILE);
	exit(0);
} 

my $h;
my $vector_variation_type;
my $bitv =  Bit::Vector->new($size_variants);
$vector_variation_type->{junctions_object} = $bitv;
$h->{bitvector} = $bitv;
$rocks3->put_batch_vector_chromosome('junctions_object' ,$bitv);


my $v1 = Bit::Vector->new($size_variants);
$v1->Fill;
$rocks3->put_batch("all",$v1);
$rocks3->write_batch();


########
	


#end fork now I have a file with all variation and  json  it's time to sort this file and store it in lmdb database
#construct intspan 	for patient I will store lmdb_id in intspan
my @patient_names = sort { $a cmp $b } map { $_->name } @{ $project->getPatients };
my $index_patients = 0;
my $hpatients;
my $categories_junctions = Cache_Commons::spec_categories_junctions();
my @categorie_patient = ( "all", "ri", "se", 'junc_ratio_10', 'junc_ratio_20', 'junc_ratio_30', 'junc_ratio_40', 'junc_ratio_50', 'junc_ratio_60', 'junc_ratio_70', 'junc_ratio_80', 'junc_ratio_90', );
foreach my $c ( keys %{ $categories_junctions->{patients} } ) {
	push(@categorie_patient, $c);
}

foreach my $pname (@patient_names) {
	$hpatients->{$pname}->{index} = $index_patients++;
	$hpatients->{$pname}->{name}  = $pname;
	foreach my $c (@categorie_patient) {
		$hpatients->{$pname}->{vector}->{$c} = Bit::Vector->new($size_variants);
	}
}


#initialisation global categorie
my $intspan_global_type;
my $categories = Cache_Commons::categories();
foreach my $g ( keys %{ $categories->{global}->{variation_type} } ) {
	$intspan_global_type->{$g} = Bit::Vector->new($size_variants);
}
foreach my $g ( keys %{ $categories_junctions->{global} } ) {
	$intspan_global_type->{$g} = Bit::Vector->new($size_variants);
}

	
foreach my $hv ( @{$all} ) {
	my $var_id = $hv->{id};
	my $junction = thaw( decompress( $hv->{obj} ) );
	$junction->{buffer} = $buffer;
	$junction->{project} = $project;
	my $index_lmdb = $h_junct_index->{$junction->id()};
	$size_variants++;
	foreach my $patient (@{ $junction->getPatients() }) {
		$hpatients->{$patient->name()}->{vector}->{all}->Bit_On($index_lmdb);
		$hpatients->{$patient->name()}->{vector}->{ri}->Bit_On($index_lmdb) if ($junction->isRI($patient));
		$hpatients->{$patient->name()}->{vector}->{se}->Bit_On($index_lmdb) if ($junction->isSE($patient));
		
		my $junction_type_description = $junction->getTypeDescription($patient);
		my $ratio = int($junction->get_percent_new_count($patient));
		
		$hpatients->{$patient->name()}->{vector}->{junc_ratio_10}->Bit_On($index_lmdb) if ($ratio >= 10);
		$hpatients->{$patient->name()}->{vector}->{junc_ratio_20}->Bit_On($index_lmdb) if ($ratio >= 20);
		$hpatients->{$patient->name()}->{vector}->{junc_ratio_30}->Bit_On($index_lmdb) if ($ratio >= 30);
		$hpatients->{$patient->name()}->{vector}->{junc_ratio_40}->Bit_On($index_lmdb) if ($ratio >= 40);
		$hpatients->{$patient->name()}->{vector}->{junc_ratio_50}->Bit_On($index_lmdb) if ($ratio >= 50);
		$hpatients->{$patient->name()}->{vector}->{junc_ratio_60}->Bit_On($index_lmdb) if ($ratio >= 60);
		$hpatients->{$patient->name()}->{vector}->{junc_ratio_70}->Bit_On($index_lmdb) if ($ratio >= 70);
		$hpatients->{$patient->name()}->{vector}->{junc_ratio_80}->Bit_On($index_lmdb) if ($ratio >= 80);
		$hpatients->{$patient->name()}->{vector}->{junc_ratio_90}->Bit_On($index_lmdb) if ($ratio >= 90);
		
		unless (exists $intspan_global_type->{$junction_type_description}) {
			$intspan_global_type->{$junction_type_description} = Bit::Vector->new($size_variants);
			$categories_junctions->{global}->{$junction_type_description} = undef;
		}
		$intspan_global_type->{$junction_type_description}->Bit_On($index_lmdb);
	}
	unless ( exists $intspan_global_type->{ $junction->type() } ) {
		warn "\n\nERROR: doesn't exists \$intspan_global_type->{".$junction->type() . "}\n\n";
		warn Dumper keys %$intspan_global_type;
		confess;
	}
	$intspan_global_type->{ $junction->type }->Bit_On($index_lmdb);
	
	my $dv = $junction->dejavu_patients("all");
	$intspan_global_type->{dejavu_5}->Bit_On($index_lmdb) if ($dv <= 5);
	$intspan_global_type->{dejavu_10}->Bit_On($index_lmdb) if ($dv <= 10);
	$intspan_global_type->{dejavu_15}->Bit_On($index_lmdb) if ($dv <= 15);
	$intspan_global_type->{dejavu_20}->Bit_On($index_lmdb) if ($dv <= 20);
	$intspan_global_type->{dejavu_25}->Bit_On($index_lmdb) if ($dv <= 25);
	$intspan_global_type->{dejavu_30}->Bit_On($index_lmdb) if ($dv <= 30);
	$intspan_global_type->{dejavu_40}->Bit_On($index_lmdb) if ($dv <= 40);
	$intspan_global_type->{dejavu_50}->Bit_On($index_lmdb) if ($dv <= 50);
	$intspan_global_type->{dejavu_60}->Bit_On($index_lmdb) if ($dv <= 60);
	$intspan_global_type->{dejavu_70}->Bit_On($index_lmdb) if ($dv <= 70);
	$intspan_global_type->{dejavu_80}->Bit_On($index_lmdb) if ($dv <= 80);
	$intspan_global_type->{dejavu_90}->Bit_On($index_lmdb) if ($dv <= 90);
	
	my $dv_r10 =  $junction->dejavu_patients("10");
	$intspan_global_type->{dejavu_5_r10}->Bit_On($index_lmdb) if ($dv_r10 <= 5);
	$intspan_global_type->{dejavu_10_r10}->Bit_On($index_lmdb) if ($dv_r10 <= 10);
	$intspan_global_type->{dejavu_15_r10}->Bit_On($index_lmdb) if ($dv_r10 <= 15);
	$intspan_global_type->{dejavu_20_r10}->Bit_On($index_lmdb) if ($dv_r10 <= 20);
	$intspan_global_type->{dejavu_25_r10}->Bit_On($index_lmdb) if ($dv_r10 <= 25);
	$intspan_global_type->{dejavu_30_r10}->Bit_On($index_lmdb) if ($dv_r10 <= 30);
	$intspan_global_type->{dejavu_40_r10}->Bit_On($index_lmdb) if ($dv_r10 <= 40);
	$intspan_global_type->{dejavu_50_r10}->Bit_On($index_lmdb) if ($dv_r10 <= 50);
	$intspan_global_type->{dejavu_60_r10}->Bit_On($index_lmdb) if ($dv_r10 <= 60);
	$intspan_global_type->{dejavu_70_r10}->Bit_On($index_lmdb) if ($dv_r10 <= 70);
	$intspan_global_type->{dejavu_80_r10}->Bit_On($index_lmdb) if ($dv_r10 <= 80);
	$intspan_global_type->{dejavu_90_r10}->Bit_On($index_lmdb) if ($dv_r10 <= 90);
}


warn 'store 2/4: rocks patients';
my $vector_patients1;
foreach my $patient (@{ $project->getPatients() }) {
	foreach my $cat (keys %{$hpatients->{$patient->name()}->{vector}}) {
		my $bitv = $hpatients->{$patient->name()}->{vector}->{$cat};
#		warn $patient->name.' -> '.$cat.': '.$bitv->Norm;
		$vector_patients1->{$patient->id."_".$cat} = $bitv if $cat eq "all" ;
		$rocks3->put_vector_patient_batch($patient,$cat,$bitv);
	}
	my $vnull =  Bit::Vector->new($size_variants);
	$rocks3->put_vector_patient_batch($patient,"he",$vnull);
	$rocks3->put_vector_patient_batch($patient,"ho",$vnull);
}
$rocks3->write_batch();

warn 'store 3/4: rocks chr categories';
foreach my $cat (keys %{$intspan_global_type}) {
	my $bitv = $intspan_global_type->{$cat};
#	warn 'chr -> '.$cat.': '.$bitv->Norm;
	$rocks3->put_batch_vector_chromosome($cat,$bitv);
}
$rocks3->write_batch();

my $vall = Bit::Vector->new($size_variants);
$vall->Fill;
$rocks3->put_batch("all",$vall);
$rocks3->write_batch();

$rocks3->write_config();
$project = undef;
$buffer = undef;

warn "\n\nEND!\n";

open (FILE, ">$file_ok");
print FILE "OK";
close(FILE);


sub get_junctions_ids {
	my ( $project, $chr ) = @_;
	my $ids = [];
	my $buffer = new GBuffer;
	$buffer->vmtouch(1);
	$project->preload_patients();
	#$project->buffer->disconnect();
	my @all;
	my @patient_names = sort { $a cmp $b } map { $_->name } @{ $project->getPatients };
	my $hpatients;
	for ( my $i = 0 ; $i < @patient_names ; $i++ ) {
		$hpatients->{ $patient_names[$i] } = $i;
	}
	my $vs = $chr->getJunctions();
	my $vs_sorted = sort_list_junctions_by_positions($vs);
	
	
	
	my $i = 0;
	my $max = scalar(@$vs_sorted)-1;
	my $h_ri_aval_amont;
	foreach my $junction ( @{$vs_sorted} ) {
		my $debug ;
		next if ($junction->getChromosome->id() ne $chr->id());
		$junction->id();
		$junction->name();
		$junction->annex();
		$junction->setPatients();
		my $h_exons_introns = $junction->get_hash_exons_introns();
		foreach my $patient (@{ $junction->getPatients() }) {
			my $is_ri_aval = $junction->is_ri_aval($patient);
			my $is_ri_amont = $junction->is_ri_amont($patient);
			if ($is_ri_aval or $is_ri_amont) {
				my $to_check;
				$to_check = 'ri_amont' if $is_ri_amont;
				$to_check = 'ri_aval' if $is_ri_aval;
				foreach my $tid (sort keys %{$h_exons_introns}) {
					my @lPos = (sort keys %{$h_exons_introns->{$tid}->{by_pos}});
					my $first_exon_intron = $h_exons_introns->{$tid}->{by_pos}->{$lPos[0]};
					my $last_exon_intron = $h_exons_introns->{$tid}->{by_pos}->{$lPos[-1]};
					$h_ri_aval_amont->{$patient->name()}->{$tid}->{$first_exon_intron}->{$to_check} = $junction->id();
				} 
			}
			
			my @lJunctions_to_check;
			my $start_search = 0;
			$start_search = $i-100 if ($i-100) >= 0;
			my $end_search = $i+100;
			$end_search = $max if ($i+100) > $max;
			while ($start_search <= $end_search) {
				push(@lJunctions_to_check, $vs_sorted->[$start_search]) if $start_search != $i;
				$start_search++;
			} 
			$junction->{hash_noise_by_patient}->{$patient->name()} = get_hash_noise($junction, $patient, \@lJunctions_to_check);
			$junction->get_noise_score($patient);
			$junction->junction_score_without_dejavu_global($patient);
			$junction->getHashSpliceAiNearStartEnd();
		}
		$i++;
	}
	
	foreach my $junction (@$vs_sorted) {
		foreach my $patient (@{ $junction->getPatients() }) {
			my $is_ri_aval = $junction->is_ri_aval($patient);
			my $is_ri_amont = $junction->is_ri_amont($patient);
			if ($is_ri_aval or $is_ri_amont) {
				my $to_check;
				$to_check = 'ri_amont' if $is_ri_aval;
				$to_check = 'ri_aval' if $is_ri_amont;
				my $h_exons_introns = $junction->get_hash_exons_introns();
				
				foreach my $tid (sort keys %{$h_exons_introns}) {
					my @lPos = (sort keys %{$h_exons_introns->{$tid}->{by_pos}});
					my $first_exon_intron = $h_exons_introns->{$tid}->{by_pos}->{$lPos[0]};
					my $last_exon_intron = $h_exons_introns->{$tid}->{by_pos}->{$lPos[-1]};
					
					if (exists $h_ri_aval_amont->{$patient->name()}->{$tid}->{$first_exon_intron}->{$to_check}) {
						my $other_junction_id = $h_ri_aval_amont->{$patient->name()}->{$tid}->{$first_exon_intron}->{$to_check};
						$junction->{get_hash_junctions_linked_to_me}->{$patient->name()}->{$other_junction_id}->{$tid}->{$first_exon_intron} = $h_ri_aval_amont->{$patient->name()}->{$tid}->{$first_exon_intron};
					}
					if (exists $h_ri_aval_amont->{$patient->name()}->{$tid}->{$last_exon_intron}->{$to_check}) {
						my $other_junction_id = $h_ri_aval_amont->{$patient->name()}->{$tid}->{$last_exon_intron}->{$to_check};
						$junction->{get_hash_junctions_linked_to_me}->{$patient->name()}->{$other_junction_id}->{$tid}->{$last_exon_intron} = $h_ri_aval_amont->{$patient->name()}->{$tid}->{$last_exon_intron};
					}
				} 
			}
		}
		
		my $hv;
		my $ref = ref($junction);
		if ($ref eq 'GenBoJunction'){
			bless $junction , 'GenBoJunctionCache';
		}
		delete $junction->{buffer};
		delete $junction->{project};
		delete $junction->{patients_object};
		delete $junction->{genes_object};
		delete $junction->{coverage_obj};
		eval {
			$hv->{obj}   = compress( freeze($junction) );
		};
		if ($@) {
			warn "\n\n";
			warn ref ($junction);
			warn $junction->{id};
			warn Dumper sort keys %$junction;
		}
		
		$hv->{start} = $junction->start;
		$hv->{end}   = $junction->end;
		$hv->{id}    = $junction->id;
		push(@$ids, $junction->id);
		push( @all, $hv );
	}
	my @sort = sort { $a->{start} <=> $b->{start} or $a->{end} <=> $b->{end} } @all;
	return (\@sort),$ids;
}

sub sort_list_junctions_by_positions {
	my ($vs) = @_;
	my (@vs_sorted, $h_junctions_list_ids);
	my $i = 0;
	foreach my $junction ( @{$vs } ) {
		my $pos = $junction->start();
		$h_junctions_list_ids->{$pos}->{$i} = undef;
		$i++;
	}
	foreach my $pos (sort {$a <=> $b} keys %$h_junctions_list_ids) {
		foreach my $id (keys %{$h_junctions_list_ids->{$pos}}) {
			push(@vs_sorted, $vs->[$id]);
		}
	}
	return \@vs_sorted;
}

sub get_hash_noise {
	my ($junction, $patient, $list_junctions_to_check) = @_;
	my $h_noise;
	my $intspan_self = $junction->getStrictGenomicSpan();
	my $h_exons_introns = $junction->get_hash_exons_introns;
	foreach my $junction2 (@$list_junctions_to_check) {
		my $has_this_patient;
		foreach my $pat (@{$junction2->getPatients()}) {
			$has_this_patient = 1 if $pat->name() eq $patient->name();
		}
		next if not $has_this_patient;
		my $h_exons_introns_2 = $junction2->get_hash_exons_introns;
		if ($h_exons_introns_2) {
			foreach my $tid (keys %{$h_exons_introns_2}) {
				next if $tid eq 'by_pos';
				next unless (exists $h_exons_introns->{$tid});
				my @lPos = (sort keys %{$junction2->get_hash_exons_introns->{$tid}->{by_pos}});
				next if (@lPos == 0);
				foreach my $pos (keys %{$junction2->get_hash_exons_introns->{$tid}->{by_pos}}) {
					if (exists $h_exons_introns->{$tid}->{by_pos}->{$pos}) {
						my $jid2 = $junction2->id();
						$jid2 =~ s/_RI//;
						$jid2 =~ s/_SE//;
						$h_noise->{$tid}->{$pos}->{$jid2} = undef;
						$h_noise->{all}->{$jid2} = undef;
					}
				}
			}
		}
		my $intspan_this = $junction2->getStrictGenomicSpan();
		my $inter1 = $intspan_self->intersection( $intspan_this );
		if (not $inter1->is_empty()) {
			my $jid2 = $junction2->id();
			$jid2 =~ s/_RI//;
			$jid2 =~ s/_SE//;
			$h_noise->{all}->{$jid2} = undef;
		}
	}
	return 	$h_noise;
}
