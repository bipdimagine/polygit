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
#use lib "/software/polyweb/poly-disk/poly-src/GenBo/lib/obj-nodb/";
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use GenBoNoSqlRocksGenome;
use GenBoNoSqlRocksVariation;
use List::Util qw(shuffle);
use Sys::Hostname;
use GenBoNoSqlRocksPolyviewerVariant;
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/packages/";
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/polyviewer/";
use PolyviewerVariant;
use Deep::Hash::Utils qw(reach slurp nest deepvalue);
use Carp;
require "$RealBin/../../../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
 my $host = hostname();


warn "*_*_*_*_*_".$host."*_*_*_*_*_";
#use Cache_Commons;



my $fork = 1;
my ($project_name, $chr_name, $no_verbose, $skip_pseudo_autosomal,$version,$annot_version);

GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'annot_version=s'    => \$annot_version,
	'chr=s'        => \$chr_name,
	'no_verbose=s' => \$no_verbose,
	'skip_pseudo_autosomal=s' => \$skip_pseudo_autosomal,
	'version=s' => \$version,
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

$project->preload_patients();
$project->buffer->disconnect();
$project->buffer->{dbh} ="-";
#system("/software/bin/vmtouch -t /data-isilon/public-data/repository/HG19/gnomad-exome/2.1/lmdb//snps/".$chr->name);
#warn "/data-isilon/public-data/repository/HG19/gnomad-exome/2.1/lmdb//snps/".$chr->name;
#warn "------------------------";
my $patients = $chr->project->getPatients();
my $regions;
my @lVarRegion;
my $new_pos       = 1;
my $old_pos       = 1;
my $nb_variants   = 0;
my $nb_var_region = 0;
#warn $fork;

my $tmp_dir = $project->rocks_pipeline_directory();

my $rg;
$rg = GenBoNoSqlRocksGenome->new(dir=>$project->rocks_pipeline_directory("genbo"),mode=>"c",index=>"genomic",chromosome=>$chr_name,genome=>"HG19",pack=>"",description=>[]);

my $rg_polyviewer = GenBoNoSqlRocksGenome->new(dir=>$project->rocks_pipeline_directory("polyviewer"),mode=>"c",index=>"genomic",chromosome=>$chr_name,genome=>"HG19",pack=>"",description=>[]);

$regions = $rg->regions;
$rg_polyviewer->regions();
my $pm = new Parallel::ForkManager($fork);
my $all;
my $process;
$pm->run_on_finish(
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;
		unless (defined($hRes) or $exit_code > 0) {
			$nbErrors++;
			print qq|No message received from child process $exit_code $pid!\n|;
			return;
		}
		delete $process->{$hRes->{end}};

	}
);
my $cpt = 1;



my $t = time;
foreach my $region (values @$regions) {
	$cpt++;
	unless ( exists $region->{chromosome} ) {
		$region->{chromosome} = $chr->id();
	}
	$t ++;
	$process->{$cpt}  ++;
	
	my $pid = $pm->start and next;
	my $time_start = time;
	my $hres = {};
	warn "start";
	my ($all,$ids)        = get_ids( $project_name, $region );
	my $time_end   = time;
	
	$hres->{ttime}    = abs( $time_start - $time_end );
	$hres->{end}    =  $cpt;
	delete $buffer->{lmdb_score};
	$pm->finish( 0, $hres );
}
$pm->wait_all_children();
if(keys %$process){
	warn Dumper $process;
	die();
}
$rg->close();
$rg_polyviewer->close();

warn "END ANNOTATION => ".abs(time-$t);
warn " START PHASE 2 ";




$rg_polyviewer = GenBoNoSqlRocksGenome->new(dir=>$project->rocks_pipeline_directory("polyviewer"),mode=>"r",index=>"genomic",chromosome=>$chr_name,genome=>"HG19",pack=>"",description=>[]);




$rg_polyviewer = GenBoNoSqlRocksGenome->new(dir=>$project->rocks_pipeline_directory("polyviewer"),mode=>"r",index=>"genomic",chromosome=>$chr_name,genome=>"HG19",pack=>"",description=>[]);
$rg = GenBoNoSqlRocksGenome->new(dir=>$project->rocks_pipeline_directory("genbo"),mode=>"r",index=>"genomic",chromosome=>$chr_name,genome=>"HG19",pack=>"",description=>[]);
my $final_dir = $project->rocks_directory();
my $finalrg = GenBoNoSqlRocksVariation->new(dir=>$project->rocks_directory("genbo"),mode=>"c",name=>$chr_name.".genbo.rocks");
my $final_index_genbo_id = GenBoNoSqlRocksVariation->new(dir=>$project->rocks_directory("index"),mode=>"c",name=>$chr_name);
my $final_polyviewer = GenBoNoSqlRocksPolyviewerVariant->new(dir=>$project->rocks_directory("polyviewer"),mode=>"c",name=>$chr_name.".polyviewer_variant");
my $fp = GenBoNoSqlRocks->new(dir=>$project->rocks_directory("polyviewer_raw"),mode=>"c",name=>$chr_name.".polyviewer_variant");
#
my $nb_regions = 0;
my $index =0;
my $array;
 $t = time;
 warn "------------------------------------------------------\n";
my $ztotal;
	foreach my $r (@{$rg->regions}) {
		my $start = $index;
		my $no = $rg->nosql_rocks($r);
		my $no_polyviewer = $rg_polyviewer->nosql_rocks($r);
		$nb_regions ++;
		my $iter = $no->rocks->new_iterator->seek_to_first;
		while (my ($var_id, $value) = $iter->each) {
			my $index = $finalrg->raw_put_batch_variation($var_id,$value);
		
			my $pv = $no_polyviewer->get($var_id);
			
			$pv->{global_vector_id} = $chr->name."!".$index;
			$final_index_genbo_id->put("g*".$pv->{id},$pv->{global_vector_id}.":".$var_id);
			$final_index_genbo_id->put("r*".$chr->name."_".$var_id,$pv->{global_vector_id}.":".$pv->{id});
			$final_index_genbo_id->put("v*".$pv->{global_vector_id},$chr->name."_".$var_id.":".$pv->{id});
			#$final_index->put_batch_raw($pv->{global_vector_id},)
			$fp->put_batch($pv->{global_vector_id},$pv);
			my $hgenes = delete $pv->{hgenes};
			my $patients = delete $pv->{hpatients};
			$pv->{global_vector_id} = $chr->name."!".$index;
			$final_polyviewer->put_batch_PolyviewerVariant($pv);
			$final_polyviewer->put_batch_patient($pv->id,$patients);
			$final_polyviewer->put_batch_gene($pv->id,$hgenes);
		
			$index++;
		}
		
		#push(@$array,[$r->{id},$start,($index-1)]);
		$no->close();
		#$finalrg->write_batch();
		$fp->write_batch;
		$final_index_genbo_id->write_batch;
	}
	$fp->close();
	$final_polyviewer->write_batch();
	$final_polyviewer->close();
	$finalrg->write_batch();
	$finalrg->close();
	$final_index_genbo_id->write_batch;
	$final_index_genbo_id->close();
	warn "------------------------------------------------------\n";
	warn "********* END   ".abs(time -$t)." :: ".$index;
	#$rg2->save_vector_index_region($array);
	#warn Dumper $array;
	$rg->close();	
	exit(0);
	warn "------------------------------------------------------\n";
	warn "+++++++++++ TEST ++++++++";
	 $final_polyviewer = GenBoNoSqlRocksPolyviewerVariant->new(dir=>$final_dir."polyviewer/",mode=>"r",name=>$chr_name.".polyviewer_variant");
	for (my $i =0 ; $i<100;$i++){
		my $pv = $final_polyviewer->getPolyviewerVariant($i);
		foreach my $gid (@{$pv->{genes_id}}){
			 $final_polyviewer->getGene($i,$gid,$pv);
		}
		foreach my $pid (@{$pv->{patients_id}}){
			 $final_polyviewer->getPatient($i,$pid,$pv);
		}
	}
	
	
  exit(0);


sub get_ids {
	my ( $project_name, $region ) = @_;
	my $ids = [];
	my $buffer = new GBuffer;
	 $buffer->vmtouch(1);
	
	
	my $project = $buffer->newProject( -name => $project_name );
	$project->preload_patients();
	$project->buffer->disconnect();
	my $no = $rg->nosql_rocks($region);
	my $no_polyviewer = $rg_polyviewer->nosql_rocks($region);
	#$project->buffer->{dbh} ="-";
	my $chr = $project->getChromosome( $region->{chromosome} );
	my $reference = $chr->getReferences( $region->{start}, $region->{end} )->[0];
	
	my @all;
	my @patient_names = sort { $a cmp $b } map { $_->name } @{ $project->getPatients };
	my $hpatients;
	for ( my $i = 0 ; $i < @patient_names ; $i++ ) {
		$hpatients->{ $patient_names[$i] } = $i;
	}
	my $t = time;
	my $vs = $reference->getStructuralVariations;#getStructuralVariations
	$t =time;
	my $nb = 0;
	foreach my $variation ( @{$vs } ) {
		my $debug ;
		$nb ++;
		warn $nb."/".scalar @{$vs } if $nb %3000 == 0; 
		die() if $variation->name() eq "11-800438-ins-116";
		$variation->gnomad("cache");
		$variation->cadd_score();
		$variation->name();
		$variation->getGenes();
		$variation->score_clinical_local();
		$variation->score_clinvar();
		$variation->text_clinvar();
		$variation->hgmd_id();		
		$variation->getChromosome();	
		$variation->annotation();
	
		my $hv;
		my $array_patients;
		my $aho = [];
		my $ap  = [];
		$variation->getPatients();
		$variation->sequencing_infos();
		$variation->split_read_infos();
		$variation->dejaVuInfosForDiag2();
		my $a	 = delete $variation->{array_dejavu};
		$variation->{cad} = pack("w".scalar(@$a),@$a);
		my $vp =  PolyviewerVariant->new();
		$vp->setLmdbVariant($variation);
		$vp->{hgenes} = {};
		$vp->{genes_id} = [];
		foreach my $g (@{$variation->getGenes}){
			my $h = $vp->set_gene($variation,$g);
			$vp->{hgenes}->{$g->{id}} = $h;
			push(@{$vp->{genes_id}},$g->{id});
		}
		$vp->{hpatients} ={};
		$vp->{patients_id} = [];
		foreach my $p (@{$variation->getPatients}){
			my $h = $vp->set_patient_cache($variation,$p);
			$vp->{patients_calling}->{$p->id} =$h; 
		}
		
		
		# line to prepare dejavu global;
		my $ref = ref($variation);
		if ($ref eq 'GenBoVariation'){
					bless $variation , 'GenBoVariationCache';
		}
		elsif  ($ref eq 'GenBoLargeDeletion'){
					bless $variation , 'GenBoLargeDeletionCache';
					
		}
		elsif  ($ref eq 'GenBoLargeInsertion'){
					bless $variation , 'GenBoLargeDuplicationCache';
		}
		elsif  ($ref eq 'GenBoDeletion'){
					bless $variation , 'GenBoDeletionCache'; 
		}
		elsif  ($ref eq 'GenBoInsertion'){
					bless $variation , 'GenBoInsertionCache';
		}
		elsif  ($ref eq 'GenBoLargeDuplication'){
					bless $variation , 'GenBoLargeDuplicationCache';
		}
		else {
			confess($ref);
		}
	
	#	my $hvariant =  update_variant_editor::construct_hash_variant_global ( $project, $variation,undef,1);
	#	delete $hvariant->{html};
		#warn Dumper $hvariant ;
		#die();
		warn Dumper $vp if $variation->id eq "7_100555983_del-51268";
		$no_polyviewer->put_batch($variation->rocksdb_id,$vp);
		delete $variation->{array_dejavu};
		delete $variation->{references_object};
		delete $variation->{dejaVuInfosForDiag2};
		delete $variation->{annex};
		delete $variation->{buffer};
		delete $variation->{project};
		
		$no->put_batch($variation->rocksdb_id,$variation);
	
		#warn Dumper $variation->annex();
		#die();
		
	}
	$no->rocks->write($no->batch);
	$no_polyviewer->write_batch();
	$rg->nosql_rocks($region)->close;
	delete $no->{batch};
	return undef;
}