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
use polyviewer_html;
use HTML::Packer;
use Compress::Zstd;

require "$RealBin/../../../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
 my $host = hostname();



#use Cache_Commons;



my $fork = 1;
my ($project_name, $chr_name, $verbose, $skip_pseudo_autosomal,$version,$annot_version);
my $ok_file;
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'annot_version=s'    => \$annot_version,
	'chr=s'        => \$chr_name,
	'verbose=s' => \$verbose,
	'skip_pseudo_autosomal=s' => \$skip_pseudo_autosomal,
	'version=s' => \$version,
	'file=s' => \$ok_file,
);


if ($ok_file && -e $ok_file) {
 	system("rm $ok_file");
}



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

if (not $verbose) {
	$project->cache_verbose(0);
}
warn "*_*_*_*_*_".$host."*_*_*_*_*_".$chr_name  if ( $project->cache_verbose() );
warn "*_*_*_*_*_ fork :".$fork."*_*_*_*_*_"  if ( $project->cache_verbose() );
warn hostname if ( $project->cache_verbose() );
$ok_file = $project->rocks_directory("logs")."/cache_store_ids.".$chr_name.".ok" if not $ok_file;

$project->preload_patients();
if ($annot_version) {
	$project->changeAnnotationVersion($annot_version);
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

#$project->buffer->disconnect();
#$project->buffer->{dbh} ="-";
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
system("vmtouch -t ".$buffer->config_path("root","public_data")."/repository/".$project->annotation_genome_version()."/*/*/rocksdb/$chr_name.rocksdb");
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

$project->exomeProjects();
$project->countExomePatients();
$project->countSimilarPatients();
$project->in_this_run_patients();
my $t = time;
foreach my $region (values @$regions) {
	$cpt++;
	unless ( exists $region->{chromosome} ) {
		$region->{chromosome} = $chr->id();
	}
	$t ++;
	$process->{$cpt}  ++;
	
	
	$project->disconnect();
	my $pid = $pm->start and next;
#	warn 'BEFORE '.$cpt.'/'.scalar(@$regions).'  -  '.$region->{start}." ".$region->{end};
	my $time_start = time;
	my $hres = {};
	get_ids( $project_name, $region );
	my $time_end   = time;
	
	$hres->{ttime}    = abs( $time_start - $time_end );
	$hres->{end}    =  $cpt;
	delete $buffer->{lmdb_score};
#	warn 'END '.$cpt.'/'.scalar(@$regions).'  -  '.$region->{start}." ".$region->{end};
	
	$project->disconnect();
	$project->buffer->disconnect();
	$project->close_rocks();
	$pm->finish( 0, $hres );
}
$pm->wait_all_children();

#warn 'AFTER WAIT CHILDREN deb';

if(keys %$process){
	warn Dumper $process;
	die();
}
$rg->close();
$rg_polyviewer->close();

$project->disconnect();
$project->buffer->disconnect();
$project->buffer->dbh_deconnect();
$project->close_rocks();

#warn 'AFTER WAIT CHILDREN';
#
#warn "END ANNOTATION => ".abs(time-$t);
#warn " START PHASE 2 ";


my $rg_polyviewer_2 = GenBoNoSqlRocksGenome->new(dir=>$project->rocks_pipeline_directory("polyviewer"),mode=>"r",index=>"genomic",chromosome=>$chr_name,genome=>"HG19",pack=>"",description=>[]);
my $rg_2 = GenBoNoSqlRocksGenome->new(dir=>$project->rocks_pipeline_directory("genbo"),mode=>"r",index=>"genomic",chromosome=>$chr_name,genome=>"HG19",pack=>"",description=>[]);

my $final_dir = $project->rocks_directory();

#warn "\n\n";
#warn $final_dir;
#warn $project->rocks_directory("genbo");
#warn $project->rocks_directory("index");
#warn $project->rocks_directory("polyviewer_raw");
#warn "\n\n";


my $finalrg = GenBoNoSqlRocksVariation->new(dir=>$project->rocks_directory("genbo"),mode=>"c",name=>$chr_name.".genbo.rocks");
my $final_index_genbo_id = GenBoNoSqlRocksVariation->new(dir=>$project->rocks_directory("index"),mode=>"c",name=>$chr_name);
#my $final_polyviewer = GenBoNoSqlRocksPolyviewerVariant->new(dir=>$project->rocks_directory("polyviewer"),mode=>"c",name=>$chr_name.".polyviewer_variant");
my $fp = GenBoNoSqlRocks->new(dir=>$project->rocks_pipeline_directory("polyviewer_raw"),mode=>"c",name=>$chr_name);
#
my $nb_regions = 0;
my $index =0;
my $array;
 $t = time;
 warn "------------------------------------------------------\n"  if ( $project->cache_verbose() );
my $ztotal;
	foreach my $r (@{$rg_2->regions}) {
		my $start = $index;
		my $no = $rg_2->nosql_rocks($r);
		my $no_polyviewer = $rg_polyviewer_2->nosql_rocks($r);
		$nb_regions ++;
		my $iter = $no->rocks->new_iterator->seek_to_first;
		while (my ($var_id, $value) = $iter->each) {
			warn $var_id if ( $project->cache_verbose() );
			my $v = $no->decode($value);
			my $index = $finalrg->put_batch_variation($var_id,$v);
			my $pv = $no_polyviewer->get($var_id);
			$pv->{global_vector_id} = $chr->name."!".$index;
			$final_index_genbo_id->put("g*".$pv->{id},$pv->{global_vector_id}.":".$var_id);
			$final_index_genbo_id->put("r*".$chr->name."_".$var_id,$pv->{global_vector_id}.":".$pv->{id});
			$final_index_genbo_id->put("v*".$pv->{global_vector_id},$chr->name."_".$var_id.":".$pv->{id});
			#$final_index->put_batch_raw($pv->{global_vector_id},)
			my $html = delete $pv->{html};
			$fp->put_batch($pv->{global_vector_id},$pv);
			#$fp2->put_batch($pv->{global_vector_id}."!",$html);
			
			my $hgenes = delete $pv->{hgenes};
			my $patients = delete $pv->{hpatients};
			$pv->{global_vector_id} = $chr->name."!".$index;
		
			$index++;
		}
		
		$no->close();
		$no_polyviewer->close();
		$finalrg->write_batch();
		$fp->write_batch;
		#$fp2->write_batch;
		$final_index_genbo_id->write_batch;
	}
	
	
#	$final_polyviewer->write_batch();
#	$final_polyviewer->close();
	
	$rg_polyviewer_2->close();
	$rg_2->close();	
	$fp->close();
	$finalrg->close();
	$final_index_genbo_id->close();
	
	warn "------------------------------------------------------\n"  if ( $project->cache_verbose() );
	warn "********* END   ".abs(time -$t)." :: ".$index  if ( $project->cache_verbose() );
	if ($ok_file) {
		my $date = system("date");
		open (LOG, ">$ok_file");
		print LOG "OK\n";
		print LOG $date;
		close (LOG);
	}
	warn 'OK FILE: '.$ok_file  if ( $project->cache_verbose() );
	$project->close_rocks();
	
	#$rg2->save_vector_index_region($array);
	#warn Dumper $array;
#	warn $project->rocks_pipeline_directory();
	exit(0);



sub get_ids {
	my ( $project_name, $region ) = @_;
	my $ids = [];
	my $buffer = new GBuffer;
	my $patient = $project->getPatients()->[0];
	my @headers;
	my $packer = HTML::Packer->init();
	my $print_html = polyviewer_html->new( project => $project, patient => $patient,header=> \@headers,bgcolor=>"background-color:#607D8B" );
	#$print_html->init();
	my $project = $buffer->newProject( -name => $project_name );
	$project->preload_patients();
	#$project->buffer->disconnect();
	
	#$project->buffer->{dbh} ="-";
	my $chr = $project->getChromosome( $region->{chromosome} );
	my $reference = $chr->getReferences( $region->{start}, $region->{end} )->[0];
	my $gene_ids = $chr->genesIntervalTree()->fetch($region->{start}, $region->{end});
#
#
	if (@$gene_ids) {
		my $z;
		for (my $i=0;$i<@$gene_ids;$i++){
			$z->[$i] = "!".$gene_ids->[$i] ;
		}
			$project->rocksGenBo->prepare($z);
	}
	
#	warn "start ". $region->{start}." ".$region->{end};
	my @all;
	my @patient_names = sort { $a cmp $b } map { $_->name } @{ $project->getPatients };
	my $hpatients;
	for ( my $i = 0 ; $i < @patient_names ; $i++ ) {
		$hpatients->{ $patient_names[$i] } = $i;
	}
	my $t = time;
#	warn '   - check variants';
	my $vs = $reference->getStructuralVariations;
	
	my @arocksid = map {$_->rocksdb_id} @$vs;
	if (scalar(@$vs) > 0 and @arocksid){
		eval{
			$chr->rocksdb("gnomad")->prepare(\@arocksid);
			$chr->rocksdb("cadd")->prepare(\@arocksid);
			$chr->rocksdb("clinvar")->prepare(\@arocksid);
			$chr->rocksdb("spliceAI")->prepare(\@arocksid);
			$chr->rocksdb("dbscSNV")->prepare(\@arocksid);
			$chr->rocksdb("revel")->prepare(\@arocksid);
			$chr->rocksdb("hgmd")->prepare(\@arocksid);
			my $db = $chr->rocks_dejavu()->get_db($vs->[0]->start);
			$db->prepare(\@arocksid);
		};
		if ($@){
			warn "\n\n";
			warn "bug3 ";
			warn "\n\n";
			die();
		}
	}
	
	
	$t =time;
	my $nb = 0;
	#warn "prepare ". $region->{start}." ".$region->{end};
	my $hvariant;
	my $hpolyviewer;
	foreach my $variation ( @{$vs } ) {
		next if ($variation->type =~ /junction/ );
		
		#TODO: a ne pas oublier d'enlever apres TEST
		#next if ($variation->isCnv() );
		
		
		my $debug ;
		$nb ++;
		warn $nb."/".scalar @{$vs } if $nb %30000 == 0; 
		
		if ($variation->type =~ /inversion/ ){
		 $debug =1;
		}
		eval {
		$variation->gnomad("cache");
 		$variation->getGenes();
 		#$variation->getTranscripts();
		$variation->cadd_score();
		$variation->name();
		$variation->score_clinical_local();
		$variation->score_clinvar();
		$variation->text_clinvar();
		$variation->hgmd_id();		
		$variation->getChromosome();	
		$variation->dejaVuInfosForDiag2();
		$variation->annotation();
		};
		if ($@){
			warn "\n\n";
			warn Dumper $@;
			warn "bug here  ";
			warn "\n\n";
			die();
		}
		my $hv;
		my $array_patients;
		my $aho = [];
		my $ap  = [];
		$variation->getPatients();
		$variation->sequencing_infos();
		#$variation->split_read_infos();
		my $a	 = delete $variation->{array_dejavu};
		$a =[] unless $a;
		$variation->{cad} = pack("w".scalar(@$a),@$a);

		my $vp =  PolyviewerVariant->new();
		$vp->setLmdbVariant($variation);
		$vp->{hgenes} = {};
		$vp->{genes_id} = [];
		my $code =0;
		foreach my $g (@{$variation->getGenes}){
			my $h = $vp->set_gene($variation,$g);
			$h->{code} = $code;
			$vp->{hgenes}->{$g->{id}} = $h;
			
			push(@{$vp->{genes_id}},$g->{id});
			$code ++;
		}
		##############
		#	next;
		##############0
		$vp->{hpatients} ={};
		$vp->{patients_id} = [];
		my $dvp;
		foreach my $pat (@{$variation->getPatients}){
			
			foreach my $p (@{$pat->getFamily()->getMembers}){
				
				next if exists $dvp->{$p->id};
				$dvp->{$p->id} ++;
				my $h = $vp->set_patient_cache($variation,$p);
				$vp->{patients_calling}->{$p->id} =$h; 
			}
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
					bless $variation , 'GenBoLargeInsertionCache';
				
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
		elsif  ($ref eq 'GenBoInversion'){
					bless $variation , 'GenBoInversionCache';
		}
		elsif  ($ref eq 'GenBoMei'){
					bless $variation , 'GenBoMeiCache';
		}
		else {
			confess($ref);
		}
		my $t = 0;
		foreach my $g (@{$variation->getGenes}){
			$t++;
			$vp->{transcripts} = $vp->{hgenes}->{$g->{id}}->{tr};
		#	$vp->{html}->{$g->id} =  $packer->minify( \$print_html->transcripts());
			delete  $vp->{transcripts} ;
		}
		#$vp->{h} = compress($vp->{h});
		#die();
		#my $hvariant =  update_variant_editor::construct_hash_variant_global ( $project, $variation,undef,1);
	#	delete $hvariant->{html};
		#warn Dumper $hvariant ;
		#die();
		my $hh;
		$hpolyviewer->{$variation->rocksdb_id} = $vp;
		$hvariant->{$variation->rocksdb_id} = $variation;
			
		delete $variation->{array_dejavu};
		delete $variation->{references_object};
		delete $variation->{dejaVuInfosForDiag2};
		delete $variation->{annex};
		delete $variation->{buffer};
		delete $variation->{project};
		
		
		#warn Dumper $variation->annex();
		#die();
		
	}
	
#	warn '   - check junctions';
	my $js = $reference->getJunctions;
	
	my $h_ri_aval_amont;
	foreach my $junction ( @{$js } ) {
		next if (not $junction->type =~ /junction/ );
		my $debug ;
		$nb ++;
#		warn $nb."/".scalar @{$js } if $nb %30000 == 0; 
		$junction->id();
		$junction->name();
		$junction->annex();
		$junction->setPatients();
		#$junction->get_hash_exons_introns();
		
		my $vp =  PolyviewerVariant->new();
		$vp->setLmdbJunction($junction);
		$vp->{hgenes} = {};
		$vp->{genes_id} = [];
		my $code =0;
		foreach my $g (@{$junction->getGenes}){
			my $h = $vp->set_gene_junction($junction,$g);
			$h->{code} = $code;
			$vp->{hgenes}->{$g->{id}} = $h;
			
			push(@{$vp->{genes_id}},$g->{id});
			$code ++;
		}
		##############
		#	next;
		##############0
		$vp->{hpatients} ={};
		$vp->{patients_id} = [];
		my $dvp;
		foreach my $pat (@{$junction->getPatients}){
			
			foreach my $p (@{$pat->getFamily()->getMembers}){
				
				next if exists $dvp->{$p->id};
				$dvp->{$p->id} ++;
				my $h = $vp->set_patient_cache($junction,$p);
				$vp->{patients_calling}->{$p->id} =$h; 
			}
		}
		bless $junction , 'GenBoJunctionCache';
		my $t = 0;
		
		$hpolyviewer->{$junction->rocksdb_id} = $vp;
#		delete $junction->{annex};
		delete $junction->{array_dejavu};
		delete $junction->{references_object};
		delete $junction->{dejaVuInfosForDiag2};
		delete $junction->{coverage_obj};
		delete $junction->{buffer};
		delete $junction->{project};
		
		$hvariant->{$junction->rocksdb_id} = $junction;
		
		warn $junction->id;
		
		



	}
#	warn "   - end objs  ". $region->{start}." ".$region->{end};
	
#	warn "   - before save  ". $region->{start}." ".$region->{end};
	$project->disconnect();
	
	
	my $no_polyviewer = $rg_polyviewer->nosql_rocks($region);
	my $no = $rg->nosql_rocks($region);
	foreach my $k (keys %$hpolyviewer){
		$no_polyviewer->put_batch($k,$hpolyviewer->{$k});
		$no->put_batch($k,$hvariant->{$k});
	}
	$no_polyviewer->close;
	$no->close;
#	warn "   - end variant  ". $region->{start}." ".$region->{end};
	return;
}


