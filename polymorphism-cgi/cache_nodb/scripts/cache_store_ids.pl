#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use Parallel::ForkManager;
use Storable qw(store retrieve freeze thaw dclone);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Cwd 'abs_path';
use Digest::MD5 qw(md5_hex);
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use lib "$RealBin/../../../GenBo/lib/obj-nodb/polyviewer";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON::XS;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use List::Util qw(shuffle);
require "$RealBin/Cache_Commons.pm";
use Sys::Hostname;
use PolyviewerVariant;
 my $host = hostname();

 my $himpact_sorted = {
	"high" => "4",
	"moderate" =>"3",
	"low" =>"1",
};
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
system("/software/bin/vmtouch -t /data-isilon/public-data/repository/HG19/cadd/1.6/lmdb/".$chr->name.".uc "."/data-isilon/public-data/repository/HG19/gnomad-genome/2.1/lmdb//snps/".$chr->name." "."/data-isilon/public-data/repository/HG19/gnomad-exome/2.1/lmdb//snps/".$chr->name);#/software/bin/vmtouch -t ".$no1->filename." "
system("/software/bin/vmtouch -t /data-isilon/public-data/repository/HG19/cadd/1.6/lmdb/".$chr->name.".uc");#/software/bin/vmtouch -t ".$no1->filename." "
system("/software/bin/vmtouch -t /data-isilon/public-data/repository/HG19/gnomad-genome/2.1/lmdb//snps/".$chr->name);

$project->preload_patients();
$project->buffer->disconnect();
$project->buffer->{dbh} ="-";
my $patients = $chr->project->getPatients();
my $regions;
my @lVarRegion;
my $new_pos       = 1;
my $old_pos       = 1;
my $nb_variants   = 0;
my $nb_var_region = 0;
#warn $fork;
$regions = Cache_Commons::get_regions( $chr, $fork);

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
#		warn "@@@@@ ".$hRes->{end}." ".$hRes->{ttime};
		delete $process->{$hRes->{end}};
#		warn "freeze ".$hRes->{freeze};
		confess($hRes->{freeze}) unless -e $hRes->{freeze};
		$nb_variants += scalar(@{$hRes->{'variants_ids'}});
		foreach my $hv (@{$hRes->{'variants_ids'}}) {
			$chr->{cache_hash_get_var_ids}->{$hv}++;
		}
	}
);
my $cpt = 1;
#warn Dumper $regions;
my $hregion;
foreach my $region (@$regions) {
	
	$region->{freeze_tmp}= File::Temp->new( TEMPLATE =>$project->name.".".$region->{chromosome}.".".$region->{start}."-".$region->{end}.'XXXXXX',
                        DIR => "/tmp",
                        SUFFIX => '.freeze');
                        
    $region->{freeze}= $region->{freeze_tmp}->filename;                
	$region->{dir_freeze}= $buffer->config->{project_pipeline}->{tmp};
	$region->{full_name} = $project->name.".".$region->{chromosome}.".".$region->{start}."-".$region->{end};
	$hregion->{$region->{chromosome}.".".$region->{start}."-".$region->{end}} = $region;
	
}
foreach my $region (values %$hregion) {
	warn "$cpt/".scalar(@$regions) if ( $project->cache_verbose() );
	$cpt++;
	unless ( exists $region->{chromosome} ) {
		$region->{chromosome} = $chr->id();
	}
	$process->{$region->{full_name}} ++;
	my $pid = $pm->start and next;
	
	my $time_start = time;
	my $hres = {};
	my ($all,$ids)        = get_ids( $project_name, $region );
	my $time_end   = time;
	
	
	$hres->{variants_ids} = $ids;
	$hres->{ttime}    = abs( $time_start - $time_end );
	$hres->{end}    =  $region->{full_name};
	$hres->{freeze}    =  $region->{freeze};
	my $freeze_file = $region->{freeze};
	unlink $freeze_file if -e $freeze_file;
	store($all, $freeze_file);
	delete $buffer->{lmdb_score};
	$pm->finish( 0, $hres );
}
$pm->wait_all_children();
if(keys %$process){
	warn Dumper $process;
	die();
}

#warn "end process" if ( $chr->project->cache_verbose() );

if ($nbErrors > 0) {
	confess("\n\nERRORS: $nbErrors errors found... confess...\n\n");
}
if (scalar keys %{$chr->{cache_hash_get_var_ids}} == 0) {
	my $cmd = "touch ".$project->getCacheBitVectorDir()."/lmdb_cache/".$chr->id().".empty";
	`touch $cmd`;
	my $no2 = $chr->get_lmdb_variations("c");
	$no2->create();
	$no2->close();
	warn "empty";
	store( {}, $project->lmdb_cache_dir . "/$chr_name.dv.freeze" ) ;    
	exit(0);
}


#end fork now I have a file with all variation and  json  it's time to sort this file and store it in lmdb database
#construct intspan 	for patient I will store lmdb_id in intspan
my $project = $chr->project;
my @patient_names = sort { $a cmp $b } map { $_->name } @{ $project->getPatients };
my $index_patients = 0;
my $hpatients;
my @categorie_patient = ( "all", "he", "ho" );
my $idpatients;
foreach my $pname (@patient_names) {
	my $op = $project->getPatient($pname);
	$idpatients->{$op->id} = $pname;
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
my $t = time;
warn 'store 1/3: lmdb variations' if ( $project->cache_verbose() );
my $dir_out = "/tmp/";


my $no2 = $chr->get_lmdb_variations("c",1);    #open lmdb database
my $no_annex = $chr->get_lmdb_annex("c");
#ok sort and read the filewarn
my $hh;
my $uniq;
my $size_variants = 0;
my $xtime =time;
my $tids;
foreach my $region (@$regions) {
	my $freeze = $region->{freeze};
	confess($freeze) unless -e $freeze;
	my $hall = retrieve $freeze;
	foreach my $hv ( @{$hall} ) {
		my $var_id = $hv->{id};
		next if exists $uniq->{$var_id};
		$uniq->{$var_id} ++;
	
		my $variation = thaw( decompress(  $hv->{obj} ) );
		$variation->sequence();
		my $toto = dclone($variation);
		
		push(@$tids,$var_id);
		my $annex =  delete $toto->{annex};
		
		$toto->{annex} = {};
		delete $toto->{annex};
		
		my $index_lmdb = $no2->put( $var_id, $toto );
		########
		#
		#########
		$variation->{buffer} = $buffer;
		$variation->{project} = $project;
		
		$no_annex->put($index_lmdb,$annex);
		$size_variants++;
		$hh->{$var_id} = $variation->{heho_string};
	 
		foreach my $pn ( keys %{$variation->patients_object() } ) {
			my $pnn = $idpatients->{$pn};
			my $vtype = $variation->sequencing_infos->{$pn}->{max}->[2];
			my $type = "he";
			$type = "ho" if $vtype eq "ho";
			die() unless $vtype;
			$hpatients->{$pnn}->{intspan}->{all}->add($index_lmdb);
			$hpatients->{$pnn}->{intspan}->{$type}->add($index_lmdb);
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
#
my $ff = $no_annex->filename();
$no_annex->close();
warn $ff;
my $zstd = $buffer->software("zstdmt");
system("$zstd $ff --rm -f --ultra ");
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
warn "time : ".abs($t-time);
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
	my $ids = [];
	my $buffer = new GBuffer;
	 $buffer->vmtouch(1);
	my $project = $buffer->newProject( -name => $project_name );
	$project->preload_patients();
	$project->buffer->disconnect();
	#$project->buffer->{dbh} ="-";
	my $chr = $project->getChromosome( $region->{chromosome} );
	my $reference = $chr->getReferences( $region->{start}, $region->{end} )->[0];
	my @all;
	my @patient_names = sort { $a cmp $b } map { $_->name } @{ $project->getPatients };
	my $hpatients;
	for ( my $i = 0 ; $i < @patient_names ; $i++ ) {
		$hpatients->{ $patient_names[$i] } = $i;
	}
	my $vs = $reference->getSVs();
	foreach my $variation ( @{$vs } ) {
		my $debug ;
	#	next unless $variation->id eq "Y_13484007_C_A";
#	next if $variation->start < 42628400;
#	next if $variation->start > 42628409;
		$debug =1 if $variation->name eq "5-42628404-cnv-del-2690";
		#warn $variation->name();
#		$variation->gnomad;
		$variation->getGnomadAC_Male;
		$variation->getGnomadAN;
		$variation->getGnomadAC;
		$variation->getGnomadHO;
		$variation->frequency;
		$variation->min_pop_freq;
		$variation->max_pop_freq;
#		delete $variation->{gnomad};
		$variation->revel_score();
		$variation->cadd_score();
		$variation->ncboost_score();
		$variation->dbscsnv_ada();
		$variation->dbscsnv_rf();
		$variation->spliceAI();
#		$variation->dejaVuInfosForDiag();
#		#$variation->dejaVuInfosForDiag();
		$variation->name();
#		$variation->cosmic();
#		$variation->getGenes();
#		$variation->score_clinical_local();
#		$variation->score_clinvar();
#		$variation->text_clinvar();
#		$variation->comment_clinical_local();
		$variation->hgmd_id();		
		$variation->annotation();
#
#		$variation->get_codon_text(1);
#		$variation->get_codon_text(-1);
		$variation->getGenes();
		$variation->getTranscripts();
		$variation->getNGSScore();
		$variation->isLargeDeletion();
		$variation->isLargeDuplication();
	
	 ####################
	#DEJAVU
	#####################
	$variation->{value}->{other_project} = $variation->other_projects();
	$variation->{value}->{other_patients} = $variation->other_patients();
	$variation->{value}->{other_patients_ho} = $variation->other_patients_ho();
	$variation->{value}->{similar_projects} = $variation->similar_projects();
	$variation->{value}->{similar_patients} = $variation->similar_patients();
	$variation->{value}->{similar_patients_ho} = $variation->similar_patients_ho();
	$variation->{value}->{this_run_patients} =  $variation->in_this_run_patients();
	#$hvariation->{value}->{this_run_patients} = $v->in_this_run_patients()."/".scalar(@{$project->getPatients});
	
	
		$variation->isCnv();
		
#
#		  	foreach my $tr ( @{$variation->getTranscripts}){
#		  		 $variation->getNomenclature($tr);
#		  	 	foreach my $p (@{$tr->getProteins}){
#		  	 		$variation->polyphenScore($p);
#		  	 		$variation->siftScore($p);
#		  	 	}
#
#		  	}
#		  
my $short_evt =1;

  $short_evt = undef if $variation->isLargeDeletion or $variation->isLargeDuplication;
	foreach my $gene (@{$variation->getGenes}){
		my $max_value;
		my $max_cat = '-';
		$variation->{annotation_polyviewer}->{$gene->id}->{spliceAI_score} = $variation->max_spliceAI_score($gene);
		$variation->{annotation_polyviewer}->{$gene->id}->{spliceAI_cat} = $variation->max_spliceAI_categorie($gene);
		my @t;
		foreach my $tr1 ( sort { ($himpact_sorted->{$variation->effectImpact($b)} <=>  $himpact_sorted->{$variation->effectImpact($a)}) or ($a->appris_level <=> $b->appris_level) }  @{$gene->getTranscripts}){
			my $htr ={};
			$htr->{enst} = $tr1->id();	
			$htr->{nm} = $tr1->external_name;
			$htr->{ccds} = $tr1->ccds_name;	
			$htr->{appris} = $tr1->appris_type;
			$htr->{consequence} = $variation->variationTypeInterface($tr1);
			$htr->{impact_score_text} = $variation->effectImpact($tr1);
			my @coding_infos = ("sift","polyphen","prot","codons","codons_AA");
		#	if ($v->isLargeDeletion or $v->isLargeDuplication){
		#	value_html_badge($htr,"prot","-");
		#	value_html_badge($htr,"codons","-");
		#	value_html_badge($htr,"codons_AA","-");
		#	}
		
			if ($variation->isCoding($tr1) && $short_evt){
				my $prot = $tr1->getProtein();
				$htr->{prot_pos} = $variation->getProteinPosition($prot);
				$htr->{codons} = $variation->getCodons($tr1);
				$htr->{codons_AA} = $variation->protein_nomenclature($prot);
				$htr->{sift} = $variation->siftScore($tr1->getProtein);
				$htr->{polyphen} = $variation->polyphenScore($tr1->getProtein);
			}
			my $te = $tr1->findExonNumber($variation->start, $variation->end);
			$htr->{exon} = $te;
		 	if ($te == -1){
		 		$htr->{exon} = $tr1->findNearestExon($variation->start, $variation->end);
		 	}
		 	if ($short_evt){
				$htr->{nomenclature} =  $variation->getNomenclature($tr1);
		 	}
		 	
		 	$htr->{impact_score} = $himpact_sorted->{$variation->effectImpact($tr1)};
		 	
		 	$htr->{main}  = 1 if $tr1->isMain();
			push(@t,$htr);
		}
		$variation->{annotation_polyviewer}->{$gene->id}->{trs} = \@t;
	}
		  
		  
		my $hv;
		my $array_patients;
		my $aho = [];
		my $ap  = [];
		my $ddebug;
		 if ($variation->md5_id eq "ad441b5a0700ab160d1016c410d160b8"){
	#$ddebug =1;
	 }
	  if ($variation->id eq "Y_59019735_G_A"){
	  		$ddebug =1;
	  }
	  $variation->sequencing_infos();
	  $variation->dp_infos();
	  $variation->split_read_infos();
		foreach my $pat ( @{ $variation->getPatients() } ) {
			$variation->isFoundBySVCaller($pat);

			my $pn = $pat->name();
			my $hp;
			$hp->{name} = $pn;
			my $patient_id = $hpatients->{$pn};
			push( @$ap, $patient_id );
			push( @$aho, $patient_id ) if ( $variation->isHomozygote($pat) );
		}
#		die() if $ddebug;

	
		$variation->{heho_string} = join( ",", sort { $a <=> $b } @$ap );
		if ( scalar(@$aho) ) {
			$variation->{heho_string} = $variation->{heho_string}." ".join( ",", sort { $a <=> $b } @$aho ) . " HO";
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
		elsif  ($ref eq 'GenBoLargeDeletion'){
					bless $variation , 'GenBoLargeDeletionCache';
		}
#		my $pv = return_PolyViewerVariant($variation);
#			warn Dumper($pv);
#		die(); 
		delete $variation->{buffer};
		delete $variation->{project};	
		#warn Dumper $variation->annex();
		#die();
		$hv->{obj}   = compress( freeze($variation) );
		$hv->{start} = $variation->start;
		$hv->{end}   = $variation->end;
		$hv->{id}    = $variation->id;
	
		push(@$ids, $variation->id);
		#die() if $debug;
		push( @all, $hv );
	}

	my @sort = sort { $a->{start} <=> $b->{start} or $a->{end} <=> $b->{end} } @all;
	return (\@sort),$ids;
}

sub return_PolyViewerVariant {
	my ($vh) = @_;
		my $pv =   PolyviewerVariant->new();  
		#$vh->{project} = $project;
		#$vh->{buffer} = $project->buffer;
		#$pv->gene($gene);
		$pv->id($vh->id);
		$pv->start($vh->start);
		$pv->end($vh->end);
		$pv->ref_allele($vh->ref_allele);	
		$pv->allele($vh->alternate_allele);	
		$pv->gnomad_id($vh->gnomad_id);
		$pv->chromosome($vh->name);
		$pv->name($vh->name);
		$pv->type($vh->type);
		#######################
		#hgmd et clinvar
		#######################
		
		
		$pv->hgmd($vh->hgmd_class);
		$pv->hgmd_id($vh->hgmd_id);
		foreach my $gid (keys %{$vh->genes_pathogenic_DM}){
		if (exists $vh->genes_pathogenic_DM->{$gid} && $vh->genes_pathogenic_DM->{$gid}->{DM} ){
				$pv->dm($vh->isDM);
				$pv->hgmd_phenotype($vh->hgmd_phenotype);
			}
			if (exists $vh->genes_pathogenic_DM->{$gid} && $vh->genes_pathogenic_DM->{$gid}->{pathogenic} ){
				$pv->clinvar_pathogenic($vh->isClinvarPathogenic);
			}
			
		}
		
		$pv->clinvar_id($vh->clinvar_id);
		$pv->clinvar($vh->clinvar_class);
	
		
		
		################
		# Calling
		##################
		$pv->isCnv($vh->isCnv);
		
		if ($vh->isCnv){
			foreach my $patient (@{$vh->getPatients}){ 
			#	$pv->setCnvValues($vh->getChromosome,$patient,$vh);
			}
		}
		else {
			foreach my $patient (@{$vh->getPatients}){ 
				foreach my $p (@{$patient->getFamily()->getMembers}) {
					$pv->patients_calling->{ $p->id }->{gt} = $vh->getSequencingGenotype($p);
					$pv->patients_calling->{ $p->id }->{pc} = $vh->getRatio($p);
					$pv->patients_calling->{ $p->id }->{dp} = $vh->getDP($p);
					#$pv->patients_calling->{ $p->id }->{model} = $vh->getTransmissionModelType($p->getFamily(),$p);
				}
			}
		}
		

		$pv->patients_calling();
		################################
		#
		# CNV
		#
		####################@
#		$pv->isCnv($vh->isCnv);
#		if ($pv->isCnv){
#		$pv->isDude(1); ### test
#		
#		$pv->patients_calling($pv->parseTrioTableCnv($vh->{html}->{trio},$project));
#			if ($pv->isDude){
#				my $chr = $project->getChromosome($pv->chromosome());
#				my $primers = $project->getPrimersByPosition($chr,$pv->start,$pv->end);
#				my $hpatients;
#				  map {$hpatients->{$_->id} = $_} @{$project->getPatients()};
#				foreach my $primer (@$primers){
#					foreach my $ph (keys %{$pv->patients_calling}){
#						my ($patient) = $hpatients->{$ph};
#							$pv->patients_calling->{ $patient->id }->{norm_depth} = $patient->cnv_region_ratio_norm($chr->name,$primer->start,$primer->end);
#							$pv->patients_calling->{ $patient->id }->{dude_score} = $primer->cnv_score($patient);
#					}
#				}
#			}
			
			
			
		
		#################
		# gnomad
		##############@
		
		$pv->gnomad_ac($vh->getGnomadAC);
		$pv->gnomad_an($vh->getGnomadAN);
		$pv->gnomad_min_pop_name($vh->min_pop_name);
		$pv->gnomad_min_pop($vh->min_pop_freq);
		$pv->gnomad_max_pop_name($vh->max_pop_name);
		$pv->gnomad_max_pop($vh->max_pop_freq);
		$pv->gnomad_ho($vh->getGnomadHO);
		$pv->gnomad_ho_male($vh->getGnomadAC_Male);
		
	#########
	# DEJAVU
	#########
	
	$pv->dejavu_other_projects($vh->other_projects());
	$pv->dejavu_other_patients($vh->other_patients());
	$pv->dejavu_other_patients_ho($vh->other_patients_ho());
	$pv->dejavu_similar_projects( $vh->similar_projects());
	$pv->dejavu_similar_patients($vh->similar_patients());
	$pv->dejavu_similar_patients_ho($vh->similar_patients_ho());
	$pv->dejavu_this_run_patients($vh->in_this_run_patients);# = '-';
	$pv->text_caller([]);
	
	return $pv;
	
	
	
}


