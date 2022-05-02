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
warn Dumper @z; 

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }

my $nbErrors = 0;
my $buffer = new GBuffer;
$buffer->vmtouch(1);
warn $buffer->config->{'public_data_annotation'}->{root};
my @pbd = ("/data-isilon/public-data","/data-isilon/public-data","/data-beegfs/public-data_nfs");#,"/data-beegfs/public-data_nfs");
@pbd = ("/data-isilon/public-data","/data-isilon/public-data");
$buffer->config->{'public_data_annotation'}->{root} = $pbd[rand @pbd];
warn $buffer->config->{'public_data_annotation'}->{root};
#my $color = $colors[ rand @colors ];
my $project = $buffer->newProject( -name => $project_name );

$project->{gencode_version} ='34';



if ($annot_version) {
	$project->changeAnnotationVersion($annot_version, 1);
}
warn 'Project annotation '.$project->annotation_version();
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
		warn "@@@@@ ".$hRes->{end}." ".$hRes->{ttime};
		delete $process->{$hRes->{end}};
		warn "freeze ".$hRes->{freeze};
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
                        DIR => "/tmp",#$buffer->config->{project_pipeline}->{tmp},
                        SUFFIX => '.freeze');
                        
    $region->{freeze}= $region->{freeze_tmp}->filename;                
	$region->{dir_freeze}= $buffer->config->{project_pipeline}->{tmp};
	$region->{full_name} = $project->name.".".$region->{chromosome}.".".$region->{start}."-".$region->{end};
	$hregion->{$region->{chromosome}.".".$region->{start}."-".$region->{end}} = $region;
	
}
#warn Dumper $regions if scalar(@$regions) == 1;
#die() if scalar(@$regions) == 1;

#$project->buffer->disconnect();
#$project->buffer->{dbh} ="-";
foreach my $region (values %$hregion) {
	warn "$cpt/".scalar(@$regions) if ( $project->cache_verbose() );
	$cpt++;
	unless ( exists $region->{chromosome} ) {
		$region->{chromosome} = $chr->id();
	}
	$process->{$region->{full_name}} ++;
	 @pbd = shuffle @pbd;
	my $pid = $pm->start and next;
#	warn "------------------------";
	#system("/software/bin/vmtouch -t /data-isilon/public-data/repository/HG19/cadd/1.6/lmdb/".$chr->name.".uc "."/data-isilon/public-data/repository/HG19/gnomad-genome/2.1/lmdb//snps/".$chr->name);#/software/bin/vmtouch -t ".$no1->filename." "
	#system("/software/bin/vmtouch -t /data-isilon/public-data/repository/HG19/gnomad-genome/2.1/lmdb//snps/".$chr->name);

	#system("/software/bin/vmtouch -t /data-isilon/public-data/repository/HG19/gnomad-exome/2.1/lmdb//snps/".$chr->name);
#warn "/data-isilon/public-data/repository/HG19/gnomad-exome/2.1/lmdb//snps/".$chr->name;
#warn "------------------------";
	my $time_start = time;
	my $hres = {};
	my ($all,$ids)        = get_ids( $project_name, $region );
	my $time_end   = time;
	
	#$hres->{variants} = $all;
	
	$hres->{variants_ids} = $ids;
	#$hres->{region}   = $region;
	$hres->{ttime}    = abs( $time_start - $time_end );
	$hres->{end}    =  $region->{full_name};
	$hres->{freeze}    =  $region->{freeze};
	warn "$cpt end :".$hres->{ttime};
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

warn "end process" if ( $chr->project->cache_verbose() );

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
my $t = time;
warn 'store 1/3: lmdb variations' if ( $project->cache_verbose() );
my $no2 = $chr->get_lmdb_variations("c");    #open lmdb database
#ok sort and read the filewarn
my $hh;
my $uniq;
my $size_variants = 0;
foreach my $region (@$regions) {
	my $freeze = $region->{freeze};
#	warn $region->{full_name};
	warn Dumper $region  unless -e $freeze;
	confess($freeze) unless -e $freeze;
	my $hall = retrieve $freeze;
	foreach my $hv ( @{$hall} ) {
		my $var_id = $hv->{id};
		next if exists $uniq->{$var_id};
		$uniq->{$var_id} ++;
		my $variation = thaw( decompress( $hv->{obj} ) );
		my $index_lmdb = $no2->put( $var_id, $variation );
		$variation->{buffer} = $buffer;
		$variation->{project} = $project;
		$size_variants++;
		$hh->{$var_id} = $variation->{heho_string};
		foreach my $patient (@{ $variation->getPatients() }) {
			my $type = $variation->getSequencingGenotype($patient);
			$hpatients->{$patient->name()}->{intspan}->{all}->add($index_lmdb);
			$hpatients->{$patient->name()}->{intspan}->{$type}->add($index_lmdb);
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

warn "close";
my $t2 = time;
$no2->close();
warn "time : ".abs($t-time)." - close ".abs($t2-time);
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
	$buffer->config->{'public_data_annotation'}->{root} = $pbd[0];
	my $project = $buffer->newProject( -name => $project_name );
	
	
	
	$project->{gencode_version} = '34';
	
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
	my $vs = $reference->getStructuralVariations;
	warn "end loading vcf";
	foreach my $variation ( @{$vs } ) {
		my $debug ;
		#$debug =1 if $variation->id eq "7_6042930_A_G";
		#warn $variation->name();
		$variation->gnomad;
#		$variation->revel_score();
		$variation->cadd_score();
#		$variation->ncboost_score();
#		$variation->dbscsnv_ada();
#		$variation->dbscsnv_rf();
#		$variation->spliceAI();
#		$variation->dejaVuInfosForDiag();
#		#$variation->dejaVuInfosForDiag();
		$variation->name();
#		$variation->cosmic();
#		$variation->getGenes();
#		$variation->score_clinical_local();
#		$variation->score_clinvar();
#		$variation->text_clinvar();
#		$variation->comment_clinical_local();
#		$variation->hgmd_id();		
		$variation->annotation();
#
#		$variation->get_codon_text(1);
#		$variation->get_codon_text(-1);
#		$variation->getGenes();
#		$variation->getTranscripts();
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
		  
		my $hv;
		my $array_patients;
		my $aho = [];
		my $ap  = [];
		foreach my $pat ( @{ $variation->getPatients() } ) {
			warn $pat->name if $debug;
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
		#warn Dumper $variation->annex();
		#die();
		$hv->{obj}   = compress( freeze($variation) );
		$hv->{start} = $variation->start;
		$hv->{end}   = $variation->end;
		$hv->{id}    = $variation->id;
		push(@$ids, $variation->id);
		push( @all, $hv );
	}
	my @sort = sort { $a->{start} <=> $b->{start} or $a->{end} <=> $b->{end} } @all;
	return (\@sort),$ids;
}