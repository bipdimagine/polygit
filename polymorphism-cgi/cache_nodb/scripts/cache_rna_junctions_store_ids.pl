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
warn $project->lmdb_cache_dir();
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
my @patient_names = sort { $a cmp $b } map { $_->name } @{ $project->getPatients };
my $index_patients = 0;
my $hpatients;
my @categorie_patient = ( "all", "RI", "SE");
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
warn 'store 1/4: lmdb variations' if ( $project->cache_verbose() );
my $no2 = $chr->get_lmdb_variations("c");    #open lmdb database
#ok sort and read the filewarn
my $uniq;
my $hh;
my $size_variants = 0;

foreach my $hv ( @{$all} ) {
	my $var_id = $hv->{id};
	next if exists $uniq->{$var_id};
	$uniq->{$var_id} ++;
	my $junction = thaw( decompress( $hv->{obj} ) );
	my $index_lmdb = $no2->put( $var_id, $junction );
	$junction->{buffer} = $buffer;
	$junction->{project} = $project;
	$size_variants++;
	$hh->{$var_id} = $junction->{heho_string};
	foreach my $patient (@{ $junction->getPatients() }) {
		$hpatients->{$patient->name()}->{intspan}->{all}->add($index_lmdb);
		$hpatients->{$patient->name()}->{intspan}->{RI}->add($index_lmdb) if ($junction->isRI($patient));
		$hpatients->{$patient->name()}->{intspan}->{SE}->add($index_lmdb) if ($junction->isSE($patient));
	}
	unless ( exists $intspan_global_type->{ $junction->type() } ) {
		warn "\n\nERROR: doesn't exists \$intspan_global_type->{".$junction->type() . "}\n\n";
		warn Dumper keys %$intspan_global_type;
		confess;
	}
	$intspan_global_type->{ $junction->type }->add($index_lmdb);
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
warn 'store 2/4: lmdb chr_name freeze' if ( $project->cache_verbose() );
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

warn 'store 3/4: lmdb patients' if ( $project->cache_verbose() );
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

#NOISE 
warn 'store 4/4: update methods calling' if ( $project->cache_verbose() );
my $buffer_cache = new GBuffer;
$buffer_cache->vmtouch(1);
my $project_cache = $buffer_cache->newProjectCache( -name => $project_name );
$project_cache->getPatients();
my $chr_cache = $project_cache->getChromosome($chr_name);
my $no5 = $chr_cache->get_lmdb_variations("w");
my $vector_junctions = $chr_cache->getJunctionsVector();
foreach my $junction (@{$chr_cache->getListVarObjects($vector_junctions)}) {
	foreach my $patient_cache (@{$junction->getPatients()}) {
		$junction->get_hash_noise($patient_cache);
	}
	my $jid = $junction->id();
	delete $junction->{buffer};
	delete $junction->{project};
	$no5->put( $jid, $junction );
}
$no5->close();

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
	my $vs = $project->getJunctions();
	
	my $h_ri_aval_amont;
	#search common ri_aval ri_amont
	foreach my $junction ( @{$vs } ) {
		my $h_exons_introns;
		foreach my $patient (@{ $junction->getPatients() }) {
			next if (not $junction->is_ri_aval($patient) and not $junction->is_ri_amont($patient));
			$h_exons_introns = $junction->get_hash_exons_introns() unless ($h_exons_introns);
			my $type_ri;
			$type_ri = 'ri_aval'  if ($junction->is_ri_aval($patient));
			$type_ri = 'ri_amont' if ($junction->is_ri_amont($patient));
			foreach my $tid (sort keys %{$h_exons_introns}) {
				my @lPos = (sort keys %{$h_exons_introns->{$tid}->{by_pos}});
				my $first_exon_intron = $h_exons_introns->{$tid}->{by_pos}->{$lPos[0]};
				my $last_exon_intron = $h_exons_introns->{$tid}->{by_pos}->{$lPos[-1]};
				$h_ri_aval_amont->{$patient->name()}->{$tid}->{$first_exon_intron}->{$type_ri} = $junction->id(); 
				$h_ri_aval_amont->{$patient->name()}->{$tid}->{$last_exon_intron}->{$type_ri} = $junction->id(); 
			} 
		}
	}
	
	foreach my $junction ( @{$vs } ) {
		next if ($junction->getChromosome->id() ne $chr->id());
		$junction->id();
		$junction->name();
		$junction->annex();
		$junction->setPatients();
		$junction->get_hash_exons_introns();
	}
	
	
	foreach my $junction ( @{$vs } ) {
		my $debug ;
		next if ($junction->getChromosome->id() ne $chr->id());
		foreach my $patient (@{ $junction->getPatients() }) {
			#check linked junctions - ri aval amont -
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
		my $array_patients;
		my $aho = [];
		my $ap  = [];
		
		# line to prepare dejavu global;
		my $ref = ref($junction);
		if ($ref eq 'GenBoJunction'){
			bless $junction , 'GenBoJunctionCache';
		}
		delete $junction->{buffer};
		delete $junction->{project};
		
		$hv->{obj}   = compress( freeze($junction) );
		$hv->{start} = $junction->start;
		$hv->{end}   = $junction->end;
		$hv->{id}    = $junction->id;
		push(@$ids, $junction->id);
		push( @all, $hv );
	}
	my @sort = sort { $a->{start} <=> $b->{start} or $a->{end} <=> $b->{end} } @all;
	return (\@sort),$ids;
}