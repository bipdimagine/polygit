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
use lib "$RealBin/../";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime); 
use JSON::XS;
use Compress::Snappy;
use Getopt::Long;
use Carp; 
use GBuffer;
use Cache_Commons;
use Sys::Hostname;
use GenBoNoSqlRocksVector;
 my $host = hostname();


warn "*_*_*_*_*_".$host."*_*_*_*_*_";


my $fork = 1;
my ($project_name, $chr_name, $annot_version);
my $ok_file;
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'annot_version=s'    => \$annot_version,
	'chr=s'        => \$chr_name,
	'file=s' => \$ok_file,
);

 if ($ok_file && -e $ok_file) {
 	system("rm $ok_file");
 }


unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }


my $buffer  = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProject( -name => $project_name );



warn "\n### CACHE: store annotations step\n";
my $nbErrors = 0;

if ($annot_version) {
	$project->changeAnnotationVersion($annot_version);
}
my $chr     = $project->getChromosome($chr_name);
my $no      = $chr->get_rocks_variations("r");

#if ($no->size == 0 or -e $project->getCacheBitVectorDir()."/lmdb_cache/".$chr->id().".empty") {
#	my $no = $chr->get_lmdb_categories("c");
#	$no->create();
#	warn $project->getCacheBitVectorDir()."/lmdb_cache/".$chr->id();
#	my $cmd = "touch ".$project->getCacheBitVectorDir()."/lmdb_cache/".$chr->id()."/genes";
#	`$cmd`;
#	my $cmd2 = "touch ".$project->getCacheBitVectorDir()."/lmdb_cache/".$chr->id()."/genes_index";
#	`$cmd2`;
#	warn "empty ". $project->getCacheBitVectorDir()."/lmdb_cache/".$chr->id().".empty";
#	warn "here";
#	warn $no->size();
#	$no->close();
#	my $no5 = $chr->get_lmdb_patients("c");
#	warn $no5->dir;
#	$no5->create();
#	my $h;
#	my $vnull =  Bit::Vector->new(0);
#	$h->{bitvector}->{'ho'} = $vnull->Shadow();
#	$h->{bitvector}->{'he'} = $vnull->Shadow();
#	$h->{bitvector}->{'all'} = $vnull->Shadow();
#	foreach my $p (@{$project->getPatients}){
#		$no5->put($p->name,$h);
#	}
#	$no5->close();
#	system("date > $ok_file") if $ok_file;
#	die();
#	exit(0);
#}
#my $no      = $chr->get_lmdb_variations("r");
my $hGenesIds;

my $ranges = $no->ranges($fork);
$no->close;
$no = undef;
my $intspan_genes_categories  = {};
my $intspan_global_categories = {};
my $categories = Cache_Commons::categories();
foreach my $c ( keys %{ $categories->{global}->{frequency} } ) {
	$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
}

foreach my $c ( keys %{ $categories->{global}->{variation_type} } ) {
	$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
}

my $scaled_score_gnomad_ho =  $project->buffer->config->{scaled_score_gnomad_ho_ac};
foreach my $c ( keys %{$scaled_score_gnomad_ho} ) {
	$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
}

my $scaled_score_gnomad =  $project->buffer->config->{scaled_score_gnomad_ac};
foreach my $c ( keys %{$scaled_score_gnomad} ) {
	$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
}
my $scaled_score_project_dejavu =  $project->buffer->config->{project_dejavu};
foreach my $c ( keys %{$scaled_score_project_dejavu} ) {
	$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
}
my $scaled_score_sample_dejavu =  $project->buffer->config->{sample_dejavu};
foreach my $c ( keys %{$scaled_score_sample_dejavu} ) {
	$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
}

my $scaled_score_sample_ho_dejavu =  $project->buffer->config->{sample_ho_dejavu};
foreach my $c ( keys %{$scaled_score_sample_ho_dejavu} ) {
	$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
}
my $functional_categorie = $project->buffer->config->{functional_annotations};
foreach my $c ( keys %{$functional_categorie} ) {
	$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
}

my $scaled_score_ncboost =  $project->buffer->config->{scaled_score_ncboost};
foreach my $c ( keys %{$scaled_score_ncboost} ) {
	$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
}
my $process;
my $pm = new Parallel::ForkManager($fork);
#for intergenic variant construct region
my $intspan_region_intergenic = Set::IntSpan::Fast::XS->new();
#and store intergenic id;
my $id_intergenic = Set::IntSpan::Fast::XS->new();
my $hintspan_patients;
### init hpatients;
my $hpatients = {};
my @categorie_patient = ( "all", "he", "ho" );
my @patient_names = sort { $a cmp $b } map { $_->name } @{ $project->getPatients };
my $index_patients = 0;
foreach my $pname (@patient_names) {
	foreach my $c (@categorie_patient) {
		$hpatients->{$pname}->{$c} = Set::IntSpan::Fast::XS->new();
	}
}

####################################
# one fork is finished i could "union"  all the intspan construct in the "sub annotation()"
#run on finish 	aggregate all result (intpsan) from child
#########################################
$pm->run_on_finish(
	sub {
		my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $hres ) = @_;
		unless ( defined($hres) or $exit_code > 0 ) {
			$nbErrors++;
			print qq|No message received from child process $exit_code $pid!\n|;
			warn Dumper $hres;
			return;
		}
		if ( exists $hres->{start} ) { print qq|1- No message received from child process $exit_code $pid!\n|;
			return;
		}
		delete $process->{ $hres->{idp} };
			
		######################
		#For Patient  Type
		######################
		foreach my $pname ( keys %{ $hres->{patients} } ) { 
			foreach my $cat ( keys %{$hres->{patients}->{$pname} } ) { 
			$hpatients->{$pname}->{$cat} = Set::IntSpan::Fast::XS->new() unless exists $hpatients->{$pname}->{$cat};
			$hpatients->{$pname}->{$cat} = $hpatients->{$pname}->{$cat}->union($hres->{patients}->{$pname}->{$cat});
			}
		}
		######################
		#For Patient  Type
		######################

		

		#################
		#For frequencies  annotation categorie (intspan on chromosome)
		#################
		$intspan_region_intergenic = $intspan_region_intergenic->union( $hres->{intergenic}->{intspan_region} );
		$id_intergenic = $id_intergenic->union( $hres->{intergenic}->{id} );
		
		foreach my $cat ( keys %{ $hres->{frequency} } ) {
			$intspan_global_categories->{$cat} = Set::IntSpan::Fast::XS->new() unless ( exists $intspan_global_categories->{$cat} );
			$intspan_global_categories->{$cat} = $intspan_global_categories->{$cat}->union( $hres->{frequency}->{$cat} );
		}
		#################
		#For genes annotation categorie intspan on gene
		#################
		foreach my $g ( keys %{ $hres->{genes} } ) {
			$hGenesIds->{$g}++;
			my $hgene = $hres->{genes}->{$g};
			unless ( exists $intspan_genes_categories->{$g} ) {
				$intspan_genes_categories->{$g} = init_genes_intspan();
				$intspan_genes_categories->{$g}->{all} = Set::IntSpan::Fast::XS->new();
			}
			foreach my $cat ( keys %{$hgene} ) {
				confess($cat) unless exists $intspan_genes_categories->{$g}->{$cat};
				$intspan_genes_categories->{$g}->{$cat} = $intspan_genes_categories->{$g}->{$cat}->union( $hgene->{$cat} );
				$intspan_genes_categories->{$g}->{all} = $intspan_genes_categories->{$g}->{all}->union( $hgene->{$cat} );
			}
		}
	}
);
#end run on finish

##############
#fork annotation
##############
 $project_name = $chr->project->name();
 $chr_name     = $chr->name();
my $true = 0;
my $idp  = 0;
while ( $true < 2 ) {   
	#check if all the region end correctly unless restart region failed
	my $cpt = 1;
	
	foreach my $r (@$ranges) {
		#warn "$cpt/".scalar(@$ranges) if ( $project->cache_verbose() );
		$cpt++;
		$process->{ $r->[0] . "-" . $r->[1] } = $r;
		my $pid;
		$pid = $pm->start and next;
		my $hres;
		warn "$cpt/".scalar(@$ranges);
		my $hres = get_annotations( $project_name, $chr_name, $r, $annot_version );
		$hres->{start} = $r->[0] . "-" . $r->[1];
		$hres->{idp} = $r->[0] . "-" . $r->[1];
		delete $hres->{start};
		$hres->{done} = 1;
		$pm->finish( 0, $hres );
	}    #end for range range
	$pm->wait_all_children();
	warn "end process" if ( $chr->project->cache_verbose() );
	
	if ($nbErrors > 0) {
		confess("\n\nERRORS: $nbErrors errors found... confess...\n\n");
	}
	
	#all region are ok
	unless ( keys %$process ) { last; }

	#at least one region failed restart only this region
	$true++;
	@$ranges = values %$process;
	#warn "second run";
}
confess("problem store process") if ( keys %$process );

my $no = $project->getChromosome($chr_name)->get_rocks_variations("r");
my $size_variants = $no->size();
$no->close();
warn $size_variants;
sleep(5);
#########################
#GLOBAL CATEGORIES
# ok now all part is finished I can store the global intspan and construct bitvector for each global categories
##########################
warn 'store 1/3: complete lmdb categories';

my $rocks3 = $chr->rocks_vector("c");
my $dir_vector = $project->rocks_directory("vector");

#my $no3 = $chr->get_lmdb_categories("c");    #open lmdb database

$rocks3->size($size_variants);
if($size_variants == 0 ){
	my $vnull =  Bit::Vector->new(0);
	foreach my $patient (@{$project->getPatients}){	
		$rocks3->put_vector_patient_batch($patient,"all",$vnull);
		$rocks3->put_vector_patient_batch($patient,"ho",$vnull);
		$rocks3->put_vector_patient_batch($patient,"he",$vnull);
	}

	$rocks3->write_batch();
	$rocks3->close();
	warn "no variants !!!!";
		system("date > $ok_file") if $ok_file;
	exit(0);
} 

foreach my $k ( keys %{$intspan_global_categories} ) {
	my $h;
	$h->{name}    = $k;
	$h->{intspan} = $intspan_global_categories->{$k};
	my $bitv = Bit::Vector->new_Enum( $size_variants, join( ',', $intspan_global_categories->{$k}->as_array ) );
	$h->{bitvector} = $bitv;
	$rocks3->put_batch_vector_chromosome($k,$bitv);
}
my $v1 = Bit::Vector->new($size_variants);
$v1->Fill;
$rocks3->put_batch("all",$v1);
$rocks3->write_batch();
#$no3->close();

#########################
#Patients  CATEGORIES
# ok now all part is finished saved patients intspan , all , he,  ho
##########################
warn 'store 2/3: lmdb patients';
#my $no5 = $chr->get_lmdb_patients("c"); #open lmdb database cretae new one or erase older

foreach my $patient (@{$project->getPatients}){	
	my $pn = $patient->name;
	my $h;
	$h->{name} = $pn;
	
	foreach my $type ( keys %{ $hpatients->{$pn} } ) {
		my $bitv = Bit::Vector->new_Enum( $size_variants, join( ',', $hpatients->{$pn}->{$type}->as_array ) );
		$rocks3->put_vector_patient_batch($patient,$type,$bitv);
	}
}
#$no5->close();
#$no5 = undef;
$rocks3->write_batch();

#################
# GENES
#################
warn 'store 3/3: lmdb genes';
my @intergenic;
my $tree = Set::IntervalTree->new;
my $iter = $intspan_region_intergenic->iterate_runs();
my $i = 0;
delete $chr->{lmdb};
delete $chr->{rocks};
my $no3 = $chr->get_rocks_variations("r");
while ( my ( $from, $to ) = $iter->() ) {
	my $id = 'intergenic_' . $chr_name . '_' . $from . '_' . $to;
	$tree->insert( $id, $from, $to + 1 );
}
my @array = $id_intergenic->as_array;
my $t = time;
foreach my $lmdb_id (@array) {
	my $v = $no3->get_index($lmdb_id);
	my $results = $tree->fetch( $v->{start}, $v->{end} + 1 );
	confess() unless @$results;
	confess() if scalar(@$results) > 1;
	my $interid = $results->[0];
	unless ( exists $intspan_genes_categories->{$interid} ) {
		$intspan_genes_categories->{$interid} = init_genes_intspan();
		$intspan_genes_categories->{$interid}->{all} = Set::IntSpan::Fast::XS->new();
		$intspan_genes_categories->{$interid}->{intergenic} = Set::IntSpan::Fast::XS->new();
	}
	$intspan_genes_categories->{$interid}->{all}->add($lmdb_id);
	$intspan_genes_categories->{$interid}->{intergenic}->add($lmdb_id);
}
construct_bitVector_for_gene( $rocks3,$intspan_genes_categories, $size_variants);
$rocks3->write_batch();
$rocks3->close();
warn 'store 4/3: lmdb fammilly';
$buffer = new GBuffer;
 $project = $buffer->newProjectCache( -name => $project_name );
 $chr = $project->getChromosome($chr_name);
 my $rocks4 = $chr->rocks_vector("w");
foreach my $family (@{$project->getFamilies}){
	foreach my $children  (@{$family->getChildren}){
		my $bitv1  = $family->getVector_individual_denovo($chr,$children,1); 
		$rocks4->put_batch_vector_transmission($children,"ind_denovo",$bitv1);
		
		my $bitv2  = $family->getVector_individual_dominant($chr,$children,1);
		$rocks4->put_batch_vector_transmission($children,"ind_dominant",$bitv2);
		
		my $bitv3  = $family->getVector_individual_recessive($chr,$children,1);
		$rocks4->put_batch_vector_transmission($children,"ind_recessive",$bitv3);
			
		my $bitv4 = $family->getVector_individual_father($chr,$children,1);
		$rocks4->put_batch_vector_transmission($children,"ind_father",$bitv4);
		
		my $bitv5 = $family->getVector_individual_mother($chr,$children,1);
		$rocks4->put_batch_vector_transmission($children,"ind_mother",$bitv5);

		my $bitv6 = $family->getVectorBothTransmission($chr,$children,1);
		$rocks4->put_batch_vector_transmission($children,"ind_both",$bitv6);

		my $bitv7 = $family->getVector_individual_uniparental_disomy($chr,$children,1);
		$rocks4->put_batch_vector_transmission($children,"ind_uniparental",$bitv7);
		
	}
	
}
$rocks4->write_batch();
$rocks4->close();
	system("date > $ok_file") if $ok_file;
exit(0);


exit(0);

sub get_annotations {
	my ( $project_name, $chr_name, $region, $annot_version ) = @_;
	my $hres;
	my $buffer = new GBuffer;
	$buffer->vmtouch(1);
	my $project = $buffer->newProject( -name => $project_name );
	$project->preload_patients();
	$project->buffer->disconnect();
	if ($annot_version) {
		$project->changeAnnotationVersion($annot_version);
	}
	#variation type 
	 $hres->{variation_type} = {};
	
	
	# prepare intspan for global categories
	my $intspan_global_categories;
	foreach my $c ( keys %{ $categories->{global}->{frequency} } ) {
		$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
	}
	foreach my $c ( keys %{ $categories->{global}->{variation_type} } ) {
		$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
	}
	my $intspan_genes_categories = {};
	#initialisattion categorie patients
	my $index_patients = 0;
	my $hpatients;
	#for intergenic variant construct region
	my $intspan_region_intergenic = Set::IntSpan::Fast::XS->new();
	#and store intergenic id;
	my $id_intergenic = Set::IntSpan::Fast::XS->new();

	my $no = $project->getChromosome($chr_name)->get_rocks_variations("r");
#	my $cursor = $no->cursor( $region->[0], $region->[1] );
	my $nb     = 0;
	my $l      = abs( $region->[0] - $region->[1] );
	my $cpt    = $region->[0];
	foreach (my $i= $region->[0];$i<=$region->[1];$i++){
	#while ( my $var_id = $cursor->next_key ) {
		my $variation;
		my $lmdb_index = $i;
		$variation = $no->get_index($i);
		
		if ( $lmdb_index < 0 ) {
			confess("index out of range <0");
		}
		#warn "************************************\n.".$var_id." ".Dumper($variation)."*--------------------*\n********" if $lmdb_index == -2;
		unless ($variation) {
			confess();
		}
		confess() unless $variation->id;
		$variation->{buffer}  = $buffer;
		$variation->{project} = $project;

		################
		# patient  Ratio
		################
		my $debug;
		$debug =1 if $variation->id eq "1_178116993_T_C";
		
		foreach my $p (@{$variation->getPatients}){
			my $type = $variation->getSequencingGenotype($p);
			warn $p->name." ".$type." ".$lmdb_index if ($debug);
			$hpatients->{$p->name}->{all} = Set::IntSpan::Fast::XS->new() unless exists $hpatients->{$p->name}->{all};
			$hpatients->{$p->name}->{$type} =Set::IntSpan::Fast::XS->new() unless $hpatients->{$p->name}->{$type};
			$hpatients->{$p->name}->{all}->add($lmdb_index);
			$hpatients->{$p->name}->{$type}->add($lmdb_index);
			my $r = $variation->getRatio($p);
			for (my $i=10;$i<=80;$i+=10) {
				$hpatients->{$p->name}->{"ratio_".$i} = Set::IntSpan::Fast::XS->new() unless exists $hpatients->{$p->name}->{"ratio_".$i};
				if ($r>=$i) {
						
						$hpatients->{$p->name}->{"ratio_".$i}->add($lmdb_index);
				}
			}
			for (my $i=1;$i<=5;$i++){
				$hpatients->{$p->name}->{"ratio_".$i} = Set::IntSpan::Fast::XS->new() unless exists $hpatients->{$p->name}->{"ratio_".$i};
				if ($r>=$i) {
						$hpatients->{$p->name}->{"ratio_".$i}->add($lmdb_index);
				}
			}
			for (my $i=10;$i<=80;$i+=10) {
				$hpatients->{$p->name}->{"lower_ratio_".$i} = Set::IntSpan::Fast::XS->new() unless exists $hpatients->{$p->name}->{"lower_ratio_".$i};
				if ($r<=$i) {
						$hpatients->{$p->name}->{"lower_ratio_".$i}->add($lmdb_index);
				}
			}
			for (my $i=1;$i<=5;$i++){
					$hpatients->{$p->name}->{"lower_ratio_".$i} = Set::IntSpan::Fast::XS->new() unless exists $hpatients->{$p->name}->{"lower_ratio_".$i};
				if ($r<=$i) {
						$hpatients->{$p->name}->{"lower_ratio_".$i}->add($lmdb_index);
				}
			}
		}
		
		##############
		# GNOMAD
		##############
		
		my $nb_ho = $variation->getGnomadHO();
	 	my $nb_ac = $variation->getGnomadAC();
	 	$nb_ac = 0 unless $nb_ac;
	 	$nb_ho = 0 unless defined $nb_ho;
	 	foreach my $c (keys %$scaled_score_gnomad_ho){
	 			$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new() unless exists $intspan_global_categories->{$c};
	 		if ($nb_ho <= $scaled_score_gnomad_ho->{$c} ){
	 			$intspan_global_categories->{$c}->add($lmdb_index);
	 		}
	 	}
	 	
	 	foreach my $c (keys %$scaled_score_gnomad){
	 			$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new() unless exists $intspan_global_categories->{$c};
	 		if ($nb_ho <= $scaled_score_gnomad->{$c} ){
	 			$intspan_global_categories->{$c}->add($lmdb_index);
	 		}
	 	}
	 	
	 	##############
		# deja_vu
		##############
		
			my $pdv = $variation->other_projects;
	 		my $sdv = $variation->other_patients ; 
	 		my $sdv_ho = $variation->other_patients_ho; 
	 		
	 		
	 		
	 			
	 		foreach my $c (sort {$a <=>$b} keys %$scaled_score_project_dejavu){
	 			$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new() unless exists $intspan_global_categories->{$c};
	 				if ($pdv <= $scaled_score_project_dejavu->{$c} ){
	 						$intspan_global_categories->{$c}->add($lmdb_index);
	 					#$vector_categories_chromosomes->{$c}->Bit_On($v->{vector_id});
	 				}

	 				
	 			}
	 			
	 			
	 			foreach my $c (sort {$a <=>$b}keys %$scaled_score_sample_dejavu){
	 			$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new() unless exists $intspan_global_categories->{$c};
	 				if ($sdv <= $scaled_score_sample_dejavu->{$c} ){
	 					$intspan_global_categories->{$c}->add($lmdb_index);
	 				}

	 				
	 			}
	 			
	 			foreach my $c (sort {$a <=>$b} keys %$scaled_score_sample_ho_dejavu){
	 				$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new() unless exists $intspan_global_categories->{$c};
	 				if ($sdv_ho <= $scaled_score_sample_ho_dejavu->{$c} ){
	 					$intspan_global_categories->{$c}->add($lmdb_index);
	 				}

	 			}	
	 		
	 	

		 	my $ncboost = $variation->ncboost_score();
	 			
	 	##############
		# ncboost
		##############			
	 			
	 	if ($ncboost and $ncboost ne '-') {
	 		$ncboost = $ncboost * 100;
		 	foreach my $c (keys %$scaled_score_ncboost){
		 		
		 		my $value_cat = $scaled_score_ncboost->{$c} * 100;
		 		$intspan_global_categories->{'ncboost_'.$c} = Set::IntSpan::Fast::XS->new() unless exists $intspan_global_categories->{'ncboost_'.$c};
		 		if ($ncboost >= $value_cat){
		 			$intspan_global_categories->{'ncboost_'.$c}->add($lmdb_index);
		 		}
		 	}
	 	}	

	 	##############
		# pathogenic
		##############	
		my $c = "dm";
		$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new() unless exists $intspan_global_categories->{$c};
		
	 		if ($variation->isDM){
	 		 		$intspan_global_categories->{$c}->add($lmdb_index);
	 		 }
	 		  $c = "hgmd";
	 		  $intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new() unless exists $intspan_global_categories->{$c};
	 		
	 		  if ($variation->hgmd_id){
	 		  		$intspan_global_categories->{$c}->add($lmdb_index);
	 		 }
	 		  $c = "clinvar_pathogenic";
	 		   $intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new() unless exists $intspan_global_categories->{$c};
	 		
	 		 if ($variation->score_clinvar == 5) {
	 		 	$intspan_global_categories->{$c}->add($lmdb_index);
	 		 }
	 		  $c = "clinvar";
	 		   $intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new() unless exists $intspan_global_categories->{$c};
	 		
		  if ($variation->clinvar_id){
	 		  		$intspan_global_categories->{$c}->add($lmdb_index);
	 		 }

	 	#####################
		# annotation_mask
		#####################
	 		 my $mask = $variation->{annotation}->{all}->{mask} ;
	 		 foreach my $c ( keys %{$project->buffer->config->{functional_annotations}}){
	 		 		$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new() unless exists $intspan_global_categories->{$c};
	 		 	 if ($mask & $project->getMaskCoding($c)){
	 		 		$intspan_global_categories->{$c}->add($lmdb_index);
	 		 	}
	 		 }
	 		 
	 		
	 		 if (scalar(@{$variation->getGenes}) >0 ) {
	 		 	my $c = "genes";
	 		 	$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new() unless exists $intspan_global_categories->{$c};
	 		 	$intspan_global_categories->{$c}->add($lmdb_index);
	 		 	
	 		 	  foreach my $g (@{$variation->getGenes}) {
	 		 	  	my $c1 = "acmg";
	 		 	  	$intspan_global_categories->{$c1} = Set::IntSpan::Fast::XS->new() unless exists $intspan_global_categories->{$c1};
	 		 		next unless exists $g->panels_name->{'ACMG-Actionable'};
	 		 		$intspan_global_categories->{$c1}->add($lmdb_index);
	 		 	}
	 		 	$project->buffer->disconnect();
	 		 }
	 		 else {
	 		 	my $c = "intergenic";
	 		 	$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new() unless exists $intspan_global_categories->{$c};
	 		 	$intspan_global_categories->{$c}->add($lmdb_index);
	 		 }
	 		 
		
	
		################
		#	get freqency categories and construct intspan for each
		###############
		my $cat_freq = $variation->categorie_frequency();
		
		$intspan_global_categories->{$cat_freq}->add($lmdb_index);
		my $cat_freq_ho = $variation->categorie_frequency_ho();
		$intspan_global_categories->{$cat_freq_ho}->add($lmdb_index);
		if ( $variation->isClinical() ) {
			$intspan_global_categories->{'pheno_snp'}->add($lmdb_index);
		}

		################
		#	get sub/ins/del/large_del categories and construct intspan for each
		###############
		if ( $variation->isVariation() ) {
			$intspan_global_categories->{'substitution'}->add($lmdb_index);
		}
		elsif ( $variation->isInsertion() ) {
			$intspan_global_categories->{'insertion'}->add($lmdb_index);
		}
		elsif ( $variation->isDeletion() ) {
			$intspan_global_categories->{'deletion'}->add($lmdb_index);
		}
		elsif ( $variation->isLargeDeletion() ) {
			$intspan_global_categories->{'large_deletion'}->add($lmdb_index);
		}
		elsif ( $variation->isLargeDuplication() ) {
			$intspan_global_categories->{'large_duplication'}->add($lmdb_index);
		}
		elsif ( $variation->isLargeInsertion() ) {
			$intspan_global_categories->{'large_insertion'}->add($lmdb_index);
		}
		else {
			confess;
		}

		################
		#	get strong/medium/low confidence categories and construct intspan for each
		###############
		my $ngs_score = $variation->getNGSScore();
		if ( $ngs_score == 2 ) {
			$intspan_global_categories->{'ngs_score2'}->add($lmdb_index);
		}
		elsif ( $ngs_score == 1 ) {
			$intspan_global_categories->{'ngs_score1'}->add($lmdb_index);
		}
		elsif ( $ngs_score == 0 ) {
			$intspan_global_categories->{'ngs_score0'}->add($lmdb_index);
		}
		
		my $cadd_score = $variation->cadd_score();
		if ($cadd_score == -1)    { $intspan_global_categories->{cadd_not}->add($lmdb_index);  }
		elsif ($cadd_score >= 40) { $intspan_global_categories->{cadd_40}->add($lmdb_index); }
		elsif ($cadd_score >= 35) { $intspan_global_categories->{cadd_35}->add($lmdb_index); }
		elsif ($cadd_score >= 30) { $intspan_global_categories->{cadd_30}->add($lmdb_index); }
		elsif ($cadd_score >= 25) { $intspan_global_categories->{cadd_25}->add($lmdb_index); }
		elsif ($cadd_score >= 20) { $intspan_global_categories->{cadd_20}->add($lmdb_index); }
		elsif ($cadd_score >= 15) { $intspan_global_categories->{cadd_15}->add($lmdb_index); }
		elsif ($cadd_score >= 10) { $intspan_global_categories->{cadd_10}->add($lmdb_index); }
		elsif ($cadd_score >= 5)  { $intspan_global_categories->{cadd_5}->add($lmdb_index);  }
		else					  { $intspan_global_categories->{cadd_0}->add($lmdb_index);  }
 
 
 
 
 
 
		################
		#	get Genes annotations prediction and functional
		###############
		my $genes = $variation->getGenes();
		my $debug;
		$debug =1 if  $variation->name eq "rs747758472";
		foreach my $g (@$genes) {
			unless ( exists $intspan_genes_categories->{ $g->id } ) {
				$intspan_genes_categories->{ $g->id } = init_genes_intspan();
			}
			
			my $h_spliceAI = $variation->spliceAI_score($g);
			if ($h_spliceAI) {
				foreach my $cat (keys %$h_spliceAI) {
					my $value = $h_spliceAI->{$cat};
					next if ($value eq '-');
					if (defined($value) and $value >= $buffer->config->{spliceAI}->{high}) {
						$intspan_global_categories->{spliceAI_high}->add($lmdb_index);
						$intspan_global_categories->{spliceAI_medium}->add($lmdb_index);
						$intspan_global_categories->{predicted_splice_site}->add($lmdb_index);
						$intspan_genes_categories->{ $g->id }->{spliceAI_high}->add($lmdb_index);
						$intspan_genes_categories->{ $g->id }->{spliceAI_medium}->add($lmdb_index);
						$intspan_genes_categories->{ $g->id }->{predicted_splice_site}->add($lmdb_index);
					}
					elsif (defined($value) and $value >= $buffer->config->{spliceAI}->{medium}) {
						$intspan_global_categories->{spliceAI_medium}->add($lmdb_index);
						$intspan_genes_categories->{ $g->id }->{spliceAI_medium}->add($lmdb_index);
						$intspan_global_categories->{predicted_splice_site}->add($lmdb_index);
						$intspan_genes_categories->{ $g->id }->{predicted_splice_site}->add($lmdb_index);
					}
				}
			}
			
			#$debug =1 if $g->id eq  "ENSG00000253317_8";
			my $cons_text = $variation->variationType($g);
			#die() if $variation->name eq "rs747758472";
			
			foreach my $c ( split( ",", $cons_text ) ) {
				confess( $cons_text . " " . $c ) unless exists $intspan_genes_categories->{ $g->id }->{$c};
				$intspan_genes_categories->{ $g->id }->{$c}->add($lmdb_index);
				my $prediction = $variation->categorie_frequency_prediction($g);
				$intspan_genes_categories->{ $g->id }->{$prediction}->add($lmdb_index);
			}
		}
			

		unless (@$genes) {
			################
			# hmm this variant seems to be intergenic
			# construct one intspan for intergenic region : $intspan_region_intergenic later I will intersect this region with real genomic intergenic region
			# and of course I also get the index in another intspan
			###############
			confess("intergenic problem") unless ( $variation->getChromosome()->intergenic_intspan()->contains( $variation->start ) );
			$intspan_region_intergenic->add_range( $variation->start - 100_000, $variation->end + 100_000 );
			$id_intergenic->add($lmdb_index);
		}
		delete $variation->{buffer};
		delete $variation->{project};
		$variation = undef;
		#end frequency
		$nb++;
	}
	
	$hres->{frequency} = $intspan_global_categories;
	$hres->{genes}     = $intspan_genes_categories;
	my $ints = $project->getChromosome($chr_name)->intergenic_intspan()->intersection($intspan_region_intergenic);
	$hres->{intergenic}->{intspan_region} = $ints;
	$hres->{intergenic}->{id} = $id_intergenic;
	$hres->{patients} = $hpatients;
	$no->close();
	my $chr     = undef;
	my $project = undef;
	my $buffer  = undef;
	return $hres;
}



sub init_genes_intspan {
	my $hintspan;
	foreach my $c ( keys %{ $categories->{genes}->{annotations} } ) {
		$hintspan->{$c} = Set::IntSpan::Fast::XS->new();
	}
	return $hintspan;
}

sub construct_bitVector_for_gene {
	my ( $no, $intspan_genes_categories, $size ) = @_;
	# construction d'une table gene du gene et surtout du bitvector pour un gene
	# table de hash du gene : {start} {end} , premier et dernier id du gene donc du intspan. {intspan} et  {size}  et bien sur j'ai rajouté {bitvector} apres tout le process et je sauvegarde
	
	
	return unless $intspan_genes_categories;
	
	foreach my $g ( keys %{$intspan_genes_categories} ) {
		my $hgene = $intspan_genes_categories->{$g}->{all};
		# let's start with all variants in genes
		my @array = $hgene->as_array;
		my $v = Bit::Vector->new_Enum($size, $array[0]."-".$array[-1]);
		$no->put_batch( $g, $v );
	}
}
