package Cache_chromosome_bitvector;
use strict;
use FindBin qw($RealBin);
use Data::Dumper;
use Parallel::ForkManager;
use Bio::DB::Sam;
use Storable qw(store retrieve freeze dclone thaw);
use Time::Duration;
use Math::Combinatorics;
use Devel::Cycle;
use Set::IntSpan::Fast::XS;
use Bit::Vector;
use Bit::Vector::Overload;
use Set::IntRange; 
use Tabix;
use JSON::XS;
use Bio::DB::Sam;
use POSIX qw(strftime);
use lib "$RealBin";
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
use validationQuery;
use GenBoNoSql;
use String::ProgressBar;
use Tabix;
use GenBoNoSqlLmdb;
use Compress::Snappy;
use Set::IntervalTree;
  use Clone 'clone';
use List::MoreUtils qw(natatime);

my $categories = {
		global =>{
				frequency =>{
						freq_none => 1,
						freq_1 => 1,
						freq_05 =>1,
						freq_01 => 1,
						freq_001 => 1,
						freq_0001 => 1,
						freq_ho_none => 1,
						freq_ho_1 => 1,
						freq_ho_05 =>1,
						freq_ho_01 => 1,
						freq_ho_001 => 1,
						freq_ho_0001 => 1,
						freq_he_none => 1,
						freq_he_1 => 1,
						freq_he_05 =>1,
						freq_he_01 => 1,
						freq_he_001 => 1,
						freq_he_0001 => 1,
						
				},
			 variation_type =>{
			 	substitution=>1,
			 	insertion=>1,
			 	deletion=>1,
			 	large_deletion=>1,
			 },
		},
		genes => {
				annotations => {
					utr=>1,
					splicing=>1,
					pseudogene=>1,
					coding=>1,
					maturemirna=>1,
					essential_splicing=>1,
					phase=>1,
					silent=>1,
					intergenic=>1,
					stop=>1,
					ncrna=>1,
					frameshift=>1,
					intronic=>1,
					"non-frameshift"=>1,
					prediction_0 => 1,
					prediction_1 => 1,
					prediction_2 => 1,
					prediction_3 => 1,
				},
			
			
		},
			patients => {
				all => 1,
				he=>1,
				ho=>1,
			}
	
};

sub getRegions {
	my ($chr,$fork) = @_;
	if ($fork == 1){
		my $regions;
		my $hregions;
		$hregions->{chromosome} = $chr->name;
		$hregions->{start} = 1;
		$hregions->{end} = $chr->length;
		push(@$regions,$hregions);
		return \@$regions;
	}
	
	my $tabix = $chr->buffer->software("tabix");
	my %hv;
	foreach my $patient (@{$chr->project->getPatients()}){
	my $calling_files = $patient->callingFiles();
	my @files;
	foreach my $m (keys %{$calling_files}){
		foreach my $v (values %{$calling_files->{$m}}){
			push(@files,$v);
		}
	}
	foreach my  $file (@files){
		my $file = $patient->getVariationsFiles->[0];
		next unless -e $file;
		open ( TOTO, "$tabix $file ".$chr->ucsc_name()." | cut -f 2 |");
		while (my $pos = <TOTO>){
			chomp($pos);
			$hv{$pos} ++;
		}
		close TOTO;
	
		}
	}
 my @snps = sort{$a <=> $b} keys %hv;
	if (scalar(@snps) == 0 ){
		my $hregions;
		$hregions->{chromosome} = $chr->name;
		$hregions->{start} = 1;
		$hregions->{end} = $chr->length;
		return [$hregions];
	}

	warn "end snp";
	my $nb;
	if ($fork == 1){
		
		$nb =scalar(@snps);
	}
	else {
		$nb = int(scalar(@snps)/($fork-1));
	}
	

	$nb = 7_000 if $nb > 7_000;
		my $regions;
	 my $iter = natatime $nb, @snps;
  while( my @tmp = $iter->() ){
   		my $hregions;
		$hregions->{chromosome} = $chr->name;
		$hregions->{start} = $tmp[0]-100;
		$hregions->{end} = $tmp[-1]+100;
		push(@$regions,$hregions);
  }
  $regions->[0]->{start} = $regions->[0]->{start} -1000;
  $regions->[-1]->{end} =  $regions->[-1]->{end}+1000;#$chr->length;

  return $regions;
}

sub getRegions1 {
	my ($chr, $fork) = @_;
		my $regions;
	if ($fork ==1 ){
		my $hregions;
		$hregions->{start} = 0;
		$hregions->{end} = $chr->length;
		$hregions->{chromosome} = $chr->name;
		push(@$regions,$hregions);
		return $regions;
	}
	

	my $len = $chr->end() + 100;
		#$this_limit = int($len / (2 * $fork));
	my $this_limit = int($len /  $fork);
	my $hRefCoord;
	my $id = 1;
	my $start = 1;
	my $next  = $this_limit;
	my $end   = $next;
	while ($end <= $chr->end()) {
		my $hregions;
		$hregions->{chromosome} = $chr->name;
		$hregions->{start} = $start;
		$hregions->{end} = $end;
			push(@$regions,$hregions);
		$start = $end;
		$end  += $next;
	}
	
	if ($end < $chr->end()){
		my $hregions;
		$hregions->{start} = $start;
		$hregions->{end} = $end;
		$hregions->{chromosome} = $chr->name;
			push(@$regions,$hregions);
	}
	
	return ($regions);
}



sub getIds {
	my ($project_name,$region) = @_;
	my $buffer  = new GBuffer;
	my $project = $buffer->newProject( -name => $project_name ); 
	my $reference = $project->getChromosome($region->{chromosome})->getReferences($region->{start}, $region->{end})->[0];
	my @all;
	my @patient_names =  sort{$a cmp $b} map{$_->name} @{$project->getPatients};
	my $hpatients;
	for (my $i=0;$i<@patient_names;$i++){
		$hpatients->{$patient_names[$i]} = $i;
	}
	


#	warn $region->{chromosome};

	foreach my $variation (@{$reference->getStructuralVariations}){
		my $hv;
		my $array_patients;
		my $aho= [];
			my $ap=[];
			foreach my $pat (@{$variation->getPatients()}) {
				my $pn = $pat->name();
				my $hp;
				$hp->{name} = $pn;
				my $patient_id = $hpatients->{$pn};
				push(@$ap,$patient_id);
				push(@$aho,$patient_id) if ($variation->isHomozygote($pat));
				$variation->{patients_details}->{$pn}->{vcf_infos} =  $variation->{check_id};
				$variation->{patients_details}->{$pn}->{he} =  $variation->{annex}->{$pat->id()}->{he};
				$variation->{patients_details}->{$pn}->{ho} =  $variation->{annex}->{$pat->id()}->{ho};
				$variation->{patients_details}->{$pn}->{he_ho_details} =  "he";
				$variation->{patients_details}->{$pn}->{he_ho_details} =  "ho"  if ($variation->{annex}->{$pat->id()}->{ho} eq '1');
				$variation->{patients_details}->{$pn}->{he_ho_details} .= ':'.$variation->{annex}->{$pat->id()}->{nb_all_ref}.':'.$variation->{annex}->{$pat->id()}->{nb_all_mut};
				$variation->{patients_details}->{$pn}->{type} = "he";
				$variation->{patients_details}->{$pn}->{type} = "ho"  if ($variation->isHomozygote($pat));
			}
			delete $variation->{buffer};
			delete $variation->{project};
			delete $variation->{annex};
			$variation->{heho_string} =  join(",",sort{$a <=> $b} @$ap);
			 if (scalar(@$aho)) {
			 	$variation->{heho_string} =  $variation->{heho_string}." ".join(",",sort{$a <=> $b} @$aho)." HO";
			}
			# line to prepare dejavu global;
			$hv->{obj} = compress(freeze ($variation));
			$hv->{start} = $variation->start;
			$hv->{end} = $variation->end;
			$hv->{id}= $variation->id;
			 push(@all,$hv);		
	}
	my @sort = sort{$a->{start} <=>$b->{start} or $a->{end} <=>$b->{end}} @all;
	return \@sort; 
}

sub getTempDir {
	my ($project,$chr) = @_;
	my $dir =  $project->getCacheBitVectorDir()."/tmp/";
	unless (-e $dir){
		mkdir $dir;
		system("mkdir -p $dir ;chmod a+rwx $dir ");
	}
	return $dir;
}

#72287

sub storeIds {
		my ($chr, $fork) = @_;
		
		my $patients = $chr->project->getPatients();
		my $regions = getRegions($chr,$fork);	
		my $pm = new Parallel::ForkManager($fork);
		
		my $all;

		
		
		#$no_obj->close();
		$pm->run_on_finish (
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;
			unless (defined($hRes) or $exit_code>0){
				print qq|No message received from child process $exit_code $pid!\n|; 
				#die();
				return;
			}
			
		
		#	
	
			#warn $no_obj->get_current_index();
			#append in the txt file all the variations found in this region and return by getIds
			
			my ($region) = grep {$_->{start} eq $hRes->{region}->{start} and $_->{end} eq $hRes->{region}->{end} } @$regions;
			die() unless $region;
			
			my $freeze_file = $chr->lmdb_cache_dir."/".$region->{chromosome}.".".$region->{start}."-".$region->{end}.'.freeze';
			$region->{freeze} = $freeze_file;
			store($hRes->{'variants'},$freeze_file);
		}
	);
	 my $cpt =1;	
		
		foreach my $region (@$regions){
				warn "$cpt/".scalar(@$regions);
				$cpt ++;
		#	next unless $region->{start} eq 12369510;
			my $pid = $pm->start and next;
			my $time_start = time;
			my $all = getIds($chr->getProject->name,$region);
			my $time_end=time;
			my $hres;
			$hres->{variants} = $all;	
			$hres->{region} = $region;	
			$hres->{ttime} = abs($time_start-$time_end);	
			$pm->finish(0, $hres);	
		}
		
	 $pm->wait_all_children();
	 warn "end process";
	#end fork now I have a file with all variation and  json  it's time to sort this file and store it in lmdb database
	#construct intspan 	for patient I will store lmdb_id in intspan 
	my $project = $chr->project;
	my @patient_names =  sort{$a cmp $b} map{$_->name} @{$project->getPatients};
	
	my $index_patients = 0;
	my $hpatients;
	
	my @categorie_patient =("all","he","ho");
	foreach my $pname (@patient_names){

			$hpatients->{$pname}->{index} = $index_patients ++;
			$hpatients->{$pname}->{name} = $pname;
			foreach my $c (@categorie_patient){
			$hpatients->{$pname}->{intspan}->{$c}   = Set::IntSpan::Fast::XS->new();
			}
	}
	
	
	#
	#initialisation global categorie
	#
	
	my @global = ("substitution","insertion","deletion","large_deletion");
	my $intspan_global_type;
	foreach my $g (keys%{$categories->{global}->{variation_type}}){
		$intspan_global_type->{$g}   = Set::IntSpan::Fast::XS->new();
	}
			
	my $no2 = $chr->get_lmdb_variations("c"); #open lmdb database
	
	#ok sort and read the file
		my $hh;
		warn "start store varid";
	#		my $hh1 = retrieve($chr->lmdb_cache_dir."/dv.freeze");
	#		warn scalar(keys %$hh1);
	#open(SORT," sort -n  -k1,1 -k2,2 $txt_file |  ");
	my $uniq;
	my $size_variants = 0;
	foreach my $region (@$regions){
		my $freeze = $region->{freeze};
		die() unless -e $freeze;
		my $hall  = retrieve $freeze;
			foreach my $hv (@{$hall}){
				my $var_id  = $hv->{id};
				next if exists $uniq->{$var_id};
			 $uniq->{$var_id} ++;
			 my $variation = thaw (decompress($hv->{obj}));;
			my $index_lmdb = $no2->put($var_id,$variation);
#			warn $no2->size()." ".$index_lmdb;
			 $size_variants ++;
			 $hh->{$var_id}  =  $variation->{heho_string};
			 foreach my $pn (keys %{$variation->{patients_details}}){
			 	my $type = $variation->{patients_details}->{$pn}->{type};
			 	 $hpatients->{$pn}->{intspan}->{all}->add($index_lmdb);
			 	$hpatients->{$pn}->{intspan}->{$type}->add($index_lmdb);
			 }
#			 	die() unless exists $intspan_global_type->{$variation->type};
			$intspan_global_type->{$variation->type}->add($index_lmdb);
			unlink $freeze;
		
	}
			 
	}

	my $chr_name = $chr->name();
	store($hh, $project->lmdb_cache_dir."/$chr_name.dv.freeze") if $hh; #store htable only fort dejavu by project purpose 	
	
	
	warn "start patients";
	my $no3 = $chr->get_lmdb_categories("c"); #open lmdb database
	foreach my $k (keys %{$intspan_global_type}){
		my $h;
		$h->{name} = $k;
		$h->{intspan} = $intspan_global_type->{$k};
		my $bitv = Bit::Vector->new_Enum( $size_variants, join(',', $intspan_global_type->{$k}->as_array) );
		$h->{bitvector} = $bitv;			
		$no3->put($k,$h);
	}
	$no3->close();	
	
	my $no4 = $chr->get_lmdb_patients("c");
	foreach my $pname (@patient_names){
		my $h;
		$h->{name} = $pname;
		foreach my $c (@categorie_patient){
			my $intspan =  $hpatients->{$pname}->{intspan}->{$c};
			$h->{intspan}->{$c}=$intspan;
			my $bitv = Bit::Vector->new_Enum( $size_variants, join(',', $intspan->as_array) );
			$h->{bitvector}->{$c}=$bitv;
			
			}
			$no4->put($pname,$h);
	}

$no4->close;
}



sub init_genes_intspan {
	
	my $hintspan;
		foreach my $c (keys %{$categories->{genes}->{annotations}}){
					$hintspan->{$c} = Set::IntSpan::Fast::XS->new();
	}
	return $hintspan;
}







###################### 
# # variations annotation
######################

sub annotations {
	my ($project_name,$chr_name,$region) = @_;
	#init genbo project 
		my $buffer = new GBuffer;
		my $project = $buffer->newProject( -name => $project_name );
		
	# prepare intspan for global categories
		my $intspan_global_categories;
		
		foreach my $c (keys %{$categories->{global}->{frequency}}){
			$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
		}
		
		my $intspan_genes_categories ={};
		
		#initialisattion categorie patients 
		
	my $index_patients = 0;
	my $hpatients;
	my @patient_names =  sort{$a cmp $b} map{$_->name} @{$project->getPatients};
	my @categorie_patient =("all","he","ho");
	foreach my $pname (@patient_names){
			foreach my $c (@categorie_patient){
			$hpatients->{$pname}->{$c}   = Set::IntSpan::Fast::XS->new();
			}
	}
	
		
		#for intergenic variant construct region 
		my $intspan_region_intergenic = Set::IntSpan::Fast::XS->new();
		#and store intergenic id;
		my $id_intergenic  = Set::IntSpan::Fast::XS->new();;
		
	
		my $no =  $project->getChromosome($chr_name)->get_lmdb_variations("r");

		#my $no_obj =  $project->getChromosome($chr_name)->get_lmdb_variations_object("r");
		my $cursor =  $no->cursor($region->[0],$region->[1]);
		my $nb =0;
		my $l = abs($region->[0] - $region->[1]);
		my $cpt = $region->[0];
		while (my $var_id = $cursor->next_key) {
				
				my $variation;
				my $lmdb_index = $cursor->current_index();
				  $variation =  $no->get($var_id);#$project->_newVariant($var_id);#$no->get($var_id);
				
				# my $vv =  $no->get($var_id);
				 if ($lmdb_index <  0) {
				 die("index out of range <0");
				 }
			
				 warn "************************************\n.". $var_id." ".Dumper($variation)."*--------------------*\n********" if $lmdb_index == -2;
				 unless ($variation){
				 		die();
				 		 $variation = $project->_newVariant($var_id);
				 }
				 die() unless $variation->id;
				# warn $variation;
				
				 $variation->{buffer} = $buffer;
				  $variation->{project} = $project;
				  
	################	
	# patient  Annotation 
	################	
#			 foreach my $pn (keys %{$variation->{patients_details}}){
#			 		my $type = $variation->{patients_details}->{$pn}->{type};
#			 	 	$hpatients->{$pn}->{all}->add($lmdb_index);
#			 		$hpatients->{$pn}->{$type}->add($lmdb_index);
#			 } 				  
				  
	################				  
	#	get freqency categories and construct intspan for each
	###############				
				my $cat_freq = $variation->categorie_frequency();
				$intspan_global_categories->{$cat_freq}->add($lmdb_index);
				my $cat_freq_ho = $variation->categorie_frequency_ho();
				$intspan_global_categories->{$cat_freq_ho}->add($lmdb_index);				
				my $cat_freq_he = $variation->categorie_frequency_he();
				$intspan_global_categories->{$cat_freq_he}->add($lmdb_index);
				
	################				  
	#	get Genes annotations prediction and functional 
	###############		
				
				my $genes = $variation->getGenes();
				foreach my $g (@$genes){
						#die("intergenic problem") if  ($variation->getChromosome()->intergenic_intspan()->contains($variation->start));
						my $cons_text = $variation->variationType($g);
						unless (exists $intspan_genes_categories->{$g->id}){
							$intspan_genes_categories->{$g->id} = init_genes_intspan();
						}
						foreach my $c (split(",",$cons_text)){
						die($cons_text." ".$c) unless exists 	$intspan_genes_categories->{$g->id}->{$c};
						$intspan_genes_categories->{$g->id}->{$c}->add($lmdb_index);
						my $prediction = $variation->categorie_frequency_predicion($g);
						$intspan_genes_categories->{$g->id}->{$prediction}->add($lmdb_index);
						}
				}
				
				unless (@$genes){
					################				  
					#	hmm this variant seems to be intergenic 
					# construct one intspan for intergenic region : $intspan_region_intergenic later I will intersect this region with real genomic intergenic region
					# and of course I also get the index in another intspan
					###############		
					die("intergenic problem") unless  ($variation->getChromosome()->intergenic_intspan()->contains($variation->start));
					$intspan_region_intergenic->add_range($variation->start-100_000,$variation->end+100_000);
					$id_intergenic->add($lmdb_index);
				}
				
				 delete $variation->{buffer} ;
				 delete  $variation->{project};
				 
				#end frequency 
#				
				
				
				$nb++;
				#warn $nb."/".$l  if $nb%10000 ==0;		
				
		}
		my $hres;
		$hres->{frequency} = $intspan_global_categories;
		$hres->{genes} = $intspan_genes_categories;
		my $ints = $project->getChromosome($chr_name)->intergenic_intspan()->intersection($intspan_region_intergenic);
		$hres->{intergenic}->{intspan_region} = $ints;
		
		$hres->{intergenic}->{id} = $id_intergenic;
		$hres->{patients} = $hpatients;
		$no->close();
		return $hres;
}



sub store_annotations {
	my ($project_name,$chr_name,$fork) = @_;
	my $buffer = new GBuffer;
	my $project = $buffer->newProject( -name => $project_name );
	my $chr= $project->getChromosome($chr_name);
	my $no = $chr->get_lmdb_variations("r");
	warn $no->dir();
	my $ranges = $no->ranges($fork);
	$no->close;
	$no = undef;
	my $intspan_genes_categories ={};
	my $intspan_global_categories = {};
	
		foreach my $c (keys %{$categories->{global}->{frequency}}){
			$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
		}
		
	my $process;	
	my $pm = new Parallel::ForkManager($fork);
	
		#for intergenic variant construct region 
	my $intspan_region_intergenic = Set::IntSpan::Fast::XS->new();
	#and store intergenic id;
	my $id_intergenic= Set::IntSpan::Fast::XS->new();
	my $hintspan_patients;
	
####################################	
# one fork is finished i could "union"  all the intspan construct in the "sub annotation()" 	
#run on finish 	aggregate all result (intpsan) from child
#########################################
		$pm->run_on_finish (
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hres) = @_;

			unless (defined($hres) or $exit_code>0 ){
				print qq|No message received from child process $exit_code $pid!\n|; 
				warn Dumper $hres;
				#die();
				return;
			}
			if (exists $hres->{start} ){
				print qq|1- No message received from child process $exit_code $pid!\n|; 
				warn Dumper $hres;
				#die();
				return;
			}
		delete $process->{$hres->{idp}};
		my $no = $chr->get_lmdb_variations("r");
		#################
		# start with patients (gloabal intspan : on chromosome)
		#################
#			 foreach my $pn (keys %{$hres->{patients}}){
#			 	foreach my $type (keys %{$hres->{patients}->{$pn}}){
#			 		unless (exists $hintspan_patients->{$pn}->{$type} ){
#			 			$hintspan_patients->{$pn}->{$type} =Set::IntSpan::Fast::XS->new();
#			 		}
#			 	 	$hintspan_patients->{$pn}->{$type} = $hintspan_patients->{$pn}->{$type}->union($hres->{patients}->{$pn}->{$type});
#			 	}
#			 } 	
		#################
		#For frequencies  annotation categorie (intspan on chromosome)
		#################
		$intspan_region_intergenic = $intspan_region_intergenic->union($hres->{intergenic}->{intspan_region});
			$id_intergenic = $id_intergenic->union($hres->{intergenic}->{id});
			#For global frequency categorie
			foreach my $cat (keys %{$hres->{frequency}}){
				warn $cat unless exists $intspan_global_categories->{$cat};
				warn $cat unless exists $intspan_global_categories->{$cat};
				$intspan_global_categories->{$cat} = $intspan_global_categories->{$cat}->union($hres->{frequency}->{$cat});
			}
			
			#################
			#For genes annotation categorie intspan on gene 
			#################
			foreach my $g (keys %{$hres->{genes}}){
					my $hgene = $hres->{genes}->{$g};
					unless (exists $intspan_genes_categories->{$g}){
							$intspan_genes_categories->{$g} = init_genes_intspan();
	
							$intspan_genes_categories->{$g}->{all} = Set::IntSpan::Fast::XS->new();
						}
						
						foreach my $cat (keys %{$hgene}){
							die($cat) unless exists  $intspan_genes_categories->{$g}->{$cat};
							$intspan_genes_categories->{$g}->{$cat} = $intspan_genes_categories->{$g}->{$cat}->union($hgene->{$cat});
							$intspan_genes_categories->{$g}->{all} = $intspan_genes_categories->{$g}->{all}->union($hgene->{$cat});
						}
			}
			$no->close();
			
		}
		);
#end run on finish 


##############
#fork annotation
##############


	my $project_name = $chr->project->name();
	my $chr_name = $chr->name();
	
	my $true = 0;
	my $idp =0;
	while ($true <2 ){ #check if all the region end correctly unless restart region failed 
	foreach my $r (@$ranges){
		$process->{$r->[0]."-".$r->[1]}  = $r;
		my $pid;
		sleep(1);
		 $pid = $pm->start and next;
			my $hres;
			$hres->{start} = $r->[0]."-".$r->[1];
			
			my $hres  = annotations($project_name,$chr_name,$r);
			$hres->{idp} = $r->[0]."-".$r->[1];
			delete $hres->{start} ;
			$hres->{done} =1;
		$pm->finish(0,$hres);
	}#end for range range
	
	$pm->wait_all_children();
	###########
	#end fork 
	###########
	unless (keys %$process){ #all region are ok
		last;
	}
	#at least one region failed restart only this region
	$true ++;
	@$ranges = values %$process;
	warn "second run";
	}
	
	
	
	
	die("problem store process") if (keys %$process);
	warn Dumper $process;
		warn "end ----->";
	my $no =  $project->getChromosome($chr_name)->get_lmdb_variations("r");
	my $size_variants = $no->size();
	$no->close();
	warn "nb variants : ".$size_variants;
	#########################
	#GLOBAL CATEGORIES
	# ok now all part is finished I can store the global intspan and construct bitvector for each global categories
	##########################
	my $no3 = $chr->get_lmdb_categories("w"); #open lmdb database
	
	foreach my $k (keys %{$intspan_global_categories}){
			my $h;
		$h->{name} = $k;
		$h->{intspan} = $intspan_global_categories->{$k};
		my $bitv = Bit::Vector->new_Enum( $size_variants, join(',', $intspan_global_categories->{$k}->as_array) );
		$h->{bitvector} = $bitv;			
		$no3->put($k,$h);

	}
	
	
	$no3->close();	
	warn "end categories !!!";
	
	#########################
	#Patients  CATEGORIES
	# ok now all part is finished saved patients intspan , all , he,  ho 
	##########################
	my $no5 = $chr->get_lmdb_patients("c"); #open lmdb database cretae new one or erase older
	
	 foreach my $pn (keys %{$hintspan_patients}){
	 		my $h;
	 		$h->{name} = $pn;
			 	foreach my $type (keys %{$hintspan_patients->{$pn}}){
			 		$h->{intspan}->{$type} = $hintspan_patients->{$pn}->{$type};
			 		my $bitv = Bit::Vector->new_Enum( $size_variants, join(',', $hintspan_patients->{$pn}->{$type}->as_array) );
					$h->{bitvector}->{$type} = $bitv;
			 	}
			 	$no5->put($pn,$h);
			 } 	
	$no5->close();
	$no5 =undef;
	
	# DONE !!!  continue with genes a bit more complex because i have to deal with subvector and create new index
	
	#################
	# GENES 
	#################
		my $no4 = $chr->get_lmdb_genes("c"); #open lmdb database
		construct_bitVector_for_gene($no4,$intspan_genes_categories);
		$no4->close();
	


my @intergenic;
my $tree = Set::IntervalTree->new;

my $iter = $intspan_region_intergenic->iterate_runs();


my $i =0;

#my $no2 = $chr->get_lmdb_variations_object("r");
my $no3 = $chr->get_lmdb_variations("r");
 while (my ( $from, $to ) = $iter->()) {
 	my $id = 'intergenic_'.$chr_name.'_'.$from.'_'.$to;
 	$tree->insert($id,$from,$to+1);
  }
	

my @array = $id_intergenic->as_array;

warn scalar(@array)." intergenic variations";	
	
	my $t = time;
	warn $array[0]." ".$array[-1];
	
	my $intspan_intergenic_categories= {};
	foreach my $lmdb_id (@array){
		
		my $v = $no3->get_index($lmdb_id);
		
		my $results = $tree->fetch($v->{start},$v->{end}+1);
		
		die() unless @$results;
		die() if scalar(@$results)>1;
		my $interid = $results->[0];
		unless (exists $intspan_intergenic_categories->{$interid}){
			$intspan_intergenic_categories->{$interid} = init_genes_intspan();
			$intspan_intergenic_categories->{$interid}->{all} = Set::IntSpan::Fast::XS->new();
			$intspan_intergenic_categories->{$interid}->{intergenic} = Set::IntSpan::Fast::XS->new();
		}
		$intspan_intergenic_categories->{$interid}->{all}->add($lmdb_id);
		$intspan_intergenic_categories->{$interid}->{intergenic}->add($lmdb_id);
		
	}
	
my $no4 = $chr->get_lmdb_genes("w"); #open lmdb database
construct_bitVector_for_gene($no4,$intspan_intergenic_categories,1);
$no4->close();
#my $no4 = $chr->get_lmdb_genes("r");
#		
#my $genes = $no4->get_values();
#my $tree1 = Set::IntervalTree->new;
#my $test = Set::IntSpan::Fast::XS->new();
#foreach my $g (@$genes){
#			my $intspan = $g->{intspan}->{all};
#			my @t = $intspan->as_array();
#			warn $t[0]." ".$t[-1]." ".$g->{name};
#			warn $g->{start}." ".$g->{end};
#			$test->add_range( $g->{start},$g->{end});
#			$tree1->insert($g->{name},$t[0],$t[-1]+1);
#		}
#	#my $tt = dclone($tree1);
#
# foreach my $g (@$genes){
# 		my $intspan = $g->{intspan}->{all};
#		my @t = $intspan->as_array();
#
#		my $results = $tree1->fetch($t[0],$t[-1]+1);
#	
#		warn $g->{name}." ".$results->[0];
#		
# }
# $no4->close();
#warn "end";
#die();	
	
		
		
		
}


sub convert_small_intspan_to_bitvector {
	my ($intspan,$start,$size) = @_;
		my $bitv =  Bit::Vector->new_Enum( $size, join(',',map{$_ - $start}  $intspan->as_array) );
		return $bitv;
}

sub convert_global_intspan_to_bitvector {
	my ($intspan,$size) = @_;
	my $bitv = Bit::Vector->new_Enum( $size, join(',', $intspan->as_array) );
	return $bitv;
}


sub construct_marc_hash {
		my ($project_name,$chr_name,$fork) = @_;
		#
		#c'est ici qu eje genere les bitvectors a partir des  intspans .
		#
		
		my $buffer = new GBuffer;
		my $project = $buffer->newProject( -name => $project_name );
		my $chr= $project->getChromosome($chr_name);
		
		my $no_patients = $chr->get_lmdb_patients("r"); #intspan pour les patients 
		
		my $no_genes = $chr->get_lmdb_genes("r"); #intspan pour les genes 
		my $no_categories = $chr->get_lmdb_categories("r"); #intspan pour les autres categories  
		
		my $genes = $no_genes->get_values();
		my $tree1 = Set::IntervalTree->new;
		my $test = Set::IntSpan::Fast::XS->new();
		
		#je ne sais plus pourquoi je fais ca mais dans mon idée c'etait de généré un interval tree pour connaitre le premier gene a regarder rapidement
		#c'etait sans doute un test qui ne sert plus 
		
		warn "start tree ";
		foreach my $g (@$genes){
				my $intspan = $g->{intspan}->{all};
				my @t = $intspan->as_array();
				$tree1->insert($g->{name},$t[0],$t[-1]+1);
		}
		
		warn "end tree ";
		
		# construct "useless"" global categorie
		# en fait c'est toute les actegories que je trouvais bizarre  qui concernent les genes mais que tu avais posé sur le chromosome
		 my $hintspan = init_genes_intspan(); #initialisation de tous les intspans 
		 
		 # $hintspan->{le nom de la categoire } et j'ai fais une boucle sur tous les genes pour cette categorie precise
		 
		 foreach my $gene (@$genes){
		 	foreach my $c (keys %{$categories->{genes}->{annotations}}){
		 		my $intspan = $gene->{intspan}->{$c}; #c'est initialisé plus haut pas besoin de regarder si ils existent
		 		$hintspan->{$c} = $hintspan->{$c}->union($intspan); #et j'ajoute a mon intspan global
		 	}
		 	
		 }
		 warn "end global";
		 
		 my $patients = $no_patients->get_values(); #recupere toutes les tav-bles (intspans de chaque patients)
		 
		  foreach my $patient (@$patients){
		  	my $intspan = $patient->{intspan}->{all}; #categorie all
		  	foreach my $c (keys %{$categories->{genes}->{annotations}}){
				my $i1 = $intspan->intersection($hintspan->{$c}); #j'intersecte ma categorie avec la liste  des variants de cette categories sur le chromosome
				my %hg;
				 my $iter = $i1->iterate_runs();
				 my $cpt =0;
				  while (my ( $from, $to ) = $iter->()) {
				  	$cpt ++;
				  	 my $results = $tree1->fetch($from,$to+1);
				  	 map{$hg{$_} ++} @$results;
				  	
				  }
				  warn join(";",keys %hg);
				warn $cpt;
			
			}
		  	
		  }
		  
		
		 
		
		

		
} 



sub construct_bitVector_for_gene{
	my ($no,$intspan_genes_categories,$debug) = @_;
	#construction d'une table gene du gene et surtout du bitvector pour un gene
	#table de hash du gene : {start} {end} , premier et dernier id du gene donc du intspan. {intspan} et  {size}  et bien sur j'ai rajouté {bitvector} apres tout le process et je sauvegarde
	return unless $intspan_genes_categories;
	foreach my $g (keys %{$intspan_genes_categories}){
		my $hgene = $intspan_genes_categories->{$g};
		# let's start with all variants in genes
			my @array = $hgene->{all}->as_array;
			my $h;
			$h->{start} = $array[0];
			$h->{end} = $array[-1];
			$h->{size} = ($h->{end} - $h->{start})+1;
			$h->{name} = $g;
			$h->{bitvector}->{all} =  Bit::Vector->new_Enum( $h->{size}, join(',',map{$_ - $h->{start}}  @array) );
			
		foreach my $cat (keys %{$hgene}){
			next if ($hgene->{$cat}->is_empty());
			$h->{intspan}->{$cat} = $hgene->{$cat};
			my @array =  map{$_ - $h->{start}} $hgene->{$cat}->as_array;
			$h->{bitvector}->{$cat} = Bit::Vector->new_Enum( $h->{size}, join(',',@array) );
			#warn $cat.' - '.$h->{bitvector}->{$cat}->Size();
		}
		#warn Dumper $hgene;
		#warn "\n\n";
		#warn Dumper $h; die;
		#$h->{intspan} = $hgene;
		$no->put($g,$h);
	}
	
}

sub cache_lite_for_dejavu {
	my ($project_name,$fork) = @_;

	#my $root_dir = "/data-isilon/dejavu/projects/";
	warn "*** CAche For Deja Vu *****";
	my $buffer1 = new GBuffer;
	my $project = $buffer1->newProject( -name => $project_name, -verbose =>1 );
	my $root_dir =$project->deja_vu_lite_dir()."/projects/";

		mkdir $root_dir unless -e $root_dir;
	unlink $root_dir."/".$project_name.".lite" if -e $root_dir."/".$project_name.".lite";
	my @chr_names = map{$_->name} @{$project->getChromosomes};
	#my $pm = new Parallel::ForkManager(1);
	my @patient_names =  sort{$a cmp $b} map{$_->name} @{$project->getPatients};
	my $dir_out = 	 $project->getCacheBitVectorDir()."/lmdb_cache";
	my $hpatients;
	for (my $i=0;$i<@patient_names;$i++){
		$hpatients->{$patient_names[$i]} = $i;
	}
	my $no = GenBoNoSql->new(dir=>$root_dir,mode=>"c");
	$no->put($project_name,"patients",$hpatients);
	my $atleast;
	foreach my $chr (@{$project->getChromosomes}) {
		my $fileout = $dir_out."/".$chr->name.".dv.freeze";; 
		warn "miss $fileout " unless -e $fileout;
		next unless -e $fileout;
		$atleast ++;
		my $h =  retrieve $fileout;
		$no->put($project_name,$chr->name,$h);
	}
	
	confess() unless $atleast;

}

sub cache_cnv {
	my ($project_name,$fork) = @_;
	my $buffer1 = new GBuffer;
	my $projectP = $buffer1->newProject( -name => $project_name, -verbose =>1 );
	#my @transcripts_cgi = @{$projectP->bundle_transcripts() } ;
	#my @transcripts = sort{$a->getGene->external_name cmp $b->getGene->external_name} map{$projectP->newTranscript($_)} @transcripts_cgi ;
	warn "###################################################\n";
	warn "Prepare Cache for  CNV  !!!\n";
	warn "###################################################\n";
	my $primers = $projectP->getPrimers();
	warn scalar(@$primers);
	
	foreach my $patient (@{$projectP->getPatients}){
	warn $patient->name;		
	foreach my $primer (@$primers){
		warn Dumper $primer->cached_cnv($patient);
		#print ".";
	}
	}
	#print "\n";
#	preload_coverage::load_cnv_score($projectP,$projectP->getPatients,\@transcripts);
	
	warn "###################################################\n";
	warn "END CNV  !!!\n";
	warn "###################################################\n";
	warn "\n";
}
1;
