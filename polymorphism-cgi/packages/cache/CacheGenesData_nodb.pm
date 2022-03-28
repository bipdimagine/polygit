package CacheGenesData_nodb;
use strict;

use FindBin qw($RealBin);
use Data::Dumper;
use Parallel::ForkManager;
 use Bio::DB::Sam;
#use PBS::Client;
use Storable qw(store retrieve freeze dclone thaw);
use Time::Duration;
use Math::Combinatorics;


sub basic_gene_infos2 {
	my ($g,$coverage) = @_;
	 my $debug ;
	$debug = 1 if $g->id eq "ENSG00000169093_X";
	my $pos = $g->position($g->getChromosome);
	warn $pos->start if $debug;
	my %gene_data1;
	$gene_data1{start}  = $pos->start();
	$gene_data1{end}    = $pos->end();
	$gene_data1{chromosome}    = $g->getChromosome->name();
	$gene_data1{strand} = $pos->strand();
	$gene_data1{name}   = $g->name();
	$gene_data1{id}   = $g->id();
	$gene_data1{cover} = 0;
	if ($g->isGene()){
		$gene_data1{gene} = 1;
		$gene_data1{reference} = join(";", map {$_->name()} @{$g->getReferences});    
		              
		$gene_data1{xref}      = $g->external_name();
		$gene_data1{description} = $g->description();
		$gene_data1{transcripts} = join (";",map {$_->name} @{$g->getTranscripts()} );
		$gene_data1{cover} = int($coverage->{total} / $coverage->{len}) if exists $coverage->{len};
	}
	return \%gene_data1;
}

sub basic_gene_infos {
	my ($gene_data,$g,$coverage) = @_;
	my $id = $g->id;
	
	return if $gene_data->{$id}->{ok} == 1;
	
		# $gene_data->{$id} = {};
	 my $pos = $g->position($g->getChromosome);
	
	$gene_data->{$id}->{start}  = $pos->start();
	$gene_data->{$id}->{start}  = 1 if $gene_data->{$id}->{start} <=1;
	$gene_data->{$id}->{end}    = $pos->end();
	$gene_data->{$id}->{chromosome}    = $g->getChromosome->name();
	$gene_data->{$id}->{strand} = $pos->strand();
	$gene_data->{$id}->{name}   = $g->name();
	$gene_data->{$id}->{id}   = $g->id();
	$gene_data->{$id}->{cover} = 0;
	if ($g->isGene()){
		$gene_data->{$id}->{gene} = 1;
		$gene_data->{$id}->{reference} = join(";", map {$_->name()} @{$g->getReferences});    
		              
		$gene_data->{$id}->{xref}      = $g->external_name();
		$gene_data->{$id}->{description} = $g->description();

   
	
		
		$gene_data->{$id}->{transcripts} = join (";",map {$_->name} @{$g->getTranscripts()} );
		$gene_data->{$id}->{cover} = int($coverage->{total} / $coverage->{len}) if exists $coverage->{len};
	}
	
	$gene_data->{$id}->{ok} = 1;

	return;
	
}
 
 
 sub basic_regions_infos {
	my ($gene_data,$chr,$start) = @_;
	my $id = $chr->name."_".$start;
	return $id if $gene_data->{$id}->{ok} == 1;
		# $gene_data->{$id} = {};
	 
	$gene_data->{$id}->{intergenic}  = 1;
	$gene_data->{$id}->{start}  = $start;
	$gene_data->{$id}->{start}  = 1 if  $start <=1;
	$gene_data->{$id}->{end}    = $start;
	$gene_data->{$id}->{chromosome}    = $chr->name;
	$gene_data->{$id}->{strand} = 1;
	$gene_data->{$id}->{name}   = "inter_".$id;
	$gene_data->{$id}->{id}   = $id;
	$gene_data->{$id}->{cover} = 0;
	
	$gene_data->{$id}->{ok} = 1;

	return $id;
	
}

sub kyoto_dejavu {
	my ($project, $chr_name) = @_;
	my $dir = $project->getDejaVuDir();
	return $dir."/".$project->name().".".$chr_name.".kct";
}

sub open_kyoto{
	my ($file) = @_;
	my $db1 = new KyotoCabinet::DB;
	if (!$db1->open($file, $db1->ONOLOCK | $db1->OCREATE |$db1->OWRITER )){
			if (-e $file){
				unlink($file);
				if (!$db1->open($file, $db1->ONOLOCK | $db1->OCREATE |$db1->OWRITER )){
					printf STDERR ("open error: %s\n", $db1->error);
					die();
				}
			}
		
			
		}
		return $db1;
} 
	
sub create_cache_genes {	
	my ( $project,$chr,$hvars,$hgenes,$fgenes,$hpatients) = @_;
	$| =1;
	my $buffer = $chr->buffer();

	my $file_kyoto = kyoto_dejavu($project, $chr->name());
	my $db1 = open_kyoto($file_kyoto);
	$db1->clear();
	
	my @array;
	

	
	
#	warn $nb_total;
	my $nb_references = 0;
	
	my $chr_data;
	my $gene_coverage;
	my $test_coding = 0;
	
	my $gene_data = {};
	
	my $test;
	
	#my @refs_ids = map {$_->name} @$refs; 
	

	my $variation_data;
	my $sc_waiting;
	
	my $scale_pourcent;
	my $data_variations;
	my $limit_pourcent = 10;
	my $references = $chr->getReferences();
	#my $vars  = $chr->getStructuralVariations();
	
	my $newgene;
	
	my @references_name = map{$_->name()} @$references;
	my $chr_name= $chr->name();
	my $project_name = $project->name();
	
	#my $vars  = $chr->getStructuralVariations();
#	my $nb_ref_total = scalar(@$vars);
	my $nb_ref_total = scalar(@$references);
foreach my $reference ( @$references) {
	
		  $nb_references ++;

	my $t =time;
	my @vars = sort {$a->start <=> $b->start } @{$reference->getStructuralVariations()};
	my $tstart = time;
	my $t1 = time;
	my $var_total = scalar(@vars);
	my $nb_var =0;
	while (my $variation = shift (@vars) ){
		#deja_vu file
			$db1->set($variation->id,$variation->json_for_kyoto().";");
				
			$nb_var ++;
			 my $pourcent = int(($nb_var / $var_total)*100 );
		     if ($pourcent >0 && ($pourcent % $limit_pourcent) == 0){
		     ##	$buffer->test_flush();
		     	my $ttime = (time - $tstart);
		     	
		     	my $rem = ($ttime * 100) / $pourcent;
		     	my $still = 1- $pourcent/100;
		     	warn "\tchr :".$chr->name()." $pourcent\% remaining time : ".concise(duration($rem*$still))."\n";
		     	#warn "\t".$chr->name()."  ".$pourcent."% \t".(time - $t1)."\t".(time - $tstart)."\n";
		     	$t1 = time;
		     	$limit_pourcent +=20;
		     }
		#next if $variation->id ne 9924659;
			my @trs;
			
		

			
		
		
		 	my $items;
		 	my $pos1 = $variation->position($chr);
		 	
		 	$items->{structural_type} = $variation->getStructuralType();
			$items->{chromosome}  = $chr->name();
			$items->{type}  = $variation->getType();
			
			$items->{id}    = $variation->id();
			$items->{start} = $pos1->start();
			$items->{end} = $pos1->end();
			$items->{polyphen_status} = $variation->polyphenStatus();
			$items->{sift_status} = $variation->siftStatus();
			if ($variation->isStop){
				$items->{polyphen_status} = 5;
				$items->{sift_status} = 5;
			}
			if ($variation->isPhase){
				$items->{polyphen_status} = 4;
				$items->{sift_status} = 4;
			}
			$items->{text}  = $variation->getChromosome->sequence($pos1->start(),$pos1->end())."/".$variation->sequence();

			#
			
			my $vdb = $variation->origin_database();
			my $db_freq = undef;
			unless ($variation->isClinical){
				$db_freq = "1p10000"  if ($variation->frequency() <= 0.0001 &&  $variation->frequency >0);
				$db_freq = "1p1000"  if ($variation->frequency() <= 0.001 &&  $variation->frequency >0.0001);
				$db_freq = "1p100"  if ($variation->frequency() <= 0.01 && $variation->frequency >0.001);
				$db_freq = "5p100"  if ($variation->frequency() <= 0.05 && $variation->frequency >0.01);
				$db_freq = "common"  if ($variation->frequency() > 0.05);
			}
			else {
				$vdb = "pheno_snp";
			}
			
			$items->{public} = $vdb;
			
			$items->{"consequence!all"} = $variation->variationType();
			$items->{reference} = $variation->getReference->name();
			$items->{name}     = $variation->name();
			my $opatients = $variation->getPatients();
			$items->{nb_patients} = $variation_data->{nb_patient} = scalar @{$opatients};
			$items->{filter} = $variation->getNGSScore();
			if (scalar @{$opatients} == 0){
				warn $chr->name()." ".$variation->id()." ".$variation->position($chr)->toString;
			}
			foreach my $p (sort { $a->name() cmp $b->name() } @{$opatients}) {	
						push(@{$variation_data->{$variation->id}->{patients}},$p->name);
				}
			$items->{allpatients} = join( ";", @{$variation_data->{$variation->id}->{patients}} );
			

				### public or not
				
			
				############
				# nb patient
				############
			
				
				
					
				my @ids_references;
				
				if ( scalar(@{$variation->getGenes()}) ) {
					$test_coding ++;
						 
					foreach my $g (@{$variation->getGenes()}) {
						push(@ids_references,$g->id);
					 
						push(@{$variation_data->{$variation->id}->{genes}},$g->name);
						my $pos = $g->position($chr);
						basic_gene_infos($gene_data,$g,$gene_coverage->{$g->name}) ;
					unless ( exists $gene_data->{$g->id}){
						warn scalar ( keys % $gene_data);
						warn join (";", keys % $gene_data);
						warn $g->id;
							confess() unless $gene_data->{$g->id}->{start};
						}
						my $cons_text = $variation->variationType($g);
						
						push( @{ $items->{genes} }, $g->name );
						
					
						$items->{$g->name()."_consequence"} = $cons_text;	
						my @gene_variation_type =  split(",",lc($cons_text));
						foreach my $typev (@gene_variation_type){
							#	warn $variation->variationType($g);
							push(@{$gene_data->{$g->id}->{data}->{$typev}},$variation->id);
						}
						
						# polyphen && sift
						my $key_polyphen = "polyphen".$variation->polyphenStatus($g);
						push(@{$gene_data->{$g->id}->{data}->{$key_polyphen}},$variation->id);
						
						my $key_sift = "sift".$variation->siftStatus($g);
						push(@{$gene_data->{$g->id}->{data}->{$key_sift}},$variation->id);
							
						my $type_score = "score".$variation->getNGSScore();
						#warn $type_score;
						push(@{$gene_data->{$g->id}->{data}->{$type_score}},$variation->id);
						( $items->{tab_consequences}) = returnConsequences($variation,$items,$gene_data);
					#cover ?
					}
				#returnConsequences($variation,$items,$gene_data);
				} #end if gene 
				else {
					my $reference = $variation->getReference();
					 my $id = basic_regions_infos($gene_data,$chr,$variation->start);
					push(@ids_references,$id);
					push( @{ $items->{genes} }, $gene_data->{$id}->{name});
					
						unless ( exists $gene_data->{$id}){
							confess() unless $gene_data->{$id}->{start};
						}
					
				
					( $items->{tab_consequences}) = returnConsequences($variation,$items,$gene_data);
					#$gene_data->{$reference->id}->{intergenic} = 1;
				( $items->{tab_consequences}) = [];
				#	push(@{$gene_data->{$id}->{positions}->{$variation->start}},$variation->id);
					push(@{$gene_data->{$id}->{data}->{intergenic}},$variation->id);
				}
				
		
				
				
				###############################
				##### infos par patient
				#####################################""
					my $cover;
					my $nb_cover;
					
	
					
					### hetero or homo
				my $homo;
				my $hetero;
				my $patients;
				my $debug;
				foreach my $patient ((@{$variation->getPatients()})){
					push( @{ $items->{patient_id} },   $patient->id() );
					push( @{ $items->{patient_name} }, $patient->name() );	
					my $annex =  $variation->annex()->{$patient->id};
					
					$patients->{$patient->name()}->{homozygote} += $variation->annex()->{$patient->id}->{ho};
					$patients->{$patient->name()}->{heterozygote} = $variation->annex()->{$patient->id}->{he};
					$patients->{$patient->name()}->{score} = $annex->{score} if $annex->{score} > $patients->{$patient->name()}->{score};
					#warn Dumper  $patient->getCallingMethods();
					my ($m) =  $patient->getCallingMethods()->[0];
					#confess() if scalar @{$patient->getCallingMethods()} >1;
					my $hash_name = "patient_unifiedgenotyper";#.$m;
					
					push( @{$items->{$hash_name} },$annex->{score} );
					push( @{$items->{$hash_name."2"} },$annex->{nb_all_ref});
					push( @{$items->{$hash_name."3"} },$annex->{nb_all_mut});
					push( @{$items->{$hash_name."4"} },$annex->{dp});
					my $m = $patient->getCallingMethods()->[0];
						push( @{$items->{$hash_name."5"} },$variation->methodCalling($patient));
					#else {
						#push( @{$items->{$hash_name."5"} },"uni");
					#}
					my $text = "he";
					$text = "ho" if $annex->{ho} >0;
					push( @{$items->{$hash_name."heho"} },$text);
					$items->{homozygote} += $annex->{ho};
					$items->{heterozygote} += $annex->{he};
					$items->{$m} = int( $annex->{score} );
					my $base;
					$base = $annex->{var_allele} unless $base;
					push( @{ $items->{patient_base} },  $base );
				}
					
					$items->{homo_hetero} += 1 if ($items->{homozygote}>0);
					$items->{homo_hetero} += 2 if ($items->{heterozygote}>0);
					
				
				#####################################
				##### infos par genes ou ids
				#####################################
			
			 my $debug;
			
				foreach my $id (@ids_references) {
					unless ( exists $gene_data->{$id}){
							confess() unless $gene_data->{$id}->{start};
						}
					push(@{$gene_data->{$id}->{all_variations}},$variation->id);
					$gene_data->{$id}->{hash_all_variations}->{$variation->id} =undef ;
					push(@{$gene_data->{$id}->{data}->{homozygote}},$variation->id) if $variation_data->{$variation->id}->{homozygote} > 0;
					push(@{$gene_data->{$id}->{data}->{heterozygote}},$variation->id) if $variation_data->{$variation->id}->{heterozygote} > 0;
					#warn "ins " if  $variation->isInsertion();
					#warn "del " if $variation->isDeletion();
					push(@{$gene_data->{$id}->{data}->{insertion}},$variation->id) if $variation->isInsertion();
					push(@{$gene_data->{$id}->{data}->{deletion}},$variation->id) if $variation->isDeletion() && !($variation->isLargeDeletion());
					push(@{$gene_data->{$id}->{data}->{large_deletion}},$variation->id) if $variation->isLargeDeletion();
					push(@{$gene_data->{$id}->{data}->{substitution}},$variation->id) if $variation->isVariation();
					push(@{$gene_data->{$id}->{data}->{cosmic}},$variation->id) if $variation->isCosmic();
					push(@{$gene_data->{$id}->{data}->{notcosmic}},$variation->id) unless $variation->isCosmic();
					push(@{$gene_data->{$id}->{data}->{$vdb}},$variation->id);
				
					push(@{$gene_data->{$id}->{data}->{$db_freq}},$variation->id) if $db_freq;
		
						push(@{$gene_data->{$id}->{data}->{homozygote}},$variation->id) if ($items->{homozygote}>0);
						push(@{$gene_data->{$id}->{data}->{heterozygote}},$variation->id) if ($items->{heterozygote}>0);
					foreach my $p (@{$variation->getPatients()}) {
						my $patient_name = $p->name();
						my $pid = $p->id;
						my $nb_ref = $variation->annex()->{$pid}->{nb_all_ref};
						my $nb_mut = $variation->annex()->{$pid}->{nb_all_mut};
						my $type_heho = "ho";
						$type_heho = "he" if $variation->annex()->{$pid}->{he} == 1;
						 $hpatients->{$patient_name}->{$variation->id} = $type_heho.":".$nb_ref.":".$nb_mut;
					
					#	warn Dumper $variation->annex()->{$pid};
						
					
						#die();
						
						my $ho_name = "homozygote_".$patient_name; 
						my $he_name = "heterozygote_".$patient_name; 
						push(@{$gene_data->{$id}->{data}->{$ho_name}},$variation->id) if 	$patients->{$p->name()}->{homozygote} >0;
						push(@{$gene_data->{$id}->{data}->{$he_name}},$variation->id)  if 	$patients->{$p->name()}->{heterozygote} >0;	
						push(@{$gene_data->{patients}->{$p->name}->{homozygote}},$variation->id)  if 	$patients->{$p->name()}->{homozygote} >0;
						push(@{$gene_data->{patients}->{$p->name}->{heterozygote}},$variation->id)  if 	$patients->{$p->name()}->{heterozygote} >0;			
						push(@{$gene_data->{$id}->{data}->{$p->name}},$variation->id); 
					}
					
					
			
			} 
		
		
			$hvars->{$variation->id} = freeze $items;
		
			$variation->getProject->purge_memory($variation->end()) if ($nb_var %1000 == 0);
			
			delete $variation->getProject->{objects}->{$variation->getType}->{$variation->id};
			$variation = undef;
		} #end var

	}#end refrences
#	warn Dumper $gene_data->{reference_chr22_1_51304566};

	 my $z =0;
	 my $zz =0;
	 my %res;
	 my $info_patient = $gene_data->{patients};
	 my @types = ("evs","1000genomes","dbsnp","intronic","1000genomes_1p","evs_1p","pseudo","dbsnp_none","intergenic");
	 my @ct;
	 my %delete_genes;
	 for (my $i =2;$i<=@types;$i++){
	 	push(@ct,combine($i,@types));
	 }
   	foreach my $gene_id (keys %{$gene_data}){
   		next if $gene_id eq "patients";
   		unless (keys %{$gene_data->{$gene_id}}){
   			confess();
   		}
   	
   		$z ++;
		$hgenes->{$gene_id} = freeze $gene_data->{$gene_id};
		
		my $nb =0;
		
		foreach my $atypes (@ct){
		my $vhash = dclone $gene_data->{$gene_id}->{hash_all_variations};
		my $cst = join(";",sort{$a cmp $b} @$atypes);
		
		my $zz=0;
		foreach my $type (@$atypes){
			
			next unless exists $gene_data->{$gene_id}->{data}->{$type};
			my $ref = $gene_data->{$gene_id}->{data}->{$type};
		
			foreach my $v (@{$ref}){
			#	warn $v;
				delete $vhash->{$v};
					last if scalar(keys %$vhash) == 0;
			}
			last if scalar(keys %$vhash) == 0;
		}
		 $delete_genes{$cst}->{$gene_id} = undef  if scalar(keys %$vhash) > 0;
		}
		#$zz++;
		$gene_data->{$gene_id} = undef;
			
		} 

	foreach my $t (keys %delete_genes){
		my $nb = scalar (keys %{$delete_genes{$t}});
		
	 	#warn $t .":".scalar (keys $delete_genes{$t});
	 	my $zx = $delete_genes{$t};
		$fgenes->{$t} = freeze $zx;
		
	}
	

	$db1->close();
	
	$gene_data->{coverage} = $chr_data->{chromosomes}->{$chr->name}->{coverage};
	
	return ($gene_data,$data_variations);

}





sub returnConsequences {
	my ($variation,$items1,$hgenes) = @_;
#	return returnConsequencesForCnv($variation) if $variation->isCNV();
	my $cpt=0;
	my $array;
	my $varpos= $items1->{start};
	foreach my $tr (@{$variation->getTranscripts()}) {
			##################
			### nomenclature 
			#################
			
			
				my $xref = $tr->name;
		
		$items1->{ "nomenclature!" . $xref } = $variation->getNomenclature($tr);
		$items1->{ "consequence!" . $xref }  = $variation->variationType($tr);
		
		###
		#
		###
#		$items1->{"polyphen!$xref"} = -10;
#
#		if ($variation->isVariation() )
#		{
#			
#			my $prot = $tr->getProtein();
#			my $stp  = $variation->polyphenStatus($prot);
#			#$stp = 4 if $variation->changeAA($prot) eq "*";
#			$items1->{"polyphen!$xref"} =  + 5;
#		}	
			
			my %items;
			my $gene = $tr->getGene();
			my $prot = $tr->getProtein();
			
#			
#			#info gene
			$items{gene} = $gene->name() . " (" .$gene->external_name. ")";
			$items{description} = $tr->getGene()->description();# = $ensGene->description();
			$items{name}        = $variation->name() . "-" . $cpt;
#			#info transcripts							
			$items{transcript} = $tr->name()."+".$tr->{external_name} ;#= getExternalName( $tr, "refseq" );	
#			
			$items{cdna_position}  = $tr->translate_position($varpos); #=$variation->position( $prot->getTranscript() )->start;	
			my $nbe = $tr->findExonNumber($varpos);
		
			$items{exon}        = $nbe; 
			$items{exon} = $tr->findNearestExon($variation->start) if $items{exon} == -1;
			if ($prot){
				
				$items{protein} = $prot->name()."+".$tr->{external_protein_name};#getExternalName( $prot, "uniprot" );
			}
#			
#			# protein external name
#
			if ($variation->isCoding($tr)){
				$items{nomenclature} = $variation->getNomenclature($tr);
				$items{cds_position} = $variation->getOrfPosition($prot);
				$items{AA_variation} = $variation->changeAA($prot);
				$items{AA_protein}   = $variation->getProteinAA($prot);
				$items{protein_position}  = $items{AA_position}  = $variation->getProteinPosition($prot);
				
				#warn $variation->name()." ".$items{AA_position};
				$items{consequence}  = $items{AA_protein}."/".$items{AA_variation}." (".$variation->getCodons($tr).")";
				if ($variation->isInsertion){
					$items{consequence}  = "ins:".$variation->sequence();
				}
				if ($variation->isDeletion){
					$items{consequence}  = "del:".$variation->delete_sequence($tr);
				}
				my $polyphen = $variation->polyphenStatus($tr->getProtein);
				my $sift =  $variation->siftStatus($tr->getProtein);
				if ($variation->isStop($tr)){ 
					$sift=$polyphen = 5;
		 		}
		 		if ($variation->isPhase($tr)){ 
					$sift=$polyphen = 4;
		 		}
				$items{polyphen_status} = $polyphen."+".$variation->polyphenScore($tr->getProtein);
				$items{sift_status} = $sift."+".$variation->siftScore($tr->getProtein);
				
			}
			else {
				$items{consequence}  = $variation->variationType($tr);
			}
#			
#			$items{exon}        = ""; #getExonName( $tr, $variation );
			$items{polyphen_html} = "-";
			push( @$array, \%items );
	}#end for transcript

	return ( $array );
}



1; 
