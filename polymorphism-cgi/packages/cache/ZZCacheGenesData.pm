package CacheGenesData;
use strict;
use FindBin qw($Bin);
use Data::Dumper;
use Parallel::ForkManager;
 use Bio::DB::Sam;
#use PBS::Client;
use Storable qw(store retrieve freeze);
use Time::Duration;


sub basic_gene_infos {
	my ($g,$coverage) = @_;

	my $pos = $g->position($g->getChromosome);
	my %gene_data;
	$gene_data{start}  = $pos->start();
	$gene_data{end}    = $pos->end();
	$gene_data{chromosome}    = $g->getChromosome->name();
	$gene_data{strand} = $pos->strand();
	$gene_data{name}   = $g->name();
	$gene_data{id}   = $g->id();
	$gene_data{cover} = 0;
	if ($g->isGene()){
		$gene_data{gene} = 1;
		$gene_data{reference} = join(";", map {$_->name()} @{$g->getReferences});    
		              
		$gene_data{xref}      = $g->getEnsemblObject->external_name();
		$gene_data{description} = $g->getEnsemblObject->description();
		# @{$tr->get_all_DBEntries()};
	
#		foreach my $link ( @{ $g->getEnsemblObject->get_all_DBLinks } ) {
#  			if ( $link->database eq "GO" ) {
#   				 my $term_id = $link->display_id;
#    			 my $term_name = '-';
#   				 my $term = $goa->fetch_by_accession($term_id);
#  				  if($term and $term->name){
#   					   $term_name = $term->name;
#   					 } 
#   			warn $g->name.": $term_id ($term_name)\n";

    #fetch complete GO hierachy
   
	
		
		$gene_data{transcripts} = join (";",map {$_->name} @{$g->getTranscripts()} );
		$gene_data{cover} = int($coverage->{total} / $coverage->{len}) if exists $coverage->{len};
	}
	
	return \%gene_data;
}
 

sub create_cache_genes {	
	my ( $project,$chr,$hvars) = @_;
	$| =1;
	my $buffer = $chr->buffer();

	
	my @array;
	

	
	
#	warn $nb_total;
	my $nb_references = 0;
	
	my $chr_data;
	my $gene_coverage;
	my $test_coding = 0;
	
	
	

	
	#my @refs_ids = map {$_->name} @$refs; 
	
	my $gene_data;
	my $variation_data;
	my $sc_waiting;
	
	my $scale_pourcent;
	my $data_variations;
	my $limit_pourcent = 10;
	my $references = $chr->getReferences();
	#my $vars  = $chr->getStructuralVariations();
	
	my $tstart = time;
	my $t1 = time;
	
	my @references_name = map{$_->name()} @$references;
	my $chr_name= $chr->name();
	my $project_name = $project->name();
	my $nn = $references->[0]->name();
	my $ref = $chr->getReference($nn);
	#my $vars  = $chr->getStructuralVariations();
#	my $nb_ref_total = scalar(@$vars);
	my $nb_ref_total = scalar(@$references);

foreach my $reference ( @$references) {

	
		  $nb_references ++;
		     my $pourcent = int(($nb_references / $nb_ref_total)*100 );
		     if ($pourcent >0 && ($pourcent % $limit_pourcent) == 0){
		     ##	$buffer->test_flush();
		     	my $ttime = (time - $tstart);
		     	
		     	my $rem = ($ttime * 100) / $pourcent;
		     	my $still = 1- $pourcent/100;
		     	warn "\tchr :".$chr->name()." $pourcent\% remaining time : ".concise(duration($rem*$still))."\n";
		     	#warn "\t".$chr->name()."  ".$pourcent."% \t".(time - $t1)."\t".(time - $tstart)."\n";
		     	$t1 = time;
		     	#die();
		     	$limit_pourcent +=20;
		     }
	
	my $vars = $reference->getStructuralVariations();
	foreach my $variation (@$vars){
		#next if $variation->id ne 9924659;
			
			my @trs;
			#warn "id ".$variation->id;
		 #@trs = sort { $a->position($reference)->start() <=> $b->position($reference)->start() } @{ $reference->getTranscripts() };
		
		   
		    # print "$nb_references/$nb_ref_total\r";
		   
		
		
		 	my $items;
		 	my $pos1 = $variation->position($chr);
		 
		 	$items->{structural_type} = $variation->getStructuralType();
			$items->{chromosome}  = $chr->name();
			$items->{type}  = $variation->type->name();
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
		
			my $vdb = $variation->origin_database();
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
			
			
			if ( $project->projectType()->name() eq "classic" ) {
				next unless $variation->isGoodBipd();

			}
				### public or not
				
			
				#die();
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
						
						$gene_data->{$g->id} = basic_gene_infos($g,$gene_coverage->{$g->name}) unless exists $gene_data->{$g->id};
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
					
					#cover ?
					}
					( $items->{tab_consequences}) = returnConsequences($variation,$items,$gene_data);
				} #end if gene 
				else {
					my $reference = $variation->getReference();
					push(@ids_references,$reference->id);
					push( @{ $items->{genes} }, $reference->name );
					$gene_data->{$reference->id} = basic_gene_infos($reference) unless exists $gene_data->{$reference->id};		
			
					push(@{$gene_data->{$reference->id}->{data}->{intergenic}},$variation->id);
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
					foreach my $trace (@{$variation->getTraces()}) {
						my $patient = $trace->getPatient();
						
						push( @{ $items->{patient_id} },   $patient->id() );
						push( @{ $items->{patient_name} }, $patient->name() );	
						my $base;
						foreach my $m (@{$project->getCallingMethods()}) {
							my $annex = $variation->getAnnex($trace,$m);
							#warn $trace->{id}." ".$m." ".$variation->id();
							#die() unless $annex;	
							#my $annex = $trace->getAnnex( $variation, $m );
							
							#next unless $annex;	
							
							
						#	warn  $patient->name()." ".$annex->{homozygote}." ".$annex->{heterozygote}." ".$annex->score();
							$patients->{$patient->name()}->{homozygote} += $annex->{ho};
							$patients->{$patient->name()}->{heterozygote} += $annex->{he};
							$patients->{$patient->name()}->{score} = $annex->{score} if $annex->{score} > $patients->{$patient->name()}->{score};
							$homo += $annex->{ho};
							$hetero += $annex->{he};
							$cover = $annex->{score} if $annex->{score} > $cover;
							
							my $hash_name = "patient_".$m;
							push( @{$items->{$hash_name} },$annex->{score} );
							push( @{$items->{$hash_name."2"} },$annex->{score2});
							push( @{$items->{$hash_name."3"} },$annex->{score3});
							push( @{$items->{$hash_name."4"} },$annex->{score4});
							my $text = "he";
							$text = "ho" if $annex->{ho} >0;
							push( @{$items->{$hash_name."heho"} },$text);
							$items->{homozygote} += $annex->{ho};
							$items->{heterozygote} += $annex->{he};
							$items->{$m} = int( $annex->{score} );
							$base = $annex->{text} unless $base;
							
						}
						push( @{ $items->{patient_base} },  $base );
					
					}
					
						$items->{homo_hetero} += 1 if ($items->{homozygote}>0);
						$items->{homo_hetero} += 2 if ($items->{heterozygote}>0);
					#	$items->{polyphen_status} = "-99";
				
				#####################################
				##### infos par genes ou ids
				#####################################
				
				foreach my $id (@ids_references) {
					push(@{$gene_data->{$id}->{all_variations}},$variation->id);
					$gene_data->{$id}->{hash_all_variations}->{$variation->id} =undef ;
					push(@{$gene_data->{$id}->{data}->{homozygote}},$variation->id) if $variation_data->{$variation->id}->{homozygote} > 0;
					push(@{$gene_data->{$id}->{data}->{heterozygote}},$variation->id) if $variation_data->{$variation->id}->{heterozygote} > 0;
					push(@{$gene_data->{$id}->{data}->{insertion}},$variation->id) if $variation->isInsertion();
					push(@{$gene_data->{$id}->{data}->{deletion}},$variation->id) if $variation->isDeletion();
					push(@{$gene_data->{$id}->{data}->{substitution}},$variation->id) if $variation->isVariation();
					
				
					push(@{$gene_data->{$id}->{data}->{$vdb}},$variation->id);
					### he/ho by patients
					
					
					
				#	die();
					
				
		
					#
					#die();	

						push(@{$gene_data->{$id}->{data}->{homozygote}},$variation->id) if ($items->{homozygote}>0);
						push(@{$gene_data->{$id}->{data}->{heterozygote}},$variation->id) if ($items->{heterozygote}>0);
					foreach my $p (@{$variation->getPatients()}) {
						my $patient_name = $p->name();
						my $ho_name = "homozygote_".$patient_name; 
						my $he_name = "heterozygote_".$patient_name; 
						push(@{$gene_data->{$id}->{data}->{$ho_name}},$variation->id) if 	$patients->{$p->name()}->{homozygote} >0;
						push(@{$gene_data->{$id}->{data}->{$he_name}},$variation->id)  if 	$patients->{$p->name()}->{heterozygote} >0;	
						push(@{$gene_data->{patients}->{$p->name}->{homozygote}},$variation->id)  if 	$patients->{$p->name()}->{homozygote} >0;
						push(@{$gene_data->{patients}->{$p->name}->{heterozygote}},$variation->id)  if 	$patients->{$p->name()}->{heterozygote} >0;			
						push(@{$gene_data->{$id}->{data}->{$p->name}},$variation->id); 
					}
					
					
				#}
			
				# $data_variations->{$variation->id} = CacheVariationsData::getVariationInfo($project, $variation, \@trs );
			
			} 
			$hvars->{$variation->id} = freeze $items;
			#$project->buffer->flushMemoryForReference();
		} #end var
		
		 $buffer->flushMemoryForReference();
	}#end refrences
		$gene_data->{coverage} = $chr_data->{chromosomes}->{$chr->name}->{coverage};
		
		#$buffer->flushMemoryForChromosome();
		#return $gene_data;
	#	die();
		#warn " ========== > ".(time - $tstart)."";
		return ($gene_data,$data_variations);

}



sub returnConsequences {
	my ($variation,$items1,$hgenes) = @_;
	
	return returnConsequencesForCnv($variation) if $variation->isCNV();
	my $cpt=0;
	my $array;
	my $varpos= $items1->{start};
		foreach my $tr (@{$variation->getTranscripts()}){
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
			$items{gene} = $gene->name() . " (" .$hgenes->{$gene->id}->{xref}. ")";
			$items{description} = $tr->{hash}->{gene_description};# = $ensGene->description();
			$items{name}        = $variation->name() . "-" . $cpt;
#			#info transcripts							
			$items{transcript} = $tr->name()."+".$tr->{hash}->{external_name} ;#= getExternalName( $tr, "refseq" );	
#			
			$items{cdna_position}  = $tr->translate_position($varpos); #=$variation->position( $prot->getTranscript() )->start;	
			my $nbe = $tr->findExonNumber($varpos);
		
			$items{exon}        = $nbe; #getExonName( $tr, $variation );
			
			if ($prot){
				
				$items{protein} = $prot->name()."+".$tr->{hash}->{external_protein_name};#getExternalName( $prot, "uniprot" );
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
			#	$items{polyphen_score} = $variation->polyphenScore($tr->getProtein);
			#	$items{sift_score} = $variation->polyphenScore($tr->getProtein);
			#	$items{polyphen_all}= $items{polyphen_status}."_".$items{polyphen_score};
			#	$items{sift_all}= $items{sift_status}."_".$items{sift_score};
				
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
