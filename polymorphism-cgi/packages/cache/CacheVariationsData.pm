package CacheVariationsData;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo";
use lib "$Bin/../../GenBo/lib/GenBoDB";
use lib "$Bin/../../GenBo/lib/obj-lite";
use lib "$Bin";

use Data::Dumper;
use Parallel::ForkManager;
use GenBoStorable;




sub getVariationInfo {
	my ( $project,$var, $trs, $debug ) = @_;
	
	my %items;
	
	my $chromosome = $var->getChromosome();
	my $contig     = $var->getContig();
	my $genes      = $var->getGenes();
	
	if (scalar(@$genes) >0){
	foreach my $g (@$genes) {
		push( @{ $items{genes} }, $g->name );
		$items{$g->name()."_consequence"} = $var->variationType($g);
		#$items{$g->name()."_polyphen"} = $var->polyphenStatus($g);
	
		#$items{genes} .= $g->name().";";
		}
	}
	else {
		push( @{ $items{genes} }, $var->getReference()->name );
	}
	$items{project} = $project->id();
	$items{structural_type} = $var->getStructuralType();
	$items{chromosome}  = $var->getChromosome()->name();
	$items{type}  = $var->type->name();
	$items{id}    = $var->id();
	$items{start} = $var->position($chromosome)->start();
	$items{end} = $var->position($chromosome)->end();
	$items{text}  = "" ;
	
	#$items{text}  = $var->getReference->sequence($var)."/".$var->sequence() unless $var->isCNV;
	
	


	
	$items{"consequence!all"} = $var->variationType();
	my $vdb = $var->origin_database();
	$items{public} = $vdb;
	
	
	if ($var->isCNV()) {
		
		if ($var->isEnsembl()) {
			$items{public} = "dgv";
		
			$items{linkdgv}=$var->{cnv_public};
		} 
	}
	
	
	$items{valid}     = $var->valid;
	$items{reference} = $var->getRootReference()->name();

	### Star score
	$items{filter} ;#= $var->filterScore();

	### Score methods
	


	



	$items{name}     = $var->name();
	$items{nb_patients} = scalar( @{$var->getPatients()});


	$items{polyphen_status} = "-99";

	$items{allpatients} = join( ";",
		map    { $_->name }
		  sort { $a->name() cmp $b->name() } @{ $var->getPatients() } );
		  
	#returnAnnexes( $project,$var, \%items );
	$items{homo_hetero} += 1 if ($items{homozygote}>0);
	$items{homo_hetero} += 2 if ($items{heterozygote}>0);
	( $items{tab_consequences}, $items{polyphen_status} ) = returnConsequences($var,\%items);

	
	return \%items;
}


#annex pour le tableaude droite celui des traces 

sub returnAnnexes {
	my ( $project,$var, $items ) = @_;
	my $data2;
	
	foreach my $trace (@{ $var->getTraces() } ){
		my $patient = $trace->getPatient();
		push( @{ $items->{patient_id} },   $patient->id() );
		push( @{ $items->{patient_name} }, $patient->name() );
		my $base;
		 foreach my $m (@{$project->getCallingMethods($var->type->name())}) {
				my $annex = $var->getAnnex( $trace, $m );
				my $score = 0;	
				$score = $annex->{score} if $annex;
				my $hash_name = "patient_".$m;
				push( @{$items->{$hash_name} },$score );
				push( @{$items->{$hash_name."2"} },$annex->{score2});
				push( @{$items->{$hash_name."3"} },$annex->{score3});
				push( @{$items->{$hash_name."4"} },$annex->{score4});
				$items->{homozygote} += $annex->{ho};
				$items->{heterozygote} += $annex->{he};
				$items->{$m} = int( $annex->{score} );
				$base = $annex->{text} unless $base;
				
		} 	
		push( @{ $items->{patient_base} },  $base );
		

	}

	return $data2;

}


sub returnAnnexes2 {
	my ( $project,$var, $items ) = @_;
	my $data2;

	foreach my $patient ( sort { $a->name() cmp $b->name() } @{ $var->getPatients() } )
	{
		my $patientName = $patient->name();
		warn $patient->name();
		my ($trace) = grep { $_->getPatient($patientName) } @{ $var->getTraces() };
		
		my $obj;
		
		$obj=$trace;
		
		push( @{ $items->{patient_id} },   $patient->id() );
		push( @{ $items->{patient_name} }, $patient->name() );
		my $base;
		 foreach my $m (@{$project->getCallingMethods($var->type->name())}) {
				my $annex = $var->getAnnex( $obj, $m );
				my $score = 0;	
				$score = $annex->{score} if $annex;
				my $hash_name = "patient_".$m;
				push( @{$items->{$hash_name} },$score );
				push( @{$items->{$hash_name."2"} },$annex->{score2});
				push( @{$items->{$hash_name."3"} },$annex->{score3});
				push( @{$items->{$hash_name."4"} },$annex->{score4});
				$base = $annex->{text} unless $base;
				
		} 
		
		push( @{ $items->{patient_base} },  $base );
		unless ($var->getProject->is_ngs){
		my %items2;


			my $alltraces = $var->getAllTracesByPatient($patient);

			next if scalar(@$alltraces) == 0;
			my $htraces;
			my $ntraces;
			foreach my $t (@$alltraces) {
				push( @$htraces, $t->id );
				push( @$ntraces, $t->name );

			}
			push( @{ $items->{traces_id} },   $htraces );
			push( @{ $items->{traces_name} }, $ntraces );

			my $annex = $var->getAnnex( $alltraces->[0] )->[0];

			$items2{nb_traces} = '-';
			if ($annex) {
				$items2{nb_traces} =	  $var->getAnnex( $alltraces->[0] )->[0]->score;
			}
			$items2{nb_all_traces} = scalar(@$alltraces);
		}

	}

	return $data2;

}



#fonction utilise pour afficher les  consquences d'une variation

sub returnConsequencesForCnv {
	my ($variation) = @_;
	my $cpt   = 77;
	my %global;
	my $array;
	foreach my $tr ( @{ $variation->getTranscripts() } ) {
		$cpt++;
	 	my $prot = $tr->getProtein();
	 	my $gene = $tr->getGene();
	 	my %items;
	 	my $ensTranscrit = $tr->getEnsemblObject();
		
		my $ensGene      = $gene->getEnsemblObject();
		$items{gene} = $gene->name() . " (" . $ensGene->external_name() . ")";
		$items{description} = $ensGene->description();
		$items{name}        = $variation->name() . "-" . $cpt;
		$items{transcript} = getExternalName( $tr, "refseq" );
		$items{exon}        = getExonName( $tr, $variation );
		
	 	if ($prot && $variation->isCoding($prot)){
	 		
			my $consequence = $variation->consequence($prot);

			

		die( $variation->name() ) if exists $global{ $variation->name() . " " . $cpt };
		$global{ $variation->name() . " " . $cpt }++;

		# transcript external name

	

		# protein external name

		$items{protein} = getExternalName( $prot, "uniprot" );

		
		$items{consequence} = ".";

		#$items{sequence}    =  $prot->sequence();
		my $sequence = ">" . $prot->name() . "\n" . $prot->sequence();
		
		
		$items{protein_position} = "-";

		$items{cdna_position}    = $variation->position( $prot->getTranscript() )->start;
		$items{AA_variation} = "XX";
		$items{AA_protein}   ="X";
		$items{AA_position}  = "-";
		
		$items{consequence}  = "- / -";   

		$items{polyphen} = "-";
		
		$items{polyphen_status} = "-";
		$items{polyphen_html} = "-";	 		
	 	}
	 	else {
		$items{cdna_position} = $variation->position($tr)->start;
		
		$items{protein}       = "";
		$items{consequence}   = "UTR";;
		$items{polyphen}      = "-";
		$items{polyphen_html} = "-";
		
	 		
	 	}
	 	push( @$array, \%items );
	}
	return ( $array, -1 );
}
sub returnConsequences {
	my ($variation,$items1) = @_;
	
	return returnConsequencesForCnv($variation) if $variation->isCNV();
	my $cpt=0;
	my $array;
	my $varpos= $variation->position($variation->getChromosome())->start;
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
		$items1->{"polyphen!$xref"} = -10;

		if ($variation->isVariation() )
		{
			
			my $prot = $tr->getProtein();
			my $stp  = $variation->polyphenStatus($prot);
			#$stp = 4 if $variation->changeAA($prot) eq "*";
			$items1->{"polyphen!$xref"} =  + 5;
		}	
			
			my %items;
			my $gene = $tr->getGene();
			my $prot = $tr->getProtein();
			
#			
			my $ensGene      = $gene->getEnsemblObject();	
#			#info gene
			$items{gene} = $gene->name() . " (" . $ensGene->external_name() . ")";
			$items{description} = $tr->{hash}->{gene_description};# = $ensGene->description();
			$items{name}        = $variation->name() . "-" . $cpt;
#			#info transcripts							
			$items{transcript} = $tr->name()."+".$tr->{hash}->{external_name} ;#= getExternalName( $tr, "refseq" );	
#			
			$items{cdna_position}  = $tr->translate_position($varpos); #=$variation->position( $prot->getTranscript() )->start;	
			if ($prot){
				#my $ensProtein   = $prot->getEnsemblObject();
				$items{protein} = $prot->name()."+".$tr->{hash}->{external_protein_name};#getExternalName( $prot, "uniprot" );
			}
#			
#			# protein external name
#
			if ($variation->isCoding($tr)){
				$items{AA_variation} = $variation->changeAA($prot);
				$items{AA_protein}   = $variation->getProteinAA($prot);
				$items{AA_position}  = $variation->getProteinPosition($prot);
				#warn $variation->name()." ".$items{AA_position};
				$items{consequence}  = $items{AA_protein}."/".$items{AA_variation};
			
				$items{polyphen} = "-";
			}
			else {
				$items{consequence}  = $variation->variationType($tr);
			}
#			
#			$items{exon}        = ""; #getExonName( $tr, $variation );
			$items{polyphen_html} = "-";
			push( @$array, \%items );
	}#end for transcript
	my $max_polyphen_status = -99;
	return ( $array, $max_polyphen_status );
}

sub returnConsequences2 {
	my ($variation) = @_;
	return returnConsequencesForCnv($variation) if $variation->isCNV();
	my %dejavu;
	my $stop = undef;
	my $text;

	my $gene;
	my $consequence;
	my $cpt   = 77;
	my $myurl ="";
	#$cgi->url( -path_info => 1 )."?type=polyphenhtml&projectName=".$project_name . "&";
	my $debug;
	my $array;
	my %global;
	my $max_polyphen_status = -99;
	#my $exons               = $variation->getExons();
	my $num_exon            = 0;
	
	foreach my $prot ( @{ $variation->getProteins() } ) {
	
		my $exon;
		my $transcrit = $prot->getTranscript();

		$cpt++;
		my %items;
		$dejavu{ $prot->getTranscript()->name() }++;

		$gene = $prot->getGene();

		my $ensTranscrit = $transcrit->getEnsemblObject();
		my $ensProtein   = $prot->getEnsemblObject();
		my $ensGene      = $gene->getEnsemblObject();

		$consequence = $variation->consequence($prot);

		$items{gene} = $gene->name() . " (" . $ensGene->external_name() . ")";
		$items{description} = $ensGene->description();
		$items{name}        = $variation->name() . "-" . $cpt;

		die( $variation->name() )
		  if exists $global{ $variation->name() . " " . $cpt };
		$global{ $variation->name() . " " . $cpt }++;

		# transcript external name

		$items{transcript} ;#= getExternalName( $transcrit, "refseq" );

		# protein external name

		$items{protein};# = getExternalName( $prot, "uniprot" );

		$items{exon}        = "";#getExonName( $transcrit, $variation );
		$items{consequence} = ".";

		#$items{sequence}    =  $prot->sequence();
		#my $sequence = ">" . $prot->name() . "\n" . $prot->sequence();
		
		
		$items{protein_position} = $variation->position($prot)->start;
		$items{cdna_position}    = $variation->position( $prot->getTranscript() )->start;
		$items{AA_variation} = $variation->changeAA($prot);
		$items{AA_protein}   = $prot->sequence($variation);
		$items{AA_position}  = $variation->position($prot)->start();
		$items{consequence}  = 
		    $prot->sequence($variation) . " / "
		  . $variation->changeAA($prot);    #$consequence;

		$items{polyphen} = "-";

	
		#my $stp = $variation->polyphenStatus($prot);
		#$stp = 4 if $variation->changeAA($prot) eq "*";

		#$items{polyphen_status} = $stp + 5;
		#$max_polyphen_status = $items{polyphen_status}
		 # if $max_polyphen_status < $items{polyphen_status};
		$items{polyphen_html} = "-";

		push( @$array, \%items );
	}
	foreach my $tr ( @{ $variation->getTranscripts() } ) {
		my %items;
		next if exists $dejavu{ $tr->name() };

		$cpt++;
		$items{name} = $variation->name() . "-" . $cpt;
		die() if exists $global{ $variation->name() . "-" . $cpt };
		$global{ $variation->name() . " " . $cpt }++;
		$consequence = "UTR";

		$gene = $tr->getGene();
		my $ensTranscrit = $tr->getEnsemblObject();
		my $ensGene      = $gene->getEnsemblObject();
		$items{gene} = $gene->name() . " (" . $ensGene->external_name() . ")";
		$items{description}   = $ensGene->description();
		$items{cdna_position} = $variation->position($tr)->start;
		$items{transcript}    = getExternalName( $tr, "refseq" );
		$items{exon}          = getExonName( $tr, $variation );
		$items{protein}       = "";
		$items{consequence}   = $consequence;
		$items{polyphen}      = "-";
		$items{polyphen_html} = "-";
		push( @$array, \%items );
	}

	return ( $array, $max_polyphen_status );

}
 
sub getExternalName {

	my ( $obj, $type ) = @_;

	my $ensObj = $obj->getEnsemblObject();
	my $name;
	my (@ext_reftr) =
	  grep { $_->dbname() =~ /$type/i } @{ $ensObj->get_all_DBEntries() };

	$name = $obj->name();
	if (@ext_reftr) {
		$name .= "+" . join( ";", map { $_->display_id() } @ext_reftr );
	}
	return $name;

}


sub getExonName {
	my ( $tr, $variation ) = @_;
	return "";
	my $exons = $variation->getExons();

	my $exon;
	my $num_exon = 0;
	foreach
	  my $e ( sort { $a->position($tr)->start <=> $b->position($tr)->start }
		@{ $tr->getExons() } )
	{
		$num_exon++;
		($exon) = grep { $e->id eq $_->id } @$exons;
		last if $exon;
	}

	return $num_exon . " " . $exon->name();
}
1;
