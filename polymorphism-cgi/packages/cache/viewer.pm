package viewer;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "..";
use strict;
use Storable qw/freeze thaw nfreeze store_fd fd_retrieve retrieve/;
use JSON::XS;

sub referenceData {
	my ($buffer,$project,$type_name) = @_; 
	my %data;
	#my $r = $project->getReference($referenceName);
	my $wid = GenBoStorable::insertWaiting($buffer->dbh,$project->id,$type_name);
	my $nb_total =0;
	warn scalar(@{$project->getStructuralVariations()});
	foreach my $ref (@{$project->getReferences()}){
		$nb_total += scalar(@{$ref->getTranscripts()});
	}
	my $is_ngs  = ($project->projectType()->name() eq "ngs");
	my $cover_span = Set::IntSpan::Fast::XS->new();
	my @positions;
	#creation d'une reference virtuelle sur la reference passée en paramètre (correspondant au gène selectionné dans l'interface) 
	 #my $virtualRef = $r->getVirtualReference($selectGene);
	my $nbv =0;
	my $step = int($nb_total * 1/75) + 1; 
	foreach my $ref (@{$project->getReferences()}) {
	#my $ref = $reference;
	my $chr = $ref->getChromosome(); 
	$data{name} = $chr->name();
	$data{length} = $chr->length();
	#permet de remplacer dans le buffer les variations connues dans ensembl par le nom ensembl
	
	#pour chacun des transcrits de la reference virtuelle
		


	foreach my $tr (@{$ref->getTranscripts()}){
			$nbv ++;
			#warn int(($nbv/$nb_total)*100);
		GenBoStorable::updateWaiting($buffer->dbh,$wid,int(($nbv/$nb_total)*75)) if $nbv%$step ==0;
		#recupération de la proteine codée par le transcrit
		#next if  $tr->position($chr)->start() <0;
		
		#next if  $tr->position($chr)->end() >$ref->length() ;
		my $prot = $tr->getProtein();
	
		#stockage des infos sur le transcrit dans la table de hash %transcript
		my %transcript;
		$transcript{name} =  $tr->getXrefName("refseq");#$tr->name();
		$transcript{chromosome} = $tr->getChromosome()->name();
		$transcript{gene_name} = $tr->getGene()->getEnsemblObject()->external_name();
		$transcript{x} = $tr->position($chr)->start();
	
		$transcript{length} = $tr->position($chr)->length();
		#pour chacun des exons des transcrits
		 foreach my $ex (@{$tr->getExons()}){
		 	my %exon;
		 	#si la proteine existe
		 	if (defined $prot){
		 	#on recupère la position des exons sur la reference virtuelle
		 	my $posprot = $prot->position($chr);
		 	#si la position de debut de la proteine sur la reference virtuelle est sup à la position de debut de l'exon sur la ref virtuelle
		 	#et si la position de debut de la proteine est inf à la position de fin de l'exon sur la ref virtuelle
		 	#alors on definit la position de l'UTR
		 	if ($posprot->start() > $ex->position($chr)->start() &&  $posprot->start() <  $ex->position($chr)->end()){
		 		$exon{utrstart} =  $ex->position($chr)->start();
		 		$exon{utrend} = $posprot->start() ;
		 	}
		 	#si la position de fin de la proteine sur la reference virtuelle est sup à la position de debut de l'exon sur la ref virtuelle
		 	#et si la position de fin de la proteine est inf à la position de fin de l'exon sur la ref virtuelle
		 	#alors on definit la position de l'utr
		 	elsif ($posprot->end() > $ex->position($chr)->start() &&  $posprot->end() <  $ex->position($ref)->end()){
		 		$exon{utrstart} =  $posprot->end() ;
		 		$exon{utrend} = $ex->position($chr)->end();
		 	}
		 	#autre cas de figure pour definir la position de l'utr
		 	elsif ($posprot->start() > $ex->position($chr)->end()){
		 		$exon{utrstart} =  $ex->position($chr)->start();
		 		$exon{utrend} = $ex->position($chr)->end();
		 	}
		 	#autre cas de figure pour definir la position de l'utr
		 	elsif ($posprot->end() < $ex->position($chr)->start()){
		 		$exon{utrstart} =  $ex->position($chr)->start();
		 		$exon{utrend} = $ex->position($chr)->end();
		 	}
		 	}
	#stockage des informations sur les xons dans la table de hash %exon			
			$exon{name} = $ex->name;
			$exon{x} = $ex->position($chr)->start();
			$exon{length} = $ex->position($chr)->length();
			$exon{strand} = $ex->position($chr)->strand();
			#on ajoute les infos sur les exons à la table %transcrit du transcrit possèdant ces exons
			push(@{$transcript{exon}},\%exon);
		}
		#on ajoute à la table %data, les infos sur les transcrits
		push(@{$data{transcript}},\%transcript);
	}
	}
	warn "\t end transcript";
	#fin de la partie permettant de tracer la partie supérieure de la vue graphique cad la partie dont les données proviennent d'ensembl

	#on recupère tous les patients du projets
	my (@patients) = sort{$a->name() cmp $b->name()} @{$project->getPatients()};
	my %variations;
	#pour chacun des patients
	$nb_total = scalar(@{$project->getStructuralVariations()});
	
	$nbv =0;
	$step = int($nb_total * 1/25) + 1; 
	foreach my $p (@patients) {
		#creation d'une table de hash %patient destinée à stocker toutes les infos pour un patient
		my %patient;
		$patient{name} = $p->name();
	
	#on recupère toutes las variations du patients
	#my $contigs = $ref->getContigs();
	
	
	foreach my $vr (@{$p->getStructuralVariations()}){
		
		#my ($true) = grep {$vr->getContig->name() eq $_->getContig->name()} @$contigs;
		$nbv++;
		GenBoStorable::updateWaiting($buffer->dbh,$wid,int(($nbv/$nb_total)*25)+75) if $nbv%$step ==0;
	
		my $ref = $vr->getChromosome();
		#next unless $true;
		next if (!$is_ngs && !($vr->isGoodBipd())) ;
		my %variation;
		#pour chaque variation  si elle existe dejà on associe la table %variation au patient et on passe à la variation suivante
		if (exists $variations{$vr->id}){
			push(@{$patient{variation}},$variations{$vr->id}); 
			next;
		}
		#si la table %variation n'existe pas encore pour cette variation alors on la crée.
		#cette table stcoke toutes les informations concernant la variation du patient
		$variation{chromosome} = $vr->getChromosome()->name();	
		$variation{name} = $vr->name();	
		$variation{type} = $vr->variationType();	
		$variation{x} = $vr->position($ref->getChromosome())->start();
	
		$variation{start} = $vr->position($ref->getChromosome())->start();;
		$variation{end} = $vr->position($ref->getChromosome())->end();
		$variation{length} =  $vr->position($ref->getChromosome())->length();
		$variation{polyphred} = $vr->isPolyphred();

		$variation{bipd} = $vr->isBipd();
		$variation{ensembl} = $vr->isEnsembl();
		
		$variation{consequence} = $vr->allConsequences();	 
		$variation{contig} = $vr->getContig()->name();
		$variation{reference} = $vr->getRootReference()->name();
		$variation{id} = $vr->id;
		$variation{validate} = $vr->valid();
		$variation{filterscore} = $vr->filterScore();
		if ( $vr->isCoding() ) {
			$variation{polyphen_status} = getPolyphenImage( $vr->polyphenStatus() );
			
		}
		else {
			$variation{polyphen_status} = qq{/icons/polyicons/cancel.png};
		}
		$variations{$vr->id} = \%variation;
		#on associe la table %variation au patient correspondant
		push(@{$patient{variation}},\%variation); 
	}
	
	#pour chacune des traces du patient
#	my $dir =  "/temporary/test-maq/";
#	my $patientName = $p->name();
#	my $dir = util_file::get_align_dir({project=>$project,method=>"maq"});

	#my $filename = $dir.$patientName."/cover-storable_20_bis.txt";
#	my $filename =  $dir."/".$p->name.".20.cover";

	if ($is_ngs){# && -e $filename) {
		#foreach my $ref (@{$project->getReferences()}){
			my $ref;
			#my $cover2 = getCoverForNGS($p,\%patient,$ref,$project);#\@positions) ;

			#$cover_span = $cover_span->union($cover2);
		#}
		
	}
	else {
	my $span_traces_fwd; 
	my $span_traces_rev;
	warn "traces";
	foreach my $tr (@{$p->getTraces()}){
		#on recupère laposition de la trace sur la reference virtuelle
		my $chrs = $tr->getChromosomes();
		next unless scalar(@$chrs);
		#warn scalar(@$chrs);
		my $chr = $tr->getChromosome();
		my $pos = $tr->position($chr);
		my $chr_name = $chr->name();
		$span_traces_fwd->{$chr_name} = Set::IntSpan::Fast::XS->new() unless exists $span_traces_fwd->{$chr_name};
		$span_traces_rev->{$chr_name} = Set::IntSpan::Fast::XS->new() unless exists $span_traces_rev->{$chr_name};
		$data{span_cover}->{$chr_name} =  Set::IntSpan::Fast::XS->new() unless exists $data{span_cover}->{$chr_name};
		$data{span_cover}->{$chr_name}->add_range( $pos->start(), $pos->end());
		
		if ($pos->strand == 1){
			$span_traces_fwd->{$chr_name}->add_range($pos->start,$pos->end);
		}
		else {
			$span_traces_rev->{$chr_name}->add_range($pos->start,$pos->end);
		}
	#	die() if $pos->start()==109954441;
#		
#		#si la trace ne s'aligne pas on passe à la suivante
#		next if $pos->start() == -1;
#		my %trace ;
#		#on stocke toutes les infos concernant la trace dans la table %trace
#		$trace{name} = $tr->name();
#		$trace{chromosome} = $tr->getChromosome()->name();
#		my $chr_name =  $tr->getChromosome()->name();
#		$trace{start} = $pos->start();
#		$trace{end} = $pos->end();
#		$trace{x} = $pos->start();
#		$trace{y} = $pos->end();
#		$trace{strand} = $pos->strand();
#		#on associe la table %trace au patient
#		push(@{$patient{trace}},\%trace); 
#		my $trace_span =  Set::IntSpan::Fast::XS->new($pos->start()."-".$pos->end());
#		$cover_span = $cover_span->union($trace_span);
#	
#		#on stocke les positions de la trace sur la reference virtuelle dans un tableau
#		push (@positions,$pos);
	}
	#on associe la table %patient à la table %data
	$patient{span_traces_fwd} = $span_traces_fwd;
	$patient{span_traces_rev} = $span_traces_rev;
	}
	
	push(@{$data{patient}},\%patient);
	}
	
#	foreach my $pos (split(",",$cover_span->as_string)){
#		my %contig;
#		my ($start,$end) = split("-",$pos);
#		warn $start." ".$end;
#		$contig{name} = "XX".$start;
#		$contig{x} = $start;
#		$contig{length} = ($end-$start)+1;
#		push(@{$data{contig}},\%contig);
#		
#	}
	
	#renvoit les données au format json

GenBoStorable::insertStorable($buffer->dbh,$project->id,$project->id,$type_name,\%data);
GenBoStorable::deleteWaiting($buffer->dbh,$wid) if $wid;
return (GenBoStorable::getStoreId($buffer->dbh,$project->id,$project->id,$type_name));
}

#fonction qui permet de renvoyer une image en fonction du statut polyphen de la variation
sub getPolyphenImage {
	my ($val) = @_;
	my $path = "/icons/myicons4/weather2/";

	#my $path = "/icons/myicons4/nuvola/16x16/actions/";
	my $img = "no.png";
	$img = "04.png" if ( $val == 2 );
	$img = "28.png" if ( $val == 1 );
	$img = "32.png" if ( $val == 0 );
	$img = "NA.png" if ( $val == -2 );
	$img = "stop_hunabkuc_software.png"  if ( $val == 4 );
	my $img = $path . "/" . $img;
	$img = "/icons/Polyicons/cancel.png" if ( $val == -5 );
	return $img;

}

 
sub getCoverForNGS {
	my ($patient,$hpatient,$ref,$project) = @_;
	
	 my $outputfile = util_file::get_cover_file(
		{
			project      => $project,
			patient_name => $patient->name,
			coverage     => "all",
		});
		return unless -e $outputfile;
		my $data = retrieve($outputfile);  
		
	my $cover2 = Set::IntSpan::Fast::XS->new();
#my $dir =  "/temporary/test-maq/";
my $patientName = $patient->name();



 	foreach my $pos (split(",",$data->as_string)){
		my %trace;
		my ($start,$end) = split("-",$pos);
		
		$trace{name} = $start;
		$trace{x} = $start;
		$trace{y} =$end;
		$trace{strand} = 1;
		push(@{$hpatient->{trace}},\%trace); 
 	}
 
 return $data;
 
}


1;