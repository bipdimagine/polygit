package CacheDiag;
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Storable qw/thaw freeze/;
use Data::Dumper;
use BioTools;
use POSIX;
use Time::Piece;
 use List::Util qw( max min);

my $himpact_sorted = {
	"high" => "4",
	"moderate" =>"3",
	"low" =>"1",
};


my $kyoto_db;
my $type_db;
sub return_db {
	my ($patient) = @_;
		return $kyoto_db->{$patient->id} if exists $kyoto_db->{$patient->id};
	if ($patient->project->noSqlPolydiag()->exists_db($patient->name)){
		$kyoto_db->{$patient->id} = "sqlite";
		$type_db->{$patient->id} = "sqlite";
		return $kyoto_db->{$patient->id};
	}
	
	confess();

}

sub get_date {
		my ($patient) = @_;
		my $file_out;
		if ($patient->project->noSqlPolydiag()->exists_db($patient->name)){
			 $file_out = $patient->project->noSqlPolydiag()->dir."/".$patient->name.".".$patient->project->noSqlPolydiag()->extension;
		}
		else {
			confess();
			# $file_out = $patient->kyoto_polydiag_cache();
		}
		
		my $t = (stat $file_out )[9];
		my $t1 = localtime();
		my $date = POSIX::strftime( 
             "%d/%m/%y", 
             localtime( 
               		$t
                 )
             );
             
          
		my $days_difference = int((time - $t) / 86400);
			
           return ($date,$days_difference);
}


sub return_list_variants {
	my ($project,$patient,$tr_id) = @_;
	my $d = $project->getCacheDir();
	my $project_name= $project->name;
	my $patient_name = $patient->name();
	my $db = return_db($patient);
	my @vars;
	my $key =$tr_id;
	my $string ="";
	if ($type_db->{$patient->id} eq "sqlite"){
		$string = $project->noSqlPolydiag()->get($patient->name,"list_$key")."";
	}
	else {
	if ($db->check("list_$key") != -1 ){
		$string = $db->get("list_$key");
		
	}
	}
	#if ($string){
	return [split(";",$string)];
	#}
#	else{
#		foreach my $v (@{$patient->getStructuralVariations()}){
#      		foreach my $t (@{$v->getTranscripts()} ){
#      			next if $tr_id ne $t->id();
#      			warn $v->id;
#      			
#      			push(@vars,$v->id);
#      		}
#      	}
#      	# $db->set("list_$key",join(";",@vars));
#	}
	
	return \@vars;
	
}
sub return_list_all_transcripts {
	my ($project,$patient,$tr) = @_;
	my $d = $project->getCacheDir();
	my $project_name= $project->name;
	my $patient_name = $patient->name();
	my $db = return_db($patient);
	my $string = "";
	
		if ($type_db->{$patient->id} eq "sqlite"){
		$string = $project->noSqlPolydiag()->get($patient->name,"transcripts");
	}
	else {
	if ($db->check("transcripts") != -1 ){
		$string = $db->get("transcripts");
		
	}
	}
	if ($string){
	return [split(";",$string)];
	}
	else {
		confess();
	}
	
}

sub return_hash_variant {
	my ($project,$vid,$tr_id,$patient,$vquery) = @_;
	my $d = $project->getCacheDir();
	my $project_name= $project->name;
	my $patient_name = $patient->name();
		my $db = return_db($patient);
	my $id = join(";",$tr_id,$vid);
	#my $id = $patient->name()."_".$tr1->id."_".$v->id;
	my $z;
	if ($type_db->{$patient->id} eq "sqlite"){
		
		$z = $project->noSqlPolydiag()->get($patient->name,$id);
	}
	else {
	
	if ($db->check($id) != -1){
		$z =  thaw $db->get($id);
	}
	}
	unless ($z) {
		return undef;
		die($id."-");
#	
	}
	
	return $z;
}
sub construct_intergenic_variant {
	my ($project,$v,$patient,$vquery) = @_;

	my $hvariation;
	$hvariation->{id} = $v->id;


		
		$hvariation->{impact_text} = 0;
		$hvariation->{impact_score} = 0;
		$hvariation->{gene} = "intergenic";
		$hvariation->{var_name} = $v->name();
		if ($v->name() =~ /rs/){
			my $vn = $v->name();
			$hvariation->{var_name} = qq{<a href='http://www.ncbi.nlm.nih.gov/snp/?term=$vn' target = '_blank'>$vn</a> };	
		}
		
			my $sequence_info = "he("; 
		my $pc ="-";		
		if ($v->annex()->{$patient->id}->{nb_all_ref} eq "?"){
			$sequence_info = "??";
		}
		else {

		$sequence_info = "ho(" if $v->annex()->{$patient->id}->{ho};
		
		my $sum = $v->annex()->{$patient->id}->{nb_all_ref}+$v->annex()->{$patient->id}->{nb_all_mut};
		if ($sum >0){
		 $pc = int ($v->annex()->{$patient->id}->{nb_all_mut} *100/($sum));
		}
		$sequence_info .= $v->annex()->{$patient->id}->{nb_all_ref}."/".$v->annex()->{$patient->id}->{nb_all_mut}.")<br>";
	
		}
		
		 if ($v->validation_method eq "sanger" ) {
		 	$sequence_info = "-";
		 }
		$hvariation->{ngs} = $sequence_info;
		$hvariation->{ratio} = $pc."%";
		$hvariation->{allele} = $v->getSequence();
		$hvariation->{ref_allele} = $v->ref_allele();
		$hvariation->{genomique} = $v->getChromosome()->name.":".$v->start;
		$hvariation->{start} = $v->start;
		$hvariation->{transcript} = "-";
		$hvariation->{chromosome} = $v->getChromosome()->name;
		$hvariation->{trans} = $v->start ;
		$hvariation->{cds} = "";
		$hvariation->{prot} ="";
		$hvariation->{codons_AA} = "";
		$hvariation->{polyphen} ="-";
		$hvariation->{sift} ="-";

		$hvariation->{codons}  =  $v->getChromosome()->sequence($v->start,$v->end)."/".$v->sequence();
	
		
		my $start = $v->start;
		my $chr = $v->getChromosome()->name();
		my $vid = $hvariation->{id};
		my $a0 = $v->ref_allele;
		my $a1 = $v->var_allele;

		my $pname = $patient->name;
		
			my $qq4 = qq{	<div  data-dojo-type="dijit/Toolbar"><button dojoType="dijit.form.Button"   iconClass="igvIcon" onClick ="displayInIGV('$chr',$start,$start);"></button></div>};
			my $qq4 = qq{	<button  class="igvIcon2" onClick ="displayInIGV('$chr',$start,$start);">&nbsp;&nbsp;&nbsp;&nbsp;</button>};
			$hvariation->{igv} = $qq4; 		
			
			my $qq5 = qq{	<div  data-dojo-type="dijit/Toolbar"><button dojoType="dijit.form.Button"   iconClass="alamutView" onClick ="displayInAlamut('$chr',$start,['$a0','$a1']);"></button></div>};
			my $qq5 = qq{	<button    class="alamutView3" onClick ="displayInAlamut('$chr',$start,['$a0','$a1']);"></button>};
			$hvariation->{alamut} = $qq5; 
						
						my $qq3 = qq{
      				<div  data-dojo-type="dijit/Toolbar">
					<button dojoType="dijit.form.Button"   iconClass="dijitEditorIcon dijitEditorIconFullScreen" onClick = viewElectro('$pname','$vid')></button>
					</div>
      					};
      					my $qq3 = qq{
      				
					<button  class="alignIcon" onClick = viewElectro('$pname','$vid')></button>
			
      					};	
			$hvariation->{align} = $qq3; 	
		
		
		
		
		#warn $hvariation->{prot};
		$hvariation->{exon} = "-";
		$hvariation->{exon} = "-";
		$hvariation->{nomenclature} =  "-";
		$hvariation->{consequence} =  "-";
		
		$hvariation->{freq}  =  $v->frequency;
		$hvariation->{freq}  =  0 if $hvariation->{freq} == -999;
		#$hvariation->{freq} = $hvariation->{freq}/100;
		$hvariation->{scaled_score} = 1;
		$hvariation->{impact_score} =1;
		$hvariation->{scaled_score} =1;
		my $debug ;
		$debug =1 if $v->name eq "4_86916249_C_T";
		warn $v->frequency if $debug;
		if ($hvariation->{freq} < 0 ){
							$hvariation->{freq_score} = 0; 
							$hvariation->{freq_level} = 1; 
							$hvariation->{freq} = 0;
						
		}
		
		if ($hvariation->{freq} <= 0.01){
							$hvariation->{freq_score} = 4; 
							$hvariation->{freq_level} = 1; 
		}
		
		elsif ($hvariation->{freq} <= 0.05){
				$hvariation->{freq_score} = 2; 
					$hvariation->{freq_level} = 3; 
		}
		else {
					$hvariation->{freq_score} = 1; 
					$hvariation->{freq_level} = 4; 
			}
				
		$hvariation->{freq} = "-" if $v->origin_public_database() eq "new";
		
		return $hvariation;
}

sub construct_variant {
	my ($project,$v,$tr1,$patient,$vquery) = @_;
	my $hvariation;
	$hvariation->{id} = $v->id;
	if ($project->isSomatic){
		
		$hvariation->{cosmic} = $v->cosmic();
		my $id  = $v->cosmic();
		if ($id){
			my ($id1,$nb) = split(":",$id);
			my $text = $id1;
			my $url = "http://cancer.sanger.ac.uk/cosmic/mutation/overview?id="; 
			 $url = "http://grch37-cancer.sanger.ac.uk/cosmic/ncv/overview?id="  if ($id1 =~/COSN/);
			  $text = "$id1 [$nb]" if $nb;
			$id1 =~s/COSM//;
			$id1 =~s/COSN//;
		
			
			
			
			
			$hvariation->{cosmic} = qq{<a href=\"$url$id1\" target="_blank">$text</a>};
		}
	}
		
		$hvariation->{impact_text} = $v->effectImpact($tr1);
		$hvariation->{impact_score} = $himpact_sorted->{$v->effectImpact($tr1)};
		$hvariation->{gene} = $tr1->getGene->external_name();
		$hvariation->{var_name} = $v->name();
		if ($v->name() =~ /rs/){
			my $vn = $v->name();
			$hvariation->{var_name} = qq{<a href='http://www.ncbi.nlm.nih.gov/snp/?term=$vn' target = '_blank'>$vn</a> };	
		}
		my @asequence_info;
		my @apc;
		my @methods;
		my $nb_methods;
		foreach my $method (@{$patient->callingMethods}){
			next unless exists $v->annex()->{$patient->id}->{method_calling}->{$method}->{nb_all_ref};
			
			my $all_annex = $v->annex()->{$patient->id}->{method_calling}->{$method};
			my $nb_ref =$all_annex->{nb_all_ref};
			my $nb_alt =  $all_annex->{nb_all_mut};
			
		my $method_name = substr $method,0,3;
		push(@methods,$method_name);
		my $sequence_info = "he("; 
		my $pc ="-";		
		if ($v->annex()->{$patient->id}->{nb_all_ref} eq "?"){
			$sequence_info = "??";
		}
		else {
		$sequence_info = "ho(" if $all_annex->{ho};
		
		my $sum = $nb_ref + $nb_alt;
		if ($sum >0){
		 $pc = int ($nb_alt *100/($sum));
		}
		$sequence_info .= $nb_ref."/".$nb_alt.")";
	
		}
		$sequence_info = $method_name.":".$sequence_info;
		$pc = $method_name.":".$pc."%";
		push(@apc,$pc);
		push(@asequence_info,$sequence_info);
		$nb_methods ++;
		}
		
		 if ($v->validation_method eq "sanger" ) {
		 	#$sequence_info = "-";
		 	push(@asequence_info,"-");
		 }
		
	
		$hvariation->{ngs} = join("<br>",@asequence_info);
	
		$hvariation->{ratio} =  join("<br>",@apc);
		$hvariation->{caller} =  join("<br>",@methods);
		$hvariation->{allele} = $v->getSequence();
		$hvariation->{ref_allele} = $v->ref_allele();
		$hvariation->{genomique} = $v->getChromosome()->name.":".$v->start;
		$hvariation->{start} = $v->start;
		$hvariation->{transcript} = $tr1->name;
		$hvariation->{chromosome} = $v->getChromosome()->name;
		$hvariation->{trans} = $v->start * $tr1->strand; 
		
		$hvariation->{cds} = "";
		$hvariation->{prot} ="";
		$hvariation->{codons_AA} = "";
		$hvariation->{polyphen} ="-";
		$hvariation->{sift} ="-";

		my $debug;
		$debug = 1   if $v->name eq "rs113993960";
		#warn  $v->delete_sequence if $debug;
	
		if ($v->isDeletion){
	
			$hvariation->{codons}  =  $v->delete_sequence."/".$v->sequence();
		}
		else {
			$hvariation->{codons}  =  $v->getChromosome()->sequence($v->start,$v->end)."/".$v->sequence();
		}
		if ($tr1->strand() == -1 ){
			$hvariation->{codons}  =  BioTools::complement_sequence($v->getChromosome()->sequence($v->start,$v->end))."/".BioTools::complement_sequence($v->sequence());
			}
			
		my $start = $v->start;
		my $chr = $v->getChromosome()->name();
		my $vid = $hvariation->{id};
		my $a0 = $v->ref_allele;
		my $a1 = $v->var_allele;
	
		my $pname = $patient->name;
		
			my $qq4 = qq{	<div  data-dojo-type="dijit/Toolbar"><button dojoType="dijit.form.Button"   iconClass="igvIcon" onClick ="displayInIGV('$chr',$start,$start);"></button></div>};
			my $qq4 = qq{	<button  class="igvIcon2" onClick ="displayInIGV('$chr',$start,$start);">&nbsp;&nbsp;&nbsp;&nbsp;</button>};
			$hvariation->{igv} = $qq4; 		
			
			my $qq5 = qq{	<div  data-dojo-type="dijit/Toolbar"><button dojoType="dijit.form.Button"   iconClass="alamutView" onClick ="displayInAlamut('$chr',$start,['$a0','$a1']);"></button></div>};
			my $qq5 = qq{	<button    class="alamutView3" onClick ="displayInAlamut('$chr',$start,['$a0','$a1']);"></button>};
			$hvariation->{alamut} = $qq5; 
						
						my $qq3 = qq{
      				<div  data-dojo-type="dijit/Toolbar">
					<button dojoType="dijit.form.Button"   iconClass="dijitEditorIcon dijitEditorIconFullScreen" onClick = viewElectro('$pname','$vid')></button>
					</div>
      					};
      					my $qq3 = qq{
      				
					<button  class="alignIcon" onClick = viewElectro('$pname','$vid')></button>
			
      					};	
			$hvariation->{align} = $qq3; 	
		
		
		if ($v->isCoding($tr1)){
			my $prot = $tr1->getProtein();
		
			$hvariation->{cds} = $v->getOrfPosition($prot);
			$hvariation->{prot} = $v->getProteinPosition($prot);
			$hvariation->{codons_AA} =   $v->protein_nomenclature($prot);#$v->getProteinAA($prot).$hvariation->{prot}.$v->changeAA($prot);
			$hvariation->{polyphen} = $v->polyphenScore($tr1->getProtein);
			$hvariation->{sift} = $v->siftScore($tr1->getProtein);
			$hvariation->{codons} =   $v->getCodons($tr1);
		}
		#warn $hvariation->{codons} if $debug;
			#			die if $debug;
		#warn $hvariation->{prot};
		$hvariation->{exon} = $tr1->findExonNumber($v->start);
		$hvariation->{exon} = $tr1->findNearestExon($v->start) if $hvariation->{exon} == -1;
		$hvariation->{nomenclature} =  $v->getNomenclature($tr1);
		$hvariation->{consequence} =  $v->variationType($tr1);
		
		$hvariation->{freq}  =  $v->frequency;
		$hvariation->{freq}  =  0 if $hvariation->{freq} == -999;
		#$hvariation->{freq} = $hvariation->{freq}/100;
		$hvariation->{scaled_score} = $v->scaledScoreVariant($tr1,$patient,$vquery);
		if($nb_methods == 1 && $hvariation->{ngs} =~/dup/){
			$hvariation->{dup} = 1;
		}
		$hvariation->{score} = $v->scoreVariant($tr1,$patient,$vquery);
		
				
		if ($hvariation->{freq} < 0 ){
							$hvariation->{freq_score} = 4; 
							$hvariation->{freq_level} = 1; 
							$hvariation->{freq} = 0;
						
		}
		
		if ($hvariation->{freq} <= 0.01){
							$hvariation->{freq_score} = 4; 
							$hvariation->{freq_level} = 1; 
		}
		
		elsif ($hvariation->{freq} <= 0.05){
				$hvariation->{freq_score} = 2; 
					$hvariation->{freq_level} = 3; 
		}
		else {
					$hvariation->{freq_score} = 1; 
					$hvariation->{freq_level} = 4; 
			}
				
		$hvariation->{freq} = "-" if $v->origin_public_database() eq "new";
		return $hvariation;
}

sub update_variation_types {
	my ($hvariation) = @_;
	return if exists $hvariation->{db_type};
	my ($chr,$pos,$ref,$alt) = split("_",$hvariation->{id});
	my $lalt = length($hvariation->{allele});
	my $lref = length($hvariation->{ref_allele});
	$hvariation->{db_allele} =$hvariation->{allele};
	if ($hvariation->{allele} eq "-"){
		#deletion 
		$hvariation->{db_type} = "deletions";
		$hvariation->{db_allele} = $hvariation->{ref_allele};
		
	}
	elsif ($hvariation->{ref_allele} eq "-"){
		#deletion 
		$hvariation->{db_type} = "insertions";
	}
	elsif ($lalt eq $lref && $lref eq 1)
	{
		$hvariation->{db_type} = "snps";
	}
	else {
		warn Dumper $hvariation;
		die();
	}
}
sub update_populations {
		 my ($project,$hvariation) = @_;
		  return if exists $hvariation->{population};
		  	update_variation_types($hvariation);
		   $hvariation->{population} = 1;
		   
			my $res = $project->buffer->get_populations_frequencies($hvariation->{chromosome},$hvariation->{db_type},$hvariation->{start},$hvariation->{db_allele});
			my $min = 999;
			my $hmin;
			my $max = -1;
			my $hmax;
		   if ($res){
		   		foreach my $pop (keys %$res){
		   			unless ($pop eq 'ALL') { 
		   				$hvariation->{$pop} = $res->{$pop}->{F};
		   				if ($res->{$pop}->{F} ne "-" && $res->{$pop}->{F} < $min){
		   					
		   					$hmin->{pop} = $pop;
		   					$min = $res->{$pop}->{F} ;
		   					$hmin->{f} = $min;
		   				}
		   				if ($res->{$pop}->{F} > $max){
		   					$hmax->{pop} = $pop;
		   					$max = $res->{$pop}->{F} ;
		   					$hmax->{f} = $max;
		   				}
		   			}
		   			else {
		   				my $value = ($res->{$pop}->{HO}*2)/$res->{$pop}->{AN};
		   				$hvariation->{freq_ho} = sprintf("%.4f", $value );
		   			}
		   		}
		   		
		   		$hvariation->{max_pop} = $hmax->{pop}.":".$hmax->{f};
		   		$hvariation->{min_pop} = $hmin->{pop}.":".$hmin->{f};
		   }
		   
}
sub update_clinvar {
	 my ($project,$hvariation) = @_;
	 return if exists $hvariation->{clinvar};
	 update_variation_types($hvariation);
	 my $buffer = $project->buffer();
	 my $pub = $buffer->get_lmdb_database("clinvar",$hvariation->{chromosome},$hvariation->{db_type})->get_with_sequence($hvariation->{start},$hvariation->{db_allele});
	  $hvariation->{clinvar}  = "" ;
	     $hvariation->{clinvar_alert}  = 0 ;
	   if ($pub){
	   	my $v = max( split(";",$pub->{sig}));
	  	if ($v == 0){
	  		 $hvariation->{clinvar}  = "Uncertain"
	  	}
	  	elsif ($v ==1 ){$hvariation->{clinvar}  = "not provided"}
	  	elsif ($v ==2 ){$hvariation->{clinvar}  = "Benign"}
	  	elsif ($v ==3 ){$hvariation->{clinvar}  = "Likely benign"}
	  	elsif ($v ==4 ){$hvariation->{clinvar}  = "Likely pathogenic"; }
	  	elsif ($v ==5 ){$hvariation->{clinvar}  = " pathogenic";	}
	  	elsif ($v ==6 ){$hvariation->{clinvar}  = " histocompatibility";}
	  	else {$hvariation->{clinvar}  = "other";}
	 	#warn $hvariation->{clinvar};
	  	if ($v>3 && $v<6){
	  		$hvariation->{scaled_score} = 4;
	  		$hvariation->{clinvar_alert}++;
	  		
	  	}
	  	$hvariation->{clinvar}." ".$hvariation->{clinvar_alert};
	   }
}

sub update_cadd {
	 my ($project,$hvariation) = @_;
	 return if exists $hvariation->{cadd};
	 update_variation_types($hvariation);
	 my $buffer = $project->buffer();
	

	  $hvariation->{cadd}  = "-";
	  return  unless $hvariation->{db_type} eq "snps";
	my $caad =  $buffer->get_lmdb_cadd($hvariation->{chromosome}, $hvariation->{db_type})->get($hvariation->{start});
	 my @values = unpack("S4",$caad);
	my $ACGT = {'A'=>0,'C'=>1,'G'=>2,'T'=>3};
	$hvariation->{cadd} = $values[$ACGT->{$hvariation->{db_allele} }] /10;
		 
}
sub update_hotspot {
	 my ($project,$tr,$hvariation) = @_;
	 update_variation_types($hvariation);
	 my $buffer = $project->buffer();
	  $hvariation->{hs}  = "-";
	  
	  return unless $hvariation->{db_type} eq "snps";
		my $hs =  $buffer->get_lmdb_database("hotspot",$hvariation->{chromosome}, $hvariation->{db_type})->get_with_sequence($hvariation->{start},$hvariation->{db_allele});
		$hvariation->{hs}  = "*" if $hs;
}

sub update_deja_vu{
	 my ($project,$tr,$hvariation) = @_;
	#my $dejavu =$tr->getChromosome->hash_kyoto_dejavu;

	my $vid = $hvariation->{id};
	my $debug;
	$debug=1 if $vid eq "14_93670213_A_AT";
	
	warn $vid  if $debug;;

	my $similar = $project->similarProjects();
	my $hres = $project->getDejaVuInfosForDiag($vid);

	my $nb_pat =0;
	my $pname = $project->name();
	my $proj_dejavu =0;
	$hvariation->{sim_deja_vu} = $hres->{similar_patient};
	$hvariation->{sim_proj_deja_vu} = $hres->{similar_project};
	
	
		my $nb_dejavu = $project->getDejaVuThisProject($vid);
		die() if $nb_dejavu eq 0;
		$hvariation->{this_deja_vu} =$nb_dejavu;
		$hvariation->{diff_project_deja_vu} =$hres->{other_project}  ; 
		$hvariation->{project_deja_vu} = $hres->{other_project} + $hres->{similar_project}; 
		$hvariation->{deja_vu} =$hres->{other_project}.":".$hres->{other_patient}; 
		$hvariation->{in_this_run} = $nb_dejavu."/". scalar(@{$project->getPatients});
		$hvariation->{similar_projects} = $hres->{similar_project}.":". $hres->{similar_patient};
		
		my $freq_score;
		my $freq_level;
		
		if ($hvariation->{freq} < 0 ){
							$hvariation->{freq_score} = 4; 
							$hvariation->{freq_level} = 1; 
							$hvariation->{freq} = 0;
						
		}
		
		if ($hvariation->{freq} <= 0.01){
							$hvariation->{freq_score} = 4; 
							$hvariation->{freq_level} = 1; 
		}
		
		elsif ($hvariation->{freq} <= 0.05){
				$hvariation->{freq_score} = 2; 
					$hvariation->{freq_level} = 3; 
		}
		else {
					$hvariation->{freq_score} = 1; 
					$hvariation->{freq_level} = 4; 
			}
	
			
		if ($hvariation->{freq_score} == 4) {
				if ($hvariation->{diff_project_deja_vu} == 0){
					$freq_score = 4; 
					$freq_level =1;
					
				}
				elsif ($hvariation->{diff_project_deja_vu} <= 20 ){
							$freq_score = 3; 
							$freq_level =2;
					
				}
				elsif ($hvariation->{diff_project_deja_vu} <= 50){
						$freq_score = 2; 
						$freq_level =3;
					
				}
				else {
						$freq_score = 1;
						$freq_level =4; 
				}
		}	
		
		$hvariation->{freq_score} = $freq_score if $freq_score < $hvariation->{freq_score} && $freq_score;
		$hvariation->{freq_level} = $freq_level if $freq_level > $hvariation->{freq_level} && $freq_level ;
		
}

sub update_edit {
	 my ($patient,$hvariation) = @_;
	 my $lists = $patient->getListValidatedVariants();
	 
	 	$hvariation->{type} = "other";
	$hvariation->{sanger} = "-";
	#$hvariation->{user_name} = "";
	my $id =$hvariation->{id};

	return unless exists $lists->{$id};
	$hvariation->{edit} =1;
	my $v =  $lists->{$id};
	if ($v->{validation_sanger} == -5 ){
				$hvariation->{type} = "rejected";
				$hvariation->{sanger} = "rejected";
				#$hvariation->{type_confirmed_ngs} = "sanger";
			}
	elsif ($v->{validation_sanger} == 3 ){
			
			
				$hvariation->{type} = "confirmed";
				$hvariation->{sanger} = "confirmed (ho)"; 
	}
	elsif ($v->{validation_sanger} == 2 ){
		$hvariation->{type} = "confirmed";
	$hvariation->{sanger} = "confirmed (he)";
	}	
	
	else {
		 if ($v->{validation_ngs} == -3) {
			$hvariation->{type} = "todo";
			$hvariation->{sanger} = "todo";
			$hvariation->{type_confirmed} = "-";
			$hvariation->{type_confirmed_ngs} = "todo";
		}
		elsif ( $v->{validation_ngs} == -1){
			$hvariation->{type} = "rejected";
			$hvariation->{type_confirmed} = "ngs";
			$hvariation->{type_confirmed_ngs} = "fp";
		}
		elsif (  $v->{validation_ngs} == 1){
			$hvariation->{type} = "validated";
			$hvariation->{type_confirmed} = "ngs";
			$hvariation->{type_confirmed_ngs} = "ho";
		}
		elsif (  $v->{validation_ngs} == 2){
				$hvariation->{type} = "validated";
				$hvariation->{type_confirmed} = "ngs";
				$hvariation->{type_confirmed_ngs} = "he";
		}
   

		}
		
	 
}

1;