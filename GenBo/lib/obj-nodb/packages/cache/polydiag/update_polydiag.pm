package update;
use strict;
use FindBin qw($Bin);
use lib "$Bin";
use Storable qw/thaw freeze/;
use Data::Dumper;
use BioTools;
use POSIX;
use Time::Piece;
 use List::Util qw( max min sum);
  use Carp qw(confess croak);
 use Data::Printer;
 
 
 
 my $himpact_sorted = {
	"high" => "4",
	"moderate" =>"3",
	"low" =>"1",
};

sub printBadge {
	my ($value,$types) = @_;
	my $color = "#4CAF50";
	 #$color = "#D62C1A" if $value > $type[0] ;
	 $color = "#FF8800" if $value > $types->[0] ;
	 $color = "#FF0025" if $value > $types->[1] ;
	 return qq{<span class="badge badge-success badge-xs" style="border-color:$color;background-color:#FFFFFF;color:$color;font-size :8px;">$value</span>} ;

}

sub printSimpleBadge {
	my ($value) = @_;
	my $color = "black";
	 return qq{<span class="badge badge-success badge-xs" style="border-color:black;background-color:#FFFFFF;color:$color;font-size :8px;">$value</span>} ;
}

sub printInvBadge {
	my ($value,$types) = @_;
	
	my $color = "#4CAF50";
	 #$color = "#D62C1A" if $value > $type[0] ;
	 $color = "#FF8800" if $value < $types->[0] ;
	 $color = "#FF0025" if $value < $types->[1] ;
	 
	 $value ="-" unless defined $value;
	 $color = "#FF0025" if $value eq "-" ;
	
	 return qq{<span class="badge badge-success badge-xs" style="border-color:$color;background-color:#FFFFFF;color:$color;font-size :8px;">$value</span>} ;
}

sub printInvButton {
	my ($value,$types,$text,$othercss) = @_;
	my $btn_class  = qq{class= "btn btn-xs btn-primary " style="background-color: #D0D0D0;font-size: 7px;font-family:  Verdana;color:black"};
	$btn_class = qq{class= "btn btn-xs btn-primary " style="background-color: #FF8800;font-size: 7px;font-family:  Verdana;;color:black"} if $value <  $types->[0] ;
	$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: #e74c3c;font-size: 7px;font-family:  Verdana;color:white"}  if $value < $types->[1] ;
	 $value ="-" unless defined $value;
	$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: #e74c3c;font-size: 7px;font-family:  Verdana;"}  if $value eq "-" ;
	$text  = $value unless $text;
	return  qq{<button type="button" $btn_class $othercss>$text</button>};
}
sub printButton {
	my ($value,$types,$text,$othercss) = @_;

	my $btn_class  = qq{class= "btn btn-xs btn-primary " style="background-color: #D0D0D0;font-size: 7px;font-family:  Verdana;color:black"};
	$btn_class = qq{class= "btn btn-xs btn-primary " style="background-color: #FF8800;font-size: 7px;font-family:  Verdana;;color:white"} if $value >=  $types->[0] ;
	$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: #e74c3c;font-size: 7px;font-family:  Verdana;color:white"}  if $value >=  $types->[1] ;
	 $value ="-" unless defined $value;
	 
	$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: red;font-size: 7px;font-family:  Verdana;"}  if $value eq "-" ;
	$text  = $value unless $text;
	return  qq{<button type="button" $btn_class $othercss>$text</button>};
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
		$debug = 1   if $v->name eq "rs151344528";
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
		$hvariation->{exon} = $tr1->findExonNumber($v->start, $v->end);
		$hvariation->{exon} = $tr1->findNearestExon($v->start, $v->end) if $hvariation->{exon} == -1;
		$hvariation->{nomenclature} =  $v->getNomenclature($tr1);
		$hvariation->{consequence} =  $v->variationType($tr1);

		if ($v->isUpstream($tr1)){
			$hvariation->{consequence} = "upstream";
		}
		if ($v->isDownstream($tr1)){
			$hvariation->{consequence} = "downstream";
		}
		
		$hvariation->{freq}  =  $v->frequency;
		$hvariation->{freq}  =  0 if $hvariation->{freq} == -999;
		#$hvariation->{freq} = $hvariation->{freq}/100;
		$hvariation->{scaled_score} = $v->scaledScoreVariant($tr1,$patient,$vquery);
		if($nb_methods == 1 && $hvariation->{ngs} =~/dup/){
			$hvariation->{dup} = 1;
		}
		#$hvariation->{score} = $v->scoreVariant($tr1,$patient,$vquery);

		if ($hvariation->{freq} < 0 ){
							$hvariation->{freq_score} = 4; 
							$hvariation->{freq_level} = 1; 
							$hvariation->{freq} = 0;
						
		}
		
		elsif ($hvariation->{freq} <= 0.01){
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
		$hvariation->{freq} = "-" unless $v->frequency();
		my $scaled =  $v->scaled_score_frequence(); 
		$hvariation->{freq_level} = $scaled->{freq_level};
		$hvariation->{freq_score} = $scaled->{freq_score};
			
		#clinical local 
	 my $v1 = $v->score_clinical_local();
	# clinvar($project,$hvariation);
	  $hvariation->{clinical_local}  = "" ;
	 if ($v1){
	 	 	$hvariation->{clinical_local}  ++ ;
	 	 #	$hvariation->{clinvar_alert}  ++ ;
	 	 	my $cm = $v->comment_clinical_local();
	 	 	my $txt =  qq{<span class="badge badge-warning" style="font-size:8px">}.$hvariation->{clinvar}."</span>".qq{<span class="badge badge-warning"  onClick='alert ($cm)' style="font-size:8px">Local Clinical</span>};
	 	 	#$hvariation->{clinvar}  = " $txt ";
	 		$hvariation->{scaled_score} = 4 ;#f $hvariation->{freq_level} <= 2 ;	
	 }
			
		return $hvariation;
}

 
 
 
 
 
 sub variation_types2 {
	my ($project,$hvariation) = @_;
	return if exists $hvariation->{db_type};
	unless ($hvariation->{obj}){
	 	 $hvariation->{obj} = $project->_newVariant($hvariation->{id});
	 }
	
	$hvariation->{db_type} =  $hvariation->{obj}->type_public_db;
	$hvariation->{db_allele} =  $hvariation->{obj}->alternate_allele;

}

 sub variation_types {
	my ($project,$hvariation) = @_;

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
		warn $hvariation->{ref_allele}." lalt ".$lalt." lref ".$lref;
		warn $hvariation->{allele};

	confess($hvariation->{id});
	}
}


sub annotations{
	 my ($project,$hvariation) = @_;
	
	 populations($project,$hvariation);
	 
	 clinvar($project,$hvariation);
	 cadd($project,$hvariation);
	
	 ratio_score($project,$hvariation);
	 my $sc2 = $hvariation->{impact_score} + $hvariation->{freq_score};
	 $hvariation->{scaled_score} =4 if $sc2 >= 7 && $hvariation->{cadd} > 20;
	 delete $hvariation->{obj} if exists $hvariation->{obj};
}



sub ratio_score {
	 my ($project,$hvariation) = @_;


	 my @p = $hvariation->{ratio} =~ /(\d+)/g;
 	unless(@p){
 	 $hvariation->{ratio_score} = 1;
 	 return;
	 }
	my $v = sum(@p)/scalar(@p);
	 if ($v>35){
	 	 $hvariation->{ratio_score} = 4;
	 }
	 elsif ($v>25){
	 	 $hvariation->{ratio_score} = 3;
	 }
	 elsif ($v>15){
	 	 $hvariation->{ratio_score} = 2;
	 }
	  else {
	 	 $hvariation->{ratio_score} = 1;
	 }
}



sub populations {
		 my ($project,$hvariation) = @_;
		  return if exists $hvariation->{population};
		  	variation_types($project,$hvariation);
		   $hvariation->{population} = 1;
		    unless ($hvariation->{obj}){
	 
	 		 $hvariation->{obj} = $project->_newVariant($hvariation->{id});
	 #	confess();
	 	}
	 	my $v = $hvariation->{obj};
	 	$hvariation->{max_pop} ="-";
	 	$hvariation->{max_pop} = $v->max_pop_name.":".sprintf("%.4f", $v->max_pop_freq ) if $v->max_pop_name;
	 	$hvariation->{max_pop} = printSimpleBadge($hvariation->{max_pop});
	 	$hvariation->{min_pop} = "-";
		$hvariation->{min_pop} = $v->min_pop_name.":".sprintf("%.4f", $v->min_pop_freq ) if $v->min_pop_name;
			$hvariation->{min_pop} = printSimpleBadge($hvariation->{min_pop});
		$hvariation->{freq_ho} = "-"; 
		$hvariation->{freq_ho} = sprintf("%.4f", $v->frequency_homozygote ) if $v->frequency_homozygote ;
	 	$hvariation->{ac} = $v->getGnomadAC;
	 	$hvariation->{an} = $v->getGnomadAN; 
	 	$hvariation->{ac} = printInvButton($hvariation->{ac},[200,10]);
	 

		$hvariation->{an} = $v->getGnomadAN; 
			
		
		$hvariation->{an} = printInvBadge($hvariation->{an},[0,0]);
	
		 
		
		$hvariation->{ac_ho} = $v->getGnomadHO;  
		$hvariation->{ac_ho} = printInvButton($hvariation->{ac_ho},[50,5]);
		if  ($v->getChromosome->name eq "X" or $v->getChromosome->name eq "Y") {
			$hvariation->{ac_ho} = printInvButton($v->getGnomadHO,[50,5], $v->getGnomadHO." -  ".qq{&nbsp<i class="fa fa-mars" > </i> &nbsp;} .$v->getGnomadAC_Male);#printInvButton($v->getGnomadAC_Male,[50,5]) unless ($v->is_in_pseudoautosomal_region );  
		}
}

sub clinvar {
	 my ($project,$hvariation) = @_;
	
	 return if exists $hvariation->{clinvar};
	 unless ($hvariation->{obj}){
	 
	 	 $hvariation->{obj} = $project->_newVariant($hvariation->{id});
	 #	confess();
	 }
	 
	 
	 my $debug;
	 
	 
	 my $v = $hvariation->{obj}->score_clinvar();
	   $hvariation->{clinvar}  = "" ;
	     $hvariation->{clinvar_alert}  = 0 ;
	 if ($v){
	 	
	 	$hvariation->{clinvar}  = $hvariation->{obj}->text_clinvar();
	 	if (($v == 4 || $v==5)    ){
	 		$hvariation->{scaled_score} = 4 if $hvariation->{freq_level} <= 2 ;	
	 		$hvariation->{clinvar_alert}++;
	 	}
	 	
	 }
}

sub hgmd {
	 my ($project,$hvariation,$return_cmd) = @_;
	
	 #return if exists $hvariation->{hgmd};
	 unless ($hvariation->{obj}){
	 
	 	 $hvariation->{obj} = $project->_newVariant($hvariation->{id});
	 #	confess();
	 }
	 
	 
	 
	 
	 if ($hvariation->{obj}->hgmd_id()){
	 	if  ($hvariation->{obj}->isDM){
	 		$hvariation->{scaled_score} = 4;
	 		$hvariation->{clinvar_alert}++;
			 		
	 	}
	 	my $txt = $hvariation->{obj}->hgmd->{phen}." - ".$hvariation->{obj}->hgmd->{class}." - ".$hvariation->{obj}->hgmd_id;
	 	$txt =~s/\"//g;
	 	my $n1 = $project->name;
	 	my $n2 = $hvariation->{obj}->hgmd_id;
	 	my $n3 = $hvariation->{obj}->id;
	 	my $cmd = qq{zoomHgmd(\'$n1\',\'$n2\',\'$n3\')};# "zoomHgmd(\'.$project->name."\',\'"." ','".$hvariation->{obj}->id."')";
	 	return ($hvariation->{obj}->hgmd->{class}, $cmd) if ($return_cmd);
 		$hvariation->{hgmd}  =qq{<span class="badge badge-warning"  onClick="$cmd" style="font-size:8px">}.$hvariation->{obj}->hgmd->{class}."</span>";#$hvariation->{obj}->hgmd()->{hgmd_id}."<br>".$hvariation->{obj}->hgmd()->{class}."<br>".$hvariation->{obj}->hgmd()->{phen};
	 }
}
sub tclinical_local {
		my ($project,$hvariation,$patient,$gene) = @_;
		
	 my $val_id = $gene->id."!".$hvariation->{id};
	 my $local_validation = $patient->project->getValidationVariation($val_id,$patient);
	
		if ($local_validation){
				my $saved = $local_validation->{validation};
				$hvariation->{"local"} = $patient->buffer->value_validation->{$saved};
				$hvariation->{"local"} = printButton($saved,[3,5],$patient->buffer->value_validation->{$saved},$saved) ;
				$hvariation->{clinical_local} ++;
		}
		
}


sub clinical_local {
	 my ($project,$hvariation,$url) = @_;
	 return if exists $hvariation->{clinical_local};
	 unless ($hvariation->{obj}){
	my ($chr_name,$start,$ref,$alt) = split('_', $hvariation->{id});
		my $chr = $project->getChromosome($chr_name);
		
		my $var_obj = $chr->cache_lmdb_variations->get($hvariation->{id});	 
		if ($var_obj){
			 $hvariation->{obj} = $project->_newVariant($hvariation->{id});
		}
		else {
			return;
		}
	 	
	 #	confess();
	 }
	 
	 
	 my $debug;
	 my $v = $hvariation->{obj}->score_clinical_local();
	   $hvariation->{clinical_local}  = "" ;
	 if ($v){
	 	 	$hvariation->{clinical_local}  ++ ;
	 	 	$hvariation->{clinvar_alert}  ++ ;
	 	 	my $cm = $hvariation->{obj}->comment_clinical_local();
	 	 	my $txt =  qq{<span class="badge badge-warning" style="font-size:8px">}.$hvariation->{clinvar}."</span>".qq{<span class="badge badge-warning"  onClick='alert (\"$cm\")' style="font-size:8px">Local Clinical</span>};
	 	 	$hvariation->{clinvar}  = " $txt ";
	 		$hvariation->{scaled_score} = 4 ;#f $hvariation->{freq_level} <= 2 ;	
	 	
	 }
}
sub cadd {
	 my ($project,$hvariation) = @_;
	 return if exists $hvariation->{cadd};
	  unless ($hvariation->{obj}){
	 	 $hvariation->{obj} = $project->_newVariant($hvariation->{id});
	 }
	 $hvariation->{cadd} = $hvariation->{obj}->cadd_score();
}

sub trio {
	 my ($project,$tr,$hvariation,$patient,$cgi,$print) = @_;
	# return $hvariation->{trio} if exists $hvariation->{trio};
	 my $fam = $patient->getFamily();
	
	 unless ($fam->isTrio){
	 	 $hvariation->{trio} = "-";
	 	 return;
	 }
	# ENST00000301067 
	 # return if exists $hvariation->{trio};
	  $hvariation->{trio} = 0 unless $project->isFamilial();
	  
	  return $hvariation->{trio} unless $project->isFamilial();

	     unless ($hvariation->{obj}){
	 	 $hvariation->{obj} = $project->_newVariant($hvariation->{id});
	 }
	  my $children;
	  my $father;
	  my $mother;
	  my $all_tab;
	  my $var = $hvariation->{obj};
	  my $fam = $patient->getFamily();
	  my $debug ;
	  $debug =1 if $var->id eq '9_334048_A_ACTTCTGG';
	  
	 
	#  if ($patient->isChild){
	  	$children = $fam->getChildren();
	  	$mother = $fam->getMother;
	  	$father = $fam->getFather;
	  	warn $mother->name.' '.$father->name if $debug;
	  	
#	  	unless ($father ){
#	  		$hvariation->{trio} ="-";
#	  		return ;
#	  	}
#	  	unless ($mother){
#	  		$hvariation->{trio} ="-";
#	  		return ;
#	  	}
	  	my $html =$cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"font-size: 8px;font-family:  Verdana;margin-bottom:0px"});
		my $html_line;
		my $type = "-";
		if ($mother){
			$type = "ho" if $var->isHomozygote($mother);
			$type = "he" if $var->isHeterozygote($mother);
			my $style = "";
		
				if ($var->getPourcentAllele($mother) eq "-"){
					my $ps = "-";
					$ps ="0" if return_coverage($mother,$var->getChromosome->ucsc_name,$var->start)>5;
					
					my $tab = [$mother->name,qq{<i class="fa fa-female fa-2x" aria-hidden="true" style="color:ligthgrey"></i>},$type,$ps."%",""];
					$html .=  $cgi->start_Tr({style=>"border: 1px solid black;color:lightgrey;$style"}).$cgi->td([$mother->name,qq{<i class="fa fa-female fa-2x" aria-hidden="true" style="color:ligthgrey"></i>},$type,$ps."%",""]).$cgi->end_Tr();
				#splice @$tab, 1, 1;
					$html_line = join(":",@$tab);
					$html_line .=" - &nbsp;";
			}
			else {
				my $tab = [$mother->name,qq{<i class="fa fa-female fa-2x" aria-hidden="true" style="color:pink"></i>},$type,$var->getPourcentAllele($mother)."%",""];
				$html .=  $cgi->start_Tr({style=>"border: 1px solid black;color:black;$style"}).$cgi->td($tab).$cgi->end_Tr();
				#splice @$tab, 1, 1;
				$html_line .= join(":",@$tab);
				$html_line .=" - &nbsp;";
			}
		}
		if ($father){
			my $style = "";
		#$style = "background-color:#F7C9C9" ;
		$type = "-";
		$type = "ho" if $var->isHomozygote($father);
		$type = "he" if $var->isHeterozygote($father);
		#$style = "background-color:#779ECB" ;
		if ($var->getPourcentAllele($father) eq "-"){
			my $ps = "-";
			$ps ="0" if return_coverage($father,$var->getChromosome->ucsc_name,$var->start)>5;
			my $tab = [$father->name,qq{<i class="fa fa-male fa-2x" aria-hidden="true" style="color:lightgrey"></i>},$type,$ps."%",""];
			$html .=  $cgi->start_Tr({style=>"border: 1px solid black;color:lightgrey;$style"}).$cgi->td($tab).$cgi->end_Tr();
			push(@$all_tab,@$tab);
			#splice @$tab, 1, 1;
			$html_line .= join(":",@$tab);
			$html_line .="- &nbsp;";
	
		}
		else {
			my $tab = [$father->name,qq{<i class="fa fa-male fa-2x" aria-hidden="true" style="color:blue"></i>},$type,$var->getPourcentAllele($father)."%",""];
			$html .=  $cgi->start_Tr({style=>"border: 1px solid black;color:black;$style"}).$cgi->td($tab).$cgi->end_Tr();
			push(@$all_tab,@$tab);
			#splice @$tab, 1, 1;
			$html_line .= join(":",@$tab);
			$html_line .=" - &nbsp;";
		}
		}
		foreach my $child (@$children) {
			$fam = $child->getFamily;
		 $type = "-";
		$type = "ho" if $var->isHomozygote($child);
		$type = "he" if $var->isHeterozygote($child);
		warn $type.' '.$child->name if $debug;
			if ($type eq "-"){
				my $tab = [$child->name,qq{<i class="fa fa-child fa-2x" aria-hidden="true" style="color:black"></i> },$type,"-","-"];
				$html .=  $cgi->start_Tr({class=>"table-info",style=>"border: 1px solid black;color:black;background-color:lightgrey"}).$cgi->td({class=>"table-primary"},$tab).$cgi->end_Tr();
					push(@$all_tab,@$tab);
					#splice @$tab, 1, 1;
				
				$html_line .= join(":",@$tab);
				$html_line .=" - &nbsp;";
				next;
			}
			
		my $model = $var->getTransmissionModel($fam,$child,$debug);
		warn $mother->name.' '.$father->name if $debug;
		warn $fam.' --> '.$child->name if $debug;
		warn $model if $debug;
		
		if ($model =~ /denovo/ && scalar (keys%{$fam->parents()}) ==1){
			$model = "denovo/?";
		}
		$hvariation->{transmission_model} = lc ($model);
		$hvariation->{transmission_model} = "xor" if $hvariation->{transmission_model} eq "mother" or $hvariation->{transmission_model} eq "father";#lc ($model);
		my $color = "white";
		my $model2 = qq{<i class="fa fa-male  fa-2x" style="color:lightgrey"></i><i class="fa fa-female  fa-2x" style="color:lightgrey"></i>};
		$color = "#779ECB" if (lc($model) eq "father");
		$model2 = qq{<i class="fa fa-male  fa-2x style="color:$color""></i>} if (lc($model) eq "father");
		
		$color = "#F7C9C9" if (lc($model) eq "mother");
		$model2 = qq{<i class="fa fa-female  fa-2x style="color:$color""></i>} if (lc($model) eq "mother");
		$color = "red" if (lc($model) =~ "denovo");
		$model2 = qq{Denovo} if (lc($model) =~ "denovo");
		$model2 = qq{Strict Denovo} if (lc($model) =~ "strict");
		$color = "#E55137" if (lc($model) eq "denovo/?");	
		$model2 = qq{Denovo/?} if (lc($model) eq "denovo/?");
		$color = "violet" if (lc($model) eq "recessive");
		$model2 = qq{Recessive} if (lc($model) eq "recessive");
		$color =  "#45B8AC" if (lc($model) eq "uniparental disomy");
		$model2 = qq{UniParental disomy} if (lc($model) eq "uniparental disomy");
		$color =  "#F9885C" if (lc($model) =~ "mosaic");
		$model2 = qq{mosaic} if (lc($model) =~ "mosaic");
		my $tab = [$child->name,qq{<i class="fa fa-child fa-2x" aria-hidden="true" style="color:black"></i> },$type,$var->getPourcentAllele($child)."%",$model2];
	
	
		
		$html .=  $cgi->start_Tr({class=>"table-info",style=>"border: 1px solid black;color:black;background-color:$color"}).$cgi->td({class=>"table-primary"},$tab).$cgi->end_Tr();
		push(@$all_tab,@$tab);
		#splice @$tab, 1, 1;
		$html_line .= join(":",@$tab);
		$html_line .=" - &nbsp;";
#		my $model = $var->getTransmissionModel($fam,$child);
#		$hvariation->{transmission_model} = lc ($model);
#		my $color = "white";
#		$color = "#779ECB" if (lc($model) eq "father");
#		$color = "#F7C9C9" if (lc($model) eq "mother");
#		$color = "red" if (lc($model) eq "denovo");
#		$color = "red" if (lc($model) eq "recessive");
#		$html .= $cgi->start_Tr({style=>"border: 1px solid black;background-color:$color;align:center"}).$cgi->th({colspan=>"4",style=>"font-size: 8px;font-family:  Verdana;margin-bottom:0px"},[$model]).$cgi->end_Tr();
		}
		$html.= $cgi->end_table();
	#	warn $html;
	#	die();
		$hvariation->{trio} = $html;
		
		if ($print ==1){
#		$hvariation->{trio} = $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"font-size: 8px;font-family:  Verdana;margin-bottom:0px"});
#		$hvariation->{trio} .=$cgi->start_Tr();
#		$hvariation->{trio} .=$cgi->td($all_tab);
#		$hvariation->{trio} .=$cgi->end_Tr();
#		$hvariation->{trio} .= $cgi->end_table();

		$hvariation->{trio} =$html_line;
		}
		# die() if $debug;	
	 # }
	
	 
}

sub return_coverage {
	my ($p,$chr,$start) = @_;
	my $tabix = $p->tabix_coverage();
	 my $res = $tabix->query_full($chr,$start,$start+1);# if $start;
	 	my($a,$p,$c) = split(" ",$tabix->read($res));
	 	return $c;
}

sub hotspot {
	 my ($project,$tr,$hvariation) = @_;
	 variation_types($project,$hvariation);
	 my $buffer = $project->buffer();
	  $hvariation->{hs}  = "-";
	  
	  return unless $hvariation->{db_type} eq "snps";
		my $hs =  $buffer->get_lmdb_database("hotspot",$hvariation->{chromosome}, $hvariation->{db_type})->get_with_sequence($hvariation->{start},$hvariation->{db_allele});
		$hvariation->{hs}  = "*" if $hs;
}


sub deja_vu{
	 my ($project,$tr,$hvariation,$debug) = @_;
	 
	my $vid = $hvariation->{id};
	my $debug;
	
	my $similar = $project->similarProjects();
	my $hres = $project->getDejaVuInfosForDiag2($vid);

	my $nb_pat =0;
	my $pname = $project->name();
	my $proj_dejavu =0;
	$hvariation->{sim_deja_vu} = $hres->{similar_patients};
	$hvariation->{sim_proj_deja_vu} = $hres->{similar_projects};
	
		my ($nb_dejavu, $nb_ho) = $project->getDejaVuThisProject($vid);
		die() if $nb_dejavu eq 0;
		$hvariation->{this_deja_vu} =$nb_dejavu;
		$hvariation->{diff_project_deja_vu} =$hres->{other_projects}  ; 
		$hvariation->{project_deja_vu} = $hres->{other_projects} + $hres->{similar_projects}; 
		$hvariation->{deja_vu} =$hres->{other_projects}.":".$hres->{other_patients}; 
		$hvariation->{in_this_run} = $hres->{in_this_run_patients}."/". $hres->{total_in_this_run_patients};#scalar(@{$project->getPatients});
		$hvariation->{similar_projects} = $hres->{similar_projects}.":". $hres->{similar_patients};
		
		my $freq_score;
		my $freq_level;
		
		if ($hvariation->{freq} < 0 ){
							$hvariation->{freq_score} = 4; 
							$hvariation->{freq_level} = 1; 
							$hvariation->{freq} = 0;
						
		}
		unless ($hvariation->{freq}) {
				$hvariation->{freq_score} = 4; 
			$hvariation->{freq_level} = 1; 
		}
		elsif ($hvariation->{freq} <= 0.01){
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
	
			
	if ($hvariation->{freq_score} == 4 or $hvariation->{freq_score} == 2) {
				if ($hvariation->{diff_project_deja_vu} == 0){
					$freq_score = 4; 
					$freq_level =1;
					$hvariation->{dejavu_score} = 4; 
					
				}
				elsif ($hvariation->{diff_project_deja_vu} <= $project->buffer->{config}->{dejavu}->{rare}  ){
							$freq_score = 3; 
							$freq_level =2;
							$hvariation->{dejavu_score} = 3; 
					
				}
				elsif ($hvariation->{diff_project_deja_vu} <=  $project->buffer->{config}->{dejavu}->{occasional}){
						$freq_score = 2; 
						$freq_level =3;
						$hvariation->{dejavu_score} = 2; 
					
				}
				else {
					$hvariation->{dejavu_score} = 1; 
						$freq_score = 1;
						$freq_level =4; 
				}
					$hvariation->{freq_score} = $freq_score;
					$hvariation->{freq_level} = $freq_level;
		}	
		#$hvariation->{freq_score} = $freq_score if $freq_score < $hvariation->{freq_score} && $freq_score;
		#$hvariation->{freq_level} = $freq_level if $freq_level > $hvariation->{freq_level} && $freq_level ;
		if ($debug){
			
		}
		
}

sub edit {
	 my ($patient,$hvariation) = @_;
	 my $lists = $patient->getListValidatedVariants();
	 	$hvariation->{type} = "other";
	$hvariation->{sanger} = "-";
	#$hvariation->{user_name} = "";
	my $id =$hvariation->{id};

	return unless exists $lists->{$id};
	
	my $v =  $lists->{$id};
	$hvariation->{always_keep} = 1;
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
 
