package update;
use strict;
use FindBin qw($Bin);
use lib "$Bin";
use Storable qw/thaw freeze/;
#use BioTools;
use POSIX;
use Time::Piece;
 use List::Util qw( max min sum);
  use Carp qw(confess croak);
 use Data::Printer;
 use JSON::XS;
use Data::Dumper;
 our $fts = "8px";
 our $print =0;
 
 my $himpact_sorted = {
	"high" => "4",
	"moderate" =>"3",
	"low" =>"1",
};

#F1DC97

sub printBadge {
	my ($value,$types) = @_;
	my $color = "#4CAF50";
	 #$color = "#D62C1A" if $value > $type[0] ;
	 $color = "#FF8800" if $value > $types->[0] ;
	 $color = "#FF0025" if $value > $types->[1] ;
	 return qq{<span class="badge badge-success badge-xs" style="border-color:$color;background-color:#FFFFFF;color:$color;font-size :$fts;">$value</span>} ;

}

sub printSimpleBadge {
	my ($value, $type) = @_;
	return $value if $print;
	my $color = "black";
	 return qq{<span class="badge badge-success badge-xs" style="font-weight:550;border-color:white;background-color:#EFF7EB;color:black;font-size :$fts;border-style: solid;border-width: 0px;font-family:verdana">$value</span>} if $type ==1;
	 return qq{<span class="badge badge-success badge-xs" style="font-weight:550;border-color:black;background-color:#EDF1F7;color:black;font-size :$fts;border-style: solid;border-width: 0px;font-family:verdana">$value</span>} if $type ==2;
	 return qq{<span class="badge badge-success badge-xs" style="font-weight:550;border-color:#2C3E50;background-color:#F7EEEE;color:black;font-size :$fts;border-style: solid;border-width: 0px;font-family:verdana">$value</span>} if $type ==3;
	 return qq{<span class="badge badge-success badge-xs" style="font-weight:550;border-color:#E8B5CE;background-color:#FFFFFF;color:$color;font-size :$fts;border-style: solid;border-width: 0px;font-family:verdana">$value</span>} ;
}

sub printInvBadge {
	my ($value,$types) = @_;
	
	my $color = "#4CAF50";
	 #$color = "#D62C1A" if $value > $type[0] ;
	 $color = "#FF8800" if $value < $types->[0] ;
	 $color = "#FF0025" if $value < $types->[1] ;
	 
	 $value ="-" unless defined $value;
	 $color = "#FF0025" if $value eq "-" ;
	
	 return qq{<span class="badge badge-success badge-xs" style="border-color:$color;background-color:#FFFFFF;color:$color;font-size :$fts;">$value</span>} ;
}

sub printInvButton {
	my ($value,$types,$text,$othercss) = @_;
	return $value if $print;
	my $btn_class  = qq{class= "btn btn-xs btn-primary " style="background-color: #D0D0D0;font-size: $fts;font-family:Verdana;color:black"};
	$btn_class = qq{class= "btn btn-xs btn-primary " style="background-color: #FF8800;font-size: $fts;font-family:Verdana;;color:black"} if $value <  $types->[0] ;
	$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: #e74c3c;font-size: $fts;font-family:Verdana;color:white"}  if $value < $types->[1] ;
	 $value ="-" unless defined $value;
	$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: #e74c3c;font-size: $fts;font-family:Verdana;"}  if $value eq "-" ;
	$text  = $value unless $text;
	return  qq{<button type="button" $btn_class $othercss>$text</button>};
}
sub printButton {
	my ($value,$types,$text,$othercss) = @_;
	return $text if $print;
	my $btn_class  = qq{class= "btn btn-sm btn-primary " style="background-color: #D0D0D0;font-size: $fts;font-family:Verdana;color:black;border-style: solid;border-width: 1px;border-color:white"};
	$btn_class = qq{class= "btn btn-sm btn-primary " style="background-color: #FA6800;font-size: $fts;font-family:Verdana;;color:white;border-style: solid;border-width: 1px;border-color:white"} if $value >=  $types->[0] ;
	$btn_class = qq{class= "btn btn-sm  btn-primary" style="background-color: #AA00FF;font-size: $fts;font-family:Verdana;color:white;border-style: solid;border-width: 1px;border-color:white"}  if $value >=  $types->[1] ;
	 $value ="-" unless defined $value;
	 
	$btn_class = qq{class= "btn btn-sm  btn-primary" style="background-color: red;font-size: $fts;font-family:Verdana;"}  if $value eq "-" ;
	$text  = $value unless $text;
	return  qq{<button type="button" $btn_class $othercss>$text</button>};
}

sub vcadd {
	my ($v,$hvariation) = @_;
	$hvariation->{value}->{cadd} = $v->cadd_score();
	$hvariation->{html}->{cadd} = printBadge($hvariation->{cadd},[20,30]);
}
sub vdbscsnv {
	my ($v,$hvariation) = @_;
	$hvariation->{value}->{dbscsnv} = "-";
	$hvariation->{value}->{dbscsnv_rf} = $v->dbscsnv_rf;
	$hvariation->{value}->{dbscsnv_ada} = $v->dbscsnv_ada;
	$hvariation->{html}->{dbscsnv} = printBadge($v->dbscsnv_rf,[0.6,0.9]).printBadge($v->dbscsnv_ada,[0.6,0.9])  if $v->dbscsnv_rf ne "-";
}
sub vrevel {
	my ($v,$hvariation) = @_;
	$hvariation->{value}->{revel} = "-";
	$hvariation->{value}->{revel} = $v->revel_score;
	$hvariation->{html}->{revel} = 	$hvariation->{revel} = printBadge($v->revel_score,[0.5,0.9]);
}

sub vrevel {
	my ($v,$hvariation) = @_;
	$hvariation->{value}->{revel} = "-";
	$hvariation->{value}->{revel} = $v->revel_score;
	$hvariation->{html}->{revel} = 	$hvariation->{revel} = printBadge($v->revel_score,[0.5,0.9]);
}

sub vcds {
		my ($v,$hvariation) = @_;
		$hvariation->{value}->{cds} = "";
		
}

sub construct_table_transcript {
	my ($v,$cgi,$header_transcripts,$level,$gene,$only_one) = @_;

	my  $transcripts = $v->getTranscripts();
	my $html= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"font-size: $fts;font-family:  Verdana;margin-bottom:0px"});
	$html.= $cgi->start_Tr();
	$html.=$cgi->th($header_transcripts);
	$html.= $cgi->end_Tr();
	my $hvariation;
	$hvariation->{cadd} = $v->cadd_score();
	$hvariation->{cadd} = printBadge($hvariation->{cadd},[20,30]);
	
		$hvariation->{dbscsnv} = "-";
		
		$hvariation->{dbscsnv} = printBadge($v->dbscsnv_rf,[0.6,0.9]).printBadge($v->dbscsnv_ada,[0.6,0.9])  if $v->dbscsnv_rf ne "-";
			$hvariation->{revel} = "-";
		
		$hvariation->{revel} = printBadge($v->revel_score,[0.5,0.9]);
		my @colors = ("#F9F6FF","#F9F9F9");
		my $nb =0;
		my $nb_skip = 0;
		my $rids =[];
	foreach my $tr1 (sort { ($himpact_sorted->{$v->effectImpact($b)} <=>  $himpact_sorted->{$v->effectImpact($a)}) or ($a->appris_level <=> $b->appris_level)} @$transcripts){
		next if $tr1->getGene->id ne $gene->id;
		$hvariation->{impact_score} = $himpact_sorted->{$v->effectImpact($tr1)};
		my $hide;
		
		$hide = "display:none;"  if ($hvariation->{impact_score}  <  $level and $nb >0) or ($only_one ==1 && $nb > 0)    ;
		$nb_skip ++ if $hide;
		$hvariation->{cds} = "";
		$hvariation->{prot} ="";
		$hvariation->{codons_AA} = "";
		$hvariation->{polyphen} ="-";
		$hvariation->{sift} ="-";
		$hvariation->{transcript} =$tr1->name;
		my $u = qq{https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=}.$tr1->name;
		$hvariation->{enst} =  qq{<a href="$u" target="_blank">}.$tr1->name."</a>"; 
		$hvariation->{nm} =$tr1->external_name;
		$hvariation->{nm} ="-" unless $hvariation->{nm};
		my $uc = qq{https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=}.$tr1->ccds_name;
		$hvariation->{ccds} = qq{<a href="$uc" target="_blank">}.$tr1->ccds_name."</a>";
	
		$hvariation->{impact_score} = $himpact_sorted->{$v->effectImpact($tr1)};
	
	# return qq{<span class="badge badge-success badge-xs" style="border-color:$color;background-color:#FFFFFF;color:$color;font-size :8px;">$value</span>} ;
		$hvariation->{ccds} ="-" unless $hvariation->{ccds};
		$hvariation->{appris} =$tr1->appris_type;
		if ($v->isCoding($tr1)){
			my $prot = $tr1->getProtein();
		
			$hvariation->{cds} = $v->getOrfPosition($prot);
			$hvariation->{prot} = $v->getProteinPosition($prot);
			$hvariation->{codons_AA} =   printBadge($v->protein_nomenclature($prot));#$v->getProteinAA($prot).$hvariation->{prot}.$v->changeAA($prot);
			$hvariation->{polyphen} = $v->polyphenScore($tr1->getProtein);
			$hvariation->{polyphen} = printBadge($hvariation->{polyphen},[0.446,0.908]);
			$hvariation->{sift} = $v->siftScore($tr1->getProtein);
			$hvariation->{sift} = printBadge($hvariation->{sift},[1,1]);
			$hvariation->{codons} =   printBadge($v->getCodons($tr1));
		}
		
		$hvariation->{exon} = $tr1->findExonNumber($v->start, $v->end);
		 if ($hvariation->{exon} == -1){
		 	$hvariation->{exon} = $tr1->findNearestExon($v->start, $v->end) ;
		 }
		 else {
		 	$hvariation->{exon} =  printBadge($hvariation->{exon});
		 }
		#$hvariation->{exon} = $tr1->findNearestExon($v->start, $v->end) if $hvariation->{exon} == -1;
		#$hvariation->{exon} = printBadge($hvariation->{exon});
		$hvariation->{nomenclature} =  printBadge($v->getNomenclature($tr1));
		my $value = $v->variationTypeInterface($tr1);
		my $x = $hvariation->{impact_score};
		$hvariation->{consequence} = printButton($x,[3,4],$value);
		#$hvariation->{consequence} = qq{<button type="button" $btn_class >$value-$x</button>};
		my $color = "black";
		$color = "#FF8800" if $hvariation->{impact_score} == 3 ;
		$color = "coral" if $hvariation->{impact_score} == 4 ; 
		my $c = $colors[$nb%2];
		$nb ++;
		#$hvariation->{consequence} =  qq{<span class="badge badge-success badge-xs" style="border-color:$color;background-color:#FFFFFF;color:$color;font-size :8px;">$value</span>} ;;
		my $rid = "row_".$v->id."_".$tr1->id;
		push(@$rids,$rid) if $hide;
		$html .=  $cgi->start_Tr({id=>$rid,style=>"border: 1px solid;background-color:$c ;".$hide});#{style=>"border: 1px solid black;color:lightgrey;$style"}).$cgi->td([$mother->name,qq{<i class="fa fa-female fa-2x" aria-hidden="true" style="color:ligthgrey"></i>},$type,$ps."%",""]);
		 foreach my $c (@$header_transcripts){
		 	$html.= $cgi->td($hvariation->{$c});
		 }

	}
	 if ($nb_skip ==0){
		$html.=$cgi->end_table();
		return $html
	 }
	my $js = encode_json $rids;
		my $za = $v->id."_tr".time;
	$html .=  $cgi->start_Tr({id=>$za});

	$html.= $cgi->td({style=>"border: 1px solid;background-color:#CECFCE;",colspan=>scalar(@$header_transcripts),onClick=>qq{showTranscripts($js,"$za");}},qq{<span class="glyphicon glyphicon-plus"></span> }."view $nb_skip Transcripts");
	$html.= $cgi->end_Tr();
	
	$html.=$cgi->end_table();
	return $html;
}

sub construct_table_transcript2 {
	my ($v,$cgi,$header_transcripts,$level,$gene) = @_;

	my  $transcripts = $v->getTranscripts();
	my $html= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"font-size: $fts;font-family:  Verdana;margin-bottom:0px"});
	$html.= $cgi->start_Tr();
	$html.=$cgi->th($header_transcripts);
	$html.= $cgi->end_Tr();
	my $hvariation = {};
	
	#$hvariation->{cadd} = $v->cadd_score();
	#$hvariation->{cadd} = printBadge($hvariation->{cadd},[20,30]);

	  vdbscsnv($v,$hvariation);
	 vrevel($v,$hvariation);
	
		my @colors = ("#F9F6FF","#F9F9F9");
		my $nb =0;
		my $nb_skip = 0;
		my $rids =[];
		my $all = {};
	foreach my $tr1 (sort { ($himpact_sorted->{$v->effectImpact($b)} <=>  $himpact_sorted->{$v->effectImpact($a)}) or ($a->appris_level <=> $b->appris_level)} @$transcripts) {
		next if $tr1->getGene->id ne $gene->id;
		my $value;
		my $html;
		
		$value->{impact_score}  =  $himpact_sorted->{$v->effectImpact($tr1)};
		$html->{impact_score}  =  $himpact_sorted->{$v->effectImpact($tr1)};
		$hvariation->{impact_score} = $himpact_sorted->{$v->effectImpact($tr1)};
		my $hide;
		
		#$hide = "display:none;"  if $hvariation->{impact_score}  <  $level and $nb >0  ;
		
		$nb_skip ++ if $hide;
		$hvariation->{cds} = "";
		$value->{cds} = "";
		$html->{cds} = "";
		
		$hvariation->{cadd} = $v->cadd_score();
		$value->{cadd} = $v->cadd_score();
		$html->{cadd} = printBadge($hvariation->{cadd},[20,30]);
		
		$hvariation->{dbscsnv} = "-";
		$value->{dbscsnv_rf} = $v->dbscsnv_rf;
		$value->{dbscsnv_ada} = $v->dbscsnv_ada;
		$hvariation->{dbscsnv} = printBadge($v->dbscsnv_rf,[0.6,0.9]).printBadge($v->dbscsnv_ada,[0.6,0.9])  if $v->dbscsnv_rf ne "-";
		$html->{dbscsnv} = printBadge($v->dbscsnv_rf,[0.6,0.9]).printBadge($v->dbscsnv_ada,[0.6,0.9])  if $v->dbscsnv_rf ne "-";
		
		$hvariation->{value}->{revel} = "-";
		$value->{revel} = $v->revel_score;
		$hvariation->{revel} = 	$hvariation->{revel} = printBadge($v->revel_score,[0.5,0.9]);
		$html->{revel} = 	$hvariation->{revel} = printBadge($v->revel_score,[0.5,0.9]);
		
		$hvariation->{prot} ="";
		$value->{prot} ="";
		$html->{prot} ="";
		
		$hvariation->{codons_AA} = "";
		$value->{codons_AA} = "";
		$html->{codons_AA} = "";
		
		$hvariation->{polyphen} ="-";
		$value->{polyphen} ="-";
		$html->{polyphen} ="-";
		
		$hvariation->{sift} ="-";
		$value->{sift} ="-";
		$html->{sift} ="-";
		
		$hvariation->{transcript} =$tr1->name;
		$value->{transcript} =$tr1->name;
		$html->{transcript} =$tr1->name;
		
		my $u = qq{https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=}.$tr1->name;
		
		
		
		$hvariation->{enst} =  qq{<a href="$u" target="_blank">}.$tr1->name."</a>"; 
		$value->{enst} =  $tr1->name; 
		$html->{enst} =  qq{<a href="$u" target="_blank">}.$tr1->name."</a>"; 
		
		$hvariation->{trid} =  $tr1->id; 
		$value->{trid} =  $tr1->id; 
		$html->{trid} =  $tr1->id; 
				
		$hvariation->{nm} =$tr1->external_name;
		$hvariation->{nm} ="-" unless $hvariation->{nm};
		$value->{nm} =$tr1->external_name;
		$html->{nm} =$tr1->external_name;
		
		 
		my $uc = qq{https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=}.$tr1->ccds_name;
		$hvariation->{ccds} = qq{<a href="$uc" target="_blank">}.$tr1->ccds_name."</a>";
		$value->{ccds} = $tr1->ccds_name;
		$html->{ccds} = qq{<a href="$uc" target="_blank">}.$tr1->ccds_name."</a>";
		$hvariation->{ccds} ="-" unless $hvariation->{ccds};
		$hvariation->{impact_score} = $himpact_sorted->{$v->effectImpact($tr1)};
	

	# return qq{<span class="badge badge-success badge-xs" style="border-color:$color;background-color:#FFFFFF;color:$color;font-size :8px;">$value</span>} ;
		
		$hvariation->{appris} =$tr1->appris_type;
		$value->{appris} =$tr1->appris_type;
		$html->{appris} =$tr1->appris_type;
		
		if ($v->isCoding($tr1)){
			my $prot = $tr1->getProtein();
		
			$hvariation->{cds} = $v->getOrfPosition($prot);
			$value->{cds} = $v->getOrfPosition($prot);
			$html->{cds} = $v->getOrfPosition($prot);
			
			$hvariation->{prot} = $v->getProteinPosition($prot);
			$value->{prot} = $v->getProteinPosition($prot);
			$html->{prot} = $v->getProteinPosition($prot);
			
			
			$hvariation->{codons_AA} =   printBadge($v->protein_nomenclature($prot));#$v->getProteinAA($prot).$hvariation->{prot}.$v->changeAA($prot);
			$value->{codons_AA} =   $v->protein_nomenclature($prot);#$v->getProteinAA($prot).$hvariation->{prot}.$v->changeAA($prot);
			$html->{codons_AA} =   printBadge($v->protein_nomenclature($prot));#$v->getProteinAA($prot).$hvariation->{prot}.$v->changeAA($prot);
			
			$hvariation->{polyphen} = $v->polyphenScore($tr1->getProtein);
			$value->{polyphen} = $v->polyphenScore($tr1->getProtein);
			$html->{polyphen} = printBadge($hvariation->{polyphen},[0.446,0.908]);
			$hvariation->{polyphen} = printBadge($hvariation->{polyphen},[0.446,0.908]);
			
			$hvariation->{sift} = $v->siftScore($tr1->getProtein);
			$value->{sift} = $v->siftScore($tr1->getProtein);
			$hvariation->{sift} = printBadge($hvariation->{sift},[1,1]);
			$html->{sift} = printBadge($hvariation->{sift},[1,1]);
			
			$hvariation->{codons} =   printBadge($v->getCodons($tr1));
			$value->{codons} = $v->getCodons($tr1);
			$html->{codons} =   printBadge($v->getCodons($tr1));
		}
		
		my $te = $tr1->findExonNumber($v->start, $v->end);
		 if ($te == -1){
		 	my $tc = $tr1->findNearestExon($v->start, $v->end);
		 	$value->{exon} = $tc;
		 	$html->{exon} = printBadge($tc);
		 	$hvariation->{exon} = $tc ;
		 }
		 else {
		 	$value->{exon} = $te;
		 	$hvariation->{exon} =  printBadge($te);
		 	 $html->{exon} =  printBadge($te);
		 	 	
		 	
		 }
		#$hvariation->{exon} = $tr1->findNearestExon($v->start, $v->end) if $hvariation->{exon} == -1;
		#$hvariation->{exon} = printBadge($hvariation->{exon});
		$hvariation->{nomenclature} =  printBadge($v->getNomenclature($tr1));
		$html->{nomenclature} =  printBadge($v->getNomenclature($tr1));
		$value ->{nomenclature} =  $v->getNomenclature($tr1);
		
		my $vvalue = $v->variationTypeInterface($tr1);
		my $x = $hvariation->{impact_score};
		$hvariation->{consequence} = printButton($x,[3,4],$vvalue);
		$html->{consequence} = printButton($x,[3,4],$vvalue);
		$value->{consequence} = $x;
		
		my $color = "black";
		$color = "#FF8800" if $hvariation->{impact_score} == 3 ;
		$color = "coral" if $hvariation->{impact_score} == 4 ; 
		my $c = $colors[$nb%2];
		$nb ++;
		my $rid = "row_".$v->id."_".$tr1->id;
		push(@{$all->{value}},$value);
		push(@{$all->{html}},$html);
		push(@$rids,$rid) if $hide;
		
		$html .=  $cgi->start_Tr({id=>$rid,style=>"border: 1px solid;background-color:$c ;".$hide});#{style=>"border: 1px solid black;color:lightgrey;$style"}).$cgi->td([$mother->name,qq{<i class="fa fa-female fa-2x" aria-hidden="true" style="color:ligthgrey"></i>},$type,$ps."%",""]);
		 foreach my $c (@$header_transcripts){
		 	$html.= $cgi->td($hvariation->{$c});
		 }

	}
	 if ($nb_skip ==0){
		$html.=$cgi->end_table();
		return $html
	 }
	my $js = encode_json $rids;
		my $za = $v->id."_tr".time;
	$html .=  $cgi->start_Tr({id=>$za});

	$html.= $cgi->td({style=>"border: 1px solid;background-color:#CECFCE;",colspan=>scalar(@$header_transcripts),onClick=>qq{showTranscripts($js,"$za");}},qq{<span class="glyphicon glyphicon-plus"></span> }."view $nb_skip Transcripts");
	$html.= $cgi->end_Tr();
	
	$html.=$cgi->end_table();
	return $html;
}

sub print_table_transcripts {
	 my ($v,$cgi,$header_transcripts,$level,$gene) = @_;
}

sub table_value_html_badge {
	my ($hv,$type,$value1,$value2) = @_;
	
	
	
}

sub vcosmic {
	my ($v,$hvariation) = @_;
	 	$hvariation->{post_bug} = 1; 
	 	$hvariation->{cosmic} = $v->cosmic();
		$hvariation->{value}->{cosmic} = $v->cosmic();
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
			$hvariation->{html}->{cosmic} = qq{<a href=\"$url$id1\" target="_blank">$text</a>};
		}
	
}

sub update_cosmic {
	my ($hvariation,$project) = @_;
	return if 	exists $hvariation->{post_bug};
	my $version = $project->buffer->get_version_database("cosmic");
	return if ($version != 98);
	
	 my $v  = $project->_newVariant($hvariation->{id});
	 	 
	 	$hvariation->{cosmic} = $v->cosmic();
		$hvariation->{value}->{cosmic} = $v->cosmic();
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
			$hvariation->{html}->{cosmic} = qq{<a href=\"$url$id1\" target="_blank">$text</a>};
		}
}

sub vname {
	my ($v,$hvariation) = @_;
	
	my $vn=$v->vcf_id;
	$vn =~ s/_/-/g;
	$vn=~ s/chr//;
	my $pp = $v->getChromosome->name."-".$v->start;
	$hvariation->{var_name} = printSimpleBadge(qq{<a href='https://gnomad.broadinstitute.org/region/$pp' target = '_blank' style="color:black">$vn</a> });;
	$hvariation->{value}->{var_name} = $vn;
	$hvariation->{html}->{var_name} = printSimpleBadge(qq{<a href='https://gnomad.broadinstitute.org/region/$pp' target = '_blank' style="color:black">$vn</a> });
	if ($v->name() =~ /rs/){
			my $vname = $v->name();
			
			$hvariation->{var_name} = printSimpleBadge(qq{<a href='https://gnomad.broadinstitute.org/variant/$vn' target = '_blank' style="color:black"><i class="fa fa-users fa-2x" style="color:coral"></i>&nbsp$vn</a> });	;#qq{<a href='http://www.ncbi.nlm.nih.gov/snp/?term=$vn' target = '_blank'>$vn</a> };	
			$hvariation->{html}->{var_name}  = printSimpleBadge(qq{<a href='https://gnomad.broadinstitute.org/variant/$vn' target = '_blank' style="color:black"><i class="fa fa-users fa-2x" style="color:coral"></i>&nbsp$vn</a> });	;
	}
}

sub vsequencing  {
	my ($v,$hvariation,$patient) = @_;
	my @asequence_info;
	my @apc;
	my @methods;
	my $nb_methods;
	my $max_pc =-1;
		my $max_dp =  -1;
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
		 $pc = int ($nb_alt *10000/($sum))/100;
		 $pc = ceil($pc) if $pc >1;
		}
		$max_pc = $pc if $pc > $max_pc;
		$max_dp = $sum if $sum > $max_dp;
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
		
		$hvariation->{max_dp} = $max_dp;
		$hvariation->{value}->{max_dp} = $max_dp;
		$hvariation->{html}->{max_dp} = $max_dp;
		
		$hvariation->{max_pc} = $max_pc;
		$hvariation->{value}->{max_pc} = $max_pc;
		$hvariation->{html}->{max_pc} = $max_pc;
		
		$hvariation->{value}->{ngs} = \@asequence_info;
		
		$hvariation->{ngs} = printSimpleBadge(join("<br>",@asequence_info));
		$hvariation->{value}->{ngs} = \@asequence_info;
		$hvariation->{html}->{ngs} = printSimpleBadge(join("<br>",@asequence_info));
		
		
		$hvariation->{ratio} =  printSimpleBadge(join("<br>",@apc));
		$hvariation->{value}->{ratio} =  \@apc;
		$hvariation->{html}->{ratio} =  printSimpleBadge(join("<br>",@apc));
		$hvariation->{caller} =  printSimpleBadge(join("<br>",@methods));
		$hvariation->{value}->{caller} =  \@methods;
		$hvariation->{html}->{caller} =  printSimpleBadge(join("<br>",@methods));
	
}

sub value_html {
	my ($hvariation,$type,$v1,$v2) = @_;
	$v2 = $v1 unless defined $v2;
	$hvariation->{value}->{$type} = $v1;
	$hvariation->{html}->{$type} = $v2 ; 
}

sub vdivers {
	my ($v,$hvariation) = @_;
		value_html($hvariation,"allele",$v->getSequence());
		value_html($hvariation,"ref_allele", $v->ref_allele());
		value_html($hvariation,"genomique", $v->getChromosome()->name.":".$v->start,printSimpleBadge($v->getChromosome()->name.":".$v->start));
		value_html($hvariation,"genomique_value", $v->getChromosome()->name.":".$v->start);
		value_html($hvariation,"start",$v->start);
		value_html($hvariation,"end",$v->end);
		value_html($hvariation,"chromosome",$v->getChromosome()->name);
	
	
}
sub valamut_igv {
	my ($v,$hvariation,$patient) = @_;
		my $start = $v->start();
		my $chr = $v->getChromosome();
		
		my $qq4 = qq{	<div  data-dojo-type="dijit/Toolbar"><button dojoType="dijit.form.Button"   iconClass="igvIcon" onClick ="displayInIGV('$chr',$start,$start);"></button></div>};
		my $qq4 = qq{	<button  class="igvIcon2" onClick ="displayInIGV('$chr',$start,$start);">&nbsp;&nbsp;&nbsp;&nbsp;</button>};
		value_html($hvariation,"igv",$chr.":".$start,$qq4);
		
		$hvariation->{igv} = $qq4; 		
				my @bams;
					my @names;
					foreach my $p (@{$patient->getFamily->getPatients()}){
						push(@bams,$p->bamUrl);
						push(@names,$p->name());
					}
					
				my $f =  join(";",@bams);#$patient->{obj}->bamUrl;;
					my $l = $v->getChromosome()->name.":".$v->start;
					my $v1 = $hvariation->{ref_allele}."/".$hvariation->{allele};	
					my $gn = $patient->project->getVersion();
					my $project_name = $patient->project->name;
					my $pnames = join(";",@names);
					
					#$text =qq{<button dojoType="dijit.form.Button"   iconClass='igvIcon' onclick='launch_web_igv("$project_name","$pnames","$f","$l","$v","$gn")' style="color:black">toto</button>};
					#$text =qq{<button onclick='alert("coucou");' style="color:black">toto</button>};
					my $text =qq{<button dojoType="dijit.form.Button"   iconClass='igvIcon' onclick='launch_web_igv("$project_name","$pnames","$f","$l","$v1","$gn")' style="color:black"></button>};
					value_html($hvariation,"igv_web",$chr->name.":".$start,$text);
						
			my $a0 = $v->ref_allele;
		my $a1 = $v->var_allele;			
			my $qq5 = qq{	<button    class="alamutView3" onClick ="displayInAlamut('$chr',$start,['$a0','$a1']);"></button>};
			value_html($hvariation,"alamut",$chr->name.":".$start.$a0."_".$a1,$qq5);
			
		
		
			
} 

sub vgnomad {
	my ($v,$hvariation) = @_;
	 
		die();
	 	my $max  ="-";
	 	$max = $v->max_pop_name.":".sprintf("%.4f", $v->max_pop_freq ) if $v->max_pop_name;
	 	value_html($hvariation,"max_pop",$max, printSimpleBadge($max));
	 	
	 	my $min = "-";
		$min = $v->min_pop_name.":".sprintf("%.4f", $v->min_pop_freq ) if $v->min_pop_name;
		value_html($hvariation,"min_pop",$min, printSimpleBadge($min));
		
		my $freq_ho = "-"; 
		$freq_ho = sprintf("%.4f", $v->frequency_homozygote ) if $v->frequency_homozygote ;
		value_html($hvariation,"freq_ho",$freq_ho, printSimpleBadge($freq_ho));
		
		value_html($hvariation,"ac",$v->getGnomadAC, printInvButton($v->getGnomadAC,[200,10]));
		value_html($hvariation,"an",$v->getGnomadAN,  printInvBadge($v->getGnomadAN,[0,0]));
		
	
		 
		value_html($hvariation,"ac_ho",$v->getGnomadHO, printInvButton($v->getGnomadHO,[50,5]));
		
		
		if  ($v->getChromosome->name eq "X" or $v->getChromosome->name eq "Y") {
			$hvariation->{ac_ho} = printInvButton($v->getGnomadHO,[50,5], $v->getGnomadHO." -  ".qq{&nbsp<i class="fa fa-mars" > </i> &nbsp;} .$v->getGnomadAC_Male);#printInvButton($v->getGnomadAC_Male,[50,5]) unless ($v->is_in_pseudoautosomal_region );  
			value_html($hvariation,"ac_ho",$v->getGnomadHO.":".$v->getGnomadAC_Male,printInvButton($v->getGnomadHO,[50,5], $v->getGnomadHO." -  ".qq{&nbsp<i class="fa fa-mars" > </i> &nbsp;} .$v->getGnomadAC_Male));
		}
}

sub vclinical_local {
		my ($v,$hvariation,$patient) = @_;
	 my $v1 = $v->score_clinical_local();
	 $hvariation->{clinical_local}  = "" ;
	 if ($v1){
	 	 	$hvariation->{clinical_local}  ++ ;
	 	 #	$hvariation->{clinvar_alert}  ++ ;
	 	 	my $cm = $v->comment_clinical_local();
	 	 	my $txt =  qq{<span class="badge badge-warning" style="font-size:8px">}.$hvariation->{clinvar}."</span>".qq{<span class="badge badge-warning"  onClick='alert ($cm)' style="font-size:8px">Local Clinical</span>};
	 	 	#$hvariation->{clinvar}  = " $txt ";
	 		$hvariation->{scaled_score} = 4 if $hvariation->{freq_level} <= 2 ;	
	 		
	 }

		
}


sub vclinvar {
		my ($v,$hvariation) = @_;
		 my $cl  = "" ;
		 my $alert = 0;
	 my $v1 = $v->score_clinvar();();
	   #  $hvariation->{clinvar_alert}  = 0 ;
	 if ($v1){
	 	
	 	my $uc = qq{https://www.ncbi.nlm.nih.gov/clinvar/?term=}.$v->clinvar->{id}."[alleleid]";
	 	my $a = qq{<a href="$uc" target="_blank" style="color:white">}.$v->text_clinvar()."</a>"; 
	 	my $oc = qq{onClick='window.open("$uc")'};
	 	value_html($hvariation,"clinvar",$v->text_clinvar(),  printButton($v,[3,4],$v->text_clinvar(),$oc));
	 	if (($v == 4 || $v==5)    ){
	 		$alert  = 4  if $hvariation->{value}->{freq_level} <= 2 ;
	 		$alert ++;
	 	}
	 }
	 	 value_html($hvariation,"clinvar_alert",$alert,$alert);
}

sub vhgmd {
		my ($v,$hvariation) = @_;
		 
		 
	 if ($v->hgmd_id()){
	 	value_html($hvariation,"hgmd",1,1);
	 	if  ($v->isDM){
	 		my $dm =1;
	 		value_html($hvariation,"dm",1,1);
	 		value_html($hvariation,"scaled_score",4,4) if ($v->getGnomadAC < 20);
	 		value_html($hvariation,"scaled_score",2,2) if ($v->getGnomadAC > 1000);
	 		$hvariation->{clinvar_alert}++;
			 		
	 	}
	 	my $txt = $v->hgmd->{phen}." - ".$v->hgmd->{class}." - ".$v->hgmd_id;
	 	$txt =~s/\"//g;
	 	my $n1 = $v->project->name;
	 	my $n2 = $v->hgmd_id;
	 	my $n3 = $v->id;
	 	my $cmd = qq{zoomHgmd(\'$n1\',\'$n2\',\'$n3\')};# "zoomHgmd(\'.$project->name."\',\'"." ','".$hvariation->{obj}->id."')";
		$hvariation->{hgmd}    = printButton(4,[3,4],$v->hgmd->{class},qq{onClick="$cmd"}); 
		value_html($hvariation,"hgmd",$v->hgmd->{class}.":".$v->hgmd->{id},printButton(4,[3,4],$v->hgmd->{class},qq{onClick="$cmd"}));
	 }
}


sub vfreq_level {
	my ($v,$hvariation) = @_;
		my $scaled =  $v->scaled_score_frequence(); 
		 value_html($hvariation,"freq_level",$scaled->{freq_level},$scaled->{freq_level});
		  value_html($hvariation,"freq_score",$scaled->{freq_score},$scaled->{freq_score});
}

my @header_transcripts = ("consequence","enst","nm","ccds","appris","exon","nomenclature","codons","codons_AA", "polyphen","sift","cadd","revel","dbscsnv");




sub construct_hash_variant {
	my ($project,$v,$vquery,$patient) = @_;
	my $hvariation;
	$hvariation->{value}->{id} =  $v->id;
	$hvariation->{html}->{id} =  $v->id;
	
	
	$hvariation->{value}->{type} = $v->type;
	$hvariation->{html}->{type} = $v->type;
	foreach my $p (@{$v->getPatients}){
		$hvariation->{html}->{patients} .= $p->name." ";
		
		
	}
	$hvariation->{html}->{infos} .= $v->{infos};
	vcosmic($v,$hvariation)	if ($project->isSomatic);
	vgnomad($v,$hvariation);
	vname($v,$hvariation);
	vsequencing($v,$hvariation,$patient);
	vdivers($v,$hvariation);
	valamut_igv($v,$hvariation,$patient);
	vclinvar($v,$hvariation);
	vhgmd($v,$hvariation);
	delete $hvariation->{obj};
	my $gs;
	my $cgi          = new CGI();
	foreach my $g (@{$v->getGenes}){
		push(@$gs,$g->id);
		$hvariation->{value}->{genes}->{$g->id};
		$hvariation->{html}->{genes}->{$g->id} = construct_table_transcript($v,$cgi,\@header_transcripts,2,$g);
	}
	return $hvariation;
	
	
	
	
}


sub vgenes {
	my ($ogene,$hgene) =@_;
		my $score = 0; 
		#$bilan->{$gene}->{max_score} +=0.5  if $ogene->pLI > 0.9; 
		$score += 0.2 if $ogene->pLI > 0.95;
		my $pheno = $ogene->phenotypes();
		$score += 1 if $pheno =~/intellectual/ or $pheno =~/mental/ or $pheno =~/retar/;
		my $hpanels = $ogene->buffer->queryPanel()->getPanelsForGeneName($ogene->external_name);
		$score += 0.4 * scalar(keys %{$hpanels});# if keys %{$bilan->{$gene}->{panels}};
#		my @t = grep{$_>5} @{$bilan->{$gene}->{scores}};
#		push(@t, grep {$_>2.5} @{$bilan->{$gene}->{score_m}});
#		push(@t, grep {$_>2.5} @{$bilan->{$gene}->{score_p}});
#		$bilan->{$gene}->{max_score} -= 0.5 if scalar(@t) >= 5 ;
#		$bilan->{$gene}->{max_score} -=  1 if scalar(@t) >= 10 ;
		value_html($hgene,"max_score",$score,$score);
		return 
}









sub construct_variant {
	my ($project,$v,$tr1,$patient,$vquery) = @_;
	my $hvariation;
	$hvariation->{id} = $v->id;
	$hvariation->{value}->{id} =  $v->id;
	$hvariation->{html}->{id} =  $v->id;
	$hvariation->{type} = $v->type;
	my $debug ;
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
		$hvariation->{max_impact_text} = $v->effectImpact($tr1->getGene);
		$hvariation->{impact_score} = $himpact_sorted->{$v->effectImpact($tr1)};
		$hvariation->{dbscsnv} = "-";
		
		$hvariation->{dbscsnv} = printBadge($v->dbscsnv_rf,[0.6,0.9]).printBadge($v->dbscsnv_ada,[0.6,0.9]);
			$hvariation->{revel} = "-";
		
		$hvariation->{revel} = printBadge($v->revel_score,[0.5,0.9]);
		#$v->dbscsnv_ada.":".$v->dbscsnv_rf if $v->dbscsnv_rf;
		$hvariation->{gene} = $tr1->getGene->external_name();
		$hvariation->{gene_id} = $tr1->getGene->id();
		#$hvariation->{var_name} = $v->name()." ".$v->{scaled_score};
			my $vn=$v->vcf_id;
			$vn =~ s/_/-/g;
			$vn=~ s/chr//;
		
		my $pp = $v->getChromosome->name."-".$v->start;
		$hvariation->{var_name} = printSimpleBadge(qq{<a href='https://gnomad.broadinstitute.org/region/$pp' target = '_blank' style="color:black">$vn</a> });;
		
		if ($v->name() =~ /rs/){
			my $vname = $v->name();
			
			$hvariation->{var_name} = printSimpleBadge(qq{<a href='https://gnomad.broadinstitute.org/variant/$vn' target = '_blank' style="color:black"><i class="fa fa-users fa-2x" style="color:coral"></i>&nbsp$vn</a> });	;#qq{<a href='http://www.ncbi.nlm.nih.gov/snp/?term=$vn' target = '_blank'>$vn</a> };	
		}
		my @asequence_info;
		my @apc;
		my @methods;
		my $nb_methods;
		my $max_pc =-1;
		my $max_dp =  -1;
		foreach my $method_name (@{$v->getMethods($patient)}) {
			push(@methods,$method_name);
			my $sequence_info = $v->getSequencingGenotype($patient,$method_name)."(".$v->getNbAlleleRef($patient,$method_name).":".$v->getNbAlleleAlt($patient,$method_name).")";
			push(@asequence_info,$sequence_info);
			push(@apc,$method_name.":".$v->getRatio($patient,$method_name)."%");
			
		}
		
		 if ($v->validation_method eq "sanger" ) {
		 	#$sequence_info = "-";
		 	push(@asequence_info,"-");
		 }
		
		$hvariation->{max_dp} = $v->getRatio($patient);
		$hvariation->{max_pc} = $v->getDP($patient);;
		$hvariation->{ngs} = printSimpleBadge(join("<br>",@asequence_info));
	
		$hvariation->{ratio} =  printSimpleBadge(join("<br>",@apc));
		$hvariation->{caller} =  printSimpleBadge(join("<br>",@methods));
		$hvariation->{allele} = $v->getSequence();
		$hvariation->{ref_allele} = $v->ref_allele();
		$hvariation->{genomique} = printSimpleBadge($v->getChromosome()->name.":".$v->start);
		$hvariation->{genomique_value} = $v->getChromosome()->name.":".$v->start;
		$hvariation->{start} = $v->start;
		$hvariation->{transcript} = printSimpleBadge($tr1->name);
		$hvariation->{enst} = $tr1->name;
		$hvariation->{transcript_external_name} = $tr1->external_name;
		$hvariation->{chromosome} = $v->getChromosome()->name;
		$hvariation->{trans} = $v->start * $tr1->strand; 
		$hvariation->{cds} = "";
		$hvariation->{prot} ="";
		$hvariation->{codons_AA} = "";
		$hvariation->{polyphen} ="-";
		$hvariation->{sift} ="-";

		my $debug;
		$debug = 1   if $v->name eq "rs151344528";
	
		if ($v->isDeletion){
	
			#$hvariation->{codons}  =  $v->delete_sequence."/".$v->sequence();
		}
		else {
			#$hvariation->{codons}  =  $v->getChromosome()->sequence($v->start,$v->end)."/".$v->sequence();
		}
		if ($tr1->strand() == -1 ){
			#$hvariation->{codons}  =  BioTools::complement_sequence($v->getChromosome()->sequence($v->start,$v->end))."/".BioTools::complement_sequence($v->sequence());
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
				my @bams;
					my @names;
					foreach my $p (@{$patient->getFamily->getPatients()}){
						push(@bams,$p->bamUrl);
						push(@names,$p->name());
					}
					
				my $f =  join(";",@bams);#$patient->{obj}->bamUrl;;
				my $l = $v->getChromosome()->name.":".$v->start;
				my $v1 = $hvariation->{ref_allele}."/".$hvariation->{allele};	
				my $gn = $patient->project->getVersion();
				my $project_name = $patient->project->name;
				my $pnames = join(";",@names);
					
					#$text =qq{<button dojoType="dijit.form.Button"   iconClass='igvIcon' onclick='launch_web_igv("$project_name","$pnames","$f","$l","$v","$gn")' style="color:black">toto</button>};
					#$text =qq{<button onclick='alert("coucou");' style="color:black">toto</button>};
	#		my $text =qq{<button dojoType="dijit.form.Button"   iconClass='igvIcon' onclick='launch_web_igv("$project_name","$pnames","$f","$l","$v1","$gn")' style="color:black"></button>};
		#	$hvariation->{igv_web} = $text; 		
					#my $text =qq{<button dojoType="dijit.form.Button"   iconClass='igvIcon' onclick='launch_web_igv("$project_name","$pnames","$f","$l","$v1","$gn")' style="color:black"></button>};
					my $text =qq{<button dojoType="dijit.form.Button"   iconClass='igvIcon' onclick='view_web_igv_bam("dialog_igv", "div_igv", "$l", "$f", "$pnames")' style="color:black"></button>};

					$hvariation->{igv_web} = $text; 		
			
			
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
			$hvariation->{codons_AA} =   printBadge($v->protein_nomenclature($prot));#$v->getProteinAA($prot).$hvariation->{prot}.$v->changeAA($prot);
			$hvariation->{polyphen} = $v->polyphenScore($tr1->getProtein);
			$hvariation->{polyphen} = printBadge($hvariation->{polyphen},[0.446,0.908]);
			$hvariation->{sift} = $v->siftScore($tr1->getProtein);
			$hvariation->{sift} = printBadge($hvariation->{sift},[1,1]);
			$hvariation->{codons} =   printBadge($v->getCodons($tr1));
		}
		$hvariation->{exon} = printBadge($tr1->findExonNumber($v->start, $v->end));
		$hvariation->{exon} = printBadge($tr1->findNearestExon($v->start, $v->end)) if $hvariation->{exon} == -1;
		$hvariation->{nomenclature} =  printBadge($v->getNomenclature($tr1));
		$hvariation->{consequence} =  $v->variationTypeInterface($tr1);
		$hvariation->{consequence_gene} =  $v->variationTypeInterface($tr1->getGene);
		
		#warn $v->variationType($tr1);
		

#		if ($v->isUpstream($tr1)){
#			#$hvariation->{consequence} = "upstream";
#		}
#		if ($v->isDownstream($tr1)){
#			$hvariation->{consequence} = "downstream";
#		}

		$hvariation->{consequence} = printBadge($hvariation->{consequence});
		$hvariation->{freq}  =  $v->frequency;
		 delete $v->{frequency};
	
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
		#p $scaled;
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
	 		$hvariation->{scaled_score} = 4 if $hvariation->{freq_level} <= 2 ;	
	 		
	 }
	 	 if ($v->hgmd_id()) {
	 	if  ($v->isDM){
	 		$hvariation->{scaled_score} = 6 ;
	 		
	 		
	 		#$hvariation->{scaled_score} = 4 if ($hvariation->{ac} < 20);
	 		$hvariation->{scaled_score} = 3 if ($v->getGnomadAC > 50);
	 		$hvariation->{scaled_score} = 2 if ($v->getGnomadAC > 300);
	 		$hvariation->{dm} = 1 ;
	 		$hvariation->{clinvar_alert} ++;
	 	}
	 	my $txt = $v->hgmd->{phen}." - ".$v->hgmd->{class}." - ".$v->hgmd_id;
	 	$txt =~s/\"//g;
	 	my $n1 = $project->name;
	 	my $n2 = $v->hgmd_id;
	 	my $n3 = $v->id;
	 	my $cmd = qq{zoomHgmd(\'$n1\',\'$n2\',\'$n3\')};# "zoomHgmd(\'.$project->name."\',\'"." ','".$hvariation->{obj}->id."')";
	 	$hvariation->{hgmd}  =qq{<span class="badge badge-warning"  onClick="$cmd" style="font-size:10px">}.$v->hgmd->{class}."</span>";#$hvariation->{obj}->hgmd()->{hgmd_id}."<br>".$hvariation->{obj}->hgmd()->{class}."<br>".$hvariation->{obj}->hgmd()->{phen};
	 	
	 	$hvariation->{hgmd}  =qq{<button type="button" class="btn btn-danger btn-xs" onClick="$cmd" style="font-size:10px;border-color:black">}.$v->hgmd->{class}."</button>";#$hvariation->{obj}->hgmd()->{hgmd_id}."<br>".$hvariation->{obj}->hgmd()->{class}."<br>".$hvariation->{obj}->hgmd()->{phen};
	 	$hvariation->{hgmd}    = printButton(4,[3,4],$v->hgmd->{class},qq{onClick="$cmd"}); 
	
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
	elsif ($hvariation->{allele} eq "ALU"){
		#deletion 
		$hvariation->{db_type} = "insertions";
	}
		elsif ($hvariation->{allele} eq "LINE1"){
		#deletion 
		$hvariation->{db_type} = "insertions";
	}
		elsif ($hvariation->{allele} eq "SVA"){
		#deletion 
		$hvariation->{db_type} = "insertions";
	}
	elsif ($lalt eq $lref && $lref eq 1)
	{
		$hvariation->{db_type} = "snps";
	}
	elsif (length($ref) < length($alt)) {
		$hvariation->{db_type} = "insertions";
	}
	elsif (length($ref) > length($alt)) {
		$hvariation->{db_type} = "deletions";
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
	 $hvariation->{scaled_score} =4 if $sc2 >= 7 && $hvariation->{obj}->cadd_score > 30;
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
	
	 #return if exists $hvariation->{clinvar};
	 unless ($hvariation->{obj}){
	 
	 	 $hvariation->{obj} = $project->_newVariant($hvariation->{id});
	 #	confess();
	 }
	 
	 
	 my $debug;
	 
	 my $v = $hvariation->{obj}->score_clinvar();
	   $hvariation->{clinvar}  = "" ;
	     $hvariation->{clinvar_alert}  = 0 ;
	 if ($v){
	 	my $uc = qq{https://www.ncbi.nlm.nih.gov/clinvar/?term=}.$hvariation->{obj}->clinvar->{id}."[alleleid]";
	 	my $a = qq{<a href="$uc" target="_blank" style="color:white">}.$hvariation->{obj}->text_clinvar()."</a>"; 
	 	my $oc = qq{onClick='window.open("$uc")'};
	 	
	 	$hvariation->{clinvar}  = printButton($v,[3,4],$hvariation->{obj}->text_clinvar(),$oc); $hvariation->{obj}->text_clinvar();;
	 
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
	 	$hvariation->{hgmd} =1;
	 	if  ($hvariation->{obj}->isDM){
	 		$hvariation->{dm} =1;
	 		
	 		$hvariation->{scaled_score} = 4 if ($hvariation->{ac} < 20);
	 		$hvariation->{scaled_score} = 2 if ($hvariation->{ac} > 1000);
	 		$hvariation->{clinvar_alert}++;
			 		
	 	}
	 	my $txt = $hvariation->{obj}->hgmd->{phen}." - ".$hvariation->{obj}->hgmd->{class}." - ".$hvariation->{obj}->hgmd_id;
	 	$txt =~s/\"//g;
	 	my $n1 = $project->name;
	 	my $n2 = $hvariation->{obj}->hgmd_id;
	 	my $n3 = $hvariation->{obj}->id;
	 	my $cmd = qq{zoomHgmd(\'$n1\',\'$n2\',\'$n3\')};# "zoomHgmd(\'.$project->name."\',\'"." ','".$hvariation->{obj}->id."')";
	 	return ($hvariation->{obj}->hgmd->{class}, $cmd) if ($return_cmd);
 		$hvariation->{hgmd}  =qq{<span class="badge badge-danger"  onClick="$cmd" style="font-size:10px">}.$hvariation->{obj}->hgmd->{class}."</span>";#$hvariation->{obj}->hgmd()->{hgmd_id}."<br>".$hvariation->{obj}->hgmd()->{class}."<br>".$hvariation->{obj}->hgmd()->{phen};
	    
	    $hvariation->{hgmd}  =qq{<button type="button" class="btn btn-primary btn-xs"  onClick="$cmd" style="font-size:10px"><bold>--}.$hvariation->{obj}->hgmd->{class}."</bold></button>";#$hvariation->{obj}->hgmd()->{hgmd_id}."<br>".$hvariation->{obj}->hgmd()->{class}."<br>".$hvariation->{obj}->hgmd()->{phen};
		$hvariation->{hgmd}    = printButton(4,[3,4],$hvariation->{obj}->hgmd->{class},qq{onClick="$cmd"}); 
	 }
}
sub clinical_local {
	 my ($project,$hvariation,$url) = @_;
	 confess();
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
	my $level = "success";
	$level = "warning" if $hvariation->{cadd} >20;
	$level = "alert" if $hvariation->{cadd} >=30;
	$hvariation->{cadd} = printBadge($hvariation->{cadd},[20,30]);
}

sub trio {
	 my ($project,$tr,$hvariation,$patient,$cgi,$print) = @_;
	 return $hvariation->{trio} if exists $hvariation->{trio};
	 my $fam = $patient->getFamily();
	
	 unless ($fam->isTrio){
	 	 $hvariation->{trio} = "-";
	 	 return;
	 }
	# ENST00000301067 
	  return if exists $hvariation->{trio};
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
	#  if ($patient->isChild){
	  	$children = $fam->getChildren();
	  	$mother = $fam->getMother;
	  	$father = $fam->getFather;
#	  	unless ($father ){
#	  		$hvariation->{trio} ="-";
#	  		return ;
#	  	7
#	  	unless ($mother){
#	  		$hvariation->{trio} ="-";
#	  		return ;
#	  	}
	  	my $html =$cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"font-size: $fts;font-family:  Verdana;margin-bottom:0px"});
		my $html_line;
		my $type = "-";
	
		if ($mother){
			$type = "ho" if $var->isHomozygote($mother);
			$type = "he" if $var->isHeterozygote($mother);
			my $style = "";
		
				if ($var->getPourcentAllele($mother) eq "-"){
					my $ps = "-";
					my $depth =  return_coverage($mother,$var->getChromosome->ucsc_name,$var->start);#$mother->depth($var->getChromosome->name,$var->start,$var->end+1);
					$ps = "$depth reads";
					
					#$ps ="0" if return_coverage($mother,$var->getChromosome->ucsc_name,$var->start)>5;
					
					my $tab = [$mother->name,qq{<i class="fa fa-female fa-2x" aria-hidden="true" style="color:ligthgrey"></i>},$type,$ps,""];
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
			#$ps ="0" if return_coverage($father,$var->getChromosome->ucsc_name,$var->start)>5;
			my $depth = "-";
			 $depth =  return_coverage($father,$var->getChromosome->ucsc_name,$var->start) if $father;#$mother->depth($var->getChromosome->name,$var->start,$var->end+1);
			$ps = "$depth reads";
			my $tab = [$father->name,qq{<i class="fa fa-male fa-2x" aria-hidden="true" style="color:lightgrey"></i>},$type,$ps,""];
			$html .=  $cgi->start_Tr({style=>"border: 1px solid black;color:lightgrey;$style"}).$cgi->td($tab).$cgi->end_Tr();
			push(@$all_tab,@$tab);
			#splice @$tab, 1, 1;
			$html_line .= join(":",@$tab);
			$html_line .="- &nbsp;";
	
		}
		else {
			my $depth = $father->depth($var->getChromosome->name,$var->start,$var->end+1);
			my $tab = [$father->name,qq{<i class="fa fa-male fa-2x" aria-hidden="true" style="color:blue"></i>},$type,$var->getPourcentAllele($father)."%",""];
			$html .=  $cgi->start_Tr({style=>"border: 1px solid black;color:black;$style"}).$cgi->td($tab).$cgi->end_Tr();
			push(@$all_tab,@$tab);
			#splice @$tab, 1, 1;
			$html_line .= join(":",@$tab);
			$html_line .=" - &nbsp;";
		}
		}
		foreach my $child (@$children) {
		 $type = "-";
		$type = "ho" if $var->isHomozygote($child);
			$type = "he" if $var->isHeterozygote($child);
			if ($type eq "-"){
				my $tab = [$child->name,qq{<i class="fa fa-child fa-2x" aria-hidden="true" style="color:black"></i> },$type,"-","-"];
				$html .=  $cgi->start_Tr({class=>"table-info",style=>"border: 1px solid black;color:black;background-color:lightgrey"}).$cgi->td({class=>"table-primary"},$tab).$cgi->end_Tr();
					push(@$all_tab,@$tab);
					#splice @$tab, 1, 1;
				
				$html_line .= join(":",@$tab);
				$html_line .=" - &nbsp;";
				next;
			}
		my $model = $var->getTransmissionModel($fam,$child);
		if ($model =~ /denovo/ && scalar (keys%{$fam->parents()}) ==1){
			$model = "denovo/?";
		}
		$hvariation->{transmission_model} = lc ($model);
		$hvariation->{transmission_model_m} = 1 if $hvariation->{transmission_model} eq "mother" ;
		$hvariation->{transmission_model_p} = 1 if $hvariation->{transmission_model} eq "father" ;#$hvariation->{transmission_model} eq "father";#lc ($model);
		$hvariation->{transmission_model} = "xor" if $hvariation->{transmission_model} eq "mother" or $hvariation->{transmission_model} eq "father";#lc ($model);
		my $color = "white";
		my $model2 = qq{<i class="fa fa-male  fa-2x" style="color:lightgrey"></i><i class="fa fa-female  fa-2x" style="color:lightgrey"></i>};
		$color = "#779ECB" if (lc($model) eq "father");
		$model2 = qq{<i class="fa fa-male  fa-2x style="color:$color""></i>} if (lc($model) eq "father");
		
		$color = "#F7C9C9" if (lc($model) eq "mother");
		$model2 = qq{<i class="fa fa-female  fa-2x style="color:$color""></i>} if (lc($model) eq "mother");
		$color = "#e74c3c" if (lc($model) =~ "denovo");
		$model2 = qq{Denovo} if (lc($model) =~ "denovo");
		$model2 = qq{Strict Denovo} if (lc($model) =~ "strict");
		$color = "#E55137" if (lc($model) eq "denovo/?");	
		$model2 = qq{Denovo/?} if (lc($model) eq "denovo/?");
		$color = "violet" if (lc($model) eq "recessive");
		$model2 = qq{Recessive} if (lc($model) eq "recessive");
		my $tab = [$child->name,qq{<i class="fa fa-child fa-2x" aria-hidden="true" style="color:black"></i> },$type,$var->getPourcentAllele($child)."%",$model2];
	
		
		$html .=  $cgi->start_Tr({class=>"table-info",style=>"border: 1px solid black;color:black;background-color:$color"}).$cgi->td({class=>"table-primary ",nowrap=>"",style=>"background-color:$color"},$tab).$cgi->end_Tr();
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
		
	 # }
	
	 
}

sub return_coverage {
	my ($p,$chr,$start) = @_;
	confess() unless $p;
#	warn $p->name;
	#confess("coucocu") unless $p;
	#confess() unless $p->depth;
	return $p->depth($chr,$start,$start+1)->[0];
	warn $chr."  $start ".$p->name;
	die();
	return $p->depth($chr,$start,$start+1);
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
sub deja_vu2 {
	my ($project,$tr,$hvariation,$hsimilar,$chr) = @_;
	my $vid = $hvariation->{id};
	my $z =0;
	my $nb_similar_project;
	my $nb_similar_patients;
	my $nb_similar_patients_ho;
	foreach my $vs (values %$hsimilar){
		if (exists $vs->{$vid}){
			$nb_similar_project ++;
			my @t = split(" ",$vs->{$vid});
			#if ($t[0]){
			$nb_similar_patients += ($t[0] =~ tr/,//);
			$nb_similar_patients ++;
			#}
		}
		
	}	
		my $no = $project->lite_deja_vu();
	my $h = $no->get($chr,$vid);
	
}



sub deja_vu{
	 my ($project,$tr,$hvariation,$debug) = @_;
	#my $dejavu =$tr->getChromosome->hash_kyoto_dejavu;
	
	 
		my $vid = $hvariation->{id};
		my $debug;
	
		my $similar = $project->similarProjects();
		my $hres = $project->getDejaVuInfosForDiag($vid);
	
		my $nb_pat =0;
		my $pname = $project->name();
		my $proj_dejavu =0;
		$hvariation->{sim_deja_vu} =$hres->{similar_patients};
		$hvariation->{sim_proj_deja_vu} = $hres->{similar_projects};
	
	
		my ($nb_dejavu, $nb_ho) = $project->getDejaVuThisProject($vid);
		die() if $nb_dejavu eq 0;
		$hvariation->{this_deja_vu} = $hres->{in_this_run_patients};
		$hvariation->{diff_project_deja_vu} =$hres->{other_projects}  ; 
		$hvariation->{diff_patient_deja_vu} =$hres->{other_patients}  ;
		$hvariation->{project_deja_vu} = $hres->{other_projects} + $hres->{similar_projects}; 
		$hvariation->{deja_vu} =  $hres->{other_projects}.":". $hres->{other_patients};#printBadge($hres->{other_project}.":".$hres->{other_patient});;
		 $hvariation->{deja_vu_value} = $hres->{other_projects}.":".$hres->{other_patients};
		$hvariation->{in_this_run} =  $hres->{in_this_run_patients}."/". $hres->{total_in_this_run_patients};
		$hvariation->{similar_projects} =  $hres->{similar_projects}.":". $hres->{similar_patients}; #printBadge($hres->{similar_project}.":". $hres->{similar_patient});
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
	 my ($patient,$hvariation,$gene_id) = @_;
	 my $hval     = $patient->validations();
	# warn Dumper $hval;
	 #die();
	# warn "----";
	
	 my $lists = $patient->getListValidatedVariants();
	 #warn Dumper $lists;
	 #die();
	 $hvariation->{type} = "other";
	$hvariation->{sanger} = "-";
	#$hvariation->{user_name} = "";
	my $id =$gene_id."!".$hvariation->{id};
	
	#warn Dumper(keys %$hval);
	#warn $hval->{$id}->{validation_sanger};
	#die($hval->{$id}->{validation_sanger}) if exists $hval->{$id} ;
	return unless exists $hval->{$id};
	my $v =  $hval->{$id}->[0];
	$v->{validation_ngs} = $v->{validation};
	#my $v =  $lists->{$id};
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
			$hvariation->{user_name} = $v->{user_name};
			$hvariation->{modification_date} = $v->{modification_date}; 
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
		elsif (  $v->{validation_ngs} >= 2){
				$hvariation->{user_name} = $v->{user_name};
				$hvariation->{modification_date} = $v->{modification_date};   
				$hvariation->{type} = "validated";
				$hvariation->{type_confirmed} = "ngs";
				$hvariation->{type_confirmed_ngs} = "he";
		}
   

		}
		warn $hvariation->{type}
	
}
 sub tclinical_local {
		my ($project,$hvariation,$patient,$gene) = @_;
		my $val_id = $gene->id."!".$hvariation->{id};
	 my $local_validation = $patient->project->getValidationVariation($val_id,$patient);

		if ($local_validation){
				my $saved = $local_validation->{validation};
				$hvariation->{local} = $patient->buffer->value_validation->{$saved};
				$hvariation->{local} = printButton($saved,[3,5],$patient->buffer->value_validation->{$saved},$saved) ;
				$hvariation->{clinical_local} ++;
		}
		
}
 1;
 
