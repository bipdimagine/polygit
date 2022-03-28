#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag";

#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../packages/cache"; 
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use draw_cnv; 
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/utility.pm";
use html; 
use infos_coverage_exons;
use JSON::XS;
#use image_coverage;
use Carp;
use export_data;
use strict;
use Set::IntSpan::Fast::XS;
use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Carp;
use Set::Intersection;
use Tabix;
use Storable qw/thaw/;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use coverage;
use Spreadsheet::WriteExcel;
use POSIX;
use validationQuery;
use Date::Tiny;
use JSON::XS;
use List::MoreUtils qw{part };

my $cgi          = new CGI();
my $buffer = GBuffer->new();
my $project_name = $cgi->param('project');
my $force_cache =  $cgi->param('force_cache');
my $patient_name = $cgi->param('patients');

my $style_score ={
		1=> {style=>"background-color:#95A5A6;color:#FFFFFF"},#E43725#95A5A6
		2 => {style=>"background-color:#EFB73E;"},#F2DEDE#F7CCC8
		3=> {style=>"background-color:#FF8800"},#E3A57E#F1CE9C#F19B92
		4=> {style=>"background-color:#CE0000;color:white"},#E43725#FF4136
		99=> {style=>"background-color:#2F2F2F;color:#FFFFFF"}
};

my @variant_columns= ("igv","alamut","var_name","genomique","allele","ngs","clinvar");
my @freq_columns = ("freq","freq_ho","max_pop","min_pop","deja_vu","similar_projects","in_this_run","freq_score","impact_score");
my @prediction_columns = ("polyphen","sift","cadd");
my @transcripts_columns= ("transcript","consequence","exon","nomenclature","codons","codons_AA",);
my $project = $buffer->newProjectCache( -name 			=> $project_name, -typeFilters=>'individual' ,-cgi_object=>1);
$project->cgi_object(1);
my $patient = $project->getPatient($patient_name);


#hashtable with name of the column 
# and method is the sub calld for the text in this column

my $format = {
	"igv" => {
		name => "igv",
		method => \&igv#
	} ,
	"alamut" =>{ 
		name => "alamut",
		method => \&alamut#
	},
	"var_name" => {
		name=> "name",
		method=> \&rslink#
	},
	"allele" => {
		name => "allele",
		method=>\&allele#
	},
	"genomique" => {
		name=> "position",
		method =>\&genomique
	},
	"ngs" => {
		name=> "sequencing",
		method =>\&ngs	#
	},
	"clinvar" => {
		name=> "clinvar",
		method =>\&clinvar#
	},
	"freq" => {
		name=> "freq",
		method =>\&freq#
	},
		
	"freq_ho" => {
		name=> "freq_ho",
		method =>\&freq_ho#
	},		
	"max_pop" => {
		name=> "max_pop",
		method =>\&max_pop#
	},			
	"max_pop" => {
		name=> "min_pop",
		method =>\&min_pop#
	},
	"deja_vu" => {
		name=> "deja_vu",
		method =>\&deja_vu#
	},
	"similar_projects" => {
		name=> "similar_projects",
		method =>\&similar_projects#
	},
	"in_this_run" => {
		name=> "in_this_run",
		method =>\&in_this_run#
	},
		"polyphen" => {
		name=> "polyphen",
		method =>\&polyphen#
	},
	"sift" => {
		name=> "sift",
		method =>\&sift#
	},
		"cadd" => {
		name=> "cadd",
		method =>\&cadd#
	},
	"transcript" => {
		name=> "transcript",
		method =>\&transcript#
	},
		"consequence" => {
		name=> "consequence",
		method =>\&consequence#
	},
	"exon" => {
		name=> "exon",
		method =>\&exon#
	},
		"nomenclature" => {
		name=> "nomenclature",
		method =>\&nomenclature
	},
	"codons" => {
		name=> "codons",
		method =>\&codons
	},
		"codons_AA" => {
		name=> "codons_AA",
		method =>\&codons_AA
	},
	
};



##
# all column for the variation all this return only text 
##
sub codons_AA {
	my ($v,$tr,$style) = @_;
		return "-" unless $tr->getProtein;
		return "-" unless $v->isCoding($tr);

	#return "-";
		return  $v->protein_nomenclature($tr->getProtein);
}
sub codons {
	my ($v,$tr,$style) = @_;
		return  $v->getCodons($tr);
}
sub nomenclature {
my ($v,$tr,$style) = @_;
return $v->getNomenclature($tr);
}

sub exon {
my ($v,$tr,$style) = @_;
return $v->nearest_exons($tr);
}

sub consequence {
	my ($v,$tr,$style) = @_;
	return  $v->variationTypeInterface($tr);
	#return "td";
} 

sub transcript {
	my ($v,$tr,$style) = @_;
	return  $tr->name;
	#return "td";
} 
sub cadd{
	my ($v,$tr,$style) = @_;
	return  $v->cadd_score();
	#return "td";
} 

sub polyphen{
	my ($v,$tr,$style) = @_;
	return  $v->polyphenScore($tr->getProtein);
	#return "td";
} 
sub sift {
	my ($v,$tr,$style) = @_;
	return  $v->siftScore($tr->getProtein);
	#return "td";
} 

sub deja_vu {
	my ($v,$style) = @_;
	return "todo";
}
sub similar_projects {
	my ($v,$style) = @_;
	return "todo";
}
sub in_this_run {
	my ($v,$style) = @_;
	return "todo";
}


sub max_pop {
	my($v,$style) = @_;
	return $v->max_pop_name().":". $v->max_pop_freq();
}
sub min_pop {
	my($v,$style) = @_;
	return $v->min_pop_name().":". $v->min_pop_freq();
}
sub freq {
	my($v,$style) = @_;
	my $text ="";
	$text  =  $v->frequency;
	$text  =  0 if $text == -999;
	$text = "-" unless $v->frequency();
	return $text; 
}
sub freq_ho {
	my($v,$style) = @_;
	return  $v->frequency_homozygote;
}
sub clinvar {
		my($variation,$style) = @_;
	 	my $v = $variation->score_clinvar();
	 	return $variation->text_clinvar();
}

sub genomique {
	my($variation,$style) = @_;
	 return $variation->getChromosome()->name.":".$variation->start;
}

sub allele {
	my($variation,,$style) = @_;
	my $a0 = $variation->{ref_allele};
	my $a1 = $variation->sequence();
	return $a1."/".$a0;
}
sub ngs {
	my($variation,$style) = @_;
	my $html = "<small>".$cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered ",style=>"vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;"});
	$html .=  $cgi->start_Tr();
	$html .=  $cgi->td(["uni","100/20","20%"]);
	$html .=  $cgi->end_Tr();
	$html .=  $cgi->start_Tr();
	$html .=  $cgi->td(["free","150/1000","60%"]);
	$html .=  $cgi->end_Tr();
	$html .=  $cgi->start_Tr();
	$html .=  $cgi->td(["sam","200/300","2"]);
	$html .=  $cgi->end_Tr();
	$html.= $cgi->end_table();
	
	$html.= "</small>";
	return $html;
}
sub rslink {
	my($variation,$style) = @_;
	my $html = $variation->name;
	if ($html =~ /rs/){
			my $vn = $variation->name;
			$html = qq{<a href='http://www.ncbi.nlm.nih.gov/snp/?term=$vn' target = '_blank'>$html</a> };	
	}
	return $html;
}
sub igv {
	my($variation,$patient_name,$style) = @_;
	
	my $f = $patient->bamUrl;
		my $l = $variation->{genomique};
		my $v = $variation->{ref_allele}."/".$variation->{allele};	
					
				my $html =qq{<button dojoType="dijit.form.Button"   iconClass="igvIcon" onclick='launch_web_igv("$project_name","$patient_name","$f","$l","$v")' style="color:black"></button>};
	return $html;
}

sub alamut {
	my($variation,,$style) = @_;
		my $start = $variation->{start};
		my $chr = $variation->{chromosome};
		my $vid = $variation->{id};
		my $a0 = $variation->{allele};
		
		my $a1 = $variation->{ref_allele};
			my $html = qq{	<button    class="alamutView3" onClick ="displayInAlamut('$chr',$start,['$a0','$a1']);"></button>};
	return $html;
}

#
#sub deja_vu{
#	 my ($project,$tr,$hvariation) = @_;
#	#my $dejavu =$tr->getChromosome->hash_kyoto_dejavu;
#	
#	 
#	my $vid = $hvariation->{id};
#	#warn Dumper keys %$hvariation;
#	#die();
#	#warn $hvariation->{chr};
#	my $debug;
#	#$debug=1 if  $hvariation->{genomique} eq "5:149776385";
#	
#	my $similar = $project->similarProjects();
#	my $hres = $project->getDejaVuInfosForDiag($vid);
#	
#	my $nb_pat =0;
#	my $pname = $project->name();
#	my $proj_dejavu =0;
#	$hvariation->{sim_deja_vu} = $hres->{similar_patient};
#	$hvariation->{sim_proj_deja_vu} = $hres->{similar_project};
#	
#	
#		my $nb_dejavu = $project->getDejaVuThisProject($vid);
#		die() if $nb_dejavu eq 0;
#		$hvariation->{this_deja_vu} =$nb_dejavu;
#		$hvariation->{diff_project_deja_vu} =$hres->{other_project}  ; 
#		$hvariation->{project_deja_vu} = $hres->{other_project} + $hres->{similar_project}; 
#		$hvariation->{deja_vu} =$hres->{other_project}.":".$hres->{other_patient}; 
#		$hvariation->{in_this_run} = $nb_dejavu."/". scalar(@{$project->getPatients});
#		$hvariation->{similar_projects} = $hres->{similar_project}.":". $hres->{similar_patient};
#		
#		my $freq_score;
#		my $freq_level;
#		
#		if ($hvariation->{freq} < 0 ){
#							$hvariation->{freq_score} = 4; 
#							$hvariation->{freq_level} = 1; 
#							$hvariation->{freq} = 0;
#						
#		}
#		unless ($hvariation->{freq}) {
#				$hvariation->{freq_score} = 4; 
#			$hvariation->{freq_level} = 1; 
#		}
#		elsif ($hvariation->{freq} <= 0.01){
#							$hvariation->{freq_score} = 4; 
#							$hvariation->{freq_level} = 1; 
#		}
#		
#		elsif ($hvariation->{freq} <= 0.05){
#				$hvariation->{freq_score} = 2; 
#					$hvariation->{freq_level} = 3; 
#		}
#		else {
#					$hvariation->{freq_score} = 1; 
#					$hvariation->{freq_level} = 4; 
#			}
#	
#			
#		if ($hvariation->{freq_score} == 4) {
#				if ($hvariation->{diff_project_deja_vu} == 0){
#					$freq_score = 4; 
#					$freq_level =1;
#					$hvariation->{dejavu_score} = 4; 
#					
#				}
#				elsif ($hvariation->{diff_project_deja_vu} <= 20 ){
#							$freq_score = 3; 
#							$freq_level =2;
#							$hvariation->{dejavu_score} = 3; 
#					
#				}
#				elsif ($hvariation->{diff_project_deja_vu} <= 50){
#						$freq_score = 2; 
#						$freq_level =3;
#						$hvariation->{dejavu_score} = 2; 
#					
#				}
#				else {
#					$hvariation->{dejavu_score} = 1; 
#						$freq_score = 1;
#						$freq_level =4; 
#				}
#		}	
#		
#		$hvariation->{freq_score} = $freq_score if $freq_score < $hvariation->{freq_score} && $freq_score;
#		$hvariation->{freq_level} = $freq_level if $freq_level > $hvariation->{freq_level} && $freq_level ;
#		
#}




#######################################
#construct_hash_variations($project,$patient);
#exit(0);
######################################
my $out;
 html::print_cgi_header($cgi,$out,0,$patient_name." - PolyDiag");
print_table_html();
exit(0);



sub print_table_html {
		$project->setBundle();
		
		my $gs = $patient->getGenes();
		print start_table();
		print header_table_variations([@variant_columns,@freq_columns,@transcripts_columns]);
	
		foreach my $g (@$gs){
		#mouais a voir si c'ets la bonne idée  
		# retourne les variant de ce gene present sur ce patient ?
		my  $vs = $g->getFilteredVariants($patient);
		print print_panel_genes($g,$vs);
		print start_table();
		print header_table_variations([@variant_columns,@freq_columns,@transcripts_columns]);
		#my $vs = $g->getVariations();
		foreach my $v (@$vs) {
			#warn $v->name();
				#mouais bis a voir si c'ets la bonne idée ?
				# retourne les transcripts du variants qui appartiennent au gene 
			my $trs = $v->getTranscripts($g);
			
				my $rowspan = scalar(@$trs);
			print print_start_tr_variation();
			print print_td_variation_annotations($v,\@variant_columns,$rowspan);
			print print_td_variation_annotations($v,\@freq_columns,$rowspan);
			#my $trs = $v->getTranscripts($g);
			foreach my $tr (@$trs){
					print print_td_transcripts_annotations($v,$tr,\@transcripts_columns,1);
						print $cgi->end_Tr();
						#last;
			}
		print $cgi->end_Tr() unless $trs;
		}
			print $cgi->end_table();
			print $cgi->end_div();
		}
}

 


exit(0);



sub print_header_edition {
	my ($name,$date,$since) = @_;
	my $out;
	 my $icon = $cgi->span({class=>"glyphicon glyphicon-eye-open  pull-left",'aria-hidden'=>"true"});
  my $icon_help = $cgi->span({class=>"glyphicon glyphicon-question-sign pull-left",'aria-hidden'=>"true"});
  my $icon_calendar = $cgi->span({class=>"glyphicon glyphicon-calendar",'aria-hidden'=>"true"});
  my $icon_export =  $cgi->span({class=>"glyphicon glyphicon-open-file pull-left",'aria-hidden'=>"true"});
  
$out  .= qq{
<div class="panel panel-primary">
  <div class="panel-heading clearfix" style="min-height=30px;max-height:30px;">
    <h3 class="panel-title pull-left" style="padding-top: 1.5px;font-size: 15px;">Edition &nbsp; $name  &nbsp; &nbsp; &nbsp; $icon_calendar  Last Update : $date ($since days ago)</h3>
   <div class="btn-group  btn-sm clearfix pull-right" style="position:relative;top:-14px;float:right">
		 		<div class=" btn btn-info btn-sm " aria-label="Left Align"  onClick='collapse_all("","");' >
		 			$icon
		 		</div>
		 		<div class=" btn btn-info btn-sm" aria-label="Left Align" onClick='dijit.byId("dialog_help_1").show();' >
		 			$icon_help
		 		</div>
		 		<div class=" btn btn-danger btn-sm " aria-label="Left Align" onClick='editor(1,2);' >
		 				<b><i class="fa fa-file-excel-o pull-left"></i></b>
		 		</div>
		</div>
  </div>
</div>  

	};
	return $out;
	
}

sub print_panel_genes{
	my ($gene,$variations) = @_;
#method print_panel_genes( GenBoGeneCache :$gene!){

	#my ($name,$clinvar_alert,$max_impact,$max_score,$max_freq,$nb_impact,$saved,$todo) = @_;
	my $out;
	my $name = $gene->external_name();
	my $clinvar_alert = 1;
	my $max_score = 1;
	my $max_freq =1;
	my $nb_impact = [1,1,1,1];
	my $max_impact =1;
	my $saved =1;
	my $todo = 1;
	
	$out .=  $cgi->start_div({class=>"panel panel-success" });
	$out .=  $cgi->start_div({class=>"panel-heading panel-warning warning ",style=>" min-height:13px;max-height:13px;padding:1px;border:1px"});
		my $glyph = "";
		$glyph = qq{<span class="glyphicon glyphicon-star-empty text-default" aria-hidden="true"></span>} if $clinvar_alert == 1 ;
		$glyph = qq{<span class="glyphicon  glyphicon-alert text-warning" aria-hidden="true"></span>}  if $clinvar_alert == 2 ;
		$out .=  $cgi->start_div({class=>" btn-group btn-xs "});
		my $label_id = $name."_label";
		my $panel_id = $name."_panel";
			$out .= qq{<div class="btn  btn-info btn-xs " style="position:relative;bottom:1px;min-width:150px;" onClick='collapse("$panel_id","$label_id")'>  <span id= "$label_id" class="glyphicon glyphicon-triangle-right  " aria-hidden="true"  style="float:left;"></span> $name &nbsp $glyph</div>};

				#$out .=$cgi->span({class=>"label label-success"},$nbc) if $nbc >0;
				my $total = scalar(@{$variations});
			  	$out .=$cgi->span({class=>"label label-success"},qq{<span class='badge badge-primary badge-xs '  >$total </span> });
				$out .=$cgi->span({class=>"label label-danger"}, $nb_impact->[2]);
				$out .=$cgi->span({class=>"label label-warning"}, $nb_impact->[1]);
				$out .=$cgi->span({class=>"label label-default"},$nb_impact->[0]);
			
				
	   		$out.= $cgi->end_div();
	   		$out .=  $cgi->start_div({class=>" btn-group btn  ",style=>'position:relative;float:right;bottom:5px;'});
	   			 #write SIF LABEL 
					my %style = %{$style_score->{$max_score}};
					$style{class} ="label  "; 
					$out .=$cgi->span(\%style,"S"); 
					
					 %style =%{$style_score->{$max_impact}};
					$style{class} ="label  "; 
					$out .=$cgi->span(\%style,"I");
					
					%style = %{$style_score->{$max_freq}};
					$style{class} ="label  "; 
					$out .=$cgi->span(\%style,"F");
				#END SIF
				#start menu hamburger edit and saved ;-)
					$out .=$cgi->span({class=>"label label-success"},qq{<span class="glyphicon glyphicon glyphicon-menu-hamburger " aria-hidden="true" "></span>});
					my $clabel = " label-default";
					$clabel = "label-info" if $todo;
					$out .=$cgi->span({class=>"label $clabel"},qq{<span class="glyphicon glyphicon-pencil " aria-hidden="true" "></span>});
					$clabel = " label-default";
					$clabel = "label-info" if$saved;
					$out .=$cgi->span({class=>"label $clabel"},qq{<span class="glyphicon glyphicon-saved" aria-hidden="true" "></span>});		
		$out.= $cgi->end_div(); # end div lavel right SIF
		$out.= $cgi->end_div(); # end panel heading 
		return $out;
}

sub start_panel_collapse {
	my ($panel_id) = @_;
	return $cgi->start_div({class=>"panel-body panel-collapse collapse ",style=>"font-size: 09px;font-family:  Verdana;",id=>$panel_id});
}
sub start_table {
	return $cgi->start_table({"data-toggle"=>"table","data-search"=>"true",class=>"table table-striped table-condensed table-bordered table-hover table-mybordered",style=>"vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;"});
}


sub header_table_variations {
	my ($header) = @_;
	my $out;
	$out.= $cgi->start_Tr();
	foreach my $c (@$header){
		$out.= qq{<th >$c</th>};
	}
	$out.= $cgi->end_Tr();
	return $out;
}

sub status_code {
	my $out;
	
	return $out;
}


sub print_start_tr_variation {
	return  $cgi->start_Tr({class=>"menurow" ,style=>"vertical-align:middle"});
}
sub print_td_variation_annotations {
	my ($variation,$infos,$rowspan) = @_;

	my $bg = "background-color:#F0F0F0";
	my $out;
	#$out.= #start TR;
	foreach my $info (@$infos){
		my $text = $variation->{$info};
		if (exists $format->{$info}){
			$text = $format->{$info}->{method}->($variation);
		}
			$out .= $cgi->td({rowspan=>$rowspan,nowrap=>1,style=>$bg} , $text);
	}
	return $out;
}

sub print_td_transcripts_annotations {
	my ($variation,$tr,$infos,$rowspan) = @_;

	my $bg = "background-color:#F0F0F0";
	my $out;
	#$out.= #start TR;
	foreach my $info (@$infos){
		my $text = $variation->{$info};
		if (exists $format->{$info}){
			$text = $format->{$info}->{method}->($variation,$tr);
		}
			$out .= $cgi->td({rowspan=>$rowspan,nowrap=>1,style=>$bg} , $text);
	}
	return $out;
}


#variation annotations

#cosmic
sub cosmic_annotations {
	my ($v) = @_;
	my $id  = $v->cosmic();
		if ($id){
			my ($id1,$nb) = split(":",$id);
			my $text = $id1;
			my $url = "http://cancer.sanger.ac.uk/cosmic/mutation/overview?id="; 
			 $url = "http://grch37-cancer.sanger.ac.uk/cosmic/ncv/overview?id="  if ($id1 =~/COSN/);
			  $text = "$id1 [$nb]" if $nb;
			$id1 =~s/COSM//;
			$id1 =~s/COSN//;
			return  qq{<a href=\"$url$id1\" target="_blank">$text</a>};
		}
		return "-";
}

#functional annotation on transcript
sub functional_annotations {
	my ($v,$tr1,$hv) = @_;
		my $hvariation;
		$hvariation->{impact_text} = $v->effectImpact($tr1);
		#return;
		#$hvariation->{impact_score} = $himpact_sorted->{$v->effectImpact($tr1)};
		$hvariation->{gene} = $tr1->getGene->external_name();
		
		$hvariation->{transcript} = $tr1->name;
		
		$hv->{trans} = $v->start * $tr1->strand; 
		
		$hvariation->{cds} = "";
		$hvariation->{prot} ="";
		$hvariation->{codons_AA} = "";
		$hvariation->{polyphen} ="-";
		$hvariation->{sift} ="-";
		
		 if ($v->isCoding($tr1)) {
			my $prot = $tr1->getProtein();
		
			$hvariation->{cds} = $v->getOrfPosition($prot);
			$hvariation->{prot} = $v->getProteinPosition($prot);
			$hvariation->{codons_AA} =   $v->protein_nomenclature($prot);#$v->getProteinAA($prot).$hvariation->{prot}.$v->changeAA($prot);
			$hvariation->{polyphen} = $v->polyphenScore($tr1->getProtein);
			$hvariation->{sift} = $v->siftScore($tr1->getProtein);
			$hvariation->{codons} =   $v->getCodons($tr1);
		}
		$hvariation->{exon} = $tr1->findExonNumber($v->start);
		$hvariation->{exon} = $tr1->findNearestExon($v->start) if $hvariation->{exon} == -1;
		$hvariation->{nomenclature} =  $v->getNomenclature($tr1);
		$hvariation->{consequence} =  $v->variationType($tr1);
		if ($tr1->strand() == -1 ){
			$hvariation->{codons}  =  BioTools::complement_sequence($v->getChromosome()->sequence($v->start,$v->end))."/".BioTools::complement_sequence($v->sequence());
		}
		
		$hv->{functional}->{$tr1->id} = $hvariation;
}



sub famillial_annotation {
	my ($v,$hvariations) = @_;
	foreach my $patient (@{$project->getPatients}){
		
	}	
}




### print code html 

sub print_header_edition {
	my ($name,$date,$since) = @_;
	my $out;
	 my $icon = $cgi->span({class=>"glyphicon glyphicon-eye-open  pull-left",'aria-hidden'=>"true"});
  my $icon_help = $cgi->span({class=>"glyphicon glyphicon-question-sign pull-left",'aria-hidden'=>"true"});
  my $icon_calendar = $cgi->span({class=>"glyphicon glyphicon-calendar",'aria-hidden'=>"true"});
  my $icon_export =  $cgi->span({class=>"glyphicon glyphicon-open-file pull-left",'aria-hidden'=>"true"});
  
$out  .= qq{
<div class="panel panel-primary">
  <div class="panel-heading clearfix" style="min-height=30px;max-height:30px;">
    <h3 class="panel-title pull-left" style="padding-top: 1.5px;font-size: 15px;">Edition &nbsp; $name  &nbsp; &nbsp; &nbsp; $icon_calendar  Last Update : $date ($since days ago)</h3>
   <div class="btn-group  btn-sm clearfix pull-right" style="position:relative;top:-14px;float:right">
		 		<div class=" btn btn-info btn-sm " aria-label="Left Align"  onClick='collapse_all("","");' >
		 			$icon
		 		</div>
		 		<div class=" btn btn-info btn-sm" aria-label="Left Align" onClick='dijit.byId("dialog_help_1").show();' >
		 			$icon_help
		 		</div>
		 		<div class=" btn btn-danger btn-sm " aria-label="Left Align" onClick='editor(1,2);' >
		 				<b><i class="fa fa-file-excel-o pull-left"></i></b>
		 		</div>
		</div>
  </div>
</div>  

	};
	return $out;
	
}
