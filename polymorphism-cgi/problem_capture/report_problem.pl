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
use image_coverage;
#use Set::;
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
use coverage;
use validationQuery;
use Date::Tiny;
use JSON::XS;
use Storable qw(store retrieve freeze thaw);
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";

my @headers = ("var_name","genomique","transcript","exon","nomenclature","consequence","codons","codons_AA","freq","freq_ho","max_pop","min_pop","clinvar","deja_vu","similar_projects","in_this_run", "polyphen","sift","cadd");
my @variant_columns= ("alamut","var_name","genomique","all","clinvar");
my @freq_columns = ("freq","freq_ho","max_pop","min_pop","deja_vu","similar_projects","in_this_run","freq_score","impact_score");
my @prediction_columns = ("polyphen","sift","cadd");
my @transcripts_columns= ("transcript","consequence","exon","nomenclature","codons","codons_AA",);
my @patients_columns= ("project","patient","sex","ngs","ratio","caller","igv");

my $style_score ={
		1=> {style=>"background-color:#F2F183"},#E43725#95A5A6rgb(255, 250, 205)##FFFACD
		2 => {style=>"background-color:#EFB73E;"},#F2DEDE#F7CCC8
		3=> {style=>"background-color:#FF8800"},#E3A57E#F1CE9C#F19B92
		4=> {style=>"background-color:#CE0000;color:white"},#E43725#FF4136
		99=> {style=>"background-color:#2F2F2F;color:#FFFFFF"}
};

my $style_td ={
		1=> {class=>"alert",style=>"background-color:#FDFDFD;"},#E43725#95A5A6
		2 => {class=>"warning2",style=>"background-color:#F7CCC8"},#F2DEDE#F7CCC8
		3=> {class=>"warning2",style=>"background-color:#F19B92"},#E3A57E#F1CE9C#F19B92
		4=> {class=>"danger",style=>"background-color:#FF4136"},#E43725#FF4136
};

my $himpact= {
	"1" => "high",
	"2" =>"moderate",
	"3" =>"low",
};
my $himpact_sorted_inv= {
	"4" => "high",
	"3" =>"moderate",
	"1" =>"low",
};

my $himpact_sorted = {
	"high" => "4",
	"moderate" =>"3",
	"low" =>"1",
};

my $himpact2= {
	"1" => "high",
	"2" =>"medium",
	"3" =>"low",
};

my $hfrequence= {
	"1" => "unique",
	"2" =>"rare",
	"3" =>"occasional",
	"4" =>"all"
};

my $hfrequence_sorted = {
	"unique" => 4 ,
	"rare" => 3,
	"occasional" =>2,
	"all"=>1,
};

my $hscore_sorted = {
	"high" => 4 ,
	"medium" => 3,
	"low" =>2,
};



my $w = "500px";

my $CSS = qq{
<style type="text/css"> 
body
{
	line-height: 1.0em;
}

.warning2{
	background-color:  #D9534F;
}
a:link { 
    text-decoration: none;
     color: red;
}
a:visited { 
    text-decoration: none;
     color: green;
}
a:active { 
    text-decoration: none;
}
#test
{
	font-family: "Lucida Sans Unicode", "Lucida Grande", Sans-Serif;
	font-size: 12px;
	background: #00;
	margin: 6px;
	width: 100px;
	border-collapse: collapse;
	text-align: left;
	vertical-align: top;
	page-break-after: always
}
#test th
{
	vertical-align: top;
	
	
}

#test td
{
	vertical-align: top;
}

#coverage
{
	font-family: "Lucida Sans Unicode", "Lucida Grande", Sans-Serif;
	font-size: 10px;
	background: #fff;
	margin: 9px;
	width: 140px;
	border-collapse: collapse;
	text-align: left;
}
#coverage th
{
	font-size: 10px;
	font-weight: normal;
	color: #039;
	padding: 4px 3px;
	border-bottom: 2px solid #6678b1;
}
#coverage th1
{
	font-size: 8px;
	font-weight: normal;
	color: #039;
	padding: 4px 3px;
	border-bottom: 2px solid #6678b1;
}
#coverage td
{
	border-bottom: 1px solid #ccc;
	color: #669;
	padding: 3px 4px;
}
#coverage tbody tr:hover td
{
	color: #009;
}


#hor-minimalist-b
{
	font-family: "Lucida Sans Unicode", "Lucida Grande", Sans-Serif;
	font-size: 10px;
	background: #fff;
	margin: 9px;
	width: 1000px;
	border-collapse: collapse;
	text-align: left;
	
}
#hor-minimalist-b th
{
	font-size: 8px;
	font-weight: normal;
	color: #039;
	padding: 4px 3px;
	border-bottom: 2px solid #6678b1;
	width : $w;
	
	
}
#hor-minimalist-b th1
{
	font-size: 8px;
	font-weight: normal;
	color: #039;
	padding: 4px 3px;
	border-bottom: 2px solid #6678b1;
	width : $w;
	
}
#hor-minimalist-b td
{
	border-bottom: 1px solid #ccc;
	color: #669;
	padding: 3px 4px;
	width : $w;
	
}
#hor-minimalist-b .td1
{
	font-size: 8px;
	border-bottom: 1px solid #ccc;
	color: #000;
	width : 200px;
	
}

#hor-minimalist-b .tdorange2
{
	border-bottom: 1px solid #ccc;
	color: #BB8200;
	padding: 3px 4px;
	width : $w;	
}

#hor-minimalist-b .tdorange1
{
	border-bottom: 1px solid #ccc;
	color: #A52A2A;
	padding: 3px 4px;
	width : $w;	
}

#hor-minimalist-b tbody tr:hover td
{
	color: #009;
	width : $w;

}

.spanpatient{
	font-family: Georgia, "Times New Roman", Times, serif;
        font-size:13px;
		text-transform: uppercase;
        font-weight: normal;
        color: #222;
        letter-spacing: 0.2em;
}

h1.patient  {
    font-family: Georgia, "Times New Roman", Times, serif;
        font-size:13px;
	margin-top: 5px; margin-bottom: 0px;
		text-transform: uppercase;
        font-weight: normal;
        color: #222;
        letter-spacing: 0.2em;
	
 }
 h5.patient  {
    font-family: Georgia, "Times New Roman", Times, serif;
        font-size:13px;
	margin-top: 5px; margin-bottom: 0px;
		text-transform: uppercase;
        font-weight: bold;
        color: red;
        letter-spacing: 0.2em;
	
 }
h2.patient  {
font-family: "Lucida Grande", Tahoma;
	font-size: 12px;
	font-weight: lighter;
	font-variant: normal;
	text-transform: uppercase;
	color: #FF3960;
    margin-top: 5px;
	text-align: center!important;
	letter-spacing: 0.2em;
 }
 h3.patient  {
font-family: "Lucida Grande", Tahoma;
	font-size: 14px;
	font-weight: lighter;
	font-variant: normal;
	text-transform: uppercase;
	color: #458B2A;
       margin-top: 10px;
	text-align: center!important;
	letter-spacing: 0.3em;
 }
 h4.patient{
 	page-break-before: always
 }
 
 

 legend {
	
    font-weight: normal;
     font-family: Georgia, "Times New Roman", Times, serif;
     font-size:19px;
     padding: 5px;
 }
fieldset {
-webkit-border-radius: 8px;
-moz-border-radius: 8px;
border-radius: 8px;
}

.page-break	{ display: block; page-break-before: always; }

.menurow:hover a{
    background-color: blue;
}
.badge-error {
  background-color: #b94a48;
}
</style>

};



my $buffer = GBuffer->new();
my $cgi          = new CGI();
my $server =  'https://'.$ENV{HTTP_HOST};
my $variation_script = $ENV{SCRIPT_NAME};
$variation_script =~s/patient_/variation_/;
$server = "darwin.bipd.fr" if $server eq "bipd";
$server = "www.polyweb.fr" if $server =~/10\.200\.27/;
my $deja_vu_url = "$server//polyweb/polydejavu/dejavu.html?input=";
#html::print_cgi_header($cgi,$CSS.$out_global,$print,$patient_name." - PolyDiag");
html::print_cgi_header($cgi,$CSS,0,"  PolyDiag");


#read all parameter
my $project_name = $cgi->param('project');
my $capture_name =  $cgi->param('capture');

#alele quality 
my $filter_quality = $cgi->param('allele_quality');

my $limit_ratio = 0;
$limit_ratio = "80" if $filter_quality == 1;
$limit_ratio = "40" if $filter_quality == 2;
$limit_ratio = "20" if $filter_quality == 3;
$limit_ratio = "10" if $filter_quality == 4;

#impact filter 
my $vimpact = $cgi->param('impact');

my $impact_score_limit = $himpact_sorted->{$himpact->{$vimpact}};
warn $impact_score_limit." ".$vimpact;
#die();
#frequence filter
my $vfreq = $cgi->param('frequence');

my $projects;
 my $buffer1 = GBuffer->new();	
my @projectAll = @{$buffer1->listProjectsByAnalyse($capture_name)};
#my @projectAll = ("NGS2016_1402","NGS2016_1108","NGS2016_1222","NGS2016_1046");
foreach my $project_name (@projectAll){
	#next if $project_name !~/NGS2016/;# and $project_name !~/NGS2012/ and $project_name !~/NGS2011/ ;
	#next if $project_name =~/NGS2010/;
	my ($n1,$n2) = split("_",$project_name);
	$n1 =~ s/NGS//;
	#next  if -e "$n1/$project_name.ok";
	my $buffer = GBuffer->new();
	my $project = $buffer->newProject(-name=>$project_name);
	push(@$projects,$project_name);
}
my $all;
my $transcripts = {};
my $variations = {};
my $patients = {};
my $genes ={};
my $genes_variations;
my $transcripts_variations;
my $project_problem = [];
my $project_variations;
my $patients_variations;
my $running_projects;
foreach my $project_name (@$projects){
	my $project = $buffer->newProject(-name=>$project_name);
	my $dir1 = $project->project_root_path()."HG19/variations/cp/";
	my @files = `ls $dir1/*.gz 2>/dev/null`;
	$running_projects->{$project_name}->{calling_poor} ++  if (scalar(@files) eq scalar(@{$project->getPatients}));
	
	
	#warn $project_name if (scalar(@files) eq scalar(@{$project->getPatients}));
	my $output = $project->getCacheDir() . "/miss/";
	#system("mkdir $output && chmod a+rwx $output") unless -e $output;
	my $file_done = "$output/$project_name.done";
	my $file_freeze = "$output/$project_name.freeze";
	$running_projects->{$project_name}->{miss} ++  if (-e $file_done);
	next unless -e $file_freeze;

	my $a = retrieve($file_freeze);
	next unless scalar(keys %{$a->{genes}});
	push (@$project_problem,$project_name);
	
	
	#genes_variations 
	
		foreach my $gid (keys %{$a->{genes_variations}}){
		
			foreach my $vid (keys %{$a->{genes_variations}->{$gid}}){
				#	$project_variations->{$project_name}->{$vid} ++;
				$genes_variations->{$gid}->{$vid} ++;
			}
			
		}
	
	foreach my $gid (keys %{$a->{transcripts_variations}}){
			foreach my $vid (keys %{$a->{transcripts_variations}->{$gid}}){
				$transcripts_variations->{$gid}->{$vid} ++;
			}
			
		}
	#	warn Dumper 	$a->{transcripts_variations}->{ENST00000497038_2};
		foreach my $vid (keys %{$a->{genes}}) {
			#warn Dumper $a->{genes}->{$vid}->{name} if  $a->{genes}->{$vid}->{name} =~ /ZN/;
		$genes->{$vid} = $a->{genes}->{$vid};
	}

	
	foreach my $vid (keys %{$a->{transcripts}}) {
		$transcripts->{$vid} = $a->{transcripts}->{$vid};
	}
	
	foreach my $vid (keys %{$a->{patients}}) {
	
		$patients->{$vid} = $a->{patients}->{$vid};
		$patients->{$vid}->{id} = $vid;
	}
	foreach my $vid (keys %{$a->{variations}}) {
		unless (exists $variations->{$vid}){
			$variations->{$vid} = $a->{variations}->{$vid};
		}
		else {
			foreach my $p (keys %{$a->{variations}->{$vid}->{patients}}){
				$variations->{$vid}->{patients}->{$p} = $a->{variations}->{$vid}->{patients}->{$p};
			}
		}
		$variations->{$vid}->{projects}->{$project_name} ++;
	}
	


	#push(@{$all->{variations}},@{$a->{variations}});
	
}

#filter variations

my $filtering_variations ={};
warn Dumper $variations->{"17_78064059_C_T"};
warn "---------------------------";
foreach my $vid (keys %$variations){
	my $variation =  $variations->{$vid};
	my $debug;
	$debug = 1 if $variation->{var_name} eq "rs1049948";
	 if ($variation->{freq_level} > $vfreq) {
	#if ($variation->{freq_level} < $vfreq ){
		$filtering_variations->{$vid} ++;
		warn "coucou" if $debug ==1;
		next;
	}

	foreach my $gid (keys %{$variation->{functional}}){
		
			foreach my $trid (keys %{$variation->{functional}->{$gid}}){
				$variation->{functional}->{$gid}->{$trid}->{impact_score} = $himpact_sorted->{ $variation->{functional}->{$gid}->{$trid}->{impact_text}};
				
					#change scaled score
					my $sc2 = $variation->{functional}->{$gid}->{$trid}->{impact_score} + $variation->{freq_score};
					$variation->{functional}->{$gid}->{$trid}->{scaled_score}  = 0;
				if ($sc2 >= 7 ){
						$variation->{functional}->{$gid}->{$trid}->{scaled_score}  =4;
				}
				elsif  ($sc2 >= 6){
			
					if ($variation->{cadd} > 20){
						$variation->{functional}->{$gid}->{$trid}->{scaled_score}  = 4;
					}
					elsif  ($variation->{cadd} > 14 or $variation->{functional}->{$gid}->{$trid}->{impact_score} >=3){
						$variation->{functional}->{$gid}->{$trid}->{scaled_score}  = 3;
					}
					else {
						$variation->{functional}->{$gid}->{$trid}->{scaled_score}  = 2;
					}
			}
			elsif  ($sc2 >= 5 && $variation->{functional}->{$gid}->{$trid}->{impact_score} >=3){
					$variation->{functional}->{$gid}->{$trid}->{scaled_score}  = 1;
			}
					
				if ($variation->{functional}->{$gid}->{$trid}->{scaled_score} >1 && scalar(keys %{$variation->{patients}}) > 3){
					$variation->{functional}->{$gid}->{$trid}->{scaled_score}  --;
				}
			if ($variation->{functional}->{$gid}->{$trid}->{scaled_score} >1 && scalar(keys %{$variation->{patients}}) > 6){
					$variation->{functional}->{$gid}->{$trid}->{scaled_score}  --;
				}
					
							
					if ($variation->{functional}->{$gid}->{$trid}->{impact_score}< 	$impact_score_limit	){
						delete  $transcripts_variations->{$trid}->{$vid};
						delete $variation->{functional}->{$gid}->{$trid};
							#$filtering_variations->{$vid} ++;
						#next;
					}
				#	warn "-->  $vid ".$gid ." ".Dumper keys %{$variation->{functional}->{$gid}} if $debug;
					unless (keys %{$variation->{functional}->{$gid}}){
					#	warn "-->  delete ". $gid if $debug;
						delete $variation->{functional}->{$gid};
					}
			}
		
			$filtering_variations->{$vid} ++  unless (keys %{$variation->{functional}});
			my @find;
		foreach my $gid (keys %{$variation->{patients}}){
			foreach my $trid (keys %{$variation->{patients}->{$gid}}){
					my @all_nums    = $variation->{patients}->{$gid}->{ratio} =~ /(\d+)/g;
					
					if (scalar(@all_nums)){
						push(@find,grep{$_>15} @all_nums);
						my @t = grep{$_>$limit_ratio} @all_nums;
						delete  $variation->{patients}->{$gid} unless @t;;
					}
			}
			$filtering_variations->{$vid} ++  unless scalar (keys %{$variation->{patients}->{$gid}});
			unless (@find) {
					foreach my $gid (keys %{$variation->{functional}}){
		
						foreach my $trid (keys %{$variation->{functional}->{$gid}}){
							if ($variation->{functional}->{$gid}->{$trid}->{scaled_score} >1){
								$variation->{functional}->{$gid}->{$trid}->{scaled_score} --;
							}
						}
					}
			}	
		}
			
	}
	
}

warn "coucou " if exists $filtering_variations->{"17_78064059_C_T"};
foreach my $trid (keys %$transcripts_variations){
	foreach my $vid (keys %$filtering_variations){
		
		delete $transcripts_variations->{$trid}->{$vid};
	}
	
}


foreach my $vid (keys %$variations){
	next if exists $filtering_variations->{$vid};
		my $variation =  $variations->{$vid};
	foreach my $gid (keys %{$variation->{functional}}){
			foreach my $trid (keys %{$variation->{functional}->{$gid}}){
				my $impact =	$variation->{functional}->{$gid}->{$trid}->{scaled_score} ;
				
					foreach my $p (keys %{$variation->{projects}}){
					
					$project_variations->{$p}->{impact}->{$impact} ++;
				}
				my $g = $transcripts->{$trid}->{gene};
				foreach my $p (keys %{$variation->{patients}}){
					
					$patients_variations->{$p}->{impact}->{$impact}->{$g} ++;
				}
				
				
				
			}
	}
}

#project Statistics 


warn Dumper $variations->{"17_78064059_C_T"};
my $out = $CSS;
#start table

my $sum =  scalar(@freq_columns) +  scalar(@prediction_columns) + scalar(@variant_columns) + scalar(@transcripts_columns) +  scalar(@patients_columns);

my $string_panel = "";#join(";",@$all_panel);
my $string_label = "";#join(";",@$all_label);
  my $icon = $cgi->span({class=>"glyphicon glyphicon-eye-open  pull-left",'aria-hidden'=>"true"});
  my $icon_help = $cgi->span({class=>"glyphicon glyphicon-question-sign pull-left",'aria-hidden'=>"true"});
  my $icon_calendar = $cgi->span({class=>"glyphicon glyphicon-calendar",'aria-hidden'=>"true"});
  my $icon_export =  $cgi->span({class=>"glyphicon glyphicon-open-file pull-left",'aria-hidden'=>"true"});
  my $text_title;
  foreach my $p (sort {$a cmp $b} keys %{$project_variations}){
  	
  		$text_title.="$p:".$project_variations->{$p}->{impact}->{4}."-".$project_variations->{$p}->{impact}->{3}."-".$project_variations->{$p}->{impact}->{2}."-".$project_variations->{$p}->{impact}->{1}."-".$project_variations->{$p}->{impact}->{0}."&nbsp";
  	
  }
  
	#my $text_title = join("&nbsp-&nbsp",@$project_problem);
	
	$out  .= qq{
<div class="panel panel-primary">
  <div class="panel-heading clearfix" style="min-height=30px;max-height:30px;">
    <h3 class="panel-title pull-left" style="padding-top: 1.5px;font-size: 15px;">$capture_name &nbsp;   &nbsp; &nbsp; &nbsp; $icon_calendar  </h3>
   <div class="btn-group  btn-sm clearfix pull-right" style="position:relative;top:-14px;float:right">
		 		<div class=" btn btn-info btn-sm " aria-label="Left Align"  onClick='collapse_all("$string_panel","$string_label");' >
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
	
		$out .= qq{<div class="btn  btn-info btn-xs " style="position:relative;bottom:1px;min-width:200px;" onClick='collapse("projects_panel","projects_label")'>  <span id= "projects_label" class="glyphicon glyphicon-triangle-right  " aria-hidden="true"  style="float:left;"></span> Projects</div>};
	
		$out .=  $cgi->start_div({class=>"panel-body panel-collapse collapse ",style=>"font-size: 09px;font-family:  Verdana;",id=>"projects_panel"});	
	$out .= $cgi->start_table({class=>"table table-responsive table-striped table-condensed table-bordered table-hover table-mybordered",style=>"vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;"});
	my $td_style ={
				0=>"background:red",
				1=>"background:green",
	};
	
	foreach my $project (keys %$running_projects){
		$out .= $cgi->start_Tr();
		$out .= $cgi->td($project);
		
		$out .= $cgi->td({style=>$td_style->{$running_projects->{$project}->{calling_poor}+0}}," ");
		$out .= $cgi->td({style=>$td_style->{$running_projects->{$project}->{miss}+0}}," ");
			$out .= $cgi->end_Tr();
	}
	
	$out.= $cgi->end_table;
	$out.=$cgi->end_div;	

	my $panel_id ="patients_panel";
	$out .= qq{<div class="btn  btn-info btn-xs " style="position:relative;bottom:1px;min-width:200px;" onClick='collapse("$panel_id","patients_label")'>  <span id= "patients_label" class="glyphicon glyphicon-triangle-right  " aria-hidden="true"  style="float:left;"></span> Patients</div>};
	
	$out .=  $cgi->start_div({class=>"panel-body panel-collapse collapse ",style=>"font-size: 09px;font-family:  Verdana;",id=>"patients_panel"});	

$out .= $cgi->start_table({class=>"table table-striped table-condensed table-bordered table-hover table-mybordered",style=>"vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;"});

my $hprojects;

map{$hprojects->{$patients->{$_}->{project}} ++} keys  %$patients_variations;

foreach my $pr (keys %$hprojects){
my $rowspan = $hprojects->{$pr}+1;
$out .= $cgi->td({rowspan=>$rowspan},$pr);
my @pid = grep {$patients->{$_}->{project} eq $pr} keys  %$patients_variations;
#sort{scalar(keys %{$patients_variations->{$b}->{impact}->{4}}) <=> scalar(keys %{$patients_variations->{$a}->{impact}->{4}})} }
foreach my $p ( @pid) {
	my $patient = $patients->{$p};

	my  $text="";
	$out .= $cgi->start_Tr();
	$out .= $cgi->td($patient->{name});
	$out .= $cgi->td($patient->{sex_icon});
		for (my$i=4;$i>=0;$i--){
  			my $bg =$style_score->{$i}->{style};
  			my @genes = keys %{$patients_variations->{$p}->{impact}->{$i}};
  			
  			$out.= $cgi->td({style=>$bg},join(";",@genes));
		}
		#$out .= $cgi->li({class=>"list-group-item list-group-item-info"}, $p." ".$text);
	$out .= $cgi->end_Tr();
}
}

$out.= $cgi->end_table();
$out.=$cgi->end_div;	
$out.="<BR>";
$out.="<HR>";

foreach my $gid (keys %{$genes}){
	
	foreach my $trid (keys %{$genes->{$gid}->{transcripts}}) {
		
		my $tvariations =[];
		foreach my $vid (keys %{$transcripts_variations->{$trid}}){
			warn "coucou $trid $gid " if  $vid eq "17_78064059_C_T";
			push(@$tvariations,$variations->{$vid});
		}
	#warn scalar(@$tvariations);	
	 unless (@$tvariations){
	 	delete $genes->{$gid};
	 	next;
	 }
		$genes->{$gid}->{$trid}->{title}	 = start_gene_line($gid,$trid,$transcripts->{$trid},$tvariations);	
}
}
	
foreach my $gid (sort {$genes->{$b}->{score} <=> $genes->{$a}->{score}} keys %{$genes}){
	my $np = $genes->{$gid}->{name};
	$out.=qq{<span class="label label-primary label-xs "><big>$np</big></span>};
	foreach my $trid (keys %{$genes->{$gid}->{transcripts}}) {
		
		my $tvariations =[];
		foreach my $vid (keys %{$transcripts_variations->{$trid}}){
			push(@$tvariations,$variations->{$vid});
		}
	#warn scalar(@$tvariations);	
	next unless @$tvariations;
	$out .= $genes->{$gid}->{$trid}->{title};
	my $panel_id = "panel_$gid".$trid;
		
	$out .=  $cgi->start_div({class=>"panel-body panel-collapse collapse ",style=>"font-size: 09px;font-family:  Verdana;",id=>$panel_id});	
	#$out .=  $cgi->start_div({class=>"panel-body panel-collapse  collapse",style=>"font-size: 09px;font-family:  Verdana;"});
	#	$out.="\n";
	$out .= $cgi->start_table({class=>"table table-striped table-condensed table-bordered table-hover table-mybordered",style=>"vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;"});
	
	#header;
		$out.= $cgi->start_Tr();
	foreach my $c (@patients_columns,@variant_columns,@freq_columns,@transcripts_columns,@prediction_columns){
			$out.= qq{<th scope="col">$c</th>};
	}
	$out.= $cgi->end_Tr();
	#end header
	
	my $nbv =0;
	foreach my $vid (keys %{$transcripts_variations->{$trid}}){
			my $var = $variations->{$vid};
				my $tr_infos =$var->{functional}->{$gid};
			#warn $var->{functional}->{$gid}->{$trid}->{impact_score};
			#$variations->{scaled_score} =0;
			
		
			$var->{freq} = sprintf("%.3f %",100*$var->{freq});
	
		#next if $tr_infos->{$trid}->{consequence} eq "intronic";
		my $bg = "background-color:#F7F7F7";
		$bg = "background-color:#FFFFFF" if $nbv%2 ==0;
#		warn $tr_infos->{$trid}->{scaled_score}." ".$vid ;
			if ($tr_infos->{$trid}->{scaled_score} > 0){	
			$bg =$style_score->{$tr_infos->{$trid}->{scaled_score}}->{style};
				#$bg =$style_score->{4}->{style};
		}
		$nbv++;
		$out.= $cgi->start_Tr({style=>"vertical-align:middle"});				

			my $pinfos = $var->{patients};
			my $rowspan =  keys %$pinfos;
		
		my $text = $var->{deja_vu} ;
		my $url = $deja_vu_url.$vid;
		my $urlp = qq{http://www.polyweb.fr/cgi-bin/polymorphism-cgi//validation_variation/patient_report.pl?edit_mode=1&report_mode=1&project_summary=1&never=1&this=6&impact=3&frequence=4&allele_quality=5&project=};
		my $nb_row=0;
		my @pids = sort {$pinfos->{$a}->{project} cmp $pinfos->{$b}->{project}} keys %$pinfos;
		my $pid = shift @pids;
		my $debug;
		
		#foreach my $pid ( ){
			$nb_row ++;
			$pinfos->{$pid}->{patient} = $patients->{$pid}->{name};
			my $url2 =$urlp.$patients->{$pid}->{project}.qq{&patients=}. $patients->{$pid}->{name}."&transcripts=".$transcripts->{$trid}->{name};
			my $n = $patients->{$pid}->{name};
			my $t = $transcripts->{$trid}->{name};
			my $p = $patients->{$pid}->{project};
			#	my ($chr,$start) = split(":",$var->{genomique});
			#my $l = $chr.":".($start-20)."-".($start+20);#$var->{genomique};
			my $l = $var->{genomique};
			my $f = "/NGS/$p/HG19/align/bwa/$n.bam";
			my $u = "http://www.polyweb.fr//NGS/$p/HG19/align/bwa/$n.bam";
			$var->{all} = $var->{ref_allele}."/".$var->{allele};	
			my $v = $var->{all};
			$pinfos->{$pid}->{igv} =qq{<button dojoType="dijit.form.Button"   iconClass="igvIcon" onclick='launch_web_igv("$p","$n","$f","$l","$v")' style="color:black"></button>};
		#	$pinfos->{$pid}->{igv} =qq{<button onclick='displayOneBAMIGV("$u","$l")' style="color:black">igv</button>};
			$pinfos->{$pid}->{patient} =qq{<button onclick='load_polydiag("$p","$n","$t")' style="color:black">$n</button>};
			#$pinfos->{$pid}->{patient} =qq{<a href ="$url2" target="_blank" style="color:black;font-weight:bold">}.$patients->{$pid}->{name}."</a>";
			$pinfos->{$pid}->{ratio} =qq{<p style="color:black;font-weight:bold">}.$pinfos->{$pid}->{ratio}."</p>";
			
			$pinfos->{$pid}->{project} = qq{<a href ="http://www.polyweb.fr/polyweb/coverage.html?project=}.$patients->{$pid}->{project}.qq{" target="_blank" style="color:black;font-weight:bold">}.$patients->{$pid}->{project}."</a>";
			$pinfos->{$pid}->{sex} = "<b>".$patients->{$pid}->{sex_icon}."</b>";
			
			foreach my $t (@patients_columns){
					$out.= $cgi->td({style=>$bg} ,  $pinfos->{$pid}->{$t});
			}
				#$out.= $cgi->end_Tr() if $nb_row < $rowspan;
		#}
			
		

		
		$var->{deja_vu} = qq{<a href="$url" target="_blank" style="color:black;font-weight:bold">$text</a>};
		my $a0 = $var->{ref_allele};
		my $a1 = $var->{allele};
		my $chr = $var->{chromosome};
		my ($chr,$start) = split(":",$var->{genomique});
		my $qq5 = qq{	<button    class="alamutView3" onClick ="displayInAlamut('$chr',$start,['$a0','$a1']);"  style="color:black"></button>};
		$var->{alamut} = $qq5; 
		foreach my $t (@variant_columns){
			$out.= $cgi->td({rowspan=>$rowspan,style=>$bg} , $var->{$t});
		}
		foreach my $t (@freq_columns){
			$out.= $cgi->td({rowspan=>$rowspan,style=>$bg} , $var->{$t});
		}
		#next unless exists $var->{functional};
		
		$tr_infos->{$trid}->{transcript} = $transcripts->{$trid}->{name};
		
		foreach my $t (@transcripts_columns){
					$out.= $cgi->td({rowspan=>$rowspan,style=>$bg} ,  $tr_infos->{$trid}->{$t});
		}
		foreach my $t (@prediction_columns){
			$tr_infos->{$trid}->{$t} = $var->{$t} if $t eq "cadd"; 
			$out.= $cgi->td({rowspan=>$rowspan,style=>$bg} ,$tr_infos->{$trid}->{$t});
		}
		$out.= $cgi->end_Tr();
		foreach my $pid (@pids ){
			$nb_row ++;
			$pinfos->{$pid}->{patient} = $patients->{$pid}->{name};
			#my $url2 =$urlp.$patients->{$pid}->{project}.qq{&patients=}. $patients->{$pid}->{name}."&transcripts=".$transcripts->{$trid}->{name};
				my $n = $patients->{$pid}->{name};
			my $t = $transcripts->{$trid}->{name};
			my $p = $patients->{$pid}->{project};
		
			#	my ($chr,$start) = split(":",$var->{genomique});
			#my $l = $chr.":".($start-20)."-".($start+20);#$var->{genomique};
			my $l =$var->{genomique};
			my $f = "/NGS/$p/HG19/align/bwa/$n.bam";
			my $u = "http://www.polyweb.fr//NGS/$p/HG19/align/bwa/$n.bam";
			my $v = $var->{all};
			$pinfos->{$pid}->{igv} =qq{<button dojoType="dijit.form.Button"   iconClass="igvIcon" onclick='launch_web_igv("$p","$n","$f","$l","$v")' style="color:black"></button>};
			$pinfos->{$pid}->{patient} =qq{<button onclick='load_polydiag("$p","$n","$t")'  style="color:black">$n</button>};
		#	$pinfos->{$pid}->{patient} =qq{<a href ="$url2" target="_blank" style="color:black;font-weight:bold">}.$patients->{$pid}->{name}."</a>";
			$pinfos->{$pid}->{ratio} =qq{<p style="color:black;font-weight:bold">}.$pinfos->{$pid}->{ratio}."</p>";
			$pinfos->{$pid}->{project} = qq{<a href ="http://www.polyweb.fr/polyweb/coverage.html?project=}.$patients->{$pid}->{project}.qq{" target="_blank" style="color:black;font-weight:bold">}.$patients->{$pid}->{project}."</a>";
			$pinfos->{$pid}->{sex} = "<b>".$patients->{$pid}->{sex_icon}."</b>";
			
			foreach my $t (@patients_columns){
					$out.= $cgi->td({style=>$bg} ,  $pinfos->{$pid}->{$t});
			}
				$out.= $cgi->end_Tr();# if $nb_row < $rowspan;
		}
		
			

	
		#$out.= $cgi->end_Tr(); #end variant line;
	
		 #end variant line;
		
	}#end variants
	$out.= $cgi->end_table();
	$out.= $cgi->end_div();
	$out.= $cgi->end_div();
	$out.= $cgi->end_div();
	}#end transcripts	

$out.= $cgi->end_div();

}#end gene
print $out;
exit(0);


sub start_gene_line {

	 my ($gid,$trid,$transcript,$avariations) =  @_;
	 my $out ="";
	  my $bilan ={};	
	  my $array_impact;
	  my $b =0;
	  my $c =0;
		foreach my $v (@$avariations){
		#	warn $v->{freq_score};
			
			my $a = $v->{functional}->{$gid}->{$trid};
			
		#	 $a->{impact_score} = $himpact_sorted->{$a->{impact_text}};
			 
			$bilan->{impact}->{$a->{scaled_score}} ++;
			$v->{impact_score} =  $a->{scaled_score} ;
			
			$bilan->{max_impact} = $a->{impact_score} if $bilan->{max_impact} < $a->{impact_score};
			$bilan->{max_freq} = $v->{freq_score} if $bilan->{max_freq} < $v->{freq_score};
			#warn  10**$a->{scaled_score}." ".$a->{scaled_score};
			$bilan->{max_score} += 10**$a->{scaled_score};# if $bilan->{max_score} < $a->{scaled_score};
		}
		$genes->{$gid}->{score} =   $bilan->{max_score} if $genes->{gid}->{score} < $bilan->{max_score};
	 my $out ;
	 $out .=  $cgi->start_div({class=>"panel panel-success" });
	 #panel heading
	$out .=  $cgi->start_div({class=>"panel-heading panel-warning warning ",style=>" min-height:13px;max-height:13px;padding:1px;border:1px"});
	$out .=  $cgi->start_div({class=>" btn-group btn-xs "});
	my $panel_id = "panel_$gid".$trid;

	my $label_id = "label_".$trid;
	my $glyph = "";
	my $name = $transcript->{name};
	
		$out .= qq{<div class="btn  btn-info btn-xs " style="position:relative;bottom:1px;min-width:200px;" onClick='collapse("$panel_id","$label_id")'>  <span id= "$label_id" class="glyphicon glyphicon-triangle-right  " aria-hidden="true"  style="float:left;"></span> $name &nbsp $glyph</div>};
	  		my $nb_var = scalar(@$avariations);
	  		$out .=$cgi->span({class=>"label label-success"},qq{<span class='badge badge-primary badge-xs '  >$nb_var</span> });
				#$out .=$cgi->span({class=>"label label-success"},$nbc) if $nbc >0;
				$out .=$cgi->span({class=>"label label-danger badge-error"}, $bilan->{impact}->{4}+0);
				
				$out .=$cgi->span({class=>"label label-warning",style=>"background-color:#FF8800"}, $bilan->{impact}->{3}+0);
				$out .=$cgi->span({class=>"label label-warning",style=>"background-color:#EFB73E"}, $bilan->{impact}->{2}+0);
				my %style;
				$style{class} ="label label-danger "; 
				#style=>"background-color:#F2F183"}
					$style{style} ="background-color:#F2F183;color:#000000"; 
				# %style =%{$style_score->{1}};
				 	$out .=$cgi->span(\%style,$bilan->{impact}->{1});
				$out .=$cgi->span({class=>"label label-default"},$bilan->{impact}->{0}+0);# if $nbc >0;
			
				
	   		$out.= $cgi->end_div();

			#div lavel right 
#				$out .=  $cgi->start_div({class=>" btn-group btn  ",style=>'position:relative;float:right;bottom:5px;'});
#			
#					my $max_impact = 0;
#					
#					my $text = $himpact_sorted_inv->{$max_impact};
#					
#					my %style = %{$style_score->{1}};
#		
#					$style{class} ="label  "; 
#					$out .=$cgi->span(\%style,"S"); 
#					
#					 %style =%{$style_score->{1}};
#					$style{class} ="label  "; 
#					$out .=$cgi->span(\%style,"I");
#					my $max_freq = 1;
#					$himpact_sorted_inv->{$max_freq};
#					%style = %{$style_score->{$max_freq}};
#					$style{class} ="label  "; 
#					$out .=$cgi->span(\%style,"F");
#					
#						
#				
#			
#					$out .=$cgi->span({class=>"label label-success"},qq{<span class="glyphicon glyphicon glyphicon-menu-hamburger " aria-hidden="true" "></span>});
#					my $clabel = " label-default";
#					$out .=$cgi->span({class=>"label $clabel"},qq{<span class="glyphicon glyphicon-pencil " aria-hidden="true" "></span>});
#					$clabel = " label-default";
#					$out .=$cgi->span({class=>"label $clabel"},qq{<span class="glyphicon glyphicon-saved" aria-hidden="true" "></span>});
#			 	$out.= $cgi->end_div(); # end div lavel right 
			
			 	
			$out.= $cgi->end_div(); # end panel heading
	   		
	 return $out;
	 	
}




