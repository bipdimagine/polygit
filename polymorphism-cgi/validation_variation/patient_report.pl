#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
#use lib "/software/polyweb/poly-src/GenBo/lib/obj-nodb/";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/packages_old/cache/polydiag";
use List::MoreUtils qw{ natatime };
#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";
use lib "$Bin/../packages/validation_variation_old"; 
use lib "$Bin/../packages/cache"; 

use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag";
use draw_cnv; 
use update;
require "$Bin/../GenBo/lib/obj-nodb/packages_old/cache/polydiag/utility.pm";
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
use Storable qw/thaw/;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use coverage;
use Spreadsheet::WriteExcel;
use POSIX;
use validationQuery;
use Date::Tiny;
use List::MoreUtils qw{part};
#use PDF::API2;
#use PDF::Table;
use constant mm => 25.4 / 72;
use constant in => 1 / 72;
use constant pt => 1;
use Time::HiRes qw ( time alarm sleep );

$update::fts = "8px";

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


</style>

};

my $style_score ={
		1=> {style=>"background-color:#95A5A6;color:#FFFFFF"},#E43725#95A5A6
		2 => {style=>"background-color:#EFB73E;"},#F2DEDE#F7CCC8
		3=> {style=>"background-color:#FF8800"},#E3A57E#F1CE9C#F19B92
		4=> {style=>"background-color:#CE0000;color:white!important"},#E43725#FF4136
		99=> {style=>"background-color:#2F2F2F;color:white!important"}
		#99=> {style=>"background-color:#DFDFDF;color:#000000; border: 3px solid black;"}
};

my $style_td ={
		1=> {class=>"alert",style=>"background-color:#FDFDFD;"},#E43725#95A5A6
		2 => {class=>"warning2",style=>"background-color:#F7CCC8"},#F2DEDE#F7CCC8
		3=> {class=>"warning2",style=>"background-color:#F19B92"},#E3A57E#F1CE9C#F19B92
		4=> {class=>"danger",style=>"background-color:#FF4136"},#E43725#FF4136
};

my $server =  $ENV{HTTP_HOST};
my $variation_script = $ENV{SCRIPT_NAME};
$variation_script =~s/patient_/variation_/;


$server = "darwin.bipd.fr" if $server eq "bipd";
$server = "www.polyweb.fr" if $server =~/10\.200\.27/;
my $deja_vu_url = "http://$server//polyweb/polydejavu/dejavu.html?input=";
#my $lcdb_url = "http://$server//polyweb/polydejavu/dejavu.html?input=";
my $deja_vu_light_url = "http://$server/$variation_script";
my $lcdb_url = $deja_vu_light_url;
$lcdb_url =~s/variation_/lcdb_/;
#"http://$server/$variation_script";
my $nb_gene_by_patient = 3;
my $nb_exon_by_genes = 10;

my $buffer = GBuffer->new();

my $cgi          = new CGI();
my $project_name = $cgi->param('project');
my $force_cache =  $cgi->param('force_cache');
#my $project = $buffer->newProject(-name=>$project_name);
my $print =  $cgi->param('print');
$update::print = $cgi->param('print');
if($print){
	$update::fts = "7px";
}
my $user = $cgi->param('user_name');
my $hgmd = $buffer->queryHgmd()->getHGMD($user);

my $project;
#if ($cgi->param('cnv_coverage') ne 1){
	$project = $buffer->newProjectCache( -name 			=> $project_name, -typeFilters=>'individual' ,-cgi_object=>1);
#}
#else {
#		$project = $buffer->newProject( -name 			=> $project_name, -typeFilters=>'individual' ,-cgi_object=>1);

#}
my @headers = ("gene","var_name","ngs","trio","genomique","transcript","exon","nomenclature","consequence","codons","codons_AA","freq","freq_ho","max_pop","min_pop","hgmd","clinvar","local","deja_vu","similar_projects","in_this_run", "polyphen","sift","cadd");

if ($project->isSomatic){
	@headers = ("gene","var_name","sanger","ngs","ratio","caller","genomique","transcript","exon","nomenclature","consequence","codons","codons_AA","cosmic","clinvar","local","freq","freq_ho","max_pop","min_pop","deja_vu","similar_projects","in_this_run", "polyphen","sift","cadd");
}
if ($print == 1){
	@headers = ("var_name","ngs","transcript","exon","nomenclature","consequence","codons","codons_AA","freq","freq_ho","max_pop","min_pop","hgmd","clinvar","local","deja_vu","similar_projects","in_this_run", "polyphen","sift","cadd");
	@headers = ("var_name","ngs","transcript","exon","nomenclature","consequence","codons","codons_AA","freq","freq_ho","max_pop","min_pop","clinvar","local","deja_vu","similar_projects","in_this_run", "polyphen","sift","cadd") unless $hgmd ==1;
	@headers = ("var_name","ngs","genomique","transcript","exon","nomenclature","consequence","codons","codons_AA","freq","freq_ho","max_pop","min_pop","hgmd","clinvar","local","deja_vu","similar_projects","in_this_run", "polyphen","sift","cadd");
	
	$style_score ={
		1=> {style=>"background-color:#95A5A6;color:#FFFFFF"},#E43725#95A5A6
		2 => {style=>"background-color:#EFB73E;"},#F2DEDE#F7CCC8
		3=> {style=>"background-color:#FF8800"},#E3A57E#F1CE9C#F19B92
		4=> {style=>"background-color:#CE0000;color:white!important"},#E43725#FF4136
		99=> {style=>"background-color:#8F8F8F;color:white!important;"}
	#	99=> {style=>"background-color:#DFDFDF;color:#000000; border: 3px solid black;"}
};
}


#my $capture = $project->getCapture();
my $similar = $project->similarProjects();

my $vquery = $project->validations_query(1); #validationQuery->new(dbh=>$buffer->dbh,capture_name=>$project->validation_db());
my $list_transmission= {
denovo=> 1,
strict_denovo =>1,
'denovo/?'=>1,
recessive =>1,
both=>1,
'uniparental disomy'=>1,
'mosaic'=>1,
'mosaic/?'=>1,
xor =>1,
};

my $filter_transmission;
$filter_transmission->{denovo} = 1 if $cgi->param('denovo');
$filter_transmission->{"strict_denovo"} = 1 if $cgi->param('denovo');
$filter_transmission->{'denovo/?'} = 1 if $cgi->param('denovo');
$filter_transmission->{recessive} = 1 if $cgi->param('recessive');
$filter_transmission->{both} = 1 if $cgi->param('both');
$filter_transmission->{xor} = 1 if $cgi->param('xor');
$filter_transmission->{'uniparental disomy'} = 1 if $filter_transmission;
$filter_transmission->{'mosaic'} = 1 if $filter_transmission;
$filter_transmission->{'mosaic/?'} = 1 if $filter_transmission;
 
my $patient_name = $cgi->param('patients');

my $edit_mode = $cgi->param('edit_mode');
my $compute_coverage = 1;
my $all;

$all = $cgi->param('all');;
if ($edit_mode){
	$compute_coverage = undef;
}
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
#define parameter 
my $padding =  $cgi->param('span');
my $cov_limit =  $cgi->param('limit');
my $vimpact = $cgi->param('impact');
my $this_run = $cgi->param('this');

$this_run=6 unless $this_run;

my $never =  $cgi->param('never');


die() unless $vimpact;

my $impact_score_limit = $himpact_sorted->{$himpact->{$vimpact}};

my $vfreq = $cgi->param('frequence');
my $filter_quality = $cgi->param('allele_quality');
my $limit_ratio = -1;
$limit_ratio = "80" if $filter_quality == 1;
$limit_ratio = "40" if $filter_quality == 2;
$limit_ratio = "20" if $filter_quality == 3;
$limit_ratio = "10" if $filter_quality == 4;
if ($filter_quality < 0){
	$limit_ratio = $filter_quality * -1;
}
my $mode_report = $cgi->param('report_mode');

my $cgi_transcript =  $cgi->param('transcripts');
my $name =  $cgi->param('name');

$patient_name ="all" unless $patient_name;
my $patients = $project->get_list_patients($patient_name,",");
my $nb_patients = scalar @{$patients->[0]->getRun->getAllPatientsInfos()};

#define limit for in this run

my $hscore_this_run = {
	"1" => 1 ,
	"2" =>int(0.1 *$nb_patients ),
	"3" => int(0.25 *$nb_patients ),
	"4" => int(0.5 *$nb_patients ),
	"5" => int(0.75 *$nb_patients ),
	"6" =>$nb_patients +1,
};

my $hthisrun ={
	1=>"uniq",
	2=>10,
	3=>25,
	4=>50,
	5=>75,
	6=>"all"
};


my $patients = $project->get_list_patients($patient_name,",");

my $freq = 9999;
my @transcripts_cgi ;
if ($cgi_transcript eq "all"){
	@transcripts_cgi = @{$project->bundle_transcripts() } ;
}
else {
	@transcripts_cgi = split(",",$cgi_transcript);
}


my $table_id = "hor-minimalist-b";




$| =1;
my $out_global ="";
my $data;
if ($edit_mode){
	my $xls =  $cgi->param('xls');
	if ($xls){
		$data = construct_data();
		$out_global .= edit_mode($data,$cgi);
		html::print_cgi($cgi,$CSS.$out_global,$print,$patient_name." - PolyDiag");
		exit(0);
	}
	my $CSS = "";

	html::print_cgi_header($cgi,$CSS.$out_global,$print,$patient_name." - PolyDiag");
	
	print qq{<div style="visibility: hidden">};
 	$data = construct_data(1);
 
 	print qq{</div>};
	$out_global .= print_hotspot($data,$cgi);
	$out_global .= edit_mode($data,$cgi);
	print $out_global;
	exit(0);
	
}
else {
	 	html::print_cgi_header($cgi,$CSS.$out_global,$print,$patient_name." - PolyDiag");
	 print qq{<div style="visibility: hidden">};
	
	 $data = construct_data(1);

	 print qq{</div>};
	if ($mode_report == 1){
		$out_global = alacarte($data,$cgi);
		exit(0);
	}
	else {
		die();
		$vimpact = 'all';
		$out_global = html($data,$cgi);
	}
}
my $CSS = "";
	#print $cgi -> header;
	#print $CSS;
	#print $out;
	#exit(0);
html::print_cgi($cgi,$CSS.$out_global,$print,$patient_name." - PolyDiag");

exit(0);


my $dude;
my %hdude;

sub construct_htranscripts {
	my ($list_transcripts,$patient) = @_;
	my $hpatient;
	my @res;

	foreach my $tr (@$list_transcripts) {
		my $utr = $cgi->param('utr')+0;
		print "+";
	my $tr_id = $tr;
			my $tr1;
				my $htranscript;

				if ($tr ne "intergenic"){
				 $tr1 = $project->newTranscript($tr);
				 $utr =1 if $tr1->getChromosome()->name eq "MT";
				$tr_id = $tr1->id;
				$htranscript->{name} = $tr1->getGene->external_name();
				$htranscript->{mean} = $tr1->mean_coding_coverage($patient);
				$htranscript->{mean} = $tr1->mean_exonic_coverage($patient) if $htranscript->{mean} == 0;
				$htranscript->{obj} = $tr1;
				$htranscript->{exons} = [];
				$htranscript->{variations} = [];
				$htranscript->{all} = [];
				$htranscript->{table} = 1;
				$htranscript->{id} =$tr1->id;
				$htranscript->{external_name} = $tr1->external_name;
				#$htranscript->{ret}  = image_coverage::image_cnv ([$project], $tr1 );
			}
			else {
				$htranscript->{name} = "intergenic";
				$htranscript->{mean} = "-";
				$htranscript->{mean} = "-";
				$htranscript->{obj} = undef;
				$htranscript->{exons} = [];
				$htranscript->{variations} = [];
				$htranscript->{all} = [];
				$htranscript->{table} = 1;
				$htranscript->{external_name} = "intergenic";
			}
			my $kvars = utility::return_list_variants($project,$patient,$tr_id);
			
			if ($tr ne "intergenic"){
		
		
		if ($compute_coverage) {
				my $exons;
			if ($cgi->param('intronic') == 1){
	 					$exons = $tr1->getAllGenomicsParts();
			}
			else {
		 			$exons = $tr1->getExons();
			}
			
			my $exons_todo = $vquery->get_exons(project_name=>$project_name,sample_name=>$patient->{name});		
			my $show_utr = $utr +0;
			my $intronic =  $cgi->param('intronic') +0;
			
			my $kyoto_id = join("_",("all",$tr,$show_utr,$intronic,$cov_limit,$padding,"data"));

			my $cdata;

			unless ($cdata){
		#	unless ($cdata && $cdata->{color2} ) {

					my $ret = image_coverage::image ([$patient], $tr1,$intronic,$show_utr, $padding, $cov_limit,1);
					$cdata = $ret->{data};
			}
			
			foreach my $exon (sort{$a->start*$a->strand <=> $b->end*$b->strand }@$exons){
				#my ($exon,$patient,$tr1,$exons_todo,$data,$limit) = @_;
				 my $pdata = $cdata->{$exon->id}->{$patient->id};
				my $hexons = infos_coverage_exons::return_hash_exons2 ($exon,$patient,$tr1,$exons_todo,$pdata,$cov_limit);
				my $show_utr = $cgi->param('utr') +0;
				
				#my $hexons = infos_coverage_exons::return_hash_exons ($exon,$patient,$exons_todo,$capture_intspan,$tr1,$show_utr,$cov_limit,$padding);

				push(@{$htranscript->{exons}},$hexons);
			}
			} #end_compute_coverage;
		}
	#	foreach my $var (@{$htr_vars->{$tr1->kyotoId}}){
		foreach my $var (@{$kvars}){   
				my $debug;
				$debug =1 if $var eq "5_11385203_C_CCGG"; 
				my $hvariation = utility::return_hash_variant($project,$var,$tr_id,$patient,$vquery);
				if ($print ==1){
					$hvariation->{min_pop} =~ s/<[^>]*>//gs;
					$hvariation->{max_pop} =~ s/<[^>]*>//gs;
					$hvariation->{cadd} =~ s/<[^>]*>//gs;
				}
					update::edit($patient,$hvariation); 
				
					update::clinvar($project,$hvariation); 
					update::hgmd($project,$hvariation); 
					update::tclinical_local($project,$hvariation,$patient,$htranscript->{obj}->getGene);
					update::deja_vu($project,$tr1,$hvariation,$debug);
		
					my $zfilter = 1;
					$zfilter = undef if $hvariation->{clinvar_alert} ;	
					$zfilter = undef if $hvariation->{clinical_local} ;
					#$zfilter = undef if $hvariation->{hgmd} ;
					$zfilter = undef if  $hvariation->{type} ne "other";	
					
					my $debug;
					$debug = 1 if $hvariation->{var_name} eq "17_1940466_T_G";
				
					if ($zfilter){
						next  if $hvariation->{this_deja_vu} > $hscore_this_run->{$this_run};
						#next if $hvariation->{impact_score} < $impact_score_limit;
						next if $hvariation->{freq_level} > $vfreq ;
					}
				$hvariation->{ratio} = $hvariation->{ratio};
				my @all_nums    = $hvariation->{ratio} =~ /([+-]?[0-9]*[.]?[0-9]+)%/g;
				#$limit_ratio = 0.2;
					if (scalar(@all_nums)){
						my @t = grep{$_>=$limit_ratio} @all_nums;
						next unless @t;
					}
				update::deja_vu($project,$tr1,$hvariation,$debug);
				if ($zfilter){
						next  if $hvariation->{this_deja_vu} > $hscore_this_run->{$this_run};
					}
				
				update::trio($project,$tr,$hvariation,$patient,$cgi,$print);
				if (exists $hvariation->{transmission_model} and  $hvariation->{transmission_model}=~ /strict_denovo/){
						 $hvariation->{transmission_model} = "strict_denovo";
				}
				
 				if ( exists $hvariation->{transmission_model}){
					my $t = $hvariation->{transmission_model};
#					#$t ='strict_denovo' unless exists $list_transmission->{$t};
					confess($t.' '.$hvariation->{id}) unless exists $list_transmission->{$t};
					if ($filter_transmission){
						
						next unless exists  $filter_transmission->{$t};
					}
				}

	
				update::annotations($project,$hvariation);
	
			
				$hvariation->{freq} = sprintf("%.5f",$hvariation->{freq}) if $hvariation->{freq} ne "-";
				
				unless  (exists $hvariation->{edit}) {
					
					if ($zfilter eq 1) {
					
					next if $hvariation->{this_deja_vu} > $hscore_this_run->{$this_run};
					next if $hvariation->{impact_score} < $impact_score_limit;
					next if $hvariation->{freq_level} > $vfreq ;
					}
					#my @vration = split("<BR>",$hvariation->{ratio});
					my @all_nums    = $hvariation->{ratio} =~ /(\d+)/g;
					if (scalar(@all_nums)){
						my @t = grep{$_>$limit_ratio} @all_nums;
						next unless @t;
					}
					
				}
				if ($force_cache==1){
				update::deja_vu($project,$tr1,$hvariation);
				
				#$db_lite->put($patient_name,$id,$h);
				}
				#unless  (exists $hvariation->{edit}){
					next if $hvariation->{freq_level} > $vfreq && $zfilter;
				#}
				if ($hvariation->{diff_project_deja_vu} == 0 && $hvariation->{freq} <= 0.01){
					$hvariation->{freq_score} = 4; 
					
				}
				elsif ($hvariation->{diff_project_deja_vu} <= 5 || $hvariation->{freq} <= 0.01){
							$hvariation->{freq_score} = 3; 
					
				}
				elsif ($hvariation->{diff_project_deja_vu} <= 20 || $hvariation->{freq} <= 0.05){
						$hvariation->{freq_score} = 2; 
					
				}
				else {
						$hvariation->{freq_score} = 1; 
				}
				
				
				
				my $href = qq{<a  href = "$deja_vu_light_url};
					 if (exists $hvariation->{dup}){
					 		$href .= qq{?project=$project_name&transcript=$tr_id&variation_id=$var" target="_blank" style="color:white;font-weight:bold">};
					 }
					 else {
					$href .= qq{?project=$project_name&transcript=$tr_id&variation_id=$var" target="_blank" style="color:black;font-weight:bold">};
					 }
				if ($edit_mode == 1){
					my $vv = $href.$hvariation->{in_this_run}."</a>" ;
					$hvariation->{in_this_run} =$vv;
				}
				
				#push(@{$htranscript->{"all_variations"}},$hvariation);
				#push(@{$hpatient->{variations}->{$hvariation->{id}}},$hvariation );
				#warn Dumper $hvariation->{obj};
				delete $hvariation->{obj};
				push(@{$htranscript->{all}},$hvariation);
				
			}
				delete $htranscript->{obj};
			push(@res,$htranscript);
			#push(@{$hpatient->{transcripts_not_sorted}},$htranscript);
	}
	return \@res;
}

sub construct_data {
	my ($print_dot) = @_; 
	my $t =time;
	my $sum_time = 0;
	my $htr_vars;
	my $cpt =0;
	my $patient = $patients->[0];

	my $hpatient;
	$hpatient->{name} = $patient->name();
	$hpatient->{obj} = $patient;
	#$hpatient->{obj} = $patient;
	
	$hpatient->{machine} = $patient->getRuns->[0]->machine();
	$hpatient->{capture_type} = $patient->getCapture()->type;
	$hpatient->{capture_description} = $patient->getCapture()->description;
	$hpatient->{capture_version}  = $patient->getCapture()->version;

	
	my $nbv = scalar(@transcripts_cgi);
	
	
	my $list_transcript = \@transcripts_cgi;

	if ($all){
		$list_transcript = utility::return_list_all_transcripts($project,$patient);
	
	}
	push(@$list_transcript,"intergenic")  if $cgi->param('all') == 1;
	
	my $fork      = 8;
	my $nb        = int( scalar(@$list_transcript) / $fork + 1 );
	my $pm        = new Parallel::ForkManager($fork);
	my $iter      = natatime( $nb, @$list_transcript );
	my $hrun;
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $h ) = @_;

			unless ( defined($h) or $exit_code > 0 ) {
				print
				  qq|No message received from child process $exit_code $pid!\n|;
				die();
				return;
			}
			#	delete $jobs->{$h->{jobs});
			foreach my $htranscript (@{$h->{htranscripts}}){
				my $hvariations = $htranscript->{all};
					foreach my $hvariation (@$hvariations){
						push(@{$hpatient->{variations}->{$hvariation->{id}}},$hvariation );
						push(@{$htranscript->{$hvariation->{type}}},$hvariation);
					}
					push(@{$hpatient->{transcripts_not_sorted}},$htranscript);
				}
			
			my $id = $h->{run_id};
			delete $hrun->{ $h->{run_id} };
		}
	);
	$project->buffer->dbh_deconnect();

	my $iv =0;
	my $id = time;
	while ( my @tmp = $iter->() ) {
		$id ++;
		$hrun->{$id} ++;
		my $pid = $pm->start and next;
		my $t   = time;
		my $res;
		 $res->{htranscripts} = construct_htranscripts( \@tmp, $patient );
		 # $res->{htranscripts} = [];
		$res->{run_id} = $id;
		$res->{ttime}  = time;
		$pm->finish( 0, $res );
		$cpt++;
			print ". " if $print_dot && $cpt%5 ==0;
			
	}
	$pm->wait_all_children();
	
	confess() if keys %$hrun;
	
	@{$hpatient->{transcripts}}  =();
	@{$hpatient->{transcripts}} = sort {$a->{name} cmp $b->{name}}  @{$hpatient->{transcripts_not_sorted}} if $hpatient->{transcripts_not_sorted};
	delete 	$hpatient->{transcripts_not_sorted};
	push(@$data,$hpatient);
	return $data;
	
	exit(0);
	
}





sub construct_data1 {
	my ($print_dot) = @_; 
	my $t =time;
	my $sum_time = 0;
	my $htr_vars;
	my $cpt =0;
	my $patient = $patients->[0];

#	 ($dude) = grep {$_ eq "dude" } @{$patient->getCallingMethods};
#	if ($dude){
#		my $file = $patient->getVariationsFile("$dude");
#		open (CNV,"zcat $file | ");
#		while(<CNV>){
#			chomp();
#			my ($chr,$s,$e,$st,$type,$sc1,$id,@all) = split(" ",$_);
#			$hdude->{$chr}->{$id}->{nb} ++;
#			$hdude->{$chr}->{$id}->{status} = $s;
#			$hdude->{$chr}->{$id}->{type} = $type;
#		}
#		warn "coucou $file";
#		 my $c = $patient->getCnvs();
#		warn "end ".scalar(@$c);
#	}
#	die();
	my $hpatient;
	$hpatient->{name} = $patient->name();
	$hpatient->{obj} = $patient;
	#$hpatient->{obj} = $patient;
	
	$hpatient->{machine} = $patient->getRuns->[0]->machine();
	$hpatient->{capture_type} = $patient->getCapture()->type;
	$hpatient->{capture_description} = $patient->getCapture()->description;
	$hpatient->{capture_version}  = $patient->getCapture()->version;

	
	my $nbv = scalar(@transcripts_cgi);
	
	
	my $list_transcript = \@transcripts_cgi;

	if ($all){
		$list_transcript = utility::return_list_all_transcripts($project,$patient);
	
	}
	push(@$list_transcript,"intergenic")  if $cgi->param('all') == 1;

	my $iv =0;
	foreach my $tr (@$list_transcript) {
		$cpt++;
			print ". " if $print_dot && $cpt%5 ==0;
			my $tr_id = $tr;
			my $tr1;
				my $htranscript;

				if ($tr ne "intergenic"){
				 $tr1 = $project->newTranscript($tr);
				$tr_id = $tr1->id;
				$htranscript->{name} = $tr1->getGene->external_name();
				$htranscript->{mean} = $tr1->mean_coding_coverage($patient);
				$htranscript->{mean} = $tr1->mean_exonic_coverage($patient) if $htranscript->{mean} == 0;
				$htranscript->{obj} = $tr1;
				$htranscript->{exons} = [];
				$htranscript->{variations} = [];
				$htranscript->{all} = [];
				$htranscript->{table} = 1;
				$htranscript->{external_name} = $tr1->external_name;
			}
			else {
				$htranscript->{name} = "intergenic";
				$htranscript->{mean} = "-";
				$htranscript->{mean} = "-";
				$htranscript->{obj} = undef;
				$htranscript->{exons} = [];
				$htranscript->{variations} = [];
				$htranscript->{all} = [];
				$htranscript->{table} = 1;
				$htranscript->{external_name} = "intergenic";
			}
			my $kvars = utility::return_list_variants($project,$patient,$tr_id);
			
			if ($tr ne "intergenic"){
		
		
		if ($compute_coverage) {
				my $exons;
			if ($cgi->param('intronic') == 1){
	 					$exons = $tr1->getAllGenomicsParts();
			}
			else {
		 			$exons = $tr1->getExons();
			}
			
			my $exons_todo = $vquery->get_exons(project_name=>$project_name,sample_name=>$patient->{name});		
			my $show_utr = $cgi->param('utr') +0;
			my $intronic =  $cgi->param('intronic') +0;
			
			my $kyoto_id = join("_",("all",$tr,$show_utr,$intronic,$cov_limit,$padding,"data"));

			my $cdata;

			unless ($cdata){
		#	unless ($cdata && $cdata->{color2} ) {

					my $ret = image_coverage::image ([$patient], $tr1,$intronic,$show_utr, $padding, $cov_limit,1);
					$cdata = $ret->{data};
			}
			
			foreach my $exon (sort{$a->start*$a->strand <=> $b->end*$b->strand }@$exons){
				#my ($exon,$patient,$tr1,$exons_todo,$data,$limit) = @_;
				 my $pdata = $cdata->{$exon->id}->{$patient->id};
				my $hexons = infos_coverage_exons::return_hash_exons2 ($exon,$patient,$tr1,$exons_todo,$pdata,$cov_limit);
				my $show_utr = $cgi->param('utr') +0;
				
				#my $hexons = infos_coverage_exons::return_hash_exons ($exon,$patient,$exons_todo,$capture_intspan,$tr1,$show_utr,$cov_limit,$padding);

				push(@{$htranscript->{exons}},$hexons);
			}
		} #end_compute_coverage;
		}
	#	foreach my $var (@{$htr_vars->{$tr1->kyotoId}}){
	
		foreach my $var (@{$kvars}){   
				my $debug;
				$debug =1 if $var eq "5_11385203_C_CCGG"; 
				my $hvariation = utility::return_hash_variant($project,$var,$tr_id,$patient,$vquery);
				if ($print ==1){
					$hvariation->{min_pop} =~ s/<[^>]*>//gs;
					$hvariation->{max_pop} =~ s/<[^>]*>//gs;
					$hvariation->{cadd} =~ s/<[^>]*>//gs;
				}
					update::edit($patient,$hvariation); 
				
					update::clinvar($project,$hvariation); 
					update::hgmd($project,$hvariation); 
					update::tclinical_local($project,$hvariation,$patient,$htranscript->{obj}->getGene);
						
				#$debug = 1 if $hvariation->{genomique} eq "1:62740264";
				#next unless $debug;
					my $zfilter = 1;
				#	$zfilter = undef if $hvariation->{consequence} =~/essential/;
				#	$zfilter = undef if $hvariation->{consequence} =~/phase/;
				#	$zfilter = undef if $hvariation->{consequence} =~/stop/;	
				#	 $zfilter = undef  if $hvariation->{consequence} !~/non/ && $hvariation->{consequence} =~/frameshift/ ;
					$zfilter = undef if $hvariation->{clinvar_alert} ;	
					$zfilter = undef if $hvariation->{clinical_local} ;
					#$zfilter = undef if $hvariation->{hgmd} ;
					$zfilter = undef if  $hvariation->{type} ne "other";	
					
					my $debug;
					$debug = 1 if $hvariation->{var_name} eq "17_1940466_T_G";
				
					
					if ($zfilter){
						next  if $hvariation->{this_deja_vu} > $hscore_this_run->{$this_run};
						#next if $hvariation->{impact_score} < $impact_score_limit;
						next if $hvariation->{freq_level} > $vfreq ;
					}
				$hvariation->{ratio} = $hvariation->{ratio};
				my @all_nums    = $hvariation->{ratio} =~ /([+-]?[0-9]*[.]?[0-9]+)%/g;
				#$limit_ratio = 0.2;
					if (scalar(@all_nums)){
						my @t = grep{$_>=$limit_ratio} @all_nums;
						next unless @t;
					}
				
				update::trio($project,$tr,$hvariation,$patient,$cgi,$print);
				if (exists $hvariation->{transmission_model} and  $hvariation->{transmission_model}=~ /strict_denovo/){
						 $hvariation->{transmission_model} = "strict_denovo";
				}
 				if ( exists $hvariation->{transmission_model}){
					my $t = $hvariation->{transmission_model};
#					#$t ='strict_denovo' unless exists $list_transmission->{$t};
					confess($t.' '.$hvariation->{id}) unless exists $list_transmission->{$t};
					if ($filter_transmission){
						
						next unless exists  $filter_transmission->{$t};
					}
				}

				#confess("problem with transmission name" ) if ( exists $hvariation->{transmission_model} and ! exists $list_transmission->{$hvariation->{transmission_model}});
				#next if ($filter_transmission && exists $hvariation->{transmission_model} && ! exists $filter_transmission->{$hvariation->{transmission_model}});
	
				update::annotations($project,$hvariation);
	
				#my $t1 =time;
				update::deja_vu($project,$tr1,$hvariation,$debug);
				#$sum_time += abs(time-$t1);
				$hvariation->{freq} = sprintf("%.5f",$hvariation->{freq}) if $hvariation->{freq} ne "-";
				
				unless  (exists $hvariation->{edit}) {
					
					if ($zfilter eq 1) {
					
					next if $hvariation->{this_deja_vu} > $hscore_this_run->{$this_run};
					next if $hvariation->{impact_score} < $impact_score_limit;
					next if $hvariation->{freq_level} > $vfreq ;
					}
					#my @vration = split("<BR>",$hvariation->{ratio});
					my @all_nums    = $hvariation->{ratio} =~ /(\d+)/g;
					if (scalar(@all_nums)){
						my @t = grep{$_>$limit_ratio} @all_nums;
						next unless @t;
					}
					
				}
				if ($force_cache==1){
				update::deja_vu($project,$tr1,$hvariation);
				
				#$db_lite->put($patient_name,$id,$h);
				}
				#unless  (exists $hvariation->{edit}){
					next if $hvariation->{freq_level} > $vfreq && $zfilter;
				#}
				if ($hvariation->{diff_project_deja_vu} == 0 && $hvariation->{freq} <= 0.01){
					$hvariation->{freq_score} = 4; 
					
				}
				elsif ($hvariation->{diff_project_deja_vu} <= 5 || $hvariation->{freq} <= 0.01){
							$hvariation->{freq_score} = 3; 
					
				}
				elsif ($hvariation->{diff_project_deja_vu} <= 20 || $hvariation->{freq} <= 0.05){
						$hvariation->{freq_score} = 2; 
					
				}
				else {
						$hvariation->{freq_score} = 1; 
				}
				
				
				
				my $href = qq{<a  href = "$deja_vu_light_url};
					 if (exists $hvariation->{dup}){
					 		$href .= qq{?project=$project_name&transcript=$tr_id&variation_id=$var" target="_blank" style="color:white;font-weight:bold">};
					 }
					 else {
					$href .= qq{?project=$project_name&transcript=$tr_id&variation_id=$var" target="_blank" style="color:black;font-weight:bold">};
					 }
				if ($edit_mode == 1){
					my $vv = $href.$hvariation->{in_this_run}."</a>" ;
					$hvariation->{in_this_run} =$vv;
				}
				
				push(@{$htranscript->{$hvariation->{type}}},$hvariation);
				push(@{$hpatient->{variations}->{$hvariation->{id}}},$hvariation );
				
				push(@{$htranscript->{all}},$hvariation);
				
			}
		
			push(@{$hpatient->{transcripts_not_sorted}},$htranscript);
	}
	
	
	@{$hpatient->{transcripts}}  =();
	@{$hpatient->{transcripts}} = sort {$a->{name} cmp $b->{name}}  @{$hpatient->{transcripts_not_sorted}} if $hpatient->{transcripts_not_sorted};
	delete 	$hpatient->{transcripts_not_sorted};
	push(@$data,$hpatient);
	return $data;
	
	exit(0);
	
}







####
#
####

sub alacarte{
	my ($data,$cgi) = @_;
	$| =1;
	my $out ="";
	my $nb_case;
	my $types = {
							he=>2,
							ho=>3,
							td=>-3,
							nc=>-5
							
						};
my $hvalidation = {
							2=>"he",
							3=>"ho",
							-3=>"td",
							-5=>"nc",
							
						};
	#$out.= $CSS;
if ($print){
$out.="<body onload='this.print();'>";
}
else {
	$out.="<body '>";
}


$out.= $cgi->start_div({id=>"all_html",name=>"all_html",jsId=>"all_html",class=>"patient"});
print $out;
$out ="";

#print $out;

		foreach my $patient (@$data) {
				print  print_project_summary($patient)  if $cgi->param('project_summary') == 1;;
				print  print_graph_coverage($patient)  if $cgi->param('graphic_coverage') == 1;
				
				print printDetailCoverage($patient)   if $cgi->param('all_coverage') == 1;
				
				$out .= print_hotspot($data,$cgi,1);
				
				$out.=print_table_variants($patient) if $cgi->param('table_variations') == 1;
				$out.= printSanger($patient) if $cgi->param('sanger_variations') == 1;	
				$out.= printRejected($patient) if $cgi->param('sanger_variations') == 1;		
				$out.= printValidated($patient) if $cgi->param('validated_variations') == 1;	
				$out.= printTodo($patient) if $cgi->param('sanger_variations') == 1;			
				$out.= printOther($patient) if $cgi->param('all_variations') == 1;	
				print $out;
				$out ="";
				$out.=printCnv($patient)  if $cgi->param('cnv_coverage') == 1;
			
		}
print $out;
	return $out;
}


sub edit_mode {
	
	my ($data,$cgi) = @_;

my $out ="" ;# = html::print_cadre($cgi,"Edition ");
my $xls =  $cgi->param('xls');
	foreach my $patient (@$data){
		# $out .= html::print_cadre($cgi,"Genes ");
	#	die();
	
		unless ($xls){
		$out.=printTableGenes($patient,"All","all",2) ;
		}
		else {
			printTableGenesXls($patient,"All","all",2) ;
		}

	}
	return $out;
}

sub print_hotspot {
	
	my ($data,$cgi,$print) = @_;
	my $out1 .=printTableHotspots( $patients->[0],$print);
		return $out1;
my $out1 ="" ;# = html::print_cadre($cgi,"Edition ");
my $hotspots = $patients->[0]->getCapture()->hotspots;
return $out1 unless scalar(@$hotspots);

foreach my $hotspot (@$hotspots){
$patients->[0]->hotspot($hotspot);
}
$out1 .=printTableHotspots( $patients->[0],$hotspots,$print);


}



sub printCnv {
	my ($patient) = @_;
	my $project = $patient->{obj}->getProject();
	 my $out_temp = html::print_cadre($cgi," Dup/Del");
	 $| = 1;
	
	 	$out_temp.= $cgi->start_table({  class=>"table table-condensed table-striped table-bordered ",style=>"font-size:9px;width:1200px "});
	 #{class=>"table table-condensed table-striped table-bordered ",style=>"font-size:5px;font-family:  Verdana;"}
	 #print qq{<div style="visibility: hidden">};
	
	 my $nb_col = 30;
	 my $buffer2 = GBuffer->new();
	 my $projecto = $buffer->newProjectCache( -name 			=> $patient->{obj}->project->name);
	 my $po = $projecto->getPatient($patient->{obj}->name);
	 
	 
	 foreach my $tr (@{$patient->{transcripts}}){
	 	#print "*";
	 	#my $po =  $patient->{obj};
	 	
	 	my $tr1 = $projecto->newTranscript($tr->{id});
	 	$tr->{obj} = $tr1;
	 	my $ret;
		
		unless (exists $ret->{data}){
			$ret  = image_coverage::image_cache_cnv ([$po], $tr1 );
			
		};
	 	#my $ret  = image_coverage::image_cnv ([$po], $tr->{obj} );
	 	my $primers = $tr1->getPrimers();
	 	$out_temp.= $cgi->start_Tr({style=>"background-color:#EAF1F1;width:100%"});
	 	my $obj_gene =  $tr1->getGene;
		my $minus = qq{<span class="glyphicon  glyphicon-minus" aria-hidden="true"></span>};
		my $name = "<i class='fa fa-2x fa-caret-right'></i>&nbsp;".$obj_gene->external_name."( ".$obj_gene->name." )";
		
		$name.= "&nbsp; $minus &nbsp;".$obj_gene->omim_inheritance if $obj_gene->omim_inheritance ;
		#$name.= "&nbsp; $minus &nbsp;".$obj_gene->short_phenotypes if $obj_gene->short_phenotypes ;
		$out_temp.= $cgi->td({style=>"background-color:#EAF1F1;width:100%",colspan=>$nb_col},$name."&nbsp; $minus &nbsp;".$tr->{name}." ".$tr1->name." ".$tr1->external_name." primers : ".scalar(@$primers) );
		$out_temp.= $cgi->end_Tr();
	 	my $z=0;
	 	my $th;
	 	 
			my $ths=[];
			my $tds=[];
			my $hprimers;
			my $nprimers;
			
			foreach my $exon (sort{$a->end*$a->strand <=> $b->end*$b->strand } @{$tr1->getExons()}){
				my $primers = $exon->getPrimers();
				 unless (scalar(@$primers)){
				 	push(@$ths,$cgi->th({class=>"infos"},"NC<br>".$exon->name));
				 	push(@$tds,$cgi->td({class=>"infos",style=>"background-color:rgb(255,255,255)",nowrap=>1},"NC" ) ) ;
				 	$z++;
				 }
				 
	 		foreach my $primer (sort{$a->end*$exon->strand <=> $b->end*$exon->strand } @$primers) {
	 			unless (exists $hprimers->{$primer->id}){
	 				$nprimers ++;
	 				$hprimers->{$primer->id} = $nprimers;
	 			}
	 				if ($z%$nb_col ==0 && $z>0){
	 					$out_temp.=$cgi->Tr(join("\n",@$ths)).$cgi->Tr(join("\n",@$tds));
	 					$ths=[];
	 					$tds=[];
	 				}
	 				$z++;
				push(@$ths,$cgi->th({class=>"infos"},$primer->multiplex."-".$hprimers->{$primer->id}."<br>".$exon->name));
				
	 				my $colors = $ret->{data}->{$primer->id}->{$po->name}->{colors};
	 				my $st_color = join(",",@$colors);
 				 push(@$tds,$cgi->td({style=>"background-color:rgb($st_color)"},$primer->cnv_score($po) ) ) ;
	 	}
	 }
	 $out_temp.=$cgi->Tr(join("\n",@$ths)).$cgi->Tr(join("\n",@$tds));
	 print $out_temp;
	 $out_temp ="";
	 }
	 
		  $out_temp.= $cgi->end_table();
	   $out_temp.= html::end_cadre($cgi);
	
	   print $out_temp;
	return "";
}



sub printSanger {
	my ($patient) = @_;
	my $out = html::print_cadre($cgi,"Sanger confirmed Variations");
	my $n = 0;
	my $type = "confirmed";
	foreach my $transcript (@{$patient->{transcripts}}){
		$n+=  scalar  @{$transcript->{$type}} if  $transcript->{$type} ;
	}

	return if $n ==0 ;
	$out.= printVariations2($patient,"Sanger Confirmed  Variations","confirmed");
	$out.= html::end_cadre($cgi);
	return $out;
}

sub printRejected {
	my ($patient) = @_;
	my $out = html::print_cadre($cgi,"Sanger  Unconfirmed Variations");
	my $n = 0;
	my $type = "rejected";
	foreach my $transcript (@{$patient->{transcripts}}){
		$n+=  scalar  @{$transcript->{$type}} if  $transcript->{$type} ;
	}

	return if $n ==0 ;
	$out.= printVariations2($patient,"Sanger Unconfirmed  Variations","rejected");
	$out.= html::end_cadre($cgi);
	return $out;
}


sub printValidated {
	my ($patient) = @_;
	my $out = html::print_cadre($cgi,"NGS validated Variations");
	my $type = "validated";
	my $n = 0;
	foreach my $transcript (@{$patient->{transcripts}}){
		$n+=  scalar  @{$transcript->{$type}} if  $transcript->{$type} ;
	}
	return if $n ==0 ;
	$out.= printVariations2($patient,"Validated Variations","validated");
	$out.=html::end_cadre($cgi);
	return $out;
}

sub printTodo {
	my ($patient) = @_;
	my $out = html::print_cadre($cgi,"Todo  Variations");
	my $type = "todo";
	my $n = 0;
	foreach my $transcript (@{$patient->{transcripts}}){
		$n+=  scalar  @{$transcript->{$type}} if  $transcript->{$type} ;
	}
	return if $n ==0 ;
	$out.= printVariations2($patient,"Todo Variations","todo");
	$out.=html::end_cadre($cgi);
	return $out;
}
sub printOther {
	my ($patient) = @_;
	warn "other 1";
	my $car = $project->maskImpact;
	my $st_impact;
	foreach my $v (sort {$car->{$a} <=> $car->{$b}} keys %$car ){
		if (  $project->getMaskCoding( $himpact2->{$vimpact}) & $project->getMaskCoding( $v))
		{
			$st_impact .= $project->maskImpactTextForLegend->{$v}." - ";
		}
	}
	$st_impact = "All" if $himpact2->{$vimpact} eq "low";
	my $st_freq = "All";

	if ($hfrequence->{$vfreq} eq "unique"){
		$st_freq = "Only novel";
	}
	if ($hfrequence->{$vfreq} eq "rare"){
		$st_freq = "public <=1% ";
		my $dejavu = $project->buffer->{config}->{dejavu}->{rare};
			$st_freq .= " dejavu <=$dejavu projects";
		
	}
	if ($hfrequence->{$vfreq} eq "occasional"){
		$st_freq = "public <=5%";
			my $dejavu = $project->buffer->{config}->{dejavu}->{occasional};
			$st_freq .= " dejavu <=$dejavu projects";
	}
	my $st = "<br><small>impact: [$st_impact] &nbsp; frequency: [$st_freq]</small>";
	my $out = html::print_cadre($cgi,"Other detected Variations".$st."");
	
	$out.= printVariationsOther($patient,"Other detected Variations ".$st." ","other");
	$out.= html::end_cadre($cgi);

	return $out;
}






sub print_graph_coverage {
	my ($patient) = @_;
	my $out;
my $nbt =0;
my $nb_p = int(scalar(@{$patient->{transcripts}}) % 15)+1;
my $cpp =1;;
#$out .= qq{<div class="container-fluid" >};
	#$out.= $cgi->legend("coverage limit=$cov_limit padding=$padding $cpp/$nb_p");
$out .= html::print_cadre($cgi,"coverage limit=$cov_limit padding=$padding  <span style='float:right'> $cpp/$nb_p</span>");
	$out.= $cgi->start_table({class=>"table table-condensed table-striped table-bordered",style=>"font-size:9px;width:1200px "});
foreach my $transcript (sort {$a->{name} cmp $b->{name}} @{$patient->{transcripts}}) {
	my $debug;
		$debug =1 if $transcript->{name} eq "NF1";
	$nbt++;
	if ($nbt%15 ==0 && $print ==1){
		$out.= $cgi->end_table();
		$out .= html::end_cadre($cgi);
		#$out .=qq{<div class="page-break">.</div>};
			$cpp++;
		  	$out .= html::print_cadre($cgi,"coverage limit=$cov_limit padding=$padding  <span style='float:right'> $cpp/$nb_p</span>");
			$out.= $cgi->start_table({class=>"table table-condensed table-striped table-bordered",style=>"font-size:9px;width:1200px "});
		
	}
			my $nb_exons = scalar @{$transcript->{exons}};
		my $name = "<h2 >".$transcript->{name}."  Cov = ".$transcript->{mean}."</h2> exons :".$nb_exons;
	
	
		#$out.= $cgi->start_table({class=>"bordered"});

		$out.= $cgi->start_Tr({style=>"background-color:#EAF1F1"});
		$out.= $cgi->td({colspan=>50},$transcript->{name}." ".$transcript->{obj}->{name}." ".$transcript->{external_name}."  Mean Coverage = ".$transcript->{mean} ." exons:".$nb_exons );
	
		
		$out.= $cgi->end_Tr();
		my @tds; 
		my @ths; 
		my$url = "//www.polyweb.fr/polyweb/images/polyicons/";
		my $nbe =0;
		foreach my $exon ( @{$transcript->{exons}}){
		
			$nbe++;
			
			if ($nbe%51 ==0 && $print ==1){
			#$out .= html::end_cadre($cgi);
		  	$out.= $cgi->Tr($cgi->th({class=>"info",nowrap=>1},\@ths));
			$out.= $cgi->Tr(join("\n",@tds));
			
		  		@tds=();
		  		@ths=();
		}
			my ($r,$g,$b) = @{$exon->{color2}};
			my $hash_font = {};
			 $hash_font = {nowrap=>1, style=>"background-color:rgb($r,$g,$b);"};# if ($g == 255 );


			
			my $name = $exon->{name};
			$name =~s/ex//;
			$name =~s/intron//;
			push(@ths,$name);
			my $text = $exon->{min};
			$text = $exon->{min} ."/".$exon->{raw_min}  if exists $exon->{raw_min};
			push(@tds,$cgi->td($hash_font,$text));
		
		}

 	$out.= $cgi->Tr($cgi->th({class=>"info",nowrap=>1},\@ths));
	$out.= $cgi->Tr(join("\n",@tds));
	
} #end transcripts
$out.=$cgi->end_table();
$out.=print_legend();
$out.= html::end_cadre($cgi);
#$out .= qq{</div>};
return $out;
}


sub print_legend {
	 my $out1;
	$out1.= $cgi->start_table({class=>"bordered"});
	 $out1.= $cgi->start_table({class=>"table table-striped  table-condensed table-bordered ",style=>"font-size: 8px;font-family:  Verdana;"});
	my $hash_font = {style=>"background-color:#00FF00",nowrap=>1};
	$out1.= $cgi->start_Tr();
	$out1.= $cgi->td("legend padding  +/- ".$padding );
	$out1.= $cgi->td($hash_font,"coverage >=".$cov_limit);
	$out1.= $cgi->td({style=>"background-color:#C8C8C8",nowrap=>1}," Not captured");
	$out1.= $cgi->td({style=>"background-color:yellow",nowrap=>1}," < $cov_limit not confirmed yet");
	$out1.= $cgi->td({style=>"background-color:blue;color:white",nowrap=>1},"Sanger confirmed");
	$out1.= $cgi->td({class=>"info",nowrap=>1},"To do");
	$out1.= $cgi->td({style=>"background-color:red",nowrap=>1},"Deleted homozygote confirmed");
	$out1.= $cgi->td({class=>"background-color:orange",nowrap=>1},"Deleted heterozygote confirmed");

	$out1.= $cgi->end_Tr();
	$out1.= $cgi->end_table();
	return $out1;
}




sub html2_bis {
	my ($data,$cgi) = @_;
my $out;
my $nb_case;


my $types = {
							he=>2,
							ho=>3,
							td=>-3,
							nc=>-5
							
						};
my $hvalidation = {
							2=>"he",
							3=>"ho",
							-3=>"td",
							-5=>"nc",
							
						};

$out.= $CSS;
if ($print){
$out.="<body onload='this.print();'>";
}
else {
	$out.="<body '>";
}
$out.= $cgi->start_div({id=>"all_html",name=>"all_html",jsId=>"all_html",class=>"patient"});


foreach my $patient (@$data) {
	$out .= print_project_summary($patient);
	$out.= $cgi->start_table();
	$out.= $cgi->start_Tr().$cgi->start_td();
	$out.= $cgi->start_fieldset({class=>"box"});
	$out.= $cgi->legend("coverage");
#$out .= qq{<div class="box" style="width:80%">};
#$out.= $cgi->h2("coverage");
foreach my $transcript (@{$patient->{transcripts}}) {
		my $name = "<h2 >".$transcript->{name}."  Cov = ".$transcript->{mean}."</h2>";
		$out.= $cgi->start_table({class=>"bordered"});
		my $nb = scalar @{$transcript->{exons}}+1;
		$out.= $cgi->start_Tr();
			$out.= $cgi->th({class=>"th1",colspan=>$nb},$transcript->{name}."  Mean Coverage = ".$transcript->{mean});
			$out.= $cgi->end_Tr();
		$out.= $cgi->start_Tr();
		$out.= $cgi->th("exons");
		my @tds; 
		push(@tds,$cgi->td("&nbsp"));
		foreach my $exon ( @{$transcript->{exons}}){
			my $hash_font = {class=>'success',nowrap=>1};;
			if($exon->{type} == -2){
				$hash_font = {bgcolor=>"red",style=>'color:black;',nowrap=>1};
			}
			if($exon->{type} == -1){
				$hash_font = {bgcolor=>"orange",style=>'color:black;',nowrap=>1};
			}
			if($exon->{type} == 0){
				$hash_font = {bgcolor=>"blue",style=>'color:black;',nowrap=>1};
			}
				if($exon->{type} == 1){
				$hash_font = {bgcolor=>"yellow",style=>'color:black;',nowrap=>1};
			}
			
			my $hash_font1 = {style=>'color:black;',nowrap=>1};
			
			my $name = $exon->{name};
			$name =~s/ex//;
			
			$out.= $cgi->th($hash_font1,$name);
			push(@tds,$cgi->td($hash_font,"&nbsp"));
			#push(@tds,$cgi->td($hash_font,$exon->{mean}));
#			$out.= $cgi->td($hash_font,$exon->{mean});
#			$out.= $cgi->td($hash_font,$exon->{min});
#			$out.= $cgi->td($hash_font,$exon->{type_string});

		
		}
			$out.= $cgi->end_Tr();
		$out.= $cgi->start_Tr();
		$out .= join("\n",@tds);
		$out.= $cgi->end_Tr();
		$out.=$cgi->end_table();
		$out.="<br>";
} #end transcripts
$out.=print_legend();
$out.= $cgi->end_fieldset();

$out.= qq{<div><br></div>};
$out.= $cgi->end_td().$cgi->end_Tr().$cgi->start_Tr().$cgi->start_td();
$out.= $cgi->start_fieldset({class=>"box"});
$out.= $cgi->legend("Variations");
$out.=print_table_variants($patient);
$out.= qq{<div><br></div>};
$out.= printVariations2($patient,"Sanger Confirmed  Variations","confirmed",$out);
$out.= qq{<div><br></div>};
$out.= printVariations2($patient,"Validated Variations","validated",$out);
$out.= qq{<div><br></div>};
$out.= $cgi->end_fieldset();
} #for patient	



$out.= $cgi->end_fieldset();
$out.= $cgi->end_td().$cgi->end_Tr().$cgi->end_table();
 return $out;
}




sub print_table_variants {
	my ($patient) = @_;
	 my $out1 = html::print_cadre($cgi,"Variations Summary");
	 $out1.="<br>";
	
	 $out1 .= $cgi->start_table({class=>"table table-striped  table-condensed table-bordered ",style=>"font-size: 8px;font-family:  Verdana;"});
	 
	#$out1.= $cgi->start_table({class=>"bordered"});
	
	my @th1 = ("Genes","transcripts", "Sanger confirmed","ngs validated","other variations","rejected variations");
	my @types = ( "confirmed","validated","other","rejected");

	$out1.= $cgi->start_Tr({class=>"sucsess"});
	#$out1.= $cgi->th({colspan=>6,class=>"th2"}," Summary");
	$out1.= $cgi->th({colspan=>6}," Summary");
	$out1.= $cgi->end_Tr();
	$out1.= $cgi->start_Tr({class=>"warning"});
	$out1.= $cgi->th(\@th1);
	$out1.= $cgi->end_Tr();
	
	foreach my $transcript (@{$patient->{transcripts}}) {
			$out1.= $cgi->start_Tr();
		
			my $z =0;
			 $z  = scalar  @{$transcript->{confirmed}} if  $transcript->{confirmed} ;
			 my $hash_font={};
			 $hash_font = {class=>"danger"} if $z >0;
			 	$out1.= $cgi->td($hash_font,$transcript->{name} );
			$out1.= $cgi->td($hash_font,$transcript->{external_name} );
			foreach my $type (@types) {
				my $n =0;
				$n+=  scalar  @{$transcript->{$type}} if  $transcript->{$type} ;
				$out1.= $cgi->td($hash_font,$n );
			}
			$out1.= $cgi->end_Tr();
	}
	$out1.= $cgi->end_table();
	$out1.= html::end_cadre($cgi);
	return $out1;
}


sub print_project_summary_pdf {
	my ($patient) = @_;
	my $out;
	
	my $pn = $patient->{name};
	my @genes_name = sort {$a cmp $b} map{$_->{name} } @{$patient->{transcripts}};
	my $nb_genes = scalar(@genes_name);
	my $nb = int(($nb_genes / 10) +0.5);
	my $r =0;
	my $z =0;
	my $part_genes = [];
	my $limit = 10;
	foreach my $gene (@genes_name){
		$r++ if $z %($limit)== 0 && $z >0;
		$z++;
		push(@{$part_genes->[$r]},$gene);
	}
	$nb = $limit if (scalar(@$part_genes) > 1);
	my $form_id = "form_".$pn;
	my $machine = $patient->{machine};
	my ($d1,$h1) = split(" ",$project->creation_date) ;
	my $pdate = Date::Tiny->from_string($d1);
	my $today = Date::Tiny->now;
	
	
	my $pdf = PDF::API2->new();
	
	my %font = (
    Helvetica => {
        Bold   => $pdf->corefont( 'Helvetica-Bold',    -encoding => 'latin1' ),
        Roman  => $pdf->corefont( 'Helvetica',         -encoding => 'latin1' ),
        Italic => $pdf->corefont( 'Helvetica-Oblique', -encoding => 'latin1' ),
    },
    Times => {
        Bold   => $pdf->corefont( 'Times-Bold',   -encoding => 'latin1' ),
        Roman  => $pdf->corefont( 'Times',        -encoding => 'latin1' ),
        Italic => $pdf->corefont( 'Times-Italic', -encoding => 'latin1' ),
    },
);
     my $page = $pdf->page;
     my $width     = 842;
	my $height    = 595; 
     $page->mediabox("A4");
     my $blue_box = $page->gfx;
	$blue_box->fillcolor('#2C3E50');
	$blue_box->rect( 45, 800, 500, 10 / mm );
	$blue_box->fill;
	my $headline_text = $page->text;
	 $headline_text->font( $font{'Helvetica'}{'Bold'}, 18/pt );
	 $headline_text->fillcolor('yellow');
	 $headline_text->translate( 50, 807);
	# $headline_text->text_right('USING PDF::API2');
	  $headline_text->text("Patient Summary ".$today->day."/".$today->month."/".$today->year);
	my $pdftable = new PDF::Table;
	
    my $some_data =[
    ["Name",
    ""],
    ["identifiant","XXX"],
    [
    "Validation",
    " -"],
    #... and so on
 ];

my $left_edge_of_table = 50;
 # build the table layout
  my $cell_props = [];
    $cell_props->[0][0] = {
        #Row 2 cell 1
        background_color => '#CCCC00',
        font_color       => 'blue',
    };
 $pdftable->table(
     # required params
     $pdf,
     $page,
     $some_data,
     x => $left_edge_of_table,
     w => 495,
     start_y => 799,
     start_h => 100,
     # some optional params
     next_y  => 750,
     next_h  => 500,
     padding => 5,
     padding_right => 10,
     border             => 0.1,
     background_color_odd  => "#FFFFFF",
     background_color_even => "#F0F0F0", #cell background color for even rows
     cell_props => $cell_props,
  );
#	warn "end";
        $pdf->saveas('/data-xfs/file.pdf');
        
        exit(0);
	my $pn = $patient->{name};
	my @genes_name = sort {$a cmp $b} map{$_->{name} } @{$patient->{transcripts}};
	my $nb_genes = scalar(@genes_name);
	my $nb = int(($nb_genes / 10) +0.5);
	my $r =0;
	my $z =0;
	my $part_genes = [];
	my $limit = 10;
	foreach my $gene (@genes_name){
		$r++ if $z %($limit)== 0 && $z >0;
		$z++;
		push(@{$part_genes->[$r]},$gene);
	}
	$nb = $limit if (scalar(@$part_genes) > 1);
	my $form_id = "form_".$pn;
	my $machine = $patient->{machine};
	my ($d1,$h1) = split(" ",$project->creation_date) ;
	my $pdate = Date::Tiny->from_string($d1);
	my $today = Date::Tiny->now;
	unless($print){
	$out .= qq{<div  data-dojo-type="dijit/Toolbar">};
			$out.= qq{<span class="spanpatient">Report   </span>};
			$out .= qq{<button dojoType="dijit.form.Button"   iconClass="dijitEditorIcon dijitEditorIconSave" onClick ="save_report('$pn','$form_id');"><span class="spanpatient">Save Report </span></button>};
			$out .= qq{</div>};
	}
	#$out.= $cgi->start_table();
	#$out.= $cgi->start_Tr().$cgi->start_td();
	#$out.= $cgi->start_fieldset({class=>"box"});
		$out.= html::print_cadre($cgi,"Patient Summary ".$today->day."/".$today->month."/".$today->year);
	#$out.= $cgi->legend
	$out.= $cgi->start_table({class=>"table table-striped table-condensed table-bordered table-hover",style=>"font-size: 8px;font-family:  Verdana;"});;
	$out.= $cgi->start_Tr();
	$out.= $cgi->th({class=>"th1"},["name","<input type ='text' size ='50'>"]);
	$out.= $cgi->end_Tr();
	$out.= $cgi->start_Tr();
	$out.= $cgi->th({class=>"th1"},["identifiant",$patient->{name}]);
	$out.= $cgi->end_Tr();
	$out.= $cgi->start_Tr();
	$out.= $cgi->th({class=>"th1"},"Validation");
	if (exists $patient->{user_validated}) {
		my ($d2,$h2) = split(" ",$patient->{date_validated}) ;
		my $pdate2 = Date::Tiny->from_string($d2);
		$out.= $cgi->td($pdate2->day."/".$pdate2->month."/".$pdate2->year);
	} 
	else {
		$out.= $cgi->th({class=>"th1"},"-");
	}
	
	$out.= $cgi->end_table();
	$out.=html::end_cadre($cgi);
	####
	## table 2
	#####
	$out.= $cgi->end_td().$cgi->end_Tr().$cgi->start_Tr().$cgi->start_td();
	
	
	$out.= html::print_cadre($cgi,"Run Summary");
	$out.= $cgi->start_table({class=>"table table-striped table-condensed table-bordered table-hover",style=>"font-size: 8px;font-family:  Verdana;width:500px"});;
	
#	$out.= $cgi->start_Tr({style=>"width:100%"});
#
#	
#	$out.= $cgi->th({class=>"default",colspan=>($nb+1),nowrap=>1},"Amplification");
#	$out.= $cgi->end_Tr()
	$out.=$cgi->start_Tr({class=>"default"});
	$out.= $cgi->td({class=>"default",rowspan=>scalar(@$part_genes)},"Genes");
	foreach my $genes_line (@$part_genes){
		$out.= $cgi->td({class=>"default"},$genes_line);
		$out.= $cgi->end_Tr().$cgi->start_Tr();
	}
	#$out.= $cgi->end_Tr().$cgi->start_Tr();
	$out.= $cgi->td({class=>"default"},"Methods");
	
	#$out.= $cgi->end_Tr().$cgi->start_Tr();
	$out.= $cgi->td({class=>"default",colspan=>$nb},uc($patient->{capture_type}) );
	$out.= $cgi->end_Tr().$cgi->start_Tr({class=>"default"});
	$out.= $cgi->td({class=>"default"},"Description").$cgi->td({colspan=>$nb},$patient->{capture_description});
	$out.= $cgi->end_Tr().$cgi->start_Tr({class=>"default"});
	$out.= $cgi->td({class=>"default"},"Version").$cgi->td({colspan=>$nb},$patient->{capture_version});
	$out.= $cgi->end_Tr().$cgi->start_Tr({class=>"default"});
	$out.= $cgi->td({colspan=>($nb+1)},"Sequencing");
	$out.= $cgi->end_Tr().$cgi->start_Tr({class=>"default"});
	$out.= $cgi->td({class=>"default"},"Machine Type");

	$out.= $cgi->td({colspan=>$nb},$machine);
	
	$out.= $cgi->end_Tr().$cgi->start_Tr({class=>"default"});
	
	$out.= $cgi->th({class=>"default"},"Date");
	
	#$out.= $cgi->end_Tr().$cgi->start_Tr();
	$out.= $cgi->td({colspan=>$nb},$pdate->day."/".$pdate->month."/".$pdate->year);
	$out.= $cgi->end_Tr().$cgi->start_Tr({class=>"default"});
	$out.= $cgi->th({class=>"default"},"Reference");
	$out.= $cgi->td({colspan=>$nb,class=>"th2"},$project->name());
	
	$out.= $cgi->end_table();
	

	
	$out.= html::end_cadre($cgi);
	
	# $out.= $cgi->end_td().$cgi->end_Tr().$cgi->end_table();
	return $out;
}


sub print_project_summary {
	my ($patient) = @_;
	my $out;
	
	my $pn = $patient->{name};
	
	my @genes_name = sort {$a cmp $b} map{$_->{name} } @{$patient->{transcripts}};
	my $nb_genes = scalar(@genes_name);
	my $nb = int(($nb_genes / 10) +0.5);
	my $r =0;
	my $z =0;
	my $part_genes = [];
	my $limit = 10;
	foreach my $gene (@genes_name){
		$r++ if $z %($limit)== 0 && $z >0;
		$z++;
		push(@{$part_genes->[$r]},$gene);
	}
	$nb = $limit if (scalar(@$part_genes) > 1);
	my $form_id = "form_".$pn;
	my $machine = $patient->{machine};
	my ($d1,$h1) = split(" ",$project->creation_date) ;
	my $pdate = Date::Tiny->from_string($d1);
	my $today = Date::Tiny->now;

		$out.= html::print_cadre($cgi,"Patient Summary ".$today->day."/".$today->month."/".$today->year);
	#$out.= $cgi->legend
	$out.= $cgi->start_table({class=>"table table-striped table-condensed table-bordered table-hover",style=>"font-size: 10px;font-family:  Verdana;width:1000px"});;
	$out.= $cgi->start_Tr();
	$name = uc($name);
	$out.= $cgi->th(["name","<input type ='text' size ='50' value='$name'></input>"]);
	$out.= $cgi->end_Tr();
	$out.= $cgi->start_Tr();
	$out.= $cgi->th(["identifiant",$patient->{name}]);
	$out.= $cgi->end_Tr();
	$out.= $cgi->start_Tr();
	$out.= $cgi->th("Validation");
	if (exists $patient->{user_validated}) {
		my ($d2,$h2) = split(" ",$patient->{date_validated}) ;
		my $pdate2 = Date::Tiny->from_string($d2);
		$out.= $cgi->td($pdate2->day."/".$pdate2->month."/".$pdate2->year);
	} 
	else {
		$out.= $cgi->th("-");
	}
	
	$out.= $cgi->end_Tr().$cgi->start_Tr();
	my $sex_eval = $patient->{obj}->compute_sex(); 
	
	$out.= $cgi->th("Sex : ");	
	 if ($sex_eval == -1){
	 	$out.= $cgi->td(qq{<i class="fa fa-question-circle fa-2x"></i>});
	 }
	elsif ($sex_eval == 1){
		$out.= $cgi->td( qq{<i class="fa fa-mars fa-2x" > </i>});
	}	
	elsif ($sex_eval == 2){
		$out.= $cgi->td(qq{  <b><i class="fa fa-venus fa-2x" > </i></b>});
	}	
			
	$out.= $cgi->end_table();
	$out.=html::end_cadre($cgi);
	####
	## table 2
	#####
	$out.= $cgi->end_td().$cgi->end_Tr().$cgi->start_Tr().$cgi->start_td();
	
	
	$out.= html::print_cadre($cgi,"Run Summary");
		$out.="<h4> Amplification</h4>";
	$out.= $cgi->start_table({class=>"table table-striped table-condensed table-bordered", style=>"font-size: 11px;font-family:  Verdana;width:1000px"});;
	
	$out.= $cgi->start_Tr(style=>"width:100%");


	#$out.= $cgi->th({colspan=>($nb+1)},"Amplification");
	$out.= $cgi->end_Tr().$cgi->start_Tr({class=>"default"});
	$out.= $cgi->th({rowspan=>scalar(@$part_genes)},"Genes");
	foreach my $genes_line (@$part_genes){
		$out.= $cgi->td({class=>"default"},$genes_line);
		$out.= $cgi->end_Tr().$cgi->start_Tr();
	}
	$out.=$cgi->end_table();
	$out.= $cgi->start_table({class=>"table table-striped table-condensed table-bordered", style=>"font-size: 11px;font-family:  Verdana;width:600px"});;
		$out.= $cgi->start_Tr();
	#$out.= $cgi->end_Tr().$cgi->start_Tr();
	$out.= $cgi->th({nowrap=>1},"Methods");
	
	#$out.= $cgi->end_Tr().$cgi->start_Tr();
	$out.= $cgi->td(uc($patient->{capture_type}) );
	$out.= $cgi->end_Tr().$cgi->start_Tr();
	$out.= $cgi->th("Description").$cgi->td($patient->{capture_description});
	$out.= $cgi->end_Tr().$cgi->start_Tr();
	$out.= $cgi->th("Version").$cgi->td($patient->{capture_version});
	$out.= $cgi->end_Tr().$cgi->start_Tr();
	
	$out.= $cgi->td("<h4>Sequencing</h4>");
	$out.= $cgi->end_Tr().$cgi->start_Tr();
	$out.= $cgi->th("Machine Type");

	$out.= $cgi->td($machine);
	
	$out.= $cgi->end_Tr().$cgi->start_Tr();
	
	$out.= $cgi->th("Date");
	
	#$out.= $cgi->end_Tr().$cgi->start_Tr();
	$out.= $cgi->td($pdate->day."/".$pdate->month."/".$pdate->year);
	$out.= $cgi->end_Tr().$cgi->start_Tr();
	$out.= $cgi->th("Reference");
	$out.= $cgi->td($project->name());
	$out.= $cgi->end_Tr().$cgi->start_Tr();
	$out.= $cgi->th("Genome Version");
	$out.= $cgi->td($project->getVersion());
	$out.= $cgi->end_Tr().$cgi->start_Tr();
	$out.= $cgi->th("Coverage : ");
	my $p = $patient->{obj};
	my $cov = $p->coverage();
	$out.= $cgi->td([$cov->{mean}." (30X : ".$cov->{'30x'}."%)"]);
	$out.= $cgi->end_Tr().$cgi->start_Tr();
	$out.= $cgi->td("<h4>Pipeline</h4>");
	$out.= $cgi->end_Tr().$cgi->start_Tr();
	$out.= $cgi->th("Variant Calling ");
	
	$out.= $cgi->td(join(":",@{$patient->{obj}->getCallingMethods}));
	$out.= $cgi->end_Tr().$cgi->start_Tr();
	$out.= $cgi->th("Alignment ");
	$out.= $cgi->td($patient->{obj}->alignmentMethod);
	$out.= $cgi->end_Tr().$cgi->start_Tr();
	my $buffer = $patient->{obj}->buffer();
	$out.= $cgi->th("Gencode ");
	$out.= $cgi->td($project->get_gencode_description()->{version});
	$out.= $cgi->end_Tr().$cgi->start_Tr();
	$out.= $cgi->th("gnomad ");
	$out.= $cgi->td($buffer->description_public_lmdb_database("gnomad-exome")->{version});
	$out.= $cgi->end_Tr().$cgi->start_Tr();
	$out.= $cgi->th("hgmd");
	
	$out.= $cgi->td($buffer->description_public_lmdb_database("hgmd")->{version});
	$out.= $cgi->end_Tr().$cgi->start_Tr();
	$out.= $cgi->th("clinvar ");
	$out.= $cgi->td($buffer->description_public_lmdb_database("clinvar")->{version});
	$out.= $cgi->end_Tr().$cgi->start_Tr();
	
		
	

	$out.= $cgi->end_table();
	
	$out.= html::end_cadre($cgi);
	
	# $out.= $cgi->end_td().$cgi->end_Tr().$cgi->end_table();
	return $out;
}


sub printDetailCoverage {
	my ($patient) = @_;
	my $out = html::print_cadre($cgi,"Coverage");
	my $nb_tr= scalar(@{$patient->{transcripts}});
	my $limit =5;
	for (my $i=0;$i<$nb_tr;$i+=$limit){
	my $final = $i+$limit;
	$final = $nb_tr if $final > $nb_tr;
	my @temp_transcript;
	for (my $j=$i;$j<$final;$j++){
		push(@temp_transcript ,$patient->{transcripts}->[$j]);
	}
	#my @temp_transcript = $patient->{transcripts}->[$i..$final];
	
	$out.= $cgi->start_table({class =>"hor-minimalist-b",id=>"test"});
	$out.= $cgi->start_Tr();
	my $nb_tr=0;
	
	foreach my $transcript (@temp_transcript){
		
		my $name = "<center><h2 class='patient'>".$transcript->{name}." ".$transcript->{id}." ".$transcript->{external_name}." <br> Cov = ".$transcript->{mean}."</h2></center>";
		#my $name = $transcript->{name}." : Coverage = ".$transcript->{mean};
		$out.= qq{<th scope="col">$name</th>};
	}
	$out.= $cgi->end_Tr();
	$out.= $cgi->start_Tr();
	foreach my $transcript (@temp_transcript){
		
		
		$out.= $cgi->start_td();
		$out.= $cgi->start_table({id=>"coverage"});
		$out.= $cgi->start_Tr();
		$out.= qq{<th scope="col">exons</th>};
		$out.= qq{<th scope="col">mean</th>};
		$out.= qq{<th scope="col">min</th>};
		$out.= qq{<th scope="col">Sanger</th>};
		$out.= $cgi->end_Tr();
		
		$out.= "<tbody>";
		
		#foreach my $exon (sort {$a->{type} <=> $b->{type}} @{$transcript->{exons}}){
		foreach my $exon ( @{$transcript->{exons}}){
			my $hash_font = {nowrap=>1};
			if($exon->{type} == -2){
				$hash_font = {bgcolor=>"red",style=>'color:black;',nowrap=>1};
			}
			if($exon->{type} == -1){
				$hash_font = {bgcolor=>"orange",style=>'color:black;',nowrap=>1};
			}
			if($exon->{type} == 0){
				$hash_font = {bgcolor=>"green",style=>'color:black;',nowrap=>1};
			}
				if($exon->{type} == 1){
				$hash_font = {bgcolor=>"grey",style=>'color:black;',nowrap=>1};
			}
			$out.= $cgi->start_Tr();
			$out.= $cgi->td($hash_font,$exon->{name});
			$out.= $cgi->td($hash_font,$exon->{mean});
			$out.= $cgi->td($hash_font,$exon->{min});
			$out.= $cgi->td($hash_font,$exon->{type_string});
#			if($exon->{type} == 2){
#				$out.= "<font color='990000'>".$cgi->td($exon->{type_string})."</font>";
#			}
#			else {
#			$out.= $cgi->td($exon->{type_string});
#			}
			$out.= $cgi->end_Tr();
		}
		$out.= "</tbody>";
	$out.= $cgi->end_table();
	$out.= $cgi->end_td();
	}
	$out.= $cgi->end_Tr();
	$out.= $cgi->end_table();
	}
	$out.= html::end_cadre($cgi);
	return $out;
}

sub html {	
my ($data,$cgi) = @_;
my $out;
my $nb_case;


my $types = {
							he=>2,
							ho=>3,
							td=>-3,
							nc=>-5
							
						};
my $hvalidation = {
							2=>"he",
							3=>"ho",
							-3=>"td",
							-5=>"nc",
							
						};

$out.= $CSS;
if ($print){
$out.="<body onload='this.print();'>";
}
else {
	$out.="<body '>";
}
$out.= $cgi->start_div({id=>"all_html",name=>"all_html",jsId=>"all_html",class=>"patient"});


foreach my $patient (@$data){
 	#
 	$out .=printDetailCoverage($patient);
 	$out.= print_project_summary($patient);
	$out.= printVariations($patient,"Sanger Confirmed  Variations","confirmed",$out);
	$out.= printVariations($patient,"Validated Variations","validated",$out);
	$out.= printVariations($patient,"Other Variations","other",$out);
	$out.= $cgi->h4({class=>"patient"},"---");
	$out.=printVariations($patient,"Rejected Variations","rejected");
	$out.=printVariations($patient,"todo Variations ","todo",1);

	my $pname = $patient->{name};
	my $prname = $project->name();
	$out.= $cgi->h1({class=>"patient"},"Comment :");
	#my $sttr = $cgi->param("transcripts");
	my $user_name = $cgi->param("user_name");
	my $conclusion = $patient->{conclusion};

#$out.= qq{<form id = "$form_id">};
#$out.= qq{
#	<input type="hidden"  name="save"  value="1">
#	<input type="hidden"  name="patients"  value="$pname">	
#	<input type="hidden"  name="user_name"  value="$user_name">	
#	<input type="hidden"  name="project"  value="$prname">
#	<input type="hidden"  name="transcripts"  value="all">
#	<textarea id="textarea2" name="textarea2" data-dojo-type="dijit/form/SimpleTextarea" rows="4" cols="50" style="width:auto;">$conclusion</textarea>
#};
#$out.= qq{</form>};
$out.= $cgi->h4({class=>"patient"},"***");
}




$out.= $cgi->end_div();
# print $cgi->textarea(-name=>"conclusion",
# 					  #-id =>"textarea",
# 					  -rows=>10,
# 					   -columns=>80 );
 return $out;
exit(0);
	
}
	my $uid= 0;


sub print_variation_td{
	my ($variation,$cat,$rowspan,$nbv,$nobutton) = @_;
	my $score = "impact_score";
	my $out;
	my $category=0;		
	my $text = $variation->{$cat};
	my $bg = "background-color:#F7F7F7";
	$bg = "background-color:#FFFFFF" if $nbv%2 ==0;

	if  ($cat ne 'igv' && $cat ne 'alamut' && $cat ne 'align' && $cat ne  'hgmd'){
	if($text =~ /,/   ){
			$text =~ s/,/<BR>/;
	 }
#	elsif (length($text)>30 && $cat ne "var_name" && $cat ne "ngs"  && $cat ne "ratio" && $cat ne "in_this_run" && $cat !~/impact/ && $text !~ /span/){
#					my $a = substr $text,0,15;
#					my $b = substr $text,-8;
#					$text = $a."...".$b;
#				
#					
#				}
	}
				if ($cat eq "var_name"){
					
						$text = qq{$text};
						
				
				}
				if ($cat eq "consequence"){
					my @terms;
					foreach my $term (split("<BR>",$text)){
						push(@terms,$buffer->get_annotation_terms($term));
					}
					$text = join("<BR>",@terms);
				}
			#	warn $cat;
			#	if ($cat eq "consequence"){
			#		$text = $buffer->get_annotation_terms($text);
			#	}
				
				if ($variation->{sanger} =~ /todo/ && $nobutton ne 2 ) {
					
					$out.= $cgi->td({rowspan=>$rowspan,bgcolor=>"beige",style=>'color:black'},$text);
				}
				elsif ($variation->{sanger} =~ /confirm/ && $nobutton ne 2){
					#if ($nobutton ==2){
						#$out.= $cgi->td({rowspan=>$rowspan,style=>'color:black'},$text);
					#}
					$out.= $cgi->td({rowspan=>$rowspan,bgcolor=>"yellow",style=>'color:blue'},$text);
				}	elsif ($variation->{type} =~ /valid/ && $nobutton ne 2){
					$out.= $cgi->td({rowspan=>$rowspan,bgcolor=>"orange",style=>'color:blue;'},$text);
				}
				elsif ($variation->{sanger} =~ /rejected/ && $nobutton ne 2){
					$out.= $cgi->td({rowspan=>$rowspan,bgcolor=>"#6B0000",style=>'color:white;'},$text);
						$out.= $cgi->td({rowspan=>$rowspan,bgcolor=>"black",style=>'color:white;'},$text);
				}
				elsif ($variation->{ngs} =~ /dup/ && $nobutton ne 2){
					$out.= $cgi->td({rowspan=>$rowspan,bgcolor=>"#DDDDDD",style=>'color:white;'},$text);
				}
				else {
					my %style =%{$style_td->{$variation->{$score}}};
					$style{rowspan} = $rowspan;
					$out.= $cgi->td(\%style,$text);
				}
	return $out;
}


sub print_variation_td_edit{
	my ($variation,$cat,$rowspan,$nbv,$nobutton,$patient) = @_;
#	die();
	my $score = "impact_score";
	my $out;
	my $category=0;		
	my $text = $variation->{$cat};
	my $bg = "background-color:#F0F0F0";
	$bg = "background-color:#FFFFFF" if $nbv%2 ==0;
	if ($cat eq 'status'){
		$out= $cgi->td({rowspan=>$rowspan,nowrap=>1,style=>$bg} , $text);
		return $out;
	}
	if  ($cat ne 'igv' && $cat ne 'alamut' && $cat ne 'align' && $cat ne 'status' && $cat ne "hgmd"){
	if($text =~ /,/   ){
			$text =~ s/,/<BR>/;
	 }
#	elsif (length($text)>30 && $cat ne "var_name" && $cat ne "ngs"&& $cat ne "hgmd" && $cat ne "ratio" && $cat ne "cosmic" && $cat !~/impact/ && $cat ne "in_this_run" && $cat ne "trio" && $text !~/span/){
#				
#					my $a = substr $text,0,15;
#					my $b = substr $text,-8;
#					$text = $a."...".$b;
#				
#					
#				}
	}
				if ($cat eq "var_name"){
					if (length($text)>50 && $text !~/http/) {
					#	$text = qq{<a href="#" data-toggle="tooltip" data-placement="top" title="$text">coucou</a>};#.$variation->{nomenclature};
						my ($chr,$start,$ref,$alt) = split("_",$text);
						my $len = abs(length($ref) - length($alt));
						my $end = $start+$len;
						my $text2 = $chr.".Del $start..$end";
						$text = qq{<a href="#" data-toggle="tooltip" data-placement="top" title="$text">$text2</a>};
					}
				
					else {
						$text = update::printSimpleBadge(qq{$text});
					}
					
				}
				if($cat eq "ngs"){
					my @t = split("<br>",$text);
				 	$text = "<table style='background-color:white;color:black;border-color:black;font-size:".$update::fts.";width:auto'> <col style='width: 20%;' /> <col style='width: 10%;' /> <col style='width: 40%;' /> <col style='width: 30%;' /><tr>";
					if ($print){
					 $text = "<table  style='border-width:0px!important;padding:1px!important;border-style:none!important;background-color:white;color:black;border-color:black;font-size:".($update::fts-1).";width:auto!important '> <col style='width: 15%;' /> <col style='width: 15%!important;' /> <col style='width: 40%;' /> <col style='width: 30%;' /><tr>";
					}	
					my @a;
					my $r;
					my @color= ("#ECF7F8","#F9F9F9");
					my $z=0;
					foreach my $l (@t){
						$z++;
						$text .= "<tr>";
						($a[0],$r) = split(":",$l);
						($a[1],$r) = split(/\(/,$r);
						$r =~ s/\)//;
						$a[2] = $r;
						
						my ($a1,$a2,$r3) = split("/",$r);
						$a1 =1 if ($a1+$a2) ==0;
						$a[3] = int(($a2/($a1+$a2))*100)."%";
						#
							
		
						#				}
						if ($print){
							$text .=$cgi->td({style=>"border-top-right:0px!;border-top-left:0px!;border-top-width:1px!;border-bottom-width:1px!;border-style:solid!important;text-align: center;background-color:".$color[$z%2].";padding:1px!important;border-color:#DAEFF1"},\@a)."</tr>";
		
								}
								else {
									$text .=$cgi->td({style=>"text-align: center;background-color:".$color[$z%2].";border-style:solid;border-width:1px;padding:2px;border-color:#DAEFF1"},\@a)."</tr>";
								}
					} 
					 $text .="</table>";
				}
				if  ($cat eq "hgmd") {
					$text = update::printSimpleBadge(qq{$text});
				}
				if  ($cat eq "gene"  ) {
					$text = update::printSimpleBadge(qq{$text});
				}
				if  ($cat eq "freq" or $cat eq "freq_ho" ) {
					$text = update::printSimpleBadge(qq{$text},1);
				}
				if ($cat eq "genomique" ){
					$text = update::printSimpleBadge(qq{$text});
				}
				if ($cat eq "transcript" or $cat eq "exon" or $cat eq "nomenclature" or $cat eq "consequence" or $cat eq "codons_AA" or $cat eq "codons"){
					$text = update::printSimpleBadge(qq{$text},2);
				}
				if ($cat eq "similar_projects" ){
					$text  =~ s/<[^>]*>//gs;
					$text = update::printSimpleBadge(qq{$text},3);
				}
				if ($cat eq "polyphen"  or  $cat eq "sift" ){
					$text = update::printSimpleBadge(qq{$text});
				}
				if ($cat eq "min_pop" or $cat eq "max_pop" ){
					#($plain_text = $html_text) =~ s/<[^>]*>//gs; 
					$text  =~ s/<[^>]*>//gs;
					$text = update::printSimpleBadge(qq{$text},1);
					#$text = update::printSimpleBadge(qq{$text});
				}
				if ($cat eq "igv"){
					my $fam = $patient->{obj}->getFamily();
					my @bams;
					my @names;
					foreach my $p (@{$fam->getPatients()}){
						push(@bams,$p->bamUrl);
						push(@names,$p->name());
					}
					
					my $f =  join(";",@bams);#$patient->{obj}->bamUrl;;
					my $l = $variation->{genomique};
					my $v = $variation->{ref_allele}."/".$variation->{allele};	
					my $gn = $patient->{obj}->project->getVersion();
					my $pnames = join(";",@names);
					$text =qq{<button dojoType="dijit.form.Button"   iconClass="igvIcon" onclick='launch_web_igv_js("$project_name","$pnames","$f","$l")' style="color:black"></button>};
				}
				if ($cat eq "consequence"){
					my @terms;
					foreach my $term (split("<BR>",$text)){
						$term = "intergenic" if $term eq "-";
						push(@terms,$buffer->get_annotation_terms($term));
					}
					$text = join("<BR>",@terms);
				}
				
				elsif ($cat eq "deja_vu" or $cat eq "similar_projects" ) {
					my $url = $deja_vu_url.$variation->{id};
					 if (exists $variation->{dup}){
					 	$text = qq{<a href="$url" target="_blank" style="color:black;font-weight:bold">$text</a>};
					 }
					 else {
					$text = qq{<a href="$url" target="_blank" style="color:black;font-weight:bold">$text</a>};
					 }
					 $text = update::printSimpleBadge(qq{$text},3);
					
				}	
				if ( $cat eq "in_this_run"){

					$text =~ s/white/black/; 

					 $text = update::printSimpleBadge(qq{$text},3);
				}
		if (exists $variation->{dup}){
			
			$variation->{scaled_score} = 99;
#			warn $bg;
		}
		
		if ($variation->{scaled_score} > 1){	
			my %style;
			if (exists $style_score->{$variation->{scaled_score}}) { %style = %{$style_score->{$variation->{scaled_score}}}; }
			elsif ($variation->{scaled_score} < 1) { %style = %{$style_score->{1}}; }
			elsif ($variation->{scaled_score} < 2) { %style = %{$style_score->{2}}; }
			elsif ($variation->{scaled_score} < 3) { %style = %{$style_score->{3}}; }
			elsif ($variation->{scaled_score} < 4) { %style = %{$style_score->{4}}; }
			elsif ($variation->{scaled_score} > 4) { %style = %{$style_score->{4}}; }
			else {
				warn Dumper $variation;
				die;
			}
			$style{rowspan} = $rowspan;
			$out.= $cgi->td(\%style,$text);
		}
#		elsif ($cat eq "freq"){	
#		my %style =%{$style_score->{$variation->{freq_score}}};
#		$style{rowspan} = $rowspan;
#		$out.= $cgi->td(\%style,$text);
#		}
#		elsif ($cat eq "score"){	
#		my %style =%{$style_score->{$variation->{scaled_score}}};
#		$style{rowspan} = $rowspan;
#		$out.= $cgi->td(\%style,$text);
#		}
		else{
#		#	$style{rowspan} = $rowspan;
		$out.= $cgi->td({rowspan=>$rowspan,style=>$bg} , $text);
		}
	return $out;
}

sub printVariationsOther {
	my ($patient,$title,$type,$nobutton) = @_;
	$nobutton =1 if $print;
	my $out;
	my @buttons = ("igv","alamut","align");
 	@buttons = () unless $edit_mode;
	my @infos = @headers;
	
	my @edit_buttons=("Validation");
	@infos = (@buttons,@infos);

	my $s_id=$patient->{name};
	my %genes;
	my %nm;
	foreach my $tr (@{$patient->{transcripts}}){
		#die() if $tr->{name} eq "RPGRIP1L";
		push(@{$genes{$tr->{name} } },$tr);
		$tr->{obj} = $project->newTranscript($tr->{id});# unless $tr->{obj} ;		
		#warn $tr->{id};
		my $n =  $tr->{obj}->name;
		my $t = $tr->{external_name};
		my @te = split(";",$t);
		my $ext = $t;
		if (scalar(@te) > 1 ){
			@te = sort{$a cmp $b} @te;
			$ext = $te[0];
			$ext .= ";".$te[1] if $te[1] =~/NM_/;
		}
		$nm{$n} =  $ext;
		#warn  $tr->{external_name};
	}
	my $gene_var;
	foreach my $k (keys %genes){
		my $g = $genes{$k};
		my $find;
		foreach my $tr (@$g){
			foreach my $v (@{$tr->{$type}}){
				$find++;
				push (@{$gene_var->{$k}->{$v->{id}}},$v);
			}
		}
		delete $genes{$k} unless $find;
	}
	#warn Dumper keys %genes;
	my $nb_line = 0;
	my $first_line = 0;
	$out .= $cgi->start_table({class=>"table table-striped table-condensed table-bordered table-hover table-mybordered",style=>"vertical-align:middle;text-align: left;font-size: 6px;font-family:  Verdana;"});
	
	foreach my $k (sort {$a cmp $b} keys  %genes) {
		warn $k;
		my $debug ;
		$debug =1 if $k eq "SOX10";
		my $nb_line_next_gene = $nb_line + scalar (keys %{$gene_var->{$k}});
		
		if ( $first_line>0 && $print ==1){
	#	$out .= html::end_cadre($cgi);
		#$out .=qq{<div class="page-break">.</div>};
			#$cpp++;
	#	  	$out .= html::print_cadre($cgi,"$title ");
			$nb_line = 0;
	}
	$first_line = 1;
		my $gene = $genes{$k};
		
		my $n1 = join("-",map{$_->{external_name}} @$gene);
		my $obj_gene =  $gene->[0]->{obj}->getGenes->[0];
		#my $name = "<big>".$obj_gene->name." ".$obj_gene->external_name."[".$obj_gene->omim_inheritance."]&nbsp"." ".$obj_gene->phenotypes."</big> ".$n1;
		#my $name = "<table  class='table-bordered  table-mybordered' style='background:#D9C7ED;font-size: 10px;padding: 10px;width: 50%'> <tr><td><big>".$obj_gene->external_name."</big></td><td> ".$obj_gene->name."</td><td>[".$obj_gene->omim_inheritance."]</td><td>"." ".$obj_gene->phenotypes."</td></tr></table>";
		my $minus = qq{<span class="glyphicon  glyphicon-minus" aria-hidden="true"></span>};
		my $name = "<i class='fa fa-2x fa-caret-right'></i>&nbsp;".$obj_gene->external_name."( ".$obj_gene->name." )";
			$name.= "&nbsp; $minus &nbsp;".$obj_gene->omim_inheritance if $obj_gene->omim_inheritance ;
			$name.= "&nbsp; $minus &nbsp;".$obj_gene->short_phenotypes if $obj_gene->short_phenotypes ;
			$out .= $cgi->end_table();
			$out .= $cgi->start_table({class=>"table table-striped table-condensed table-bordered table-hover table-mybordered",style=>"vertical-align:middle;text-align: center;font-size:".$update::fts.";font-family:  Verdana;width:auto"});
			$out.= $cgi->start_Tr();
			$out.= $cgi->th({style=>"vertical-align:middle;background:#D9C7ED;border: 2px solid black;",colspan=>scalar(@infos)},"$name");
			$out.= $cgi->end_Tr();
			
		
	my $t = scalar(@infos);
	my $pc = int(100/$t);
		foreach my $c (@infos){
			if($c eq "var_name" or $c eq "ngs" or $c eq "transcript" or $c eq "clinvar" or $c eq "genomique"  ){
				$out.= qq{<col width="5%"\>};
			}
			elsif ($c eq "nomenclature" or $c eq "codons" or $c eq "consequence" or $c eq "max_pop" or $c eq "min_pop" or $c eq "hgmd" or $c eq "local"){
				$out.= qq{<col width="3%"\>};
			}
			else {
				$out.= qq{<col width="2%"\>};
			}
			
	}
my $nbv =0;
$out.= $cgi->start_Tr();
	foreach my $c (@infos){
			$out.= qq{<th scope="col">$c</th>};
	}
	$out.= $cgi->end_Tr();
foreach my $tvariations ( sort {$a->[0]->{trans} <=> $b->[0]->{trans}} values %{$gene_var->{$k}}){	
		my @variations = sort {$a->{trans} <=> $b->{trans}} @$tvariations;
		my $variation = $variations[0];
		$out.= $cgi->start_Tr({style=>"vertical-align:middle"});
		my $rowspan =scalar(@variations); 
		warn $rowspan if $debug;
		$nb_line+=$rowspan;
		my $cat_start;
		$nbv ++;
		for ( my $ind = 0; $ind< @infos;$ind++){
			
			my $cat = $infos[$ind];
			
			$out.=print_variation_td_edit($variation,$cat,$rowspan,$nbv,$nobutton);
			 $cat_start = $ind;
			last if $cat eq "genomique";
			
		}
		$cat_start ++;
		warn scalar(@variations);
		foreach my $variation (sort {$a->{trans} <=> $b->{trans}} @variations){
		warn "++".$cat_start;
		warn @infos;
			for (my $ind = $cat_start;$ind< @infos;$ind++){
				
				my $cat = $infos[$ind];
				if ($cat eq "transcript"){
					if (exists $nm{$variation->{transcript}}){
					#$variation->{transcript} = "<table><tr><td>".$variation->{transcript}."</td><tr></tr><td>".$nm{$variation->{transcript}}."</td></table>";
					$variation->{transcript} = $variation->{transcript}."++<br>".$nm{$variation->{transcript}};
					}
					
				}
				$out.=print_variation_td_edit($variation,$cat,1,$nbv,$nobutton);
			}
				$uid++;
				$out.= $cgi->end_Tr();
				$out.= $cgi->start_Tr({style=>"vertical-align:middle"});
		}
		if ( $nobutton ne 2 ){
			
			$out.= $cgi->start_Tr();
			$out.= $cgi->td({style=>"vertical-align:middle",colspan=>scalar(@infos)},[$variation->{trio}]);
			$out.= $cgi->end_Tr();
		}
		#warn $nobutton;
	}
	}
	
$out.= $cgi->end_table();
	$out.="</body>";
	return $out;
}



sub printVariations2 {
	my ($patient,$title,$type,$nobutton) = @_;
	$nobutton =1 if $print;
	my $out;
	my @buttons = ("igv","alamut","align");
 	@buttons = () unless $edit_mode;
	my @infos = @headers;#("gene","var_name","sanger","ngs","ratio","caller","genomique","transcript","exon","nomenclature","consequence","codons","codons_AA","freq","deja_vu","similar_projects","in_this_run", "polyphen","sift");
	#my @infos = ("var_name","sanger","ngs","ratio","genomique","transcript","exon","nomenclature","consequence","codons","codons_AA","freq","deja_vu","similar_projects","in_this_project", "polyphen","sift");
	
	my @edit_buttons=("Validation");
	@infos = (@buttons,@infos);
	# my @infos = ("gene","var_name","sanger","ngs","genomique","transcript","exon","nomenclature","consequence","codons","codons_AA","freq","deja_vu","in_this_project", "polyphen","sift");

	#	if ($type eq "all" && $nobutton ne 2){
	#		@infos = ("impact","type","var_name","sanger","ngs","genomique","exon","nomenclature","consequence","codons","codons_AA","freq","deja_vu","in_this_project", "polyphen","sift");
	#	}
	my $s_id=$patient->{name};
	my %genes;
	my %nm;
	foreach my $tr (@{$patient->{transcripts}}){
		push(@{$genes{$tr->{name} } },$tr);
		my $n =  $tr->{obj}->name;
		my $t = $tr->{external_name};
		my @te = split(";",$t);
		my $ext = $t;
		if (scalar(@te) > 1 ){
			@te = sort{$a cmp $b} @te;
			$ext = $te[0];
			$ext .= ";".$te[1] if $te[1] =~/NM_/;
		}
		$nm{$n} =  $ext;
	}
	my $gene_var;
	foreach my $k (keys %genes){
		my $g = $genes{$k};
		my $find;
		foreach my $tr (@$g){
			
			foreach my $v (@{$tr->{$type}}){
				$find++;
				push (@{$gene_var->{$k}->{$v->{id}}},$v);
			}
		}
		delete $genes{$k} unless $find;
	}
	my $nb_line = 0;
	my $first_line = 0;
	
	foreach my $k (sort {$a cmp $b} keys  %genes){
		my $nb_line_next_gene = $nb_line + scalar (keys %{$gene_var->{$k}});
		if ( $first_line>0 && $print ==1){
		$out .= html::end_cadre($cgi);
		#$out .=qq{<div class="page-break">.</div>};
			#$cpp++;
		  	$out .= html::print_cadre($cgi,"$title ");
		  
	
			$nb_line = 0;
			
		
	}
	$first_line = 1;
		my $gene = $genes{$k};
		my $n1 = join("-",map{$_->{external_name}} @$gene);
		#my $name = $gene->[0]->{name}." ".$n1;
		my $obj_gene =  $gene->[0]->{obj}->getGenes->[0];
		my $minus = qq{<span class="glyphicon  glyphicon-minus" aria-hidden="true"></span>};
		my $name = "<i class='fa fa-2x fa-caret-right'></i>&nbsp;".$obj_gene->external_name."( ".$obj_gene->name." )";
		
		$name.= "&nbsp; $minus &nbsp;".$obj_gene->omim_inheritance if $obj_gene->omim_inheritance ;
		$name.= "&nbsp; $minus &nbsp;".$obj_gene->short_phenotypes if $obj_gene->short_phenotypes ;
		$out .=  $cgi->start_div({class=>" panel panel-default",style=>"font-size: 11px;font-family:  Verdana;"});
		#$out.= $cgi->div({class=>"panel-heading  "},$name);
		
		$out .= $cgi->start_table({class=>"table table-striped table-condensed table-bordered table-hover table-mybordered",style=>"vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;"});
		$out.= $cgi->start_Tr();
			$out.= $cgi->th({style=>"vertical-align:middle;background:#D9C7ED;border: 2px solid black;",colspan=>scalar(@infos)},"$name");
			
			$out.= $cgi->end_Tr();
		#$out.= $cgi->start_table({class=>"bordered"});
	#$out.= $cgi->start_table({id=>"hor-minimalist-b"});
	my $t = scalar(@infos);
		$out.= $cgi->start_Tr();

	#$out.= $cgi->th({scope=>"row"},@infos);
	foreach my $c (@infos){
			$out.= qq{<th scope="col">$c</th>};
	}
	if ($edit_mode ){
		foreach my $c (@edit_buttons){
			$out.= qq{<th scope="col">$c</th>};
	}
	#	$out .= $cgi->th({scope=>"col"},@edit_buttons);
	}
	$out.= $cgi->end_Tr();
my $nbv =0;
foreach my $tvariations ( sort {$a->[0]->{trans} <=> $b->[0]->{trans}} values %{$gene_var->{$k}}){	
		my @variations = sort {$a->{trans} <=> $b->{trans}} @$tvariations;
		my $variation = $variations[0];
		$out.= $cgi->start_Tr({style=>"vertical-align:middle"});
		my $rowspan =scalar(@variations); 
		$nb_line+=$rowspan;
		my $cat_start;
		$nbv ++;
		for ( my $ind = 0; $ind< @infos;$ind++){
			
			my $cat = $infos[$ind];
			$out.=print_variation_td_edit($variation,$cat,$rowspan,$nbv,$nobutton);
			 $cat_start = $ind;
			last if $cat eq "genomique";
			
		}
		$cat_start ++;
		foreach my $variation (sort {$a->{trans} <=> $b->{trans}} @variations){
			for (my $ind = $cat_start;$ind< @infos;$ind++){
				my $cat = $infos[$ind];
					if ($cat eq "transcript"){
					
					if (exists $nm{$variation->{transcript}}){
					#$variation->{transcript} = "<table><tr><td>".$variation->{transcript}."</td><tr></tr><td>".$nm{$variation->{transcript}}."</td></table>";
					$variation->{transcript} = $variation->{transcript}."<br>".$nm{$variation->{transcript}};
					}
				}
				$out.=print_variation_td_edit($variation,$cat,1,$nbv,$nobutton);
			}
			$uid++;
				$out .= print_validation_button($patient,$variation,$uid,$variation) if $nobutton == 2;
				$out.= $cgi->end_Tr();
				$out.= $cgi->start_Tr({style=>"vertical-align:middle"});
		}
		
		if ( $nobutton ne 2 ){
			
			$out.= $cgi->start_Tr();
			$out.= $cgi->td({style=>"vertical-align:middle",colspan=>scalar(@infos)},[$variation->{trio}]);
			$out.= $cgi->end_Tr();
		}
				
	}
			$out.= $cgi->end_table();
			$out.= $cgi->end_div();
	}
	

	$out.="</body>";
	return $out;
}

sub printSortedVariations2 {
	my ($patient,$title,$type,$nobutton) = @_;
	die();
	 my $sorted_score = "impact_score";
	  my $sorted_score2 = "freq_score";
	my $text_score = $himpact_sorted;
	$nobutton =1 if $print;
	my $out;
	my @buttons = ("status","igv","alamut","align");
 	@buttons = () unless $edit_mode;
	my @infos = @headers;#("gene","var_name","sanger","ngs","ratio","caller","genomique","transcript","exon","nomenclature","consequence","codons","codons_AA","freq","deja_vu","similar_projects","in_this_run", "polyphen","sift");
	
	my @edit_buttons=("Validation");
	@infos = (@buttons,@infos);
	
	my $s_id=$patient->{name};
	my $total_variations = $patient->{variations};

	my $nb_line = 0;
	my $first_line = 0;
	
	#$out .= html::end_cadre($cgi);
	#$out .= html::print_cadre($cgi,"Variations ");
#	$out .=  $cgi->start_div({class=>" panel panel-default",style=>"font-size: 11px;font-family:  Verdana;"});
#	$out.= $cgi->div({class=>"panel-heading"},"Variations : sort by :".$sorted_score);
#	$out .= $cgi->start_table({class=>"table table-striped table-condensed table-bordered table-hover table-mybordered",style=>"vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;"});
#		$out.= $cgi->start_Tr();
#
#	#$out.= $cgi->th({scope=>"row"},@infos);
#	foreach my $c (@infos){
#			$out.= qq{<th scope="col">$c</th>};
#	}
#	if ($edit_mode ){
#		foreach my $c (@edit_buttons){
#			$out.= qq{<th scope="col">$c</th>};
#	}
#	#	$out .= $cgi->th({scope=>"col"},@edit_buttons);
#	}
#	$out.= $cgi->end_Tr();

my $nbv =0;
my $new_table = undef;
foreach my $tvariations ( sort {$b->[0]->{$sorted_score} <=> $a->[0]->{$sorted_score} || $b->[0]->{$sorted_score2} <=> $a->[0]->{$sorted_score2}} values %$total_variations) {	
	
		my @variations = @$tvariations;
		my $variation = $tvariations->[0];
		if ($variation->{$sorted_score} ne $new_table) {
			if (defined $new_table){
					
					$out .= $cgi->end_table();
					$out .= $cgi->end_div();
				}
			$out .=  $cgi->start_div({class=>" panel panel-default",style=>"font-size: 11px;font-family:  Verdana;"});
			my ($text)= grep {$text_score->{$_} == $variation->{$sorted_score}} keys %$text_score;
			$out.= $cgi->div({class=>"panel-heading"},"Variations : sort by :".$text);
			$out .= $cgi->start_table({class=>"table table-striped table-condensed table-bordered table-hover table-mybordered",style=>"vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;"});
			$out.= $cgi->start_Tr();

	#$out.= $cgi->th({scope=>"row"},@infos);
				foreach my $c (@infos){
						$out.= qq{<th scope="col">$c</th>};
				}
				if ($edit_mode ){
						foreach my $c (@edit_buttons){
							$out.= qq{<th scope="col">$c</th>};
							}
				}
				$out.= $cgi->end_Tr();
				$new_table = $variation->{$sorted_score} ;
		}
	
		$out.= $cgi->start_Tr({style=>"vertical-align:middle"});
		my $rowspan =scalar(@variations); 
		$nb_line+=$rowspan;
		my $cat_start;
		$nbv ++;
		for ( my $ind = 0; $ind< @infos;$ind++){
			
			my $cat = $infos[$ind];
			$variation->{caller} = "uni" unless exists $variation->{caller};
		
			$out.=print_variation_td($variation,$cat,$rowspan,$nbv,$nobutton);
			 $cat_start = $ind;
			last if $cat eq "genomique";
			
		}
		$cat_start ++;
		foreach my $variation (sort {$a->{$sorted_score} <=> $b->{$sorted_score}} @variations){
			
			for (my $ind = $cat_start;$ind< @infos;$ind++){
				my $cat = $infos[$ind];
				$out.=print_variation_td($variation,$cat,1,$nbv,$nobutton);
			}
			$uid++;
				$out .= print_validation_button($patient,$variation,$uid,$variation) if $nobutton == 2;
				$out.= $cgi->end_Tr();
				$out.= $cgi->start_Tr({style=>"vertical-align:middle"});
		}
		
}		

			$out.= $cgi->end_table();
			$out.= $cgi->end_div();

	

	$out.="</body>";
	return $out;
}
sub printTableGenesXls {
		my ($patient,$title,$type,$nobutton) = @_;
		print "Content-type: application/msexcel\n";
	print "Content-Disposition: attachment;filename=".$project->name()."-".$patient->{name}.".xls\n\n";
	
		my $workbook  = Spreadsheet::WriteExcel->new(\*STDOUT);
	my $worksheet = $workbook->add_worksheet();
	my $bg_color;
	my @colors = ("red","orange","blue","green","cyan","gray");

	foreach my $c (@colors){
		$bg_color->{$c} = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        #bg_color => 'white',
                                        color => $c,
                                        bold=>1,
                     );                                 
	}
	$bg_color->{strike} = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        #bg_color => 'white',
                                        color => "gray",
                                        bold=>1,
                                        font_strikeout=>1,
                     );    
                
		
		
		my $s_id=$patient->{name};
		my $total_variations = $patient->{variations};
		my ($z) =shift @headers;
		my @infos = ($z,"trans","phenotypes",@headers);
		if ($project->isSomatic){
		 @infos = ("gene","trans","phenotypes","var_name","sanger","ngs","ratio","caller","genomique","transcript","exon","nomenclature","consequence","codons","codons_AA","cosmic","clinvar","freq","deja_vu","similar_projects","in_this_run", "polyphen","sift","cadd");
		}
		# @infos = ("gene","trans","phenotypes","var_name","sanger","ngs","ratio","caller","genomique","transcript","exon","nomenclature","consequence","codons","codons_AA","clinvar","freq","deja_vu","similar_projects","in_this_run", "polyphen","sift","cadd");
		
		my $desc = $buffer->description_public_lmdb_database("gnomad-exome");
		push(@infos,grep {$_ ne "ALL"} @{$desc->{array}->{populations}});
		my $col = 0; 
		my $row =0;

		for (my $i=0;$i<@infos;$i++){
			$worksheet->write($row,$i,$infos[$i]);
		}
		
	$row++;
		
	my $bilan;
	foreach my $tvariations (values %$total_variations) {
				my %htemp;
		foreach my $v (@$tvariations){
				my $gene=$v->{gene};
			push(@{$htemp{$gene}},$v);
			
		}
		foreach my $gene (keys %htemp){
			push(@{$bilan->{$gene}->{variations}},$htemp{$gene});
			my %ids;
			map{$ids{$_->{id}} ++ } @{$htemp{$gene}};
			$bilan->{$gene}->{nb} += 	scalar( keys %ids);
		}
		
#		my $variation = $tvariations->[0];
#		my $gene=$variation->{gene};
#
#		push(@{$bilan->{$gene}->{variations}},$tvariations);
#		$bilan->{$gene}->{nb} += 	scalar(@$tvariations);
		foreach my $v (@$tvariations){
			my $gene=$v->{gene};
			$bilan->{$gene}->{impact}->{$v->{impact_score}} ++;
			$bilan->{$gene}->{max_impact} = $v->{impact_score} if $bilan->{$gene}->{max_impact}+0 < $v->{impact_score};
			$bilan->{$gene}->{max_freq} = $v->{freq_score} if $bilan->{$gene}->{max_freq}+0 < $v->{freq_score};
			$bilan->{$gene}->{max_score} = $v->{scaled_score} if $bilan->{$gene}->{max_score}+0 < $v->{scaled_score};
			
	
						my $max_score = 0;
						$max_score  +=  $v->{scaled_score};
						
					my %style_score;
					if (exists $style_score->{$max_score}) { %style_score = %{$style_score->{$max_score}}; }
					elsif ($max_score < 1) { %style_score = %{$style_score->{1}}; }
					elsif ($max_score < 2) { %style_score = %{$style_score->{2}}; }
					elsif ($max_score < 3) { %style_score = %{$style_score->{3}}; }
					elsif ($max_score < 4) { %style_score = %{$style_score->{4}}; }
					elsif ($max_score > 4) { %style_score = %{$style_score->{4}}; }
					
					$style_score{class} ="badge"; 
						$style_score{style} .=";font-size:xx-small;";
					$v->{status} .= qq{<div class="row"> };
			$v->{status} = "<h6><small>".$cgi->span(\%style_score,"S");
			
			my %style_impact = %{$style_score->{$v->{impact_score}}};
			$style_impact {class} ="badge"; 
			$style_impact{style} .=";font-size:xx-small;";
			$v->{status} .= $cgi->span(\%style_impact,"I");
			my %style_freq = %{$style_score->{$v->{freq_score}}};
			$style_freq{class} ="badge ";
				$style_freq{style} .=";font-size:xx-small;"; 
				$v->{status} .= $cgi->span(\%style_freq,"F")."</small></h6>";
				$v->{status} .="</div>";
			$bilan->{$gene}->{freq}->{$v->{freq_score}} ++;
			if ($v->{sanger} =~ /todo/  ){
					$bilan->{$gene}->{todo} = 1;
			
				}
			elsif ($v->{sanger} =~ /confirm/){
				$bilan->{$gene}->{saved} = 1;
			}
		
		}
		
	}

	foreach my $gene (sort{$bilan->{$b}->{max_score} <=>$bilan->{$a}->{max_score} || $bilan->{$b}->{max_impact} <=>$bilan->{$a}->{max_impact} || $bilan->{$b}->{max_freq} <=>$bilan->{$a}->{max_freq} || $a cmp $b} keys %$bilan){
	 my $pheno = "";
	 my $in;
	 my $to;
	 eval {
				 	 my $tr_id = $bilan->{$gene}->{variations}->[0]->[0]->{transcript}.'_'.$bilan->{$gene}->{variations}->[0]->[0]->{chromosome};
				 	my $ogene = $project->newTranscript($tr_id)->getGene();
				 	 $in = 	$ogene->omim_inheritance();		
				 	 $pheno =$ogene->short_phenotypes;	
				 	 ($pheno,$to) = split(/\[/,$ogene->description) unless $pheno;
	};
			
				
	foreach my $v (sort{$a->[0]->{trans} <=> $b->[0]->{trans}}@{$bilan->{$gene}->{variations}}){
		
		my @variations = sort {$a->{trans} <=> $b->{trans}} @$v;
		foreach my $variation (@variations){
		#my $variation = $variations[0];
		$variation->{trans} = $in;
		$variation->{phenotypes} = $pheno;
		$variation->{gene} = $variation->{gene}."(".$in.")" if $in;
		my $rowspan =scalar(@variations); 
		my $cat_start;
		for ( my $ind = 0; $ind< @infos;$ind++){
				$variation->{caller} = "uni" unless exists $variation->{caller};
				if ($ind == 1){
					
				}
				my $cat = $infos[$ind];
					my $text = $variation->{$cat};
					
				if (lc($cat) eq 'var_name') {
					my @lCol = split('-', $text);
					$text = $lCol[1].$lCol[3];
				}
					
				$text =~ s|<.+?>||g;
				$worksheet->write($row,$ind,$text);
		
			#last if $cat eq "genomique";
			
		}
	
		$row ++;
		}
	}

	}	
		$workbook->close();
		exit(0);
		
}

sub printTableGenes {
	my ($patient,$title,$type,$nobutton) = @_;
	my $out ="";

	my $s_id=$patient->{name};
	my $total_variations = $patient->{variations};
	my $nb_line = 0;
	my $first_line = 0;
	my @buttons = ("status","igv","alamut","align");
	my @infos = @headers; ("gene","var_name","ngs","ratio","caller","genomique","transcript","exon","nomenclature","consequence","codons","codons_AA","cosmic","freq","freq_ho","max_pop","min_pop","clinvar","hgmd","deja_vu","similar_projects","in_this_run", ,"polyphen","sift","cadd");
	my @edit_buttons=("Validation");
	@infos = (@buttons,@infos);
	my $bilan;
	foreach my $tvariations (values %$total_variations) {
				my %htemp;
		foreach my $v (@$tvariations){
				my $gene=$v->{gene};
			push(@{$htemp{$gene}},$v);
			
		}
		
		foreach my $gene (keys %htemp){
			push(@{$bilan->{$gene}->{variations}},$htemp{$gene});
			my %ids;
			map{$ids{$_->{id}} = $_ } @{$htemp{$gene}};
			$bilan->{$gene}->{nb} += 	scalar(keys %ids);
			
			my $clinvar;
			
			for my $v (values %ids){
				$bilan->{$gene}->{nb_impact_1} ++ if $v->{scaled_score} eq "1";
				$bilan->{$gene}->{nb_impact_2} ++ if $v->{scaled_score} eq "2";
				$bilan->{$gene}->{nb_impact_3} ++ if $v->{scaled_score} > 3;
				$bilan->{$gene}->{nb_clinvar} ++ if $v->{clinvar} ne "";
	
				$bilan->{$gene}->{nb_clinvar_alert} ++ if $v->{clinvar_alert} > 0;
			}
		#	map{;$bilan->{$gene}->{nb_impact_2} ++ if $_->{impact_score} eq "2";$bilan->{$gene}->{nb_impact_3} ++ if $_->{impact_score} eq "3"} values %ids;
			#map{ } values %ids;
		}
		
#		my $variation = $tvariations->[0];
#		my $gene=$variation->{gene};
#
#		push(@{$bilan->{$gene}->{variations}},$tvariations);
#		$bilan->{$gene}->{nb} += 	scalar(@$tvariations);
		foreach my $v (@$tvariations){
			my $gene=$v->{gene};
			$bilan->{$gene}->{impact}->{$v->{impact_score}} ++;
			$bilan->{$gene}->{max_impact} = $v->{impact_score} if $bilan->{$gene}->{max_impact}+0 < $v->{impact_score};
			$bilan->{$gene}->{max_freq} = $v->{freq_score} if $bilan->{$gene}->{max_freq}+0 < $v->{freq_score};
			my $scaled_score = $v->{scaled_score};
			#$scaled_score ++   if $v->{clinvar_alert} > 0;
			$bilan->{$gene}->{max_score} = $scaled_score if $bilan->{$gene}->{max_score}+0 < $scaled_score;
	
						my $max_score = 0;
						$max_score  +=  $v->{scaled_score};
						
					my %style_score;
					if (exists $style_score->{$max_score}) { %style_score = %{$style_score->{$max_score}}; }
					elsif ($max_score < 1) { %style_score = %{$style_score->{1}}; }
					elsif ($max_score < 2) { %style_score = %{$style_score->{2}}; }
					elsif ($max_score < 3) { %style_score = %{$style_score->{3}}; }
					elsif ($max_score < 4) { %style_score = %{$style_score->{4}}; }
					elsif ($max_score > 4) { %style_score = %{$style_score->{4}}; }
						
					$style_score{class} ="badge"; 
						$style_score{style} .=";font-size:xx-small;";
					$v->{status} .= qq{<div class="row"> };
			$v->{status} = "<h6><small>".$cgi->span(\%style_score,"S");
			
			my %style_impact = %{$style_score->{$v->{impact_score}}};
			$style_impact {class} ="badge"; 
			$style_impact{style} .=";font-size:xx-small;";
			$v->{status} .= $cgi->span(\%style_impact,"I");
			my %style_freq = %{$style_score->{$v->{freq_score}}};
			$style_freq{class} ="badge ";
				$style_freq{style} .=";font-size:xx-small;"; 
				$v->{status} .= $cgi->span(\%style_freq,"F")."</small></h6>";
				$v->{status} .="</div>";
			$bilan->{$gene}->{freq}->{$v->{freq_score}} ++;
			if ($v->{sanger} =~ /todo/  ){
					$bilan->{$gene}->{todo} = 1;
			
				}
			elsif ($v->{sanger} =~ /confirm/){
				$bilan->{$gene}->{saved} = 1;
			}
		
		}
		
	}



#$out .=  $cgi->start_div({class=>"panel-group", id=>"accordion",style=>"padding:2px"});
	
my $div_alert;	
$div_alert = "EDITION ";
my $all_panel=[];
my $all_label=[];
foreach my $gene ( keys %$bilan){
	push(@$all_panel,"panel_".$s_id."_".$gene);
	push(@$all_label,"label_".$s_id."_".$gene);
}
my $string_panel = join(";",@$all_panel);
my $string_label = join(";",@$all_label);
  my $icon = $cgi->span({class=>"glyphicon glyphicon-eye-open  pull-left",'aria-hidden'=>"true"});
  my $icon_help = $cgi->span({class=>"glyphicon glyphicon-question-sign pull-left",'aria-hidden'=>"true"});
  my $icon_alamut = $cgi->span({class=>"alamutView3 pull-left",'aria-hidden'=>"true"});
  my $icon_igv = $cgi->span({class=>"igvIcon2 pull-left",'aria-hidden'=>"true"});
  my $icon_calendar = $cgi->span({class=>"glyphicon glyphicon-calendar",'aria-hidden'=>"true"});
  my $icon_export =  $cgi->span({class=>"glyphicon glyphicon-open-file pull-left",'aria-hidden'=>"true"});
	my ($date,$since) = utility::get_date($patient->{obj});
	my $fam = $patient->{obj}->getFamily();
	
		my $car = $project->maskImpact;
	my $st_impact;
	foreach my $v (sort {$car->{$a} <=> $car->{$b}} keys %$car ){
		if (  $project->getMaskCoding( $himpact2->{$vimpact}) & $project->getMaskCoding( $v))
		{
			$st_impact .= $project->maskImpactTextForLegend->{$v}." - ";
		}
	}
	$st_impact = "All" if $himpact2->{$vimpact} eq "low";
	my $st_freq = "All";

	if ($hfrequence->{$vfreq} eq "unique"){
		$st_freq = "Only novel";
	}
	if ($hfrequence->{$vfreq} eq "rare"){
		$st_freq = "<=1%";
	}
	if ($hfrequence->{$vfreq} eq "occasional"){
		$st_freq = "<=5%";
	}
	my $st = " <small> <i class='fa fa-plus-square ' ></i> impact: [$st_impact] &nbsp;  <i class='fa fa-plus-square fa-1x'></i> frequence:[$st_freq]</small>&nbsp;";

	
	
	my $text = $st." <i class='fa fa-plus-square'></i> <small> ratio[ >=".$limit_ratio."%] &nbsp;<i class='fa fa-dot-circle' fa-2x></i> in this run [<=".$hthisrun->{$this_run}."%]</small>";
	if ($fam->isTrio){
#		die();
		my $label  = "label-default";
		$label  = "label-danger"  if $cgi->param('denovo');
		$text .=  qq{&nbsp;<h4 style="position:relative;top:-10px;float:right"><span class="label $label">denovo</span>};
		$label  = "label-default";
		$label  = "label-danger"   if $cgi->param('recessive');
		$text .=  qq{&nbsp;<span class="label $label">recessive</span>};
		$label  = "label-default";
		$label  = "label-danger"  if $cgi->param('xor');
		$text .=  qq{&nbsp; <span class="label $label">exclusive (M or F)</span>};
		$label  = "label-default";
		$label  = "label-danger"  if $cgi->param('both');
			$text .=  qq{&nbsp;<span class="label $label">both</span></h4>};

		
	}
	my $icon_sex =  qq{<i class="fa fa-mars" aria-hidden="true" style="color:cyan"></i>&nbsp};
	 $icon_sex =  qq{<i class="fa fa-venus" aria-hidden="true" style="color:pink"></i>&nbsp} unless $patient->{obj}->isMale();
	my $titlep = $icon_sex.$patient->{obj}->name;
	my $mname = "-";
	$mname = $fam->getMother->name() if  $fam->getMother;
	my $pname ="-";
	$pname = $fam->getFather->name() if  $fam->getFather;
	$titlep .= " &nbsp ".qq{<i class="fa fa-female" aria-hidden="true" style="color:pink"></i>&nbsp}.$mname."&nbsp ".qq{<i class="fa fa-male " aria-hidden="true" style="color:cyan"></i>&nbsp}.$pname if $fam->isTrio;
	
	my (@lbam_alamut, @lPatientsNames);
	foreach my $p (@{$fam->getPatients()}) {
		push(@lbam_alamut, $p->bamUrl());
		push(@lPatientsNames, $p->name());
	}
	my $string_url_bam = join(',', @lbam_alamut);
	my $string_url_names = join(',', @lPatientsNames);
	
	$out  .= qq{

<div class="panel panel-primary" style="padding:0px">
  <div class="panel-heading clearfix" style="min-height=30px;max-height:30px;padding-bottom:1px">
    <h3 class="panel-title pull-left" style="padding-top: 1.5px;font-size: 15px;"> &nbsp; $titlep  &nbsp; &nbsp; &nbsp; $icon_calendar  Last Update : $date ($since days ago)</h3>
   <div class="btn-group  btn-sm clearfix pull-right" style="position:relative;top:-14px;float:right">
		 		<div class=" btn btn-sm" aria-label="Left Align" style="background-color:#DEDFDE;padding-bottom:0px;padding-top:2px;top:0px;" onClick='LoadIGVPatient_editor("$string_url_names", "$string_url_bam");' >
		 			 $icon_igv 
		 		</div>
		 		
					
		 		<div class=" btn btn-danger btn-sm" aria-label="Left Align" style="padding-bottom:0px;padding-top:2px;top:0px;" onClick='httpGetLoadListBam("$string_url_bam");' >
		 			 $icon_alamut 
		 		</div>
		 		<div class=" btn btn-info btn-sm " aria-label="Left Align"  onClick='collapse_all("$string_panel","$string_label");' >
		 			$icon
		 		</div>
		 		<div class=" btn btn-info btn-sm" aria-label="Left Align" onClick='dijit.byId("dialog_help_1").show();' >
		 			$icon_help
		 		</div>
		 		<div class=" btn btn-success btn-sm " aria-label="Left Align" onClick='editor(1,2);' >
		 				<b><i class="fa fa-file-excel-o pull-left"></i></b>
		 		</div>
		</div>
		
	
  </div>
<div class="alert alert-success" role="alert" style="padding:0px">
			$text
	</div>	
</div>  

	};

	if ($dude){
		print_dude_infos();
	}
									

foreach my $gene (sort{$bilan->{$b}->{max_score} <=>$bilan->{$a}->{max_score} || $bilan->{$b}->{nb_clinvar_alert} <=>$bilan->{$a}->{nb_clinvar_alert} || $bilan->{$b}->{max_impact} <=>$bilan->{$a}->{max_impact} || $bilan->{$b}->{max_freq} <=>$bilan->{$a}->{max_freq} || $bilan->{$b}->{nb_clinvar} <=>$bilan->{$a}->{nb_clinvar} || $a cmp $b} keys %$bilan){
	my $panel_id = "panel_".$s_id."_".$gene;
	my $label_id = "label_".$s_id."_".$gene;
	my $max_score = 0;
	$max_score  += $bilan->{$gene}->{max_score};
	$out .=  $cgi->start_div({class=>"panel panel-success" });
	 #panel heading
		$out .=  $cgi->start_div({class=>"panel-heading panel-warning warning ",style=>" min-height:13px;max-height:13px;padding:1px;border:1px"});
			my $glyph = "";
				
				 $glyph = qq{<span class="glyphicon glyphicon-star-empty text-default" aria-hidden="true"></span>} if $bilan->{$gene}->{nb_clinvar} > 0;
				 $glyph = qq{<span class="glyphicon  glyphicon-alert text-alert" aria-hidden="true" style="color:red"></span>} if $bilan->{$gene}->{nb_clinvar_alert} > 0;
				 my ($ogene, $in, $pheno, $to);
				 eval {
				 	 my $tr_id = $bilan->{$gene}->{variations}->[0]->[0]->{transcript}.'_'.$bilan->{$gene}->{variations}->[0]->[0]->{chromosome};
				 	 $ogene = $project->newTranscript($tr_id)->getGene();
				 	 $in = 	$ogene->omim_inheritance();		
				 	 $pheno =$ogene->short_phenotypes;	
				 	 ($pheno,$to) = split(/\[/,$ogene->description) unless $pheno;
				 };
				$out .=  $cgi->start_div({class=>" btn-group btn-xs "});
				$out .= qq{<div class="btn  btn-info btn-xs " style="position:relative;bottom:1px;min-width:150px;" onClick='collapse("$panel_id","$label_id")'>  <span id= "$label_id" class="glyphicon glyphicon-triangle-right  " aria-hidden="true"  style="float:left;"></span> $gene &nbsp<sup><b>&nbsp;$in</b></sup>&nbsp $glyph</div>};
	   			my $nbv = $bilan->{$gene}->{nb};
	   		
	   			$out .=$cgi->span({class=>"label label-success"},qq{<span class='badge badge-primary badge-xs '  >$nbv </span> });
				#$out .=$cgi->span({class=>"label label-success"},$nbc) if $nbc >0;
				$out .=$cgi->span({class=>"label label-danger"}, ($bilan->{$gene}->{nb_impact_3}+0));
				$out .=$cgi->span({class=>"label label-warning"}, ($bilan->{$gene}->{nb_impact_2}+0));
				$out .=$cgi->span({class=>"label label-default"},($bilan->{$gene}->{nb_impact_1}+0));# if $nbc >0;
			
				
	   		$out.= $cgi->end_div();

			#div lavel right 
				$out .=  $cgi->start_div({class=>" btn-group btn  ",style=>'position:relative;float:right;bottom:5px;'});
			
					my $max_impact = 0;
					$max_impact  += $bilan->{$gene}->{max_impact};
					
					my $text = $himpact_sorted_inv->{$max_impact};
					
					
					my %style;
					
					
					if (exists $style_score->{$max_score}) { %style = %{$style_score->{$max_score}}; }
					elsif ($max_score < 1) { %style = %{$style_score->{1}}; }
					elsif ($max_score < 2) { %style = %{$style_score->{2}}; }
					elsif ($max_score < 3) { %style = %{$style_score->{3}}; }
					elsif ($max_score < 4) { %style = %{$style_score->{4}}; }
					elsif ($max_score > 4) { %style = %{$style_score->{4}}; }
					
					$style{class} ="label  "; 
					$out .=$cgi->span(\%style,"S"); 
					
					 %style =%{$style_score->{$max_impact}};
					$style{class} ="label  "; 
					$out .=$cgi->span(\%style,"I");
					my $max_freq = 0;
					$max_freq  += $bilan->{$gene}->{max_freq};
					$himpact_sorted_inv->{$max_freq};
					%style = %{$style_score->{$max_freq}};
					$style{class} ="label  "; 
					$out .=$cgi->span(\%style,"F");
					
						
				
			
					$out .=$cgi->span({class=>"label label-success"},qq{<span class="glyphicon glyphicon glyphicon-menu-hamburger " aria-hidden="true" "></span>});
					my $clabel = " label-default";
					$clabel = "label-info" if exists $bilan->{$gene}->{todo};
					$out .=$cgi->span({class=>"label $clabel"},qq{<span class="glyphicon glyphicon-pencil " aria-hidden="true" "></span>});
					$clabel = " label-default";
					$clabel = "label-info" if exists $bilan->{$gene}->{saved};
					$out .=$cgi->span({class=>"label $clabel"},qq{<span class="glyphicon glyphicon-saved" aria-hidden="true" "></span>});
			 	$out.= $cgi->end_div(); # end div lavel right 
			 
	   		
	   		my @t = split(";",$pheno);
	   		if (scalar(@t)>3){
	   			$pheno = $t[0]."-".$t[1]."-".$t[-1]." etc.".qq{<span class="glyphicon  glyphicon-option-horizontal" aria-hidden="true"></span>};
	   		}
	   		
	   		if ($pheno){
	   			my $color ;
	   			
	   			
	   			$pheno =  qq{<div class="label label-basic  " style="text-shadow:0px 1px 1px #000;font-size:1em;position:relative;bottom:8px;min-width:200px;font-size: 1em;color:white;font-weight: normal;"> <i class="fa fa-circle fa-xs" $color ></i> $pheno </div>};
	   		}
	   		
			 	
				$out.=$pheno;
			 	
			$out.= $cgi->end_div(); # end panel heading
		
		 $out.= "<br>";
	#  panel table
	$out .=  $cgi->start_div({class=>"panel-body panel-collapse  collapse",style=>"font-size: 09px;font-family:  Verdana;",id=>$panel_id});
		$out.="\n";
		$out .= $cgi->start_table({"data-toggle"=>"table","data-search"=>"true",class=>"table table-striped table-condensed table-bordered table-hover table-mybordered",style=>"vertical-align:middle;text-align: center;font-size: 7px;font-family:  Verdana;"});
		#$out .= $cgi->start_table({"data-toggle"=>"table","data-search"=>"true",style=>"vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;"});
	
		$out.= $cgi->start_Tr();
		
			foreach my $c (@infos){
						$out.= qq{<th >$c</th>};
				}
				if ($edit_mode ){
						foreach my $c (@edit_buttons){
							$out.= qq{<th >$c</th>};
							}
				}
				$out.= $cgi->end_Tr();
				
				my $cat_start =0;
				my $nbv=1;
				
	foreach my $v (sort{$a->[0]->{trans} <=> $b->[0]->{trans}}@{$bilan->{$gene}->{variations}}){
		
		my @variations = sort {$a->{trans} <=> $b->{trans}} @$v;
		my $variation = $variations[0];
		$out.= $cgi->start_Tr({class=>"menurow" ,style=>"vertical-align:middle"});
		my $rowspan =scalar(@variations); 
		$nb_line+=$rowspan;
		my $cat_start;
		$nbv ++;
		for ( my $ind = 0; $ind< @infos;$ind++){
				$variation->{caller} = "uni" unless exists $variation->{caller};
			my $cat = $infos[$ind];
			$out.=print_variation_td_edit($variation,$cat,$rowspan,$nbv,$nobutton,$patient);
			 $cat_start = $ind;
			last if $cat eq "genomique";
			
		}
		$cat_start ++;
	
			foreach my $variation (sort {$a->{trans}*1.0 <=> $b->{trans}*1.0} @variations) {
			
				for (my $ind = $cat_start;$ind< @infos;$ind++){
					my $cat = $infos[$ind];
					$out.=print_variation_td_edit($variation,$cat,1,$nbv,$nobutton,$patient);
				}
				$uid++;
				#warn Dumper $gene;
				#warn $ogene;
				#die();
				$out .= validation_select($patient,$variation,$ogene,$variation);
				$out .= print_validation_button($patient,$variation,$uid,$variation);
				$out.= $cgi->end_Tr();
			}
	}

	$out.= $cgi->end_table();	
	$out.="<!-- 1 -->";
	$out.= $cgi->end_div();		
	$out.="<!-- 2 -->";
	#  end panel table
	$out.= $cgi->end_div();	#$out.="<!-- 3 -->";	
	
	#$out.= $cgi->end_div();		
	#last;	
#  end panel gene
}

#$out.= $cgi->end_div();		

		#	$out.= $cgi->end_div();

	#	$out.= "<br>";
	$out.= q{
		<script type="text/javascript">
$(document).ready(function(){
    $('[data-toggle="tooltip"]').tooltip({
        placement : 'top'
    });
});
</script>
	};

	return $out;
}	

my $tdid =0;

#sub printTableHotspots2 {
#	my ($patient,$hotspots,$print) = @_;
##$out .=  $cgi->start_div({class=>"panel-group", id=>"accordion",style=>"padding:2px"});
#	my $out ="";
#	my $div_alert;	
#my $s_id = $patient->{name};
#
#	$out  .= qq{
#<div class="panel panel-default">
#  <div class="panel-heading clearfix" style="min-height=30px;max-height:30px;">
#    <h3 class="panel-title pull-left" style="padding-top: 1.5px;font-size: 15px;">Hotspots &nbsp; $s_id</h3> 
#  </div>
#</div>  
#	};
#
#foreach my $hotspot (@$hotspots){
#	my $panel_id = "panel_".$s_id."_".$hotspot->{sequence};
#	my $label_id = "label_".$s_id."_".$hotspot->{sequence};
#	my $text = $hotspot->{name}.":".$hotspot->{sequence};
#	$out .=  $cgi->start_div({class=>"panel panel-info" });
#	 #panel heading
#	$out .=  $cgi->start_div({class=>"panel-heading panel-warning warning ",style=>" min-height:13px;max-height:13px;padding:1px;border:1px"});
#	$out .= qq{<div class="btn  btn-success btn-xs " style="position:relative;bottom:1px;min-width:150px;" onClick='collapse("$panel_id","$label_id")'>  <span id= "$label_id" class="glyphicon glyphicon-triangle-right  "   style="float:left;"></span> $text &nbsp</div>};
#	   		#	$out .=$cgi->span({class=>"label label-success"},qq{<span class='badge badge-primary badge-xs'  >$nbv</span>});
#		my $nbv = scalar (keys %{$hotspot->{results}->{$s_id}});
#		$out .=$cgi->span({class=>"label label-success"},qq{<span class='badge badge-primary badge-xs'  >$nbv</span>});	
#				
#	   		$out.= $cgi->end_div();
#	#	$out.= $cgi->end_div();
#
#	
#			
#		 $out.= "<br>";
#		 
#	#  panel table
#	$out .=  $cgi->start_div({class=>"panel-body panel-collapse  ",style=>"font-size: 09px;font-family:  Verdana;",id=>$panel_id});
#	$out .= $cgi->start_table({class=>"table table-striped table-condensed table-bordered table-hover table-mybordered",style=>"font-size: 8px;font-family:  Verdana;"});
#	$out.= $cgi->start_Tr();
#	$out.=$cgi->th(["sequence","Forward",'Reverse','%']);
#	$out.= $cgi->end_Tr();
#	my $res = $hotspot->{results}->{$s_id};
#	
#	foreach my $motif ( sort{$res->{$b}->{pourcent} <=> $res->{$a}->{pourcent} } keys %$res){	
#	$out.= $cgi->start_Tr();
#	my @td;
#	my $pout = $cgi->start_table({class=>"table table-striped table-condensed table-bordered table-hover table-mybordered",style=>"font-size: 8px;font-family:  Verdana;"});
#	foreach my $l (@{$res->{$motif}->{table_align}}){
#		$pout.= $cgi->start_Tr();
#
#			$pout.=join("\n",@$l);	
#		$pout.= $cgi->end_Tr();
#	}
#	$pout.= $cgi->end_table();
#	
#	push(@td, $pout);
#	push(@td,$res->{$motif}->{p});
#	push(@td,$res->{$motif}->{m});
#	push(@td,$res->{$motif}->{pourcent}."%");
#	$out.=$cgi->td(\@td);	
#	$out.= $cgi->end_Tr();
#	}
#	$out.= $cgi->end_table();	
#	$out.="<!-- 1 -->";
#	$out.= $cgi->end_div();		
#	$out.="<!-- 2 -->";
#	#  end panel table
#	$out.= $cgi->end_div();	#$out.="<!-- 3 -->";	
#	
#	#$out.= $cgi->end_div();		
#	#last;	
##  end panel gene
#}
#
#
#
#
#
#
#
#	$out.= q{
#		<script type="text/javascript">
#$(document).ready(function(){
#    $('[data-toggle="tooltip"]').tooltip({
#        placement : 'top'
#    });
#});
#</script>
#	};
#	
#	return $out;
#}	
#

sub printTableHotspots {
	my ($patient,$print) = @_;
	my $hotspots = $patient->hotspot;
	return "" unless $hotspots;

	my $out ="";
	
	#$out .=  $cgi->start_div({class=>"panel-heading panel-alert alert ",style=>" min-height:13px;max-height:13px;padding:1px;border:1px"});
	my $label_id = "hs_".$patient->id;
	my $panel_id = "pa_".$patient->id;
	$out .= qq{<div class="btn  btn-warning btn-xs " style="position:relative;bottom:1px;min-width:150px;" onClick='collapse("$panel_id","$label_id")'>  <span id= "$label_id" class="glyphicon glyphicon-triangle-right  "   style="float:left;"></span> HOTSPOT &nbsp</div>};
	#$out .=$cgi->span({class=>"label label-success"},qq{<span class='badge badge-primary badge-xs'  >-</span>});
	$out .=  $cgi->start_div({class=>"panel-body panel-collapse  collapse",style=>"width:50%;font-size: 09px;font-family:  Verdana;",id=>$panel_id});
	my $div_alert;	
	my $s_id = $patient->{name};

my $t = time;


#$out .=  $cgi->start_div({class=>"panel-heading panel-warning warning ",style=>" min-height:13px;max-height:13px;padding:1px;border:1px"});
#	$out .= qq{<div class="btn  btn-success btn-xs " style="position:relative;bottom:1px;min-width:150px;" onClick='collapse("$panel_id","$label_id")'>  <span id= "$label_id" class="glyphicon glyphicon-triangle-right  "   style="float:left;"></span> $text &nbsp</div>};
	   		#	$out .=$cgi->span({class=>"label label-success"},qq{<span class='badge badge-primary badge-xs'  >$nbv</span>});
#		my $nbv = scalar (keys %{$hotspot->{results}->{$s_id}});
#		$out .=$cgi->span({class=>"label label-success"},qq{<span class='badge badge-primary badge-xs'  >$nbv</span>});	
	my @header = ("ID","NAME","PROT","A","C","G","T","DEL","COV");	 
	$out .= $cgi->start_table({class=>"table table-striped table-condensed table-bordered table-hover table-mybordered",style=>"font-size: 9px;font-family:  Verdana;"});
foreach my $g (keys %$hotspots){
	#$out .=  $cgi->start_div({class=>"panel panel-info" });
	 #panel heading
	 
	 # $out.= $cgi->end_div();
		#REF	POS	COV	A	C	G	T	DEL	REFSKIP	SAMPLE
	#my $var_obj = $self->cache_lmdb_variations->get($vid);
	#  panel table
	$out.= $cgi->start_Tr();
	$out.=$cgi->th({colspan=>(scalar(@header)+1),style=>"background-color:#217DBB;color:white;font-size:12px"},$g);
	$out.= $cgi->end_Tr();
	$out.= $cgi->start_Tr();
	$out.=$cgi->th({style=>"background:#E0E0FF"},["igv",@header]);
	$out.= $cgi->end_Tr();
	
	my @bams;
	my @names;
	foreach my $p (@{$patient->getFamily->getPatients()}){
		push(@bams,$patient->bamUrl);
		push(@names,$patient->name());
	}
					
	my $f =  join(";",@bams);#$patient->{obj}->bamUrl;;
	 my $pnames = join(";",@names);
	foreach my $hotspot (@{$hotspots->{$g}}){
	
	my @td;
		my $chr = $project->getChromosome($hotspot->{REF});
		my $var_obj = $chr->cache_lmdb_variations->get($hotspot->{GENBO_ID});
		#warn $chr->cache_lmdb_variations->get(0);
		my $style ={};
		 $style = {style=>"background-color:#D2386C;color:white"} if $var_obj && $var_obj->existsPatient($patient);
		 $out.= $cgi->start_Tr($style);
		 my $nba;
		 my $chrn = $chr->fasta_name;
		 my $start = $hotspot->{POS};
		 my $l = $chr->fasta_name.":".$start;
		 my $gn = $project->getVersion();
		 my $project_name = $project->name;
		
		my $text =qq{<button dojoType="dijit.form.Button"   iconClass='igvIcon' onclick='view_web_igv_bam("dialog_igv", "div_igv", "$l", "$f", "$pnames")' style="color:black"></button>};
		 
		 
		$out.=$cgi->td($text);
		foreach my $h (@header){
			if ($h eq  $hotspot->{A_ALT}){
				my $pc = int(($hotspot->{$h}/$hotspot->{COV})*1000)/10;
				my $color = "#f2dedc";
				$color = "#F7BFB9" if $pc>2;
				$color = "#E9897E" if $pc>5;
				$out.=$cgi->td({style=>"background-color:$color;color:black"},"$pc% (".$hotspot->{$h}.")");
			}
			elsif ($h eq  $hotspot->{A_REF}){
				my $pc = int(($hotspot->{$h}/$hotspot->{COV})*1000)/10;
				$out.=$cgi->td({style=>"background-color:#c7eadd;color:black"},"$pc% (".$hotspot->{$h}.")");
			}
			elsif ($h eq  "DEL" && $hotspot->{A_ALT} eq "-"){
				my $pc = int(($hotspot->{$h}/$hotspot->{COV})*1000)/10;
				my $color = "#f2dedc";
				$color = "#F7BFB9" if $pc>2;
				$color = "#E9897E" if $pc>5;
				$out.=$cgi->td({style=>"background-color:#F7BFB9;color:black"},"$pc% (".$hotspot->{$h}.")");
			}
			else {
				my $pc = int(($hotspot->{$h}/$hotspot->{COV})*1000)/10;
				$out .= $cgi->td($hotspot->{$h});
			}
			#push(@td, $hotspot->{$h});
		}
	
		#$out.=$cgi->td(\@td);	
		$out.= $cgi->end_Tr();
	}
	
	}
	$out.= $cgi->end_table();	
	$out.= $cgi->end_div();	#$out.="<!-- 3 -->";	
	$out.= "<BR>";
	return $out;
}





my $tdid =0;


sub validation_select{
	my ($hpatient,$variation,$gene) = @_;
	my $patient = $hpatient->{obj};
	my $cgi =  new CGI;# unless $cgi;
	my $buffer = $patient->buffer;
	my $project = $patient->project;
	my $out;
	my $pname = $patient->{name};
	my $vid = $variation->{id};
	my $tt = $patient->name."_".$gene->{js_id};
	my $menu = $tt."_menu";
	my $sp = $tt."_span";
	my $tdid = $gene->{js_id}."_".rand(50);
			
#	my @tds = ( qq{<a href="#" id="$bi_todo" class="button clean" onclick="this.className='button todo';validation('$pname','$vid',this,-3,'todo','$sp','$menu',$bi);">todo</a>},
#						qq{<a href="#"  id="$bi_ho" class="button clean" onclick="this.className='button todo';validation('$pname','$vid',this,1,'todo','$sp','$menu',$bi);">Ho</a>},
#						qq{<a href="#" id="$bi_he"  class="button clean" onclick="this.className='button todo';validation('$pname','$vid',this,2,'todo','$sp','$menu',$bi);">He</a>},
#						qq{<a href="#"  id="$bi_fp" class="button clean" onclick="this.className='button todo';validation('$pname','$vid',this,-1,'todo','$sp','$menu',$bi);">FP</a>}		
#	
my $bgcolor = "info";
my $val_id = $gene->id."!".$variation->{id};

my $saved  ;
my $all_validations = $patient->validations;
my $validation_term;
my $validation_value = 0;
 if (exists $all_validations->{$val_id}){
 #	my @found = grep {$_->{sample_id} eq $patient->id} @{$all_validations->{$val_id}};
 	$saved =  $all_validations->{$val_id};
 	$validation_term = $saved->[0]->{term};
 	$validation_value = $saved->[0]->{validation};
 }

#$saved = $all_validations->{$val_id}->[0]->{validation}  if exists $all_validations->{$val_id};

my $option;
foreach my $val (sort {$b <=> $a }keys %{$buffer->value_validation}){
	my $term = $buffer->value_validation->{$val};
	my $sel ="";
	if (lc($validation_term) eq lc($term) ){
		$sel = "selected";
	}
	$option .=  qq{<option $sel value="$val">$term</option>\n};
}
unless ($validation_term){
	$option .=  qq{<option selected value="0">-</option>\n};
}

	$option .="</select></div>";
	$bgcolor = "info";
	$bgcolor = "secondary" if  $validation_value >= 3;
	$bgcolor = "warning" if  $validation_value >= 4;
 	$bgcolor = "danger" if  $validation_value >= 5;
 	
	my $uniq_id ="$pname"."+"."_$tdid";
	my $label_id= "div_$uniq_id";
	my $select_id =  "select_$uniq_id";
	my $force_text = qq{
		</select>
		<label id ="$label_id">
		<input type="checkbox" onChange="document.getElementById('select_$uniq_id').disabled = false;"> force
		</label>
		</div>
	};
	my $disabled ="";
	if ($saved){
		$disabled ="disabled";
	}
	else {
 		$force_text=qq{</select><label id ="$label_id"></label></div>};
	}
	my $gene_id = $gene->{id};
	my $select_text = qq{
		<div  style="vertical-align:middle;padding:1px;border-radius: 5px; -moz-border-radius: 5px; -webkit-border-radius: 5px; border: 2px solid #FFFFFF;" >
		<select id="select_$uniq_id" style="padding:1px" onchange ="validation_acmg('$pname','$vid',this,'$uniq_id','$gene_id');" $disabled>
	};
	my $td_text = $select_text.$option.$force_text;
	$out .= $cgi->td({class=>$bgcolor,style=>"color:#000000;",id=>"td_$uniq_id"},$td_text); 
	$variation->{html}->{validation_select}->{$gene->{id}} = $out;
	#value_html($variation,"validation_select",$tt,$out);
	return $out;
}





sub print_validation_button{
	my ($patient,$variation,$iud,$variation) = @_;
	return;
	my $out;
	my $pname = $patient->{name};
	my $vid = $variation->{id};
	my $tt = $patient->{name}."_select_"."$uid";
	my $menu = $tt."_menu";
	my $sp = $tt."_span";
	$tdid++;
			
#	my @tds = ( qq{<a href="#" id="$bi_todo" class="button clean" onclick="this.className='button todo';validation('$pname','$vid',this,-3,'todo','$sp','$menu',$bi);">todo</a>},
#						qq{<a href="#"  id="$bi_ho" class="button clean" onclick="this.className='button todo';validation('$pname','$vid',this,1,'todo','$sp','$menu',$bi);">Ho</a>},
#						qq{<a href="#" id="$bi_he"  class="button clean" onclick="this.className='button todo';validation('$pname','$vid',this,2,'todo','$sp','$menu',$bi);">He</a>},
#						qq{<a href="#"  id="$bi_fp" class="button clean" onclick="this.className='button todo';validation('$pname','$vid',this,-1,'todo','$sp','$menu',$bi);">FP</a>}		
#	
my $bgcolor = "info";

my $option = qq{

  <option selected value="0">-</option>
  <option value="-3">todo</option>
  <option value="1">he</option>
  <option value ="2">ho</option>
  <option value ="-1" >FP</option>
</select>
  </div>

};
  
  my $sp=1;
  $bgcolor = "warning" if exists $variation->{type_confirmed_ngs};
  $sp = undef  if exists $variation->{type_confirmed_ngs};
  
  my $uniq_id ="$pname"."+"."_$tdid";
my $label_id= "div_$uniq_id";
my $select_id =  "select_$uniq_id";
 
 my $force_text = qq{
 	</select>
 	<label id ="$label_id">
      <input type="checkbox" onChange="document.getElementById('select_$uniq_id').disabled = false;"> force
    </label>
  </div>
 };
 my $disabled ="";

 if (  $variation->{sanger} ne "-" &&  $variation->{sanger} ne "todo" ){
 	$sp = undef;
 		$bgcolor = "danger";
 		my $tt = $variation->{sanger} ;
 		$tt =~ s/confirmed//;
	$option =  qq{
<option selected value ="0" >Sanger $tt</option>
  <option  value ="-3">todo</option>
  <option value = "2">he</option>
  <option value= "1">ho</option>
  <option value = "-1">FP</option>
    <option value="0">-</option>
	};
 	
 }
elsif ($variation->{type_confirmed_ngs} eq "todo"){
		$sp = undef;
#	$bgcolor = "success";
$option = qq{
  <option value="0">-</option>
  <option  selected value ="-3">todo</option>
  <option value = "2">he</option>
  <option value= "1">ho</option>
  <option value = "-1">FP</option>
	};
 
 
}
if ($variation->{type_confirmed_ngs} eq "he"){
			$sp = undef;
	$option =  qq{
	
  <option >-</option>
 <option  value ="-3">todo</option>
  <option selected  value = "2">he</option>
  <option value= "1">ho</option>
  <option value = "-1">FP</option>
    <option value="0">-</option>
	};
}
if ($variation->{type_confirmed_ngs} eq "ho"){
			$sp = undef;
	$option =qq{
	  <option >-</option>
<option  value ="-3">todo</option>
  <option value = "2">he</option>
  <option selected value= "1">ho</option>
  <option value = "-1">FP</option>
    <option value="0">-</option>
};
}
if ($variation->{type_confirmed_ngs} eq "fp"){
			$sp = undef;
	$option =qq{
 <option  value ="-3">todo</option>
  <option value = "2">he</option>
  <option value= "1">ho</option>
  <option selected value = "-1">FP</option>
    <option value="0">-</option>
};
}

unless (defined $sp){
 	 $disabled ="disabled";
 	 
 }
 else {
$force_text=qq{</select><label id ="$label_id"></label></div>};
 }
 
  my $select_text = qq{
 	<div  style="vertical-align:middle;padding:1px;border-radius: 5px; -moz-border-radius: 5px; -webkit-border-radius: 5px; border: 2px solid #FFFFFF;" >
 	<select id="select_$uniq_id" style="padding:1px" onchange ="validation_2('$pname','$vid',this,'$uniq_id');" $disabled>
	
 } ;
my $td_text = $select_text.$option.$force_text;

	$out .= $cgi->td({class=>$bgcolor,style=>"color:#000000;",id=>"td_$uniq_id"},$td_text); 
	return $out;
}




sub html_cnv{
my ($patient,$tr) = @_;	
my $out;
my $jsid=0;
my $cnv_patients;
my $id;
my $tab_ids;
	return unless scalar(@{$tr->getPrimers()});
#my $tr = $project->newTranscript($tr1);
	
 
	$out.=$cgi->start_Tr();
	#$out.= $cgi->th({scope=>"col", abbr=>"Configurations", class=>"nobg"}, "Dup/Del");
	$out.= draw_cnv::table_cnv_transcripts($cgi,$patient,$tr);
	$out.= $cgi->end_Tr();
	return $out;

}