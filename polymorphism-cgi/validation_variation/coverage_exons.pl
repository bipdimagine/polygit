#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/kyoto";

#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
#use Set::;
use Storable qw/store thaw retrieve/;
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
use lib "$Bin/../packages/validation_variation"; 
use draw_cnv; 
use infos_coverage_exons;
use image_coverage;
use html;

	my $CSS_TABLE = <<END;
	<style type="text/css"> 
	body {
	font: normal 6px auto "Trebuchet MS", Verdana, Arial, Helvetica, sans-serif;
	color: #4f6b72;
	background: #E6EAE9;
}

a {
	color: #c75f3e;
}

table.coveragetable {

	padding: 0;
	margin: 0;
		margin:20px;
	border:#ccc 1px solid;

	-moz-border-radius:3px;
	-webkit-border-radius:3px;
	border-radius:3px;

	-moz-box-shadow: 0 1px 2px #d1d1d1;
	-webkit-box-shadow: 0 1px 2px #d1d1d1;
	box-shadow: 0 1px 2px #d1d1d1;
	
}

table.coveragetable caption {
	padding: 0 0 5px 0;
	width: 700px;	 
	font: italic 9px "Trebuchet MS", Verdana, Arial, Helvetica, sans-serif;
	text-align: right;
}

table.coveragetable th {
	padding:10px 15px 12px 15px;
	border-top:1px solid #fafafa;
	border-bottom:1px solid #e0e0e0;

	background: #ededed;
	background: -webkit-gradient(linear, left top, left bottom, from(#ededed), to(#ebebeb));
	background: -moz-linear-gradient(top,  #ededed,  #ebebeb);
}

table.coveragetable th.nobg {
	border-top: 0;
	border-left: 0;
	border-right: 1px solid #C1DAD7;
	background: none;
}

table.coveragetable td {
	border-right: 1px solid #C1DAD7;
	border-bottom: 1px solid #C1DAD7;
	background: #fff;
	padding: 6px 6px 6px 12px;
	color: #4f6b72;
	font: normal 6px auto "Trebuchet MS", Verdana, Arial, Helvetica, sans-serif;
	
}


table.coveragetable td.alt {
	background: #F5FAFA;
	color: #797268;
}
table.coveragetable td.red {
	font: bold 9px auto "Trebuchet MS", Verdana, Arial, Helvetica, sans-serif;
	color: #000000;
	
	
}
table.coveragetable td.black {
	font: bold 9px auto "Trebuchet MS", Verdana, Arial, Helvetica, sans-serif;
	color: #000000;
	background:  url(/icons/Polyicons/line_black.png) no-repeat;
	
}
table.coveragetable th.spec {
	border-left: 1px solid #C1DAD7;
	border-top: 0;
	background: #fff url(/icons/Polyicons/coin_green.png) no-repeat;
	font: bold 9px "Trebuchet MS", Verdana, Arial, Helvetica, sans-serif;
}

table.coveragetable th.specalt {
	border-left: 1px solid #C1DAD7;
	border-top: 0;
	background: #FFE1E1 url(/icons/Polyicons/coin_red.png) no-repeat;
	font: bold 9px "Trebuchet MS", Verdana, Arial, Helvetica, sans-serif;
	color: #797268;
}
table.coveragetable th.specalt1 {
	border-left: 1px solid #C1DAD7;
	border-top: 0;
	background: #FFE1E1 url(/icons/Polyicons/coin_orange.png) no-repeat;
	font: bold 9px "Trebuchet MS", Verdana, Arial, Helvetica, sans-serif;
	color: #797268;
}
table.coveragetable td:hover td {
	background: #f2f2f2;
	background: -webkit-gradient(linear, left top, left bottom, from(#ffff00), to(#ffA500));
	background: -moz-linear-gradient(top,  #ffff00,  #ffA500);	
}


p:hover {
    /* proprietes css */
 
	cursor:pointer;
	font: bold 11px;
	color: #FFA500;

}
 .main h2 {
    margin: 0em 0 0.5em 0;
    font-weight: normal;
    position: relative;
    text-shadow: 0 -1px rgba(0,0,0,0.6);
    font-size: 15px;
    line-height: 20px;
    background: #D74B4B;
    background: rgba(53,86,129, 0.8);
    border: 1px solid #fff;
    padding: 5px 15px;
    color: white;
    border-radius: 0 10px 0 10px;
    box-shadow: inset 0 0 5px rgba(53,86,129, 0.5);
    font-family: 'Muli', sans-serif;
}
</style>




END


my $nb_gene_by_patient = 3;
my $nb_exon_by_genes = 10;

my $buffer = GBuffer->new();

my $cgi          = new CGI();
#print $cgi -> header;
#print $CSS_TABLE."\n";
	html::print_header_polydiag($cgi);
my $project_name = $cgi->param('project');

my $project = $buffer->newProject(-name=>$project_name);
my $captures = $project->getCaptures();
my $capture = $captures->[0];
my $capture = $project->getCaptures()->[0];



my $vquery;
my $seq_type = $cgi->param('exome');
my $max = 50;


if ($project->isDiagnostic){
	$max = 75;
	$vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>$project->validation_db());
}

my $patient_name = $cgi->param('patients');
my $gene =  $cgi->param('gene');
my $only_red =1 if $cgi->param('only_red');
my $only_orange =1 if $cgi->param('only_orange');
$patient_name ="all" unless $patient_name;
my $patients = $project->get_list_patients($patient_name,",");


my $gene_id =  $cgi->param('gene');
my @transcripts_cgi ;
if ($gene_id) {
	my @real_gene_id;
	foreach my $gid (split(",",$gene_id)){
		my $ensg_id = $project->getEnsgIDs($gid);
		foreach my $geid (@$ensg_id){
			my $gene = $project->newGene($geid);
			push(@transcripts_cgi, map{$_->name} @{$gene->getTranscripts});
		}
	}
	
	
}
else {
my $cgi_transcript =  $cgi->param('transcripts');
if ($cgi_transcript eq "all"){
	@transcripts_cgi = @{$project->bundle_transcripts() } ;
}
else {
	@transcripts_cgi = split(",",$cgi_transcript);
}

}

@transcripts_cgi = splice (@transcripts_cgi ,$cgi->param('start'),$cgi->param('step')) if $cgi->param('step') ; 
#if ($cgi_transcript || $cgi_transcript ne "+" ){
#	$transcripts_cgi = [split(",",$cgi_transcript)];
#
#}
#else {
	#	$transcripts_cgi = $capture->transcripts_name() ;
#}
my $splice_5 = 10;
my $limit    = 0;
my $string_id ="e"+time;
$string_id ="or" if $only_red;
$string_id ="oo" if $only_orange;
$limit = $cgi->param('limit');
if ($project->isGenome() and $limit > 15) {
	$limit = 15;
}
$splice_5 = $cgi->param('span');
$splice_5 = 0 unless $splice_5;
my $exons_red;
my $common_red;
my $patient_red;
my $data; 
my $all_htranscripts;
my $patient_exons;

if ($cgi->param('cnv')){
	html_cnv();
	exit(0);
}


my $hcdata;
my %dejavu_transcripts;
$| = 1;
print qq{<div style="visibility: hidden">};
foreach my $patient (sort{$a->name <=> $b->name} @{$patients}){
	my $capture = $patient->getCapture();
	print $patient->name().",";
	my $hpatient;
	$hpatient->{name} = $patient->name();
	my $variations_sanger = {};
	my $exons_todo  = {};
	my $variations_todo = {};
	if ($vquery){
		 $variations_sanger = $vquery->get_variations_sanger(project_name=>$project_name,sample_name=>$patient->{name});
		$exons_todo = $vquery->get_exons(project_name=>$project_name,sample_name=>$patient->{name});
		$variations_todo = $vquery->get_variations_todo(project_name=>$project_name,sample_name=>$patient->{name});
	}
	my @variations =();
	push(@variations,keys %$variations_sanger) if $variations_sanger;
	push(@variations,keys %$variations_todo) if $variations_todo;

	#my $variations_ions = $vquery->get_variations_ions(project_name=>$project_name,sample_name=>$patient->{name});
	foreach my $tr (@transcripts_cgi){
	$dejavu_transcripts{$tr} =  $project->newTranscript($tr) unless exists $dejavu_transcripts{$tr};
	my $tr1 = $dejavu_transcripts{$tr} ;
	my $htranscript;
	#next if ($gene ne 'all' && $gene ne $tr1->getGene->external_name());
	my $chr = $tr1->getChromosome->name();
	my $start = $tr1->start;
	my $end = $tr1->end;
	my $text = $tr1->getGene->external_name();
	#my $ret  = image_coverage::image ($patients, $tr1,0,0, 0, 5 );
	#warn Dumper $ret->{data};
	#die();
	my $capture_intspan = $tr1->getChromosome->getIntSpanCapture();#$capture->genomic_span($tr1->getChromosome);
	$htranscript->{name} = qq{<a href="javascript:displayInIGV('$chr',$start,$end);">$text</a><br><small>$tr</small>};
	$htranscript->{vname} = $tr1->name();
	#$htranscript->{name} = $tr1->getGene->external_name();
	$htranscript->{gene_name} = $tr1->getGene->external_name();
	$htranscript->{lname} = $tr1->getGene->external_name()."  :".$tr1->external_name;
	$htranscript->{id} = $tr1->id;
	#$htranscript->{mean} = $tr1->mean_coding_coverage($patient);
	 #$htranscript->{mean_exonic}  = $tr1->mean_exonic_coverage($patient);
	$htranscript->{exons} = [];
	$htranscript->{variations} = [];
	my $exons;
	if ($cgi->param('intronic') == 1){
	 $exons = $tr1->getAllGenomicsParts();
	}
	else {
		 $exons = $tr1->getExons();
	}
		 my $show_utr = $cgi->param('utr')+0;
		 my $intronic = $cgi->param('intronic')+0;

			unless (exists $hcdata->{$tr} ){
				 my $res;
				 
					
					
				if ($project->isNoSqlDepth){	
					
				 $res  = image_coverage::image_depth_lmdb ($patients, $tr1,$intronic,$show_utr, $splice_5, $limit,1);
				
			
 	 		}
 	 		else {
 	 			die();
 	 			 $res  = image_coverage::image($patients, $tr1,$intronic,$show_utr, $splice_5, $limit,1);
 	 		}
					$hcdata->{$tr} = $res->{data};
					
				 
			}
			
		my $cdata= $hcdata->{$tr} ;
	#my $ret = image_coverage::image ($patients, $tr1,$intronic,$show_utr, $splice_5, $limit,1);
	
#	warn Dumper($data);
	foreach my $exon (sort{$a->start*$a->strand <=> $b->end*$b->strand }@$exons){
#		warn $exon->name();
			#my ($exon,$patient,$tr1,$exons_todo,$data,$limit) = @_;
			my $pdata = $cdata->{$exon->id}->{$patient->id};
			#warn Dumper $pdata;
			my $hexons;
		 	$hexons = infos_coverage_exons::return_hash_exons2 ($exon,$patient,$tr1,$exons_todo,$pdata,$limit);

		my $variation_inside = 0;
		#warn $patient->name;
		foreach my $v (@variations ){
			#warn $v;
			my ($chrv,$startv,$reste) = split("_",$v);
			$chrv =~s/chr//;
			next if $chrv ne $chr;
			
			$variation_inside= 1 if $startv>= $exon->{start} && $startv <= $exon->{end};
		}	

		
		push(@{$htranscript->{exons}},$hexons);
		$patient_exons->{$patient->name}->{$hexons->{id}} = $hexons;
		
		#warn $exon->name." ".$mean." ".$intspan->as_string()." ".$exon->getGenomicSpan->as_string;
	}
	$all_htranscripts->{$htranscript->{id}} = $htranscript;
	push(@{$hpatient->{transcripts}},$htranscript);
	#$tr1->purge_patient($patient);
	#$all_hpatients->{$htranscript->{id}} = $htranscript;
	}

	push(@$data,$hpatient);
}
print qq{</div>};
unless ( $cgi->param('xls')){
	html2();
	exit(0);
}


sub html_cnv{

foreach my $patient (@{$patients}){
		print qq{<section class="main">};
		print $cgi->h2($patient->name());
		print qq{</section>};
		print $cgi->start_table({class=>"table table-striped table-bordered  table-condensed "});
	
	foreach my $tr1 (@transcripts_cgi){
	my $tr = $project->newTranscript($tr1);
	print draw_cnv::table_cnv_transcripts($cgi,$patient,$tr);
	}#end tr
	print $cgi->end_table();
	print "<hr>";
}#end patient

#	my $string_connect = join(",",@$tab_ids);
#	
#	print qq{
#	<div id="tooltip_cnv" data-dojo-type="dijit/Tooltip" data-dojo-props="connectId:[$string_connect],position:['above']">
#     
#	</div>
#	
#	};

}

sub html_cnv_all{


my $xr = 255;
my $xg = 10;
my $xb = 0;

my $yr = 0;
my $yg = 255;
my $yb = 0;
my $jsid=0;



foreach my $tr1 (@transcripts_cgi){
	
	my $tr = $project->newTranscript($tr1);
	
	print $cgi->start_table({class=>"coveragetable"});
	print $cgi->start_Tr();
	print $cgi->th({scope=>"col", abbr=>"Configurations", class=>"nobg"}, "Gene ".$tr->getGene()->external_name);
	my @write_patients;
	foreach my $patient (@{$patients}){
		my $pn = $patient->{name};
		print $cgi->th({scope=>"col"},$patient->name());
		push(@write_patients,$patient);
	}
	print $cgi->end_Tr();
	print $cgi->start_Tr();
	#foreach my $exon (sort{$a->start <=> $b->start } @{$tr->getExons()} ){
		
		#next unless exists $exon->{problem};
	foreach my $primer (sort{$b->getChromosome->length <=> $a->getChromosome->length || $a->start <=> $b->start } @{$tr->getPrimers()} ){
		#next if exists $primer->{notok};
		print $cgi->start_Tr();
		my $text="";
		my $gtext;
		foreach my $e (@{$primer->getExons}){
			$text.= $e->name()." ";
			
		}
		push(@$gtext, $text." Multiplex :".$primer->multiplex);
		push(@$gtext,$primer->getChromosome->ucsc_name.":".$primer->start."-".$primer->end);
		
		
		print $cgi->th({scope=>"row", abbr=>"Model", class=>"spec"}, join("<br>",@$gtext));
		foreach my $patient (@$patients){
			my $pn = $patient->{name};
			my $exon = $primer->{$patient->name}->{exon};
			my $value = "";#$patient->{name}.":".$exon->{label};
			my $mean = int($primer->mean($patient));
			my $min = $primer->min($patient);
			my $text = "<small>".$primer->mean($patient)."</small>";
			my $val = $min>$max ? $max : $min;
			$val -= $limit;
			#$val = 100 - $val;
			my $red = int(($xr + (($val * ($yr-$xr)) / $max) ));
			my $green = int(($xg + (($val * ($yg-$xg)) / $max) ));
			my $blue = 50;#int(($xb + (( $val* ($yb-$xb)) / $max) ));
			my $label =$exon; 
			my $text = $mean."(".$min.")";
		
			my $en = $exon;
			
			my $mean = int($primer->statistics()->mean()*100)/100;
			#my $text = $exon->{mean}."<small>(".$exon->{min}.")<br>[".$exon->{intspan}."]</small>";
			$jsid++;
			my $check;
			$text = $cgi->p($text);
				#$text = qq{<div class="openGraph">$text</div>} ;
				my $pn = $patient->name();
				my $en = $primer->name();
				my $click = qq{load_graph_primer('$pn','$en',$en->start,$en->end};
				my $plex = $primer->multiplex;
				my $score = $primer->cnv_score($patient);
				
				my $sdp= int($primer->statistics()->standard_deviation()*100)/100;
				my $sdm = int($primer->statistics_multiplex($patient)->standard_deviation()*100)/100;;
				my $meanm = int($primer->statistics_multiplex($patient)->mean()*100)/100;
				my $zscore = $primer->zscore->{$patient->id};
				my $itext = $score." - ".$sdp." -  $mean <br> multi: $sdm -  $meanm <br> zscore  ".$zscore;
			
				if (!$primer->is_multiplex_ok($patient)){
			#	if($sdm+$meanm>1.5 || $meanm-$sdm <0.5){
					print $cgi->td({onClick=>$click, class=>"red",style=>"background-color:rgb(90,90,90);font-size:8px;"},$itext);
				}
				elsif (!$primer->is_primer_ok) {
			#	elsif($sdp+$mean>1.5 || $mean-$sdp < 0.5){
					print $cgi->td({onClick=>$click, class=>"red",style=>"background-color:rgb(50,50,250);font-size:8px;"},$itext);
				}
			#	elsif ($score - $sdp >= 1.4){
				elsif ($zscore>= 2.5){
					my $color =int (200*($score-0.5))+50;
					print $cgi->td({onClick=>$click, class=>"red",style=>"background-color:rgb($color,$color,0);font-size:8px;"},$itext );
				}
				#elsif ($score + $sdp  <= 0.6){
				elsif($zscore<= -2.5){
					
					my $color =int (50*(0.6-$score))+205;
					print $cgi->td({onClick=>$click, class=>"red",style=>"background-color:rgb($color,0,0);font-size:8px;"},$itext);
				}
				else {
					
						print $cgi->td({onClick=>$click, class=>"red",style=>"background-color:rgb(220,220,220);font-size:8px;"},$itext );
				}
			
	
			
		}
		print $cgi->end_Tr();
	}
	print $cgi->end_table();
	print "<HR>";
}
}

sub computeCoverage {
	my ($p,$chr,$start,$end) = @_;
	my $tabix = $p->tabix_coverage();
	my $res = $tabix->query($chr,$start,$end);
	my $min = 9999999;
	my $sum =0;
	my $nb =0;
	while(my $line = $tabix->read($res)){  
		my($a,$b,$c) = split(" ",$line);
		$min = $c if $c < $min;
		$sum += $c;
		$nb ++;
		
	}
	return (int($sum/$nb),$min);
}



sub html2{

$| =1;



my $xr = 255;
my $xg = 10;
my $xb = 0;

my $yr = 0;
my $yg = 255;
my $yb = 0;
my $id =0;
my $jsid=$cgi->param('start')*1000+1;

foreach my $tr (sort {$a->{gene_name} cmp $b->{gene_name}} values %$all_htranscripts) {
	my $tn = $tr->{id};

	#print $cgi->start_table({class=>"coveragetable"});
		print $cgi->start_table({class=>"table table-striped table-bordered  table-condensed ",style=>"text-align: center;vertical-align:middle;font-size: 09px;font-family:  Verdana;"});
	print $cgi->start_Tr({class=>"info"});
	my $click3 = qq{load_graph_gene('$tn');};
	print $cgi->th({onClick=>$click3,scope=>"col", abbr=>"Configurations",style=>"background-color:#2C3E50;color:white"}, $tr->{name});
		
		
	my @write_patients;
	foreach my $patient (@{$patients}){
		my $pn = $patient->{name};
		if ($only_orange){
			next unless exists $patient_red->{$tr->{id}}->{$patient->name};
			next if $patient_red->{$tr->{id}}->{$patient->name} == scalar (keys %{$common_red->{$tr->{id} } });
	#		warn $pn." ".$patient_red->{$tr->{id}}->{$patient->name}." ".scalar (keys %{$common_red->{$tr->{id} } } );
		}
		next if $only_orange && !(exists $patient_red->{$tr->{id}}->{$patient->name});
		my $click2 = qq{load_graph_transcript('$pn','$tn');};
	
		print $cgi->th({onClick=>$click2,scope=>"col",style=>"background-color:#014B7D;color:white"},$patient->name());
		push(@write_patients,$patient);
	}
	print $cgi->end_Tr();
	print $cgi->start_Tr(style=>"text-align: center;");
	
	foreach my $exon1 (sort{$a->{id} <=> $b->{id} }@{$tr->{exons}}){
		next if $only_orange && !(exists $exons_red->{$tr->{id}}->{$exon1->{id}});
		print $cgi->start_Tr();
		my $bubble_name = "tooltip_".  $exon1->{id};
		 my @tab_ids;
		my $en = $exon1->{vname};
		my $click4 = qq{load_graph_one_exon_all_patients('$tn','$en');};
		next if $only_red && !(exists $common_red->{$tr->{id}}->{$exon1->{id}} );
		if  (exists $common_red->{$tr->{id}}->{$exon1->{id}}){
			print $cgi->th({class=>"danger",onClick=>$click4,scope=>"row", abbr=>"Model",}, $exon1->{vname});
		}
		elsif (exists $exons_red->{$tr->{id}}->{$exon1->{id}}){
			print $cgi->th({class=>"warning",onClick=>$click4,scope=>"row", abbr=>"Model"}, $exon1->{vname});
		}
		else {
			print $cgi->th({onClick=>$click4,scope=>"row", abbr=>"Model", style=>"background-color:#014B7D;color:white;font-size:11px"}, $exon1->{vname});
		}
		foreach my $patient (@write_patients){
			my $pn = $patient->{name};
			next if $only_orange && !(exists $patient_red->{$tr->{id}}->{$patient->name});
			my $exon = $patient_exons->{$patient->name}->{$exon1->{id}};
			 
			my $value = "";#$patient->{name}.":".$exon->{label};
			my $text = "<small>".$exon->{mean}."</small>".$exon->{intspan};
			my $val = $exon->{mean}>$max ? $max : $exon->{mean} ;
			$val -= $limit;
			#$val = 100 - $val;
			my $red = int(($xr + (($val * ($yr-$xr)) / $max) ));
			my $green = int(($xg + (($val * ($yg-$xg)) / $max) ));
			my $blue = 50;#int(($xb + (( $val* ($yb-$xb)) / $max) ));
			my $value = $exon->{label}->{$patient->{name}} ;#$patient->{name}.":".$exon->{label};
			my $label =$exon->{name}; 
			my $text = $exon->{mean}."(".$exon->{min}.")"; 
			if (exists $exon->{raw_mean}){
				$text .= "[".$exon->{raw_mean}."(".$exon->{raw_min}.")]"; 
			}
		
			my $en = $exon->{vname}; 
			
			#my $text = $exon->{mean}."<small>(".$exon->{min}.")<br>[".$exon->{intspan}."]</small>";
			$jsid++;
			my ($vv) = values  %{$exon->{label}};
			my @av = split(":",$vv);
			my $ee = pop @av;
			my $es =pop @av;
			
 			my $click = qq{load_graph_exon('$pn','$tn','$en',$es,$ee,'raw');};
			my $check;
			my $text_bubble = "$pn ".$tr->{vname}." ".$en;
			$text = $cgi->p({class=>"openGraph",title=>"$text_bubble",onClick=>$click},$text);
				#$text = qq{<div class="openGraph">$text</div>} ;
			if($vquery){
			 $check = qq{
					<input id="$string_id$jsid" name ="$string_id$jsid" data-dojo-type="dijit/form/CheckBox" value="$value"  onChange="test(this)" /> 
					<label for="mycheck" class="label">$text</label>
				};
				$text = $check if $exon->{red} ==1;
			}
			else {
				 $check = $text;
			}
		
			
			$id++;
			push(@tab_ids,"cov_$id");
			my ($r,$g,$b) = @{$exon->{color}};
			print $cgi->td({class=>"red",style=>"background-color:rgb($r,$g,$b);font-size:09px;",id=>"cov_$id"},$text);
#			if ($exon->{type} == 1){
#				print $cgi->td({class=>"red",style=>"background-color:rgb(50,50,200);font-size:8px;",id=>"cov_$id"},$text);
#			}
#			elsif ($exon->{type} == 0){
#				print $cgi->td({class=>"red",style=>"background-color:rgb(150,50,150);font-size:8px;",id=>"cov_$id"},$text);
#			}
#			elsif ($exon->{type} == -1 ||  $exon->{type} == -2){
#				print $cgi->td({class=>"red",style=>"background-color:rgb(250,250,0);font-size:8px;",id=>"cov_$id"},$text);
#			}
#			elsif ($exon->{color} eq "violet"){
#				print $cgi->td({class=>"red",style=>"background-color:rgb(150,50,150);font-size:8px;",	id=>"cov_$id"},$text);
#			}
#			elsif ($exon->{type} == 2){
#				print $cgi->td({class=>"black",style=>"background-color:rgb(230,20,20);font-size:8px;",id=>"cov_$id"},	$check);
#			}
#			elsif ($exon->{type} == -3){
#				print $cgi->td({class=>"black",style=>"background-color:rgb(200,200,200);font-size:8px;",id=>"cov_$id"},	"-");
#			}
#			else {		
#					print $cgi->td({class=>"red",style=>"background-color:rgb($red,$green,$blue);font-size:8px;",id=>"cov_$id"},$text);
#			}

			
		}
		print $cgi->end_Tr();
	#	warn $bubble_name." ".scalar(@tab_ids);
	#	my $string_connect = join(",",@tab_ids);
	#	my $string_connect = $tab_ids[0];
#		print qq{
#		<div id="$bubble_name" data-dojo-type="dijit/Tooltip" data-dojo-props="connectId:[$string_connect],position:['above']">
#		</div>
#	};
	}
	print $cgi->end_table();
	print "<HR>";

}

}


xls();	
exit(0);

sub xls{

	print "Content-type: application/msexcel\n";
	print "Content-Disposition: attachment;filename=".$project->name().".xls\n\n";
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
                     
	my $row_patient = 1;
	my $col = 0;
	my $row =0;
	$worksheet->write($row,$col,"Coverage :".$limit." span : ".$splice_5);
	$row++;
foreach my $patient (@$data){
	$col =0;
	$worksheet->write($row,$col,$patient->{name});
	#$row++;
	$col++;
	my $start_row = $row;
	my $row_next_patient;
	my $bold_merge;
	$bold_merge->[0]      = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        );
	foreach my $tr (@{$patient->{transcripts}}){
		
		#$worksheet->write($start_row,$col,$tr->{name});
		$worksheet->merge_range($start_row,$col,$start_row,$col+$nb_exon_by_genes -1,$tr->{lname},$bold_merge->[0]);
		my $row_exon = $start_row + 1;
	
		
		my $nb_exons = 0;
		my $start_col = $col;
		foreach my $exon (sort{$a->{id} <=> $b->{id} }@{$tr->{exons}}){
			my $color = $exon->{color};
			$nb_exons++;
			$color ="cyan" if $color eq "violet";
			$worksheet->write($row_exon,$col++,$exon->{name},$bg_color->{$color});
			if ($nb_exons % $nb_exon_by_genes == 0) {
				$row_exon ++;
				$col = $start_col;
				$row_next_patient= $row_exon if $row_next_patient < $row_exon;
			}
		}
		
		$col = $start_col + $nb_exon_by_genes +1;
	}
	
	$row = $row_next_patient +1;
}
            
                   
} 
