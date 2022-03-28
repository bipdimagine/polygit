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
use JSON::XS;

my $nb_gene_by_patient = 3;
my $nb_exon_by_genes = 10;

my $buffer = GBuffer->new();

my $cgi          = new CGI();
my $project_name = $cgi->param('project');

my $project = $buffer->newProject(-name=>$project_name);
my $capture = $project->getCapture();
my $similar = $project->similarProjects();

my $vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>$capture->validation_db());

my $patient_name = $cgi->param('patients');
$patient_name ="all" unless $patient_name;
my $patients = $project->get_list_patients($patient_name,",");;

my $cgi_transcript =  $cgi->param('transcripts');
my @transcripts_cgi ;
if ($cgi_transcript eq "all"){
	@transcripts_cgi = @{$capture->transcripts_name()} ;
}
else {
	@transcripts_cgi = split(",",$cgi_transcript);
}


my $table_id = "hor-minimalist-b";
my $data = construct_data();

my $out = html($data,$cgi);
print_cgi($out);
exit(1);

sub construct_data {
	my $data; 

foreach my $patient (@{$patients}){
	my $test = $vquery->get_report(project=>$project->name,sample=>$patient->name);
	#$test=undef;
	if ($test){
		my $hpatient = decode_json($test->{json});
		$hpatient->{date_validated} =  $test->{creation_date};
		$hpatient->{user_validated} =  $test->{user_name};
		$hpatient->{conclusion} = $test->{conclusion};
		if ($cgi->param('save')) {
		my $conclusion = $cgi->param('textarea2')."";
		my $user_name = $cgi->param('user_name');
		my @toto = ($hpatient);
		my $json = encode_json $hpatient;
		$vquery->save_report(project=>$project->name,sample=>$patient->name(),json=>$json,user_name=>$user_name,conclusion=>$conclusion);
	}
		push(@$data,$hpatient);
		next;
	}
	my $hpatient;
	$hpatient->{name} = $patient->name();
	my $variations_sanger = $vquery->get_variations_sanger(project_name=>$project_name,sample_name=>$patient->{name});
	
	my $exons_todo = $vquery->get_exons(project_name=>$project_name,sample_name=>$patient->{name});
	my $variations_todo = $vquery->get_variations_todo(project_name=>$project_name,sample_name=>$patients->[0]->{name});
	
	my $variations_ions = $vquery->get_variations_ions(project_name=>$project_name,sample_name=>$patient->{name});
	foreach my $tr (@transcripts_cgi){
		
	my $tr1 = $project->newTranscript($tr);
	my $htranscript;
	$htranscript->{name} = $tr1->getGene->external_name();
	$htranscript->{mean} = $tr1->mean_coding_coverage($patient);
	$htranscript->{mean} = $tr1->mean_exonic_coverage($patient) if $htranscript->{mean} == 0;
	$htranscript->{exons} = [];
	$htranscript->{variations} = [];
	my $exons = $tr1->getExons();
	foreach my $exon (sort{$a->start*$a->strand <=> $b->end*$b->strand }@$exons){
		
		my ($mean,$intspan,$min) = $exon->mean_intspan_coverage_coding(patient=>$patient,padding=>20,limit=>10);
		#next if (scalar($intspan->as_array)<=1);
		my $hexons;
		$hexons->{name} = $exon->name();
		$hexons->{mean} = $mean;
		$hexons->{min} = $min;
		if (exists $exons_todo->{$exon->coverage_id}) {
			$hexons->{type_string} = "Sanger ";
	
			if ($exons_todo->{$exon->coverage_id}->{done} eq 1) {
				$hexons->{type} = 0;
				$hexons->{type_string} = " OK " ;
			} 
			else {
				$hexons->{type} = 1;
				$hexons->{type_string} = "to do" ;
			}
			
		}
		else {
			$hexons->{type} = 2;
			$hexons->{type_string} .= " - " ;
		}
		push(@{$htranscript->{exons}},$hexons);
		#warn $exon->name." ".$mean." ".$intspan->as_string()." ".$exon->getGenomicSpan->as_string;
	}
	foreach my $v (@{$tr1->getGene()->getStructuralVariations}){
		#die();
		my $prot = $tr1->getProtein();
		my $hvariation;
		next unless exists $v->annex()->{$patient->id};
		
	
		if (exists  $variations_ions->{$v->vcf_id}){
		
			if ($variations_ions->{$v->id}->{validation} == -1){
				
				$hvariation->{type} = "rejected";
			}
			elsif ($variations_ions->{$v->vcf_id}->{validation} == -3){
				$hvariation->{type} = "todo";
			}
			elsif ($variations_ions->{$v->vcf_id}->{validation} > 0){
				$hvariation->{type} = "validated";
			}
			$hvariation->{sanger} = "-";
			$hvariation->{type_confirmed} = "ion";
		}
		elsif (exists  $variations_sanger->{$v->vcf_id}){
			if ($variations_sanger->{$v->vcf_id}->{validation_sanger} < 0){
				$hvariation->{type} = "rejected";
				
				$hvariation->{sanger} = "rejected"
			}
			else{
				$hvariation->{type} = "confirmed";
				$hvariation->{sanger} = "confirmed";
			}
		}
		
		else {
			$hvariation->{type} = "other";
		}
#				die() if $debug;
		my $deja_vu_this_project = $v->deja_vu->{$v->project->name()};
		delete $v->deja_vu->{$v->project->name()};
		my $deja_vu_similar;
		$hvariation->{sim_deja_vu} = 0; 

		my $nb_pat = 0;
		while( my ($pname, $hp) = each(%{$v->deja_vu})) {
			$nb_pat += scalar(keys %{$hp});
			next unless exists $similar->{$pname};
			$hvariation->{sim_deja_vu} += scalar(keys %{$hp});
    		# do something with $key and $value
    		
			}
		$hvariation->{this_deja_vu} = scalar(@{$v->getPatients});
		$hvariation->{impact} = $v->effectImpact($tr1);
		$hvariation->{gene} = $tr1->getGene->external_name();
		$hvariation->{var_name} = $v->name();
		$hvariation->{deja_vu} = scalar(keys %{$v->deja_vu}).":".$nb_pat; 
		$hvariation->{in_this_project} = scalar(@{$v->getPatients}) ."/". scalar(@{$project->getPatients});
		my $sequence_info = "he("; 
		$sequence_info = "ho(" if $v->annex()->{$patient->id}->{ho};
		$sequence_info .= $v->annex()->{$patient->id}->{nb_all_ref}."/".$v->annex()->{$patient->id}->{nb_all_mut}.")";
		$hvariation->{ion} = $sequence_info;
		$hvariation->{genomique} = $v->getChromosome()->name.":".$v->start;
		$hvariation->{start} = $v->start;
		$hvariation->{chromosome} = $v->getChromosome()->name;
		$hvariation->{trans} = $tr1->translate_position($v->start);
		$hvariation->{cds} = "";
		$hvariation->{prot} ="";
		$hvariation->{codons_AA} = "";
		$hvariation->{polyphen} ="-";
		$hvariation->{sift} ="-";
		$hvariation->{codons}  =  $v->getChromosome()->sequence($v->start,$v->start)."/".$v->sequence();
		
		if ($prot){
		
			$hvariation->{cds} = $v->getOrfPosition($prot);
			$hvariation->{prot} = $v->getProteinPosition($prot);
			$hvariation->{codons_AA} =   $v->getProteinAA($prot).$hvariation->{prot}.$v->changeAA($prot);
			$hvariation->{polyphen} = $v->polyphenScore($tr1->getProtein);
			$hvariation->{sift} = $v->siftScore($tr1->getProtein);
			$hvariation->{codons} =   $v->getCodons($tr1);
		}
		
		#warn $hvariation->{prot};
		$hvariation->{exon} = $tr1->findExonNumber($v->start);
		$hvariation->{exon} = $tr1->findNearestExon($v->start) if $hvariation->{exon} == -1;
		$hvariation->{nomenclature} =  $v->getNomenclature($tr1);
		$hvariation->{consequence} =  $v->variationType($tr1);
		
		$hvariation->{freq}  =  $v->db_by_freq->{freq};
		$hvariation->{freq}  =  0 if $hvariation->{freq} == -999;
		$hvariation->{freq} = $hvariation->{freq}/100;
		
	
		push(@{$htranscript->{$hvariation->{type}}},$hvariation);
		
		
	}

	
	push(@{$hpatient->{transcripts}},$htranscript);
	}
	if ($cgi->param('save')) {
	my $conclusion = $cgi->param('textarea2')."";
	my $user_name = $cgi->param('user_name');
	my @toto = ($hpatient);
	my $json = encode_json $hpatient;
	$vquery->save_report(project=>$project->name,sample=>$patient->name(),json=>$json,user_name=>$user_name,conclusion=>$conclusion);
	}
	push(@$data,$hpatient);
	}
	return $data;
}

sub print_cgi {
	my ($out) = @_;
	print $cgi -> header;
	print $out;
	exit(0);
}
sub html {	
my ($data,$cgi) = @_;
	
my $out;
my $nb_case;

my $w = "500px";
my $CSS = <<END;
<style type="text/css"> 
body
{
	line-height: 1.0em;
}
a:link { 
    text-decoration: none;
}
a:visited { 
    text-decoration: none;
}
a:active { 
    text-decoration: none;
}
#test
{
	font-family: "Lucida Sans Unicode", "Lucida Grande", Sans-Serif;
	font-size: 12px;
	background: #fff;
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

#hor-minimalist-b .tdred
{
	border-bottom: 1px solid #ccc;
	color: #FF0000;
	padding: 3px 4px;
	width : $w;
	
}

#hor-minimalist-b .tdorange2
{
	border-bottom: 1px solid #ccc;
	color: #BB8200;
	padding: 3px 4px;
	width : $w;	
}

#hor-minimalist-b .tdorange
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



h1  {
    font-family: Georgia, "Times New Roman", Times, serif;
        font-size:13px;
	margin-top: 5px; margin-bottom: 0px;
		text-transform: uppercase;
        font-weight: normal;
        color: #222;
        letter-spacing: 0.2em;
	
 }
 h5  {
    font-family: Georgia, "Times New Roman", Times, serif;
        font-size:13px;
	margin-top: 5px; margin-bottom: 0px;
		text-transform: uppercase;
        font-weight: bold;
        color: red;
        letter-spacing: 0.2em;
	
 }
h2  {
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
 h3  {
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
 h4{
 	page-break-before: always
 }
 
</style>

END
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

$out.= $cgi->start_div({id=>"all_html",name=>"all_html",jsId=>"all_html"});
foreach my $patient (@$data){
	$out.= $cgi->h1("Patient  :"."Identifiant :". $patient->{name}."&nbsp name: <input type ='text' size ='60'>");
	my $today = Date::Tiny->now;
	
	my @genes_name = map{$_->{name} } @{$patient->{transcripts}};
	$out.= $cgi->h1( " Genes : ".join(";",@genes_name)); 
	$out.= $cgi->h1( " Sequenceur : "." ION TORRENT "." Technique :" .$project->getCapture()->type." ".$project->getCapture()->description); 	
	$out.= $cgi->h1( " date : ".$today->day."/".$today->month."/".$today->year); 
	
	$out.= $cgi->h1("Reference :" .$project->name());
	my ($d1,$h1) = split(" ",$project->creation_date) ;
	my $pdate = Date::Tiny->from_string($d1);
	$out.= $cgi->h1("Serie :" .$project->description() ." date :".$pdate->day."/".$pdate->month."/".$pdate->year);
	if (exists $patient->{user_validated}) {
		my ($d2,$h2) = split(" ",$patient->{date_validated}) ;
		my $pdate2 = Date::Tiny->from_string($d2);
	$out.= $cgi->h5("Validation :" .$patient->{user_validated} ." date :".$pdate2->day."/".$pdate2->month."/".$pdate2->year);
	} 
	$out.= "<hr>";
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
	
	$out.= $cgi->start_table({id=>"test"});
	$out.= $cgi->start_Tr();
	my $nb_tr=0;
	foreach my $transcript (@temp_transcript){
		
		my $name = "<center><h2>".$transcript->{name}." <br> Cov = ".$transcript->{mean}."</h2></center>";
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
			$out.= $cgi->start_Tr();
			$out.= $cgi->td($exon->{name});
			$out.= $cgi->td($exon->{mean});
			$out.= $cgi->td($exon->{min});
			$out.= $cgi->td($exon->{type_string});
			$out.= $cgi->end_Tr();
		}
		$out.= "</tbody>";
	$out.= $cgi->end_table();
	$out.= $cgi->end_td();
	}
	$out.= $cgi->end_Tr();
	$out.= $cgi->end_table();
	}

	$out.= printVariations($patient,"Sanger Confirmed  Variations","confirmed",$out);
	$out.= printVariations($patient,"Validated Variations","validated",$out);
	$out.=printVariations($patient,"Other Variations","other",$out);
	$out.= $cgi->h4("---");
	$out.=printVariations($patient,"Rejected Variations","rejected");
	$out.=printVariations($patient,"todo Variations ","todo");
	
	

my $pname = $patient->{name};
my $prname = $project->name();
$out.= $cgi->h1("Comment :");
#my $sttr = $cgi->param("transcripts");
my $user_name = $cgi->param("user_name");
my $conclusion = $patient->{conclusion};
$out.= qq{<form id = "myForm">};
$out.= qq{
	<input type="hidden"  name="save"  value="1">
	<input type="hidden"  name="patients"  value="$pname">	
	<input type="hidden"  name="user_name"  value="$user_name">	
	<input type="hidden"  name="project"  value="$prname">
	<textarea id="textarea2" name="textarea2" data-dojo-type="dijit/form/SimpleTextarea" rows="4" cols="50" style="width:auto;">$conclusion</textarea>
};
$out.= qq{</form>};
$out.= $cgi->h4("***");
}




$out.= $cgi->end_div();
# print $cgi->textarea(-name=>"conclusion",
# 					  #-id =>"textarea",
# 					  -rows=>10,
# 					   -columns=>80 );
 return $out;
exit(0);
	
}

sub printVariations {
	my ($patient,$title,$type) = @_;
	my $out;
	my @infos = ("var_name","sanger","ion","genomique","exon","nomenclature","consequence","codons","codons_AA","freq","deja_vu","in_this_project", "polyphen","sift");
	
	$out.= $cgi->start_table({id=>"hor-minimalist-b"});
	$out.= "<thead>";
	$out.= $cgi->start_Tr();
	my $t = scalar(@infos);
	$out.= qq{<th scope="col" colspan=$t><h2>$title<br></h2></th>};
	$out.= $cgi->end_Tr();
	
	$out.= "</thead>";
	$out.= "<tbody>";
	foreach my $transcript (@{$patient->{transcripts}}){
		#next unless $transcript->{$type};
		my $name = $transcript->{name};
		$out.= $cgi->start_Tr();
		my $nb =0;
			next unless $transcript->{$type};
		$nb += scalar  @{$transcript->{$type}} if $transcript->{$type}  ;
		$out.= qq{<th scope="col" colspan=$t><h3>$name - $nb<h3></th>};
		$out.= $cgi->end_Tr();
		$out.= $cgi->start_Tr();

	#$out.= $cgi->th({scope=>"row"},@infos);
	foreach my $c (@infos){
			$out.= qq{<th scope="col">$c</th>};
	}
	
	$out.= $cgi->end_Tr();
	
		$transcript->{$type} = [] unless exists $transcript->{$type};
		
		foreach my $variation (sort {$a->{trans} <=> $b->{trans}} @{$transcript->{$type}}){
		
			$out.= $cgi->start_Tr();
			
			foreach my $i (@infos){
			
				my $text = $variation->{$i};
				if($text =~ /,/){
					$text =~ s/,/<BR>/;
				}
				elsif (length($text)>30){
					warn $text;
					my $a = substr $text,0,15;
					my $b = substr $text,-8;
					$text = $a."...".$b;
				
					
				}
				if ($i eq "var_name"){
					my $chr = $variation->{chromosome};
					my $start = $variation->{start};
						$text = qq{<a href="javascript:displayInIGV('$chr',$start,$start);">$text</a>};
				}
				if ($variation->{impact} eq "high" ){
					$out.= $cgi->td({class=>"tdred",style=>'width:50px'},$text);
				}
				elsif  ($variation->{impact} eq "moderate" ){
					$out.= $cgi->td({class=>"tdorange",style=>'width:50px'},$text);
				}
				elsif  ($variation->{impact} eq "low" ){
					$out.= $cgi->td({class=>"tdorange2",style=>'width:50px'},$text);
				}
				else {
					$out.= $cgi->td({style=>'width:50px'},$text);
				}
			}
			
			$out.= $cgi->end_Tr();
			
		}
			
	}
	$out.= "</tbody>";
	
	$out.= $cgi->end_table();
	return $out;
}
#
#	
	