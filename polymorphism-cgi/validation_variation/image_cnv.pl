#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../GenBo/lib/kyoto";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";
use lib "$Bin/../packages/validation_variation"; 
use html; 

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
use Storable qw/store thaw retrieve freeze/;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use coverage;
use Spreadsheet::WriteExcel;
use POSIX;
use validationQuery;
use Date::Tiny;
use lib "$Bin/../packages/validation_variation"; 
use draw_cnv; 
use infos_coverage_exons;
use GD;
#use KyotoCabinet;
use image_coverage;

my $buffer = GBuffer->new();
my $cgi          = new CGI();
binmode STDOUT;

#all parameter from cgi 
# project project name
# Transcript transcirpt name
# UTR  (1,0)?
# intronic ?
# limit ?
# span ?
#patients= patient name  "," or directly all 
# le runId si plusieurs run dans un projet (on a 2 images differentes)
# llarg= (0,1) large image means the table , if not you only want the png image
my $project_name = $cgi->param('project');
my $transcript_name = $cgi->param('transcript');
$transcript_name = $cgi->param('transcripts') unless $transcript_name;
my $run_id =  $cgi->param('run_id');
my $utr =$cgi->param('utr')+0;
my $intronic = $cgi->param('intronic')+0;
my $limit = $cgi->param('limit');
my $padding = $cgi->param('span');

my $patient_names =  $cgi->param('patients');

my $large =  $cgi->param('large');
# project constructor
my $project = $buffer->newProjectCache(-name=>$project_name);
my $panel = $cgi->param('panel');
$project->setPanel("$panel") if $panel;
if (-e $project->noSqlCnvsDir){
		bless $project , 'GenBoProjectCache';
}

#my $project = $transcript->project;

my $patients;

#get patient list , patient in the param patient AND  patient in the runID

if ($run_id && $patient_names eq "all"){
my $runs = $project->getRuns();
foreach my $run (@$runs){
my $infos = $run->getAllPatientsInfos();
next if $run->id ne $run_id;
my @p = map{$_->{patient}} grep {$_->{project} eq $project_name} @$infos;
$patient_names = join(",",@p);
last;

}
}

#specific method from GenBo return only the patients object with a list of name

my $patients =  $project->get_list_patients($patient_names,",");


#construct the object trancript from the transript name

my $tr1 = $project->newTranscript($transcript_name);


my $image;

#ok the first one is not useful , this one is the real call , this method is for eliminated in the project object all the reference of the patients object not listed in the args
# after this call if you call $project->getPatients() return only the patient lists , so it's fastest if you don't have a complete list

my $patients =  $project->get_only_list_patients($patient_names,",");

	my $kyoto_id = "cnv_".$tr1->name;			
#	yes nothing to say exoect it's stupid
my $tr1 = $project->newTranscript($transcript_name);
my $image;
my $ret;
#method in the image coverage package located in GenBo/obj-nodb/lib the heart of all the process , in this method i get all primer for the transcripts and get information of the cnv from the cache 
# when I have all the data I draw the image , and also return data of primer * patient value of the cnv
#this method return hash table with two keys 
# {image} and {data}
#structure of {data} :
#{$primer->id}->{$patient->name}->{color} directly code the color of the case 
# {$primer->id}->{$patient->name}->{cnv_score} the cnv_score ;-) 
#my $ret;
#if (-e $project->noSqlCnvsDir){
	$ret = image_coverage::image_cnv_no_cache($patients,$tr1,$cgi);
 #$ret  = image_coverage::image_cnv ($patients,$tr1 );
#}
#else {
#	$ret = image_coverage::image_cache_cnv ($patients,$tr1 );
#}

$ret->{image_png} =  $ret->{image}->png;



#here it's if you want the image
unless ($large){
print $cgi->header("image/png");
print $ret->{image_png};
exit(0);
#finish for the image
}


my $out;
$out.= html::print_cadre($cgi,$tr1->getGene->external_name());
$out.= $cgi->start_table({class=>"table table-striped table-condensed table-bordered  table-mybordered",style=>"font-size: 8px;font-family:  Verdana;width:600px"});


#strucure


my @tds;
#$tr1->getChromosome->setPrimersForPatient($patients);

#GenBO methods get all primer object for this transcripts
my $primers = $tr1->getPrimers();
my $primers2 = $tr1->getChromosome()->getPrimers;
my @ths;
$out.=$cgi->start_Tr();
	$out.=$cgi->th({class=>'info'},"");

	$out.=$cgi->th({class=>'info'},"exons");
	$out.=$cgi->th({class=>'info'},"multi");
	
	
#the column list of all the patients select previously in the $patient array;
	
foreach my $patient (sort{$a->name cmp $b->name} @$patients){
	my $sex;
	if ($patient->sex == 1){
		$sex.= qq{<i class="fa fa-mars fa-2x" > </i>};
	}	
	elsif ($patient->sex == 2){
		$sex.= qq{  <i class="fa fa-venus fa-2x"  > };
	}	
	$out.=$cgi->th({class=>'info'},$patient->name()."&nbsp;$sex")
	}
$out.=$cgi->end_Tr();
my @tds;

#first loop for all primers one line/primer in the table 
foreach my $primer (sort{$a->end*$a->strand <=> $b->end*$b->strand } @$primers) {
$out.=$cgi->start_Tr();
#die();
my $exons = $primer->getExons();
my $et;
#my $et = join(";",map{$_->name} grep{$_->getTranscript()->name eq $tr1->name()}@$exons);
# col of the table , the primer ,  exon on this primer ($et) , the multiplex
$out.=$cgi->td({class=>'primary'},$primer->name());
$out.=$cgi->td({class=>'primary'},$et);
$out.=$cgi->td({class=>'primary'},$primer->multiplex);
#and after value and color for each patient

foreach my $patient (sort{$a->name cmp $b->name} @$patients){
	my $colors = $ret->{data}->{$primer->id}->{$patient->name}->{colors};
	my $st_color = join(",",@$colors);
	my $score;

		my @lTmp = split('\+\+', $ret->{data}->{$primer->id}->{$patient->name}->{cnv_score});
		$score = $lTmp[0];
		
		$out.=$cgi->td({style=>"background-color:rgb($st_color)"},$score);
	}
			
			$out.=$cgi->end_Tr();
}#endprimer

	$out.=$cgi->end_table();
$out.= html::end_cadre($cgi,"CNV");
html::print_cgi($cgi,$out);
exit(0);

#that's all 
