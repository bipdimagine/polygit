#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/kyoto";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";
use lib "$Bin/../packages/validation_variation"; 
use html; 

#use Set::;
use Storable qw/store thaw retrieve freeze/;
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
use preload_coverage;

$| =1;
my $buffer = GBuffer->new();

my $cgi          = new CGI();
html::print_cgi_header($cgi);


my $prg =  $cgi->url(-relative=>1);
my $url_image = url(-absolute=>1);
$url_image =~ s/$prg/image_cnv.pl/;
my $url_xls = url(-absolute=>1);
 $url_xls  =~ s/$prg/table_xls_cnv.pl/;

my $project_name = $cgi->param('project');
my $order = $cgi->param('order');
$url_xls .="?project=$project_name";

my $project = $buffer->newProject(-name=>$project_name);
my $out;
my $utr =$cgi->param('utr')+0;
my $intronic = $cgi->param('intronic')+0;

my $limit = $cgi->param('limit');
my $padding = $cgi->param('span');

my $panel_name =  $cgi->param('panel');

my @transcripts_cgi;
if ($panel_name and $panel_name ne 'all'){
 	my $panel = $project->getPanel("$panel_name");
 	@transcripts_cgi = @{$panel->main_transcripts() } ;
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

my $capture = $project->getCaptures()->[0];
my $vquery;

my $res;
my $buffer1 = GBuffer->new();
my $project1 = $buffer->newProject(-name=>$project_name);
my $runs = $project1->getRuns();
foreach my $run (@$runs){
my $infos = $run->getAllPatientsInfos();
my @p = map{$_->{patient}} grep {$_->{project} eq $project_name} @$infos;

$res->{$run->id} = \@p;
#warn Dumper $res;
#die();

}
$project1 = undef;
$buffer1 = undef;


my @transcripts_by_name = sort{$a->getGene->external_name cmp $b->getGene->external_name} map{$project->newTranscript($_)} @transcripts_cgi ;


my @transcripts_by_ordered_chrom ;
foreach my $chrom (1..22, "X", "Y", "MT"){
	my @list_tr_chrom ;
	foreach my $tr (@transcripts_cgi){
		if ($project->newTranscript($tr)->getChromosome()->name() eq $chrom) {
			push(@list_tr_chrom, $tr )
			}
	}

	unless (scalar(@list_tr_chrom) eq 0) {
		my @list_tr_chrom_trie = sort{$a->getChromosome()->name() <=> $b->getChromosome()->name() || $a->start <=> $b->start} map{$project->newTranscript($_)} @list_tr_chrom ;
		foreach my $trans (@list_tr_chrom_trie) {
			push(@transcripts_by_ordered_chrom,  $trans) ;	
	}
	} 
	
} 


my @transcripts ;
if ($order eq "name") {
		@transcripts = @transcripts_by_name;
}
else {@transcripts =@transcripts_by_ordered_chrom ;}
	



my $patient_names = $cgi->param('patients');

#$url_image.="?project=".$project_name."&limit=$limit&span=$padding&utr=$utr&intronic=$intronic&patients=$patient_names&transcript=";
foreach my $r (keys %$res){
	
my $buffer2 = GBuffer->new();
my $project2 = $buffer->newProject(-name=>$project_name);
#my $patient_names = join(",",@{$res->{$r}});
#$url_image.="?project=".$project_name."&limit=$limit&span=$padding&utr=$utr&intronic=$intronic&patients=$patient_names&transcript=";
my @transcripts = sort{$a->getGene->external_name cmp $b->getGene->external_name} map{$project2->newTranscript($_)} @transcripts_cgi ;

my $patients =  $project2->get_only_list_patients($patient_names,",");
print qq{<div style="visibility: hidden">};
preload_coverage::load_cnv_score($project2,$patients,\@transcripts,1);
print "</div>";
$project2 = undef;
$buffer2 = undef;
}

my $z = qq{<a class="btn btn-danger btn-sm" href="$url_xls" role="button" target="_blank" style="float:right" ><b><i class="fa fa-file-excel-o pull-left"></i></b></a>};
$out.= html::print_cadre($cgi,"CNV $z");
$out.=$cgi->start_table({class=>"table table-condensed table-bordered table-mybordered",style=>"font-size: 08px;font-family:  Verdana;"});
#$out.= $cgi->start_div({class=>"container"});
#$out.= $cgi->start_div({class=>"row"});

my $col = 12;
$col = 10 if scalar(@{$project->getPatients()}>20);
$col = 7 if scalar(@{$project->getPatients()}>30);
$col = 5 if scalar(@{$project->getPatients()}>40);
my $nb=0;
"&transcript=";
	my @ths;
	my @tds;
#for (my $i=0;$i<@transcripts_by_ordered_chrom;$i++){
for (my $i=0;$i<@transcripts;$i++){
	if ($i%$col == 0 && $i>0){
		$out.=print_lines(\@ths,\@tds);
		@ths=();
		@tds=();
	
	}
	my $tname = 	$transcripts[$i]->name;
 	my $gene_name = $transcripts[$i]->getGene->external_name;
 	my $chr_name= $transcripts[$i]->getChromosome->name;
 	my $tr_id = $transcripts[$i]->id;
 	#last unless $transcripts[$z]; 	
 	my $td = $cgi->start_table({class=>"table table-condensed table-bordered table-mybordered",style=>"font-size: 08px;font-family:  Verdana;"});
 	
 	foreach my $r (keys %$res){
 		
 			#my $patient_names = join(",",@{$res->{$r}});
 			my $url_image_tmp = $url_image."?project=".$project_name."&limit=$limit&span=$padding&utr=$utr&intronic=$intronic&patients=$patient_names&run_id=$r&transcript=";
 	my $url2 = $url_image_tmp.$transcripts[$i]->id();


 
 		$td .= $cgi->td({style=>'background-color:#FCF8E3;background-image:../../images/polyicons/layer_his_add.png',onClick=>"zoomCnv('$tr_id','$r')"},qq{<img class="load" src="$url2"  align="top" style="box-shadow: 2px 2px 3px #aaaaaa;" lowsrc="../../images/polyicons/layer_his_add.png"></img>}."<br>$r");
  #push(@tds,$cgi->td({style=>'background-color:#FCF8E3;background-image:../../images/polyicons/layer_his_add.png',onClick=>"zoomCnv('$tr_id')"},qq{<img class="load" src="$url2"  align="top" style="box-shadow: 2px 2px 3px #aaaaaa;" lowsrc="../../images/polyicons/layer_his_add.png"></img>}));
 	
  #	 push(@ths,$cgi->th({style=>'background-color:#F5F5F5;'},qq{<div class="caption" ><h3>$tr_name $gene_name</h3></div>}));
 	}
 	$td .=$cgi->end_table();
 		  push(@tds,$cgi->td($td));
 		  
 	 	 push(@ths,$cgi->th({class=>'warning'},qq{ $gene_name &nbsp; chr:$chr_name<br>$tname}));
	
 }
 	$out .=print_lines(\@ths,\@tds);
 	$out.=$cgi->end_table();
$out.= html::end_cadre($cgi,"CNV");
print $out;
#html::print_cgi($cgi,$out);
exit(0);


sub print_lines{
	my ($ths,$tds) = @_;
	my $out = $cgi->start_Tr().join("\n",@$ths).$cgi->end_Tr();
	 $out .= $cgi->start_Tr().join("\n",@$tds).$cgi->end_Tr();
	$ths=[];
	$tds=[];
	return $out;
}