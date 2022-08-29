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
use lib "$Bin/../packages/validation_variation"; 
use html; 

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
use preload_coverage;
$|=1;
my $buffer = GBuffer->new();

        
                
my $cgi          = new CGI();
my $out;
my $prg =  $cgi->url(-relative=>1);
my $url_image = url(-absolute=>1);
$url_image =~ s/$prg/image_coverage.pl/;
$url_image =~ s/polydiag_old/validation_variation/;

my $order = $cgi->param('order');

my $project_name = $cgi->param('project');
my $only_low;
 $only_low = 1 if $cgi->param('only_low') ==1;
my $project = $buffer->newProject(-name=>$project_name);
my $utr =$cgi->param('utr')+0;
my $intronic = $cgi->param('intronic')+0;

my $limit = $cgi->param('limit');
my $padding = $cgi->param('span');

my @transcripts_cgi;
my $gene_id =  $cgi->param('gene');
if ($gene_id) {
	my @real_gene_id;
	foreach my $gid (split(",",$gene_id)){
			my $gene = $project->newGene($gid);
			
			push(@transcripts_cgi, map{$_->id} @{$gene->getTranscripts});
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

my $capture = $project->getCaptures()->[0];
my $vquery;
if ($project->isDiagnostic){
	$vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>$project->validation_db());
}
my $transcripts_exons_todo ={};
my $transcripts_exons_validated = {};

foreach my $patient (@{$project->getPatients}){
	my $exons_todo = {};
	 $exons_todo = $vquery->get_exons(project_name=>$project_name,sample_name=>$patient->{name}) if $vquery;
	foreach my $v (values %$exons_todo){
	
		if  ($v->{done} == 0){
			$transcripts_exons_todo->{$v->{transcript}}  ++ ;
		}
		else {
				$transcripts_exons_validated->{$v->{transcript}}  ++ ;
		}
		
	}
}
my $patient_names = $cgi->param('patients');

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





#my @transcripts = sort{$a->getGene->external_name cmp $b->getGene->external_name} map{$project->newTranscript($_)} @transcripts_cgi ;
my $nbt = scalar(@transcripts);
@transcripts = splice (@transcripts ,$cgi->param('start'),$cgi->param('step')) if $cgi->param('step') && !$only_low ; 
$out.= $cgi->start_table({class=>"table table-striped table-condensed table-bordered  table-mybordered",style=>"font-size: 8px;font-family:  Verdana;"});
my $col = 12;

$col = 10 if scalar(@{$project->getPatients()}>15);
$col = 6 if scalar(@{$project->getPatients()}>30);
$col = 5 if scalar(@{$project->getPatients()}>40);

my $nb=0;

$url_image.="?project=".$project_name."&limit=$limit&span=$padding&utr=$utr&intronic=$intronic&patients=$patient_names&transcript=";
my @ths;
my@tds;
html::print_cgi_header($cgi,$out,undef,"Coverage View [$project_name]");
my $no =  $project->noSqlCoverage();

my $patients = $project->getPatients();
my $nb =0;
my @problem_transcripts;
my $list_transcripts2 = $no->get($project->name."-min","minimum-".$padding."-".$utr) ;
#warn Dumper $list_transcripts2;
#die();
#$list_transcripts2 = undef;

#my $list_transcripts = $no->get($project->name,"list-".$padding."-".$utr);

 if ($only_low) {
 	my $list_transcripts = getLowTranscripts();
 	@transcripts = grep {exists $list_transcripts->{$_->id}} @transcripts;
 }



#warn Dumper $list_transcripts;
#warn scalar(@problem_transcripts);
#warn scalar(@transcripts);
#@transcripts = @problem_transcripts;


my $nbtp = scalar(@transcripts);
$out.= html::print_cadre($cgi,"Exons Coverage : Transcripts ($nbtp/$nbt) ");
for (my $i=0;$i<@transcripts;$i++){
#	my $spanTr = $no->get($patients->[0]->name,$transcripts[$i]->id);	
#	warn Dumper  $spanTr;
if ($i%$col == 0 && $i>0){
		$out.=print_lines(\@ths,\@tds);
		@ths=();
		@tds=();
	
	}


 	 my $tname= $transcripts[$i]->name;
# 	my $class = "nobg";
# $class ="th1" if $z%2 ==0;
 my $text = $transcripts[$i]->getGene()->external_name()."<br>".$tname;
 if (exists $transcripts_exons_todo->{$tname} ){
 			$text .= "<br>";
 			$text .= " validated :".$transcripts_exons_validated->{$tname} if exists $transcripts_exons_validated->{$tname} ;
 			$text .= " todo:".$transcripts_exons_todo->{$tname} if (exists $transcripts_exons_todo->{$tname} ) ;
 
 	}
 	elsif (exists $transcripts_exons_validated->{$tname} ){
 			$text .= "<br> validated :".$transcripts_exons_validated->{$tname};
 			$text .= " todo:".$transcripts_exons_todo->{$tname} if (exists $transcripts_exons_todo->{$tname} ) ;
 
 	}
	push(@ths,$cgi->th({class=>"success"},$text));
	my $tr_id = $transcripts[$i]->id;
 	my $url2 = $url_image.$tr_id;
 	
	 my $img = qq{<div style="padding:0px;min-width:36px;min-height:36px" ><img src="$url2"  align="top" style="box-shadow: 2px 2px 3px #aaaaaa;"  ></img></div>}; 
	 push(@tds,$cgi->td({onClick=>"zoomCoverage('$tr_id')",style=>"padding:0px;background-color:#DFF0D8"},$img));
	
	
 }
 $out.=print_lines(\@ths,\@tds);
$out.=$cgi->end_table();
$out.= html::end_cadre($cgi);
print $out;
#html::print_cgi($cgi,$out);
exit(0);


sub getLowTranscripts {
	my $list_transcripts2 = $no->get($project->name,"minimum-".$padding."-".$utr) ;
	unless ($list_transcripts2){
		$list_transcripts2 = preload_coverage::computeLowTranscripts($project,\@transcripts,$no,$intronic,$utr,$padding,1);
		$no->put($project->name,"minimum-".$padding."-".$utr,$list_transcripts2) ;
	}
	my $patients =  $project->get_list_patients($patient_names,",");
my $list_transcripts;
foreach my $t (@transcripts){
	my @z;
	foreach my $p (@$patients){
		
		push(@z,$list_transcripts2->{$t->id}->{$p->{id}})
	}
	my $min = min(@z);
	#warn $min;
	next if $min >= $limit;
	$list_transcripts->{$t->id} ++; 
}
return $list_transcripts;
}

sub print_lines{
	my ($ths,$tds) = @_;
	my $out = $cgi->start_Tr().join("\n",@$ths).$cgi->end_Tr();
	 $out .= $cgi->start_Tr().join("\n",@$tds).$cgi->end_Tr();
	$ths=[];
	$tds=[];
	return $out;
}