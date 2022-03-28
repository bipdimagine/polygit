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
use image_primer_coverage;

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
print $cgi -> header;
print $CSS_TABLE."\n";
my $project_name = $cgi->param('project');

my $project = $buffer->newProject(-name=>$project_name);
my $capture_name = $cgi->param('capture');

my $capture = $project->getCapture($capture_name);

my $multiplex_name = $cgi->param('multiplex');



my $patient_name = $cgi->param('patients');
my $patients = $project->get_list_patients($patient_name,",");


my $splice_5 = 10;
my $limit    = 0;
$limit = $cgi->param('limit');
$splice_5 = $cgi->param('span');
$splice_5 = 0 unless $splice_5;
	html_amplicons();

sub html_amplicons{


my $xr = 255;
my $xg = 10;
my $xb = 0;

my $yr = 0;
my $yg = 255;
my $yb = 0;

my $jsid=0;


foreach my $multi (sort {$a <=> $b} @{$capture->getMultiplex}){
		#my $ret  = image_primer_coverage::image($project->getPatients, $capture,$multiplex_name, $limit );
	#	my $ret  = image_coverage::image_cnv ($patients, $tr1 );
		#warn Dumper($ret);
		#die();
		warn $multi." :: ".$multiplex_name;
	next if $multi ne $multiplex_name  ; 
	print $cgi->start_table({class=>"coveragetable"});
	print $cgi->start_Tr();
	print $cgi->th({scope=>"col", abbr=>"Configurations", class=>"nobg"}, "multiplex ".$multi);
	my @write_patients;
	foreach my $patient (@{$patients}){
		my $pn = $patient->{name};
		print $cgi->th({scope=>"col"},$patient->name());
		push(@write_patients,$patient);
	}
	print $cgi->end_Tr();
	print $cgi->start_Tr();
	my $ret  = image_primer_coverage::image($project->getPatients, $capture,$multiplex_name, $limit );
	foreach my $primer (sort{$b->getChromosome->length <=> $a->getChromosome->length || $a->start <=> $b->start } @{$capture->getPrimersByMultiplex($multi)} ){
		print $cgi->start_Tr();
		my $text;
		my $gtext;
		my $gene_name;
		my $exon_name;
		foreach my $e (@{$primer->getExons}){
			push(@$gtext, $e->getTranscript->getGene->external_name."-".$e->getTranscript()->name."-".$e->name);
			$gene_name = $e->getTranscript->getGene->external_name;
			$exon_name = $e->name;
		}
		push(@$gtext,$primer->getChromosome->ucsc_name.":".$primer->start."-".$primer->end);
		
		
		print $cgi->th({scope=>"row", abbr=>"Model", class=>"spec"}, join("<br>",@$gtext));
		
		foreach my $patient (@$patients){
					
			my $pn = $patient->{name};
			my $stat = $ret->{data}->{$primer->id}->{$patient->id};	
			my $exon = $primer->{$patient->name}->{exon};
			my $value = "";#$patient->{name}.":".$exon->{label};
			my $mean = int($stat->{mean});
			my $min = $stat->{min};
			my $text = "<small>".$mean."</small>";
	
			my $label =$exon; 
			my $text = $mean."(".$min.")";
		
			my $en = $exon;
			
			$jsid++;
			my $check;
			my $title =$pn."-multi:".$multi." gene:".$gene_name." exon:".$exon_name;
			my ($red,$green,$blue) = @{$stat->{color}};
			$text = $cgi->p($text);
			my $pn = $patient->name();
			my $en = $primer->id();
			my $click = qq{load_graph_primer('$pn','$en');};
			print $cgi->td({class=>"red",style=>"background-color:rgb($red,$green,$blue);font-size:8px;"},$text);
	
			
		}
		print $cgi->end_Tr();
	}
	print $cgi->end_table();
	print "<HR>";
}
}

