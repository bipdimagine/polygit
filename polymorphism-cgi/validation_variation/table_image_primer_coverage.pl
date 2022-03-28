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
$| =1;

my $buffer = GBuffer->new();

my $cgi          = new CGI();

my $prg =  $cgi->url(-relative=>1);
my $url_image = url(-absolute=>1);
$url_image =~ s/$prg/image_primer_coverage.pl/;

my $out;
my $project_name = $cgi->param('project');
$out.= html::print_cadre($cgi,"Multiplex Coverage ");


my $project = $buffer->newProject(-name=>$project_name);
my $utr =$cgi->param('utr')+0;
my $intronic = $cgi->param('intronic')+0;

my $limit = $cgi->param('limit');
my $padding = $cgi->param('span');

my @transcripts_cgi;

my $captures = $project->getCaptures();
my $vquery;
if ($project->isDiagnostic){
	$vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>$project->validation_db());
}


$url_image.="?project=".$project_name."&limit=$limit&span=$padding&utr=$utr&intronic=$intronic&capture=";

#$out .= qq{<p id="p1"> TEST </p>};
$out.= $cgi->start_table({class=>"table table-striped table-condensed table-bordered  table-mybordered",style=>"font-size: 8px;font-family:  Verdana;"});
my $col = 10;
$col = 5 if scalar(@{$project->getPatients()}>40);
$col = 15 if scalar(@{$project->getPatients()}<20);
my $nb=0;
my @multi_capt;
foreach my $capture (@$captures){
	foreach my $multi (sort{$a<=>$b}@{$capture->getMultiplex()}){
		my $toto;
		$toto->{multi} = $multi;
		$toto->{capt} = $capture;
		push(@multi_capt,$toto);
	}
}
my @ths;
my@tds;

for (my $i=0;$i<@multi_capt;$i++){
#	$out .=  qq{<script type='text/javascript'> document.getElementById('p1').innerHtml= "prout".$i};
 	 	my $cname= $multi_capt[$i]->{capt}->name;
		 my $mname = $multi_capt[$i]->{multi};
 		my $text = "capture : ".$cname."<br> multiplex ".$mname;
if ($i%$col == 0 && $i>0){
		$out.=print_lines(\@ths,\@tds);
		@ths=();
		@tds=();
	
	}
	push(@ths,$cgi->th({class=>"info"},$text));
 	
 	my $url2 = $url_image.$cname."&multiplex=".$mname;
	 my $img = qq{<div><img src="$url2" style="box-shadow: 2px 2px 3px #aaaaaa;" ></img></div>}; 
	 
	 push(@tds, $cgi->td({onClick=>"zoomPrimer('$cname','$mname')",style=>"background-color:#D9EDF7"},$img) );
 }
$out.=print_lines(\@ths,\@tds);
$out.=$cgi->end_table();
$out.= html::end_cadre($cgi);
html::print_cgi($cgi,$out);
exit(0);

	



sub print_lines{
	my ($ths,$tds) = @_;
	my $out = $cgi->start_Tr().join("\n",@$ths).$cgi->end_Tr();
	 $out .= $cgi->start_Tr().join("\n",@$tds).$cgi->end_Tr();
	$ths=[];
	$tds=[];
	return $out;
}
