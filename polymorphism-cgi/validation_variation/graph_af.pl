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
use coverage;
my $nb_gene_by_patient = 3;
my $nb_exon_by_genes = 10;

my $buffer = GBuffer->new();

my $cgi          = new CGI();
my $project_name = $cgi->param('project');
#$project_name="NGS2021_3667";
my $project = $buffer->newProjectCache(-name=>$project_name);
my $chr_name = $cgi->param('chromosome');
 #$chr_name = "9";
my $chr = $project->getChromosome($chr_name);
my $chr_fasta_name = $chr->fasta_name;
my $res;
my $patient_name =  $cgi->param('patient');
my $patient;
if($patient_name){
 $patient = $project->getPatient("$patient_name");
}
else {
	my $fam = $project->getFamily();
	$patient = $fam->getChild();
	
}
my $start =  $cgi->param('start');
my $end =  $cgi->param('end');

#my($chr_name,$other) =  split(":",$region_name);
#my($start,$end) = split("-",$other);
if ($start){
$end += 1_000_000 ;
$start -= 1_000_000 ;
}
else {
	$start= 1;
	$end = $chr->length();
}

#$end = $chr->length - 1 ;

my $vector;
my $nb = $chr->getVariantsVector->Size();
my $no = $chr->cache_lmdb_variations();#->get_varid($vid);

warn $start;
#$start -= 1_000_000;
my $vector_pos = $chr->getVectorByPosition($start,$end);


my $mother = $patient->getFamily()->getMother();
my $father = $patient->getFamily()->getFather();

my $v1 = $mother->getVectorOriginHo($chr);
 $v1 -= $father->getVectorOrigin($chr);
$v1 &= $patient->getVectorOrigin($chr);
$v1 &=  $vector_pos;




my $v2 = $father->getVectorOriginHo($chr);
 $v2 -= $mother->getVectorOrigin($chr);
$v2 &= $patient->getVectorOrigin($chr);
$v2 &=  $vector_pos;



my $list2 = to_array($v2,$chr->name);





my $v3 = $v1+$v2;
my $list3 = to_array($v3,$chr->name);
$project->setListVariants($list3);
my $array;
my $x;
while(my $v = $project->nextVariant){
	$x++;
	#next if $x%3 > 0;
	my $v1 = $v->getPourcentAllele($mother);
	my $v2 = $v->getPourcentAllele($father);
	my $value = $v->getPourcentAllele($patient);
	next if $value eq "-";
	if ($v1 ne "-"){
		$array .="[".$v->start.",$value,0],";
	}
	
	$v1= "" if $v1 eq "-";
	
	
	if ($v2 ne "-"){
		$array .="[".$v->start.",0,$value]\n,";
	}
	next if  $v2 eq "-";
	$v2= "" if $v2 eq "-";
	
}

$| =1;
# die();  
print html::print_cgi_header($cgi);
my $html = qq{
<html>
  <head>
    <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>
    <script type="text/javascript">
      google.charts.load('current', {'packages':['corechart']});
      google.charts.setOnLoadCallback(drawChart);

      function drawChart() {
        var data = google.visualization.arrayToDataTable([
          ['Position', 'Mother','Father'],
			$array
			
        ]);

        var options = {
          title: 'Variant Allele Frequency Chromosome $chr_name:$start-$end',
          hAxis: {title: 'Position', minValue: $start, maxValue: $end, format: 'decimal'},
          vAxis: {title: 'VAF ', minValue: 0, maxValue: 100},
       
          explorer: { actions: ['dragToZoom', 'rightClickToReset'],maxZoomIn: 10 },
          legend: 'none',
          pointSize: 2,
        };

        var chart = new google.visualization.ScatterChart(document.getElementById('chart_div'));

        chart.draw(data, options);
      }
    </script>
  </head>
  <body>
    <div id="chart_div" style="width: 2000px; height: 500px;"></div>
  </body>
</html>
};
print $html;


sub to_array {
	my ( $v, $name ) = @_;
	my $set  = Set::IntSpan::Fast::XS->new( $v->to_Enum );
	my $iter = $set->iterate_runs();
	my @t;
	while ( my ( $from, $to ) = $iter->() ) {
		for my $member ( $from .. $to ) {
			if ($name) {
				push( @t, $name . "!" . $member );
			}
			else {
				push( @t, $member );
			}
		}
	}
	return \@t;
}
