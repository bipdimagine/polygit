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
my $project = $buffer->newProject(-name=>$project_name);
my $patient = $cgi->param('patient');
my $chr_name = $cgi->param('chromosome');
 #$chr_name = "9";
my $chr = $project->getChromosome($chr_name);
my $chr_fasta_name = $chr->fasta_name;
my $res;
my $patient_name =  $cgi->param('patient');
foreach my  $patient (@{$project->getPatients}){
	warn $patient->name;
	next if $patient->name ne $patient_name;
	my $dir = $project->getVariationsDir("wisecondor");
	my $file = $dir."/".$patient->name."_bins.bed.gz";
	open (BED,"/software/bin/tabix $file $chr_fasta_name | cut -f 2,5 | grep -v NaN |" );
	while (<BED>){
		my($a,$b) =split(" ",$_);
		chomp($b);
		
		push(@{$res->{$a}},$b);
	}
	close BED;
	 last;
}
my $array;
my $last ;
my @t = sort {$a <=> $b} keys %$res;

   foreach my $c (@t){  
   	my $b = $res->{$c}->[0] ;
      $array .= " [$c ,$b]," ;
      $last = $c;
    #  warn $c." ".$b;
   }
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
          ['Age', 'Weight'],
			$array
			
        ]);

        var options = {
          title: 'WiseCondor Raw Value chromosome: $chr_name',
          hAxis: {title: 'Position', minValue: 1, maxValue: $last, format: 'decimal'},
          vAxis: {title: 'WiseCondor', minValue: -1.5, maxValue: 1.5},
          explorer: { actions: ['dragToZoom', 'rightClickToReset'],maxZoomIn: 10 },
          legend: 'none'
        };

        var chart = new google.visualization.ScatterChart(document.getElementById('chart_div'));

        chart.draw(data, options);
      }
    </script>
  </head>
  <body>
    <div id="chart_div" style="width: 1200px; height: 500px;"></div>
  </body>
</html>
};
print $html;
