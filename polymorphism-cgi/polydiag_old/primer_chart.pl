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
use List::MoreUtils qw/pairwise/;

my $nb_gene_by_patient = 3;
my $nb_exon_by_genes = 10;

my $buffer = GBuffer->new();

my $cgi          = new CGI();
my $project_name = $cgi->param('project');

my $project = $buffer->newProject(-name=>$project_name);
my $runs = $project->getRuns();



my $patient_name = $cgi->param('patient');
my $cgi_transcript =  $cgi->param('transcript');
my $cgi_primer = $cgi->param('primer');

$patient_name ="all" unless $patient_name;
my $patients = $project->get_list_patients($patient_name," ");


my $captures = $project->getCaptures();


my $patient_name = $cgi->param('patients');
my $cgi_transcript =  $cgi->param('transcripts');

my $vp =  $project->getPatient($patient_name);
my $capture = $project->getCapture();

my $run = $vp->getRun();

my $splice_5 = 0;
my $limit    = 0;

$limit = $cgi->param('limit');
$splice_5 = $cgi->param('span');
$splice_5 = 1 unless $splice_5;
my $exons_red;
my $common_red;
my $data; 
my $all_htranscripts;
my $patient_exons;

my $cgi    = new CGI();
my @res;
my $seq;
my $pdata;
foreach my $patient (@{$run->getPatients}){
	foreach my $primer (@{$capture->getPrimers}){
		#warn $primer->id();
		next if $primer->id() ne $cgi_primer;
		unless($seq){
				my $chr = $primer->getChromosome();
			   $seq = $chr->getSequence($primer->start-$splice_5,$primer->end+$splice_5);
		}
		my @data = $primer->getStatisticsCoverage($patient)->get_data();
		if ($patient->name eq $patient_name){
			$pdata = \@data;
		}
		else {
		 @res = pairwise{$a +$b} @res,@data;
		}
	}
}

  my @seq = split("",$seq);
   my $tab;
my $nb = scalar(@{$run->getPatients});
$nb --;
for (my $i = 0;$i<@res;$i++){
	$res[$i] = int($res[$i]/$nb);
   	my $tt;
   	$tt->{seq} = $seq[$i];
   	$tt->{cov} = $pdata->[$i];

   	 $tt->{covm} = $res[$i];
   	push(@$tab,$tt);
   	
   }
 export_data::print($project,$cgi,$tab);

exit(0);




