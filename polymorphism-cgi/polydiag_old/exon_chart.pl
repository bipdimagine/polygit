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
my $cgi_exon = $cgi->param('exon');

$patient_name ="all" unless $patient_name;
my $patients = $project->get_list_patients($patient_name," ");


my $captures = $project->getCaptures();


my $patient_name = $cgi->param('patients');
my $cgi_transcript =  $cgi->param('transcripts');

my $vp =  $project->getPatient($patient_name);
my $run = $vp->getRun();

my $splice_5 = 50;
my $limit    = 0;

$limit = $cgi->param('limit');
$splice_5 = $cgi->param('span');
$splice_5 = 1 unless $splice_5;
#$splice_5= 10;
my $exons_red;
my $common_red;
my $data; 
my $all_htranscripts;
my $patient_exons;
warn $cgi_transcript;
my $tr = $project->newTranscript($cgi_transcript);
my $cgi    = new CGI();
my @res;
my $seq;
my $pdata;
my $estart;
my $eend;
my $is_exon;
foreach my $patient (@{$run->getPatients}){
	foreach my $exon (@{$tr->getAllGenomicsParts}){
		next if $exon->name() ne $cgi_exon;
	
		$is_exon=1 if $exon->isExon();
		$splice_5 =0  unless ($is_exon);
		unless($seq){
				my $chr = $exon->getChromosome();
			   $seq = $chr->getSequence($exon->start-$splice_5,$exon->end+$splice_5);
		}
		my $res = $exon->getTranscript()->get_coverage($patient)->coverage($exon->start-$splice_5,$exon->end+$splice_5);
		$data = $res->{array};
		
	
	$eend = abs($exon->end-$exon->start)+$splice_5;
	if ($patient->name eq $patient_name){
			$pdata = $data;
		}
		else {
		 @res = pairwise{$a +$b} @res,@$data;
		}
	}
}

  my @seq = split("",$seq);
   my $tab;
my $nb = scalar(@{$run->getPatients});
$nb --;
my $step = int(scalar(@res) / 500);
$step = 1 if $step == 0;
#$step = 1 if $is_exon == 1;

for (my $i = 0;$i<@res;$i+=$step){
	my $covm=0;
	my $cov = 0;
	my $seq='';
	for (my $z=0;$z<$step;$z++){
		$covm += int($res[$i+$z]/$nb);
		$cov +=  $pdata->[$i+$z];
		#$seq.=lc($tt->{seq} );
		}
	$covm = int($covm/$step);
	$cov = int($cov/$step);
   	my $tt;
   	$tt->{seq} = $seq[$i];
	$tt->{posm} = ($i-$splice_5);
   	if ($i< $splice_5){
   		$tt->{seq} = "-".lc($tt->{seq});
   		$tt->{seq} = $tt->{posm};
   	}
   	if ($i> $eend){
   		$tt->{seq} = "+".($i-$eend);
   		
   	}

   	$tt->{cov} = $pdata->[$i];
	
   	 $tt->{covm} = $covm;
   	push(@$tab,$tt);
   	
   }

 export_data::print($project,$cgi,$tab);

exit(0);




