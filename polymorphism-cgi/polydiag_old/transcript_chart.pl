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
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
my $nb_gene_by_patient = 3;
my $nb_exon_by_genes = 10;

my $buffer = GBuffer->new();

my $cgi          = new CGI();
my $project_name = $cgi->param('project');

my $project = $buffer->newProject(-name=>$project_name);
my $runs = $project->getRuns();

my $captures = $project->getCaptures();

my $vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>$captures->[0]->validation_db());

my $patient_name = $cgi->param('patient');
my $cgi_transcript =  $cgi->param('transcript');


$patient_name ="all" unless $patient_name;
my $patients = $project->get_list_patients($patient_name," ");




my $patient_name = $cgi->param('patients');
my $cgi_transcript =  $cgi->param('transcripts');

my $vp =  $project->getPatient($patient_name);
my $run = $vp->getRun();

my $splice_5 = 10;
my $limit    = 0;

$limit = $cgi->param('limit');
$splice_5 = $cgi->param('span');
$splice_5 = 1 unless $splice_5;
my $exons_red;
my $common_red;
my $data; 
my $all_htranscripts;
my $patient_exons;
my $tr = $project->newTranscript($cgi_transcript);
my $cgi    = new CGI();
my @res;
my $seq;
my $pdata;
my @texon;
my @pexon;
my $nb_patient = scalar(@{$run->getPatients}) -1;
my @exon_name;
push(@exon_name,"");
push(@pexon,0);
push(@texon,0);
my @exons = sort{$a->start*$a->strand <=> $b->start*$b->strand} @{$tr->getExons};
my $i=0;
my @hexons;
foreach my $exon (@exons){
	 
	my $h;
	$h->{name} = $exon->name;
foreach my $patient (@{$run->getPatients}){
	#next unless $patient->name eq $patient_name;
	
	
		
		my ($mean,$intspan,$min) =$exon->mean_intspan_coverage_coding($patient,0,1);
	 #  my $data  = $exon->getTranscript()->get_coverage($patient)->coverage($exon->start,$exon->end);
		#my $nn = $data->{nb};
		#my $s = $data->{sum};
		#my $mean =$data->{mean};
		
	
		if ($patient->name eq $patient_name){
				push(@pexon,int($mean));
				$h->{meanp} = int($mean);	
				push(@exon_name,$exon->name);
		}
		else {
			#push(@texon,int($mean/$nb_patient));
			$h->{meanr} += ($mean);	
			$texon[$i] += int($mean/$nb_patient);
			#push (@texon , $mean/$nb_patient);
		}
	}
	$h->{meanr} = int($h->{meanr}/$nb_patient);	
	push(@hexons,$h);
	$i++;
}
 my @seq = split("",$seq);
 my $tab;
my $nb = scalar(@{$run->getPatients});
$nb --;
foreach my $hexon (@hexons){
   	my $tt;
 	$tt->{seq} = $hexon->{name};
   	$tt->{cov} = $hexon->{meanp};
	$tt->{title} = $tr->getGene->external_name." ".$tr->name." ".$tr->external_name."  patient: ".$patient_name;;
   	 $tt->{covm} = $hexon->{meanr};
   	push(@$tab,$tt);
   	
   }
 export_data::print($project,$cgi,$tab);

exit(0);




