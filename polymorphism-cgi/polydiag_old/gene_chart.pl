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
my $cgi_exon = $cgi->param('exon');
my $cgi_transcript =  $cgi->param('transcripts');
my $tr = $project->newTranscript($cgi_transcript);


my $patients = [sort{$a->name cmp $b->name} @{$project->getPatients()}]; 

my $capture = $project->getCaptures()->[0];

my $vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>$capture->validation_db());




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
my $cgi    = new CGI();

my @pgenes;


my @exons = sort{$a->start <=> $b->start} @{$tr->getExons};
my $sum_gene;
foreach my $patient ( @{$patients}){
	my $i =0;
	my $sum_patient;
	
	
	foreach my $exon (@exons){
		if ($cgi_exon){
			next if  $cgi_exon ne $exon->name();
			
		}
		my $data = $exon->getTranscript()->get_coverage($patient)->coverage($exon->start,$exon->end);
		my $nn = $data->{nb};
		if ($nn >0){
		my $s = $data->{sum};
		$sum_patient += $s /$nn;
		$i++;
		}
	
	}
	
	my $mean = int ($sum_patient/$i);
	push(@pgenes,$mean);
	$sum_gene +=  $mean;
}

my $meang = int($sum_gene/scalar(@$patients));


 my $tab;

for (my $i = 0;$i<@$patients;$i++){
	
   	my $tt;
 	$tt->{seq} = $patients->[$i]->name;
   	$tt->{cov} = $pgenes[$i];
	$tt->{title} = $tr->getGene->external_name." ".$tr->name." ".$tr->external_name;
	$tt->{title} .= "  exon: ".$cgi_exon if $cgi_exon; 
   	 $tt->{covm} = $meang;
   	push(@$tab,$tt);
   	
   }
 export_data::print($project,$cgi,$tab);

exit(0);




