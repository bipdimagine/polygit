#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";
#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/validation_variation"; 
use html; 
use GBuffer;
#use Text::CSV qw( csv );
#use CGI::Carp qw(fatalsToBrowser);;
use Spreadsheet::WriteExcel;

$| =1;
my $buffer = GBuffer->new();

my $cgi          = new CGI();


my $prg =  $cgi->url(-relative=>1);
my $url_image = url(-absolute=>1);
$url_image =~ s/$prg/image_cnv.pl/;
my $project_name = $cgi->param('project');
my $patient_name = $cgi->param('patient');
my $project = $buffer->newProject(-name=>$project_name);
print "Content-type: application/msexcel\n";
print "Content-Disposition: attachment;filename=".$project->name().".xls\n\n";


		my $workbook  = Spreadsheet::WriteExcel->new(\*STDOUT);
	my $worksheet = $workbook->add_worksheet();
	
my $patients = $project->getPatients();
if ($patient_name){
	$patients = [$project->getPatient($patient_name)];
}
my $col =0 ;
my $row =0;
my $hscore;
my @h;
push(@h,["chromosome",'gene',"transcript","start","end","exon",map{$project_name."-".$_->name } sort{$a->name cmp $b->name} @$patients]);
foreach my $v (@h){
	$worksheet->write($row, $col,  $v);
	$col++
}
$row++;
$col =0;
#print join(',',@{$t[0]});
#print "\n";
	my @transcripts_cgi = @{$project->bundle_transcripts() } ;
	
	foreach my $tname (@transcripts_cgi){
		my $transcript = $project->newTranscript($tname);
		my $exons = $transcript->getExons();
		 my $tree = Set::IntervalTree->new;
		 foreach my $exon (@$exons){
		 	$tree->insert($exon->name,$exon->start,$exon->end);
		 }
		my @primers = sort{$a->end*$a->strand <=> $b->end*$b->strand } @{$transcript->getPrimers()};
		
		for (my $i=0;$i<@primers;$i++){
		
		my $primer = $primers[$i];
		my @l;
		push(@l, $primer->getChromosome->name);
		push(@l, $transcript->getGene->external_name);
		push(@l,$transcript->name);
		push(@l, $primer->start,$primer->end);
		my $ex = $tree->fetch($primer->start, $primer->end);
		#my $exons = $primer->getExons();
		#my $et = join(";",map{$_->name} grep{$_->getTranscript()->name eq $transcript->name} @$exons);
		push(@l, join(";",@$ex));
		foreach my $patient (sort{$a->name cmp $b->name} @$patients){
			push(@l, $primer->cnv_score($patient));
			 #warn $primer->cnv_score($patient)." ".$patient->meanDepth($primer->getChromosome->name,$primer->start,$primer->end);
			 
			 #$primer->level($patient);
		}
		foreach my $v (@l){
			$worksheet->write($row, $col,  $v);
			$col++
		}
		$row++;
		$col =0;
	
	}
	}
	$workbook->close();
	warn "ok";
	exit(0);
	#csv (in => \@t, out => \*STDOUT);