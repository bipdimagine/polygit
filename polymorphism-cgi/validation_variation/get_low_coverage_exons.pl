#!/usr/bin/perl
use CGI qw/:standard :html3/;
use FindBin qw($Bin); 
use strict;

use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";

use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use Set::IntSpan::Fast::XS; # module for span class objects
use Spreadsheet::WriteExcel; 
use Spreadsheet::WriteExcel::Utility; 
use preload_coverage;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils ':all';
use lib "$Bin/../packages/validation_variation";
use html;
#use Parallel::ForkManager;

my $cgi          = new CGI();
$|=1;

my $projectName = $cgi->param('project');
my $cov_limit =$cgi->param('covlimit'); 
my $intronic;# = $cgi->param('intronic');
my $padding =$cgi->param('padding'); 
my $utr = $cgi->param('utr');;
html::print_cgi_header($cgi);
#my $fork; #lignes inactives


getLowExons_by_patient_and_transcript($projectName, $cov_limit, $padding);

###export fichier excell (résultats complets et light), export fichier txt ou html


sub getLowExons_by_patient_and_transcript {
	my ($projectName, $limit, $padding) = @_;
	my $buffer = GBuffer->new();
	my $project = $buffer->newProjectCache( -name => $projectName );
	my $no =  $project->noSqlCoverage();
	
	#attribuer une valeur aux variables utr et padding sinon, ça plante dans la fonction statistic_coverage_coding
	if ($utr eq undef) {$utr = 0;}
 	if ($padding eq undef) {$padding = 0;}	
	
	my $fileout_base= "Exons_inf_".$limit."_padding_".$padding."_".$projectName ;
	#création dossier dans dossier du cache : /data-beegfs/polycache/HG19/NGS2016_1223
	

###Objets transcripts
my @transcripts_cgi = @{$project->bundle_transcripts() } ;
my @transcripts_by_name = sort{$a->getGene->external_name cmp $b->getGene->external_name} map{$project->newTranscript($_)} @transcripts_cgi ;

	print qq{<div style="display: none">};
	my $lPatients = $project->getPatients() ;
	my $cpat=0;
my $t = time;
	print "--";
	my @transcripts_tries = sort { $a->getGene->external_name cmp $b->getGene->external_name }	  map  { $project->newTranscript($_) } @transcripts_cgi;
	my $fork = 10;
	my $hrun;
	my $hh;
	my $pm   = new Parallel::ForkManager($fork);	
		$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $h ) = @_;

			unless ( defined($h) or $exit_code > 0 ) {
				print
				  qq|No message received from child process $exit_code $pid!\n|;
				die();
				return;
			}
			my $id = $h->{run_id};
			delete $hrun->{ $h->{run_id} };
			print "OO\n";
			 $hh->{ $h->{patient}} = $h->{line};
			
		}
	);
	
	my $id= time ;
	print "--";
	
	foreach my $patient (sort {$a->{name} cmp $b->{name}} @$lPatients){
		#last if $patient->name() ne "GUE_Gil";
			$hrun->{$id} ++;
			$id ++;
			print ".>".$patient->name;
			my $pid = $pm->start and next;
			my $line;
			$line =  "<br><b>".$patient->name()."</b><br>";
		
		
		foreach my $t (@transcripts_tries) {
				print "*";
		#	next unless exists $list_transcripts_min_min->{$t->{id}};
			
		#	next if $list_transcripts_min_min->{$t->{id}}->{$patient->id} > $cov_limit;
			#my $t = $project->newTranscript($tid);
			my $capture_intspan = $t->getChromosome->getIntSpanCapture();
			my $debug;
			$debug =1 if $t->getGene->external_name eq "CROCC";
			my $exons = $t->getExons();
			my @t;
			my $nb = 0;
			my $intspan = Set::IntSpan::Fast::XS->new();
			print "=";
			foreach my $exon (sort{$a->end*$a->strand <=> $b->end*$b->strand } @$exons){
				$nb ++;
				my $s1 = $exon->getGenomicSpan()->intersection($capture_intspan);
			 	next  if $s1->is_empty;
				my ($pos, $intspan_tmp) = $exon->return_start_end_no_utr(padding=>$padding);
				unless ($pos){ next ;} #bizarre 
				my $min = $patient->minDepth($t->getChromosome->name,$pos->{start},$pos->{end});
					if ($min <= $cov_limit) {
						push(@t,$nb);
					}
			}
			print "=";
			$line .=  $t->getGene->external_name." ".": exon ".join("," ,@t)." (".$t->name.");" if @t;
			}
		$pm->finish( 0, {id=>$id,line=>$line,patient=>$patient->name} );
		#die();
		#$no->lmdb($chr->name);
	}
$pm->wait_all_children();
	print qq{</div>\n};
	foreach my $patient (sort {$a->{name} cmp $b->{name}} @$lPatients){
		print $hh->{$patient->name};
	}
	
exit(0);
	


}





