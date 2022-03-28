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

my $buffer = GBuffer->new();
my $cgi          = new CGI();


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
	my $project = $buffer->newProject( -name => $projectName );
	my $no =  $project->noSqlCoverage();
	
	#attribuer une valeur aux variables utr et padding sinon, ça plante dans la fonction statistic_coverage_coding
	if ($utr eq undef) {$utr = 0;}
 	if ($padding eq undef) {$padding = 0;}	
	
	my $fileout_base= "Exons_inf_".$limit."_padding_".$padding."_".$projectName ;
	#création dossier dans dossier du cache : /data-beegfs/polycache/HG19/NGS2016_1223
	

###Objets transcripts
my @transcripts_cgi = @{$project->bundle_transcripts() } ;
my @transcripts_by_name = sort{$a->getGene->external_name cmp $b->getGene->external_name} map{$project->newTranscript($_)} @transcripts_cgi ;
#warn Dumper @transcripts_by_name ;
#die();

#my @transcripts_by_ordered_chrom ;
#
##foreach my $chrom (1..22, "X", "Y", "MT"){
#	my @list_tr_chrom ;
#	foreach my $tr (@transcripts_cgi){
#	#	if ($project->newTranscript($tr)->getChromosome()->name() eq $chrom) {
#			push(@list_tr_chrom, $tr );
#		#	}
#	}
#
#	unless (scalar(@list_tr_chrom) eq 0) {
#		my @list_tr_chrom_trie =  ;
#		foreach my $trans (@list_tr_chrom_trie) {
#			push(@transcripts_by_ordered_chrom, $trans) ;	
#	}
#	} 
	
#} 


#my @transcripts ; #les noms, pas des objets ? 

	

 	#liste des min des min des exons pour tous les patients  et tous les transcrits 
	#my $list_transcripts_min_min = preload_coverage::computeLowTranscripts($project,\@transcripts_cgi,$no,$intronic,$utr,$padding,1);  #boucle sur les transcrits, pas de fork dans la fonction
	#warn Dumper $list_transcripts_min_min;
	#die();
	#$no->put($project->name,"minimum-".$padding."-".$utr,$list_transcripts_min_min) ;
#	warn Dumper $list_transcripts_min_min ;
	
	#patients
	my $lPatients = $project->getPatients() ;

#my $list_low_transcripts;
#
#
#foreach my $t (@transcripts_cgi){
#	my @z;
#
#	foreach my $p (@$lPatients){
#		#récupère dans @z la liste des couvertures minimum pour les exons pour chaque transcrit, afin de faire le test du min
#		push(@z,$list_transcripts_min_min->{$t->id}->{$p->{id}}) ####problem, on a rien comme valeur =>résolu
#	}
#	#on prend le min des min
#	my $min = min(@z); 
#	next if $min >= $limit;
#
#	$list_low_transcripts->{$t->id} ++; 
#}


	#liste réduite aux transcripts dont au moins un exon sur au moins un patient présente un min de couv <limite
#	warn Dumper $list_transcripts ;
	#warn "nb transcripts avec exons inf cov_limit : ".keys(%$list_low_transcripts) ; 

#écriture des résultats dans des fichiers xls, txt et html
#	my $workbook  = Spreadsheet::WriteExcel->new("/data-xfs/dev/cfourrag/GenBo/script/ngs_exome/cecile/Low_exons_".$projectName.".xls");
#	my $workbookdigest  = Spreadsheet::WriteExcel->new("/data-xfs/dev/cfourrag/GenBo/script/ngs_exome/cecile/Low_exons_".$projectName."_rendu.xls");
#	my $worksheet2 = $workbookdigest->add_worksheet($projectName);
#	my $format2 = $workbookdigest->add_format();
#    $format2->set_bold();
#	my $formatit = $workbookdigest->add_format();
#    $formatit->set_italic();

	#faire une double boucle sur cette liste pour extraire les min <limit pour les patients concernés
	my $cpat=0;
my $t = time;
	
	warn abs(time -$t);
	my @transcripts_tries = sort { $a->getGene->external_name cmp $b->getGene->external_name }	  map  { $project->newTranscript($_) } @transcripts_cgi;
	#my $list_transcripts_min_min = preload_coverage::computeLowTranscripts($project,\@transcripts_tries,$no,$intronic,$utr,$padding,1);
	warn abs(time -$t);
		#warn Dumper $list_transcripts_min_min->{ENST00000329235_X};
		#die();
		#warn abs(time - $t);
		#warn scalar(keys %$list_transcripts_min_min);
	my $fork = 6;
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
			 $hh->{ $h->{patient}} = $h->{line};
			#print $h->{line} if $h->{line};
			
		}
	);
	
	my $id= time ;
	print qq{<div style="display: none">};
	foreach my $patient (sort {$a->{name} cmp $b->{name}} @$lPatients){
		#last if $patient->name() ne "GUE_Gil";
			$hrun->{$id} ++;
			$id ++;
			print ".";
			my $pid = $pm->start and next;
			my $line;
			$line =  "<br><b>".$patient->name()."</b><br>";
		
		warn $patient->name();
		
		foreach my $t (@transcripts_tries) {
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
			$line .=  $t->getGene->external_name." ".": exon ".join("," ,@t)." (".$t->name.");" if @t;
			}
		warn "end";
		$pm->finish( 0, {id=>$id,line=>$line,patient=>$patient->name} );
		#die();
		#$no->lmdb($chr->name);
	}
$pm->wait_all_children();
	print qq{</div>};
	foreach my $patient (sort {$a->{name} cmp $b->{name}} @$lPatients){
		print $hh->{$patient->name};
	}
	
exit(0);
	
#		
##		my $gene_name_precedent;
##		my $ENST_precedent;
##		my $flag =0 ;
##		my $debut = 0 ;
#		foreach my $t (@transcripts_tries) {
###		warn Dumper $t ;
#		my $transcript_chrom = $t->getChromosome()->name;
#		my $gene_name =$t->gene_external_name();
#		if ($gene_name eq "SRY") {next;}
##		#pour n'afficher que le résultat d'un des ENST du gène si il y a plusieurs ENST par gène
#		if ($gene_name eq $gene_name_precedent) {if ($flag==0) {print FILEOUT " (".$ENST_precedent.")" ;  $flag = 1 ; next;} else {next ;} }
#		else {$flag =0;}
#		my $exons = $t->getExons();
#		my $c =1;
#		my $chaine_num_exons ="exons: ";
#		my $hGene_lowExons ;
#		my @list_low_exons ;
#		foreach my $exon (sort{$a->end*$a->strand <=> $b->end*$b->strand } @$exons){
#			my $s1 = $exon->getGenomicSpan()->intersection($capture_intspan);
#			
#			my $sstart = $exon->{start};
#			my $send = $exon->{end};
#			my $intspan_sans_utr = $exon->intspan_no_utr();
###			warn $transcript_obj->name."-exon".$c."   ".$sstart."---".$send;
##			
#			my ($pos, $intspan_tmp) = $exon->return_start_end_no_utr(padding=>$padding);
###			warn Dumper $pos ;
#			unless ($pos){$c +=1 ; next ;} #bizarre 
#			my $sstart_sans_utr = $pos->{start};
#			my $send_sans_utr = $pos->{end};
#
##			
#				my $res2 = $exon->getTranscript->getGene->get_coverage($patient)->coverage($sstart_sans_utr,$send_sans_utr);
#				my $mean = $res2->{mean} ;
###				warn $mean ;
#				my $min = $res2->{min} ;
###				warn $min ;
##				
##				
#				if ($min <= $cov_limit) {
###					warn "exon".$c."-".$min ;
#					push(@list_low_exons, $c);
###					warn "liste ".join("," ,@list_low_exons);
#					$hGene_lowExons->{$gene_name}=@list_low_exons;
###					warn  Dumper "coucou".$hGene_lowExons ;
##					
#					my $nbcol = 0;
#					$chaine_num_exons .= $c.",";
#
#				}
#			$c +=1 ;
#			}
#
#if ($hGene_lowExons){ 
##			warn Dumper $hGene_lowExons;
##			warn "liste1 ".join("," ,@list_low_exons);
#			
#			if (scalar(@list_low_exons) eq 1) {
#				if ($debut ==0){ print  $gene_name.": exon ".join("," ,@list_low_exons); }
#				else {
#				print  "<i>".$gene_name."</i>: exon ".join("," ,@list_low_exons)." ; ";
##				$worksheet2->write($cpat+1,$nbcol2,$gene_name, $formatit);
##				$worksheet2->write($cpat+1,$nbcol2+1,": exon : ".join("," ,@list_low_exons)." ; ");
#				}
#			}
#			else {if ($debut ==0){ print  $gene_name.": exon ".join("," ,@list_low_exons); }
#				else {
#			#	print FILEOUT " ; ".$gene_name.": exons ".join("," ,@list_low_exons);
#				print  "<i>".$gene_name."</i>: exons ".join("," ,@list_low_exons)." ; ";
##				$worksheet2->write($cpat+1,$nbcol2,$gene_name, $formatit);
##				$worksheet2->write($cpat+1,$nbcol2+1,": exons ".join("," ,@list_low_exons)." ; ");
#				}
#			}
#		$debut = 1;	
##		$nbcol2++;
##		$nbcol2++;
#		$gene_name_precedent=$gene_name ;
#		$ENST_precedent=$t->name;
#		}
#		
#		}
#	$cpat++;
#	$cpat++;
#	}

#retourne cette liste de LowTranscripts (par leur id)
#return $list_low_transcripts;

}





