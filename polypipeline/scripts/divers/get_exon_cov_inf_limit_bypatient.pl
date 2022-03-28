#!/usr/bin/perl
use CGI qw/:standard :html3/;
use FindBin qw($Bin); 
use strict;

use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";

use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use Set::IntSpan::Fast::XS; # module for span class objects
use Spreadsheet::WriteExcel; 
use Spreadsheet::WriteExcel::Utility; 
use preload_coverage;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
#use Parallel::ForkManager;



my $projectName;
my $cov_limit;
my $intronic;
my $padding;
my $utr;
my $order;
#my $fork; #lignes inactives

GetOptions(
	'project=s'		=> \$projectName,
	'covlimit=s'		=> \$cov_limit,
	'intronic=s'		=> \$intronic,
	'utr=s'				=> \$utr,
	'padding=s'		=> \$padding
#	"fork=s" =>\$fork,
	);


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
	
	my $dir_output = $project->getCacheDir() . "/low_exons/"; #mais manque le mkdir.....
	system("mkdir $dir_output");
	system("chmod a+rwx $dir_output");
	warn $fileout_base ;
	warn $dir_output ;
	my $dir_file= "$dir_output$fileout_base.txt";
	open(FILEOUT, ">$dir_output$fileout_base.txt") or die("can't open $dir_output$fileout_base.txt") ;
	open(FILEOUT2, ">".$dir_output.$fileout_base.".html") or die("can't open ".$dir_output.$fileout_base.".html") ;
	print FILEOUT $projectName.": exons mal couverts avec minimum de couverture < ".$limit." et padding de ".$padding;
	print FILEOUT2 qq{<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">};
	print FILEOUT2 qq{<head>$projectName : Low_exons coverage \< $limit </head><body><title>Low exons coverage \< $limit for $projectName</title> <br>};
	
	
###Objets transcripts
warn Dumper $project->bundle_transcripts();
my @transcripts_cgi = @{$project->bundle_transcripts() } ;
my @transcripts_by_name = sort{$a->getGene->external_name cmp $b->getGene->external_name} map{$project->newTranscript($_)} @transcripts_cgi ;
#warn Dumper @transcripts_by_name ;
#die();

my @transcripts_by_ordered_chrom ;
foreach my $chrom (1..22, "X", "Y", "MT"){
	my @list_tr_chrom ;
	foreach my $tr (@transcripts_cgi){
		if ($project->newTranscript($tr)->getChromosome()->name() eq $chrom) {
			push(@list_tr_chrom, $tr );
			}
	}

	unless (scalar(@list_tr_chrom) eq 0) {
		my @list_tr_chrom_trie = sort{$a->getChromosome()->name() <=> $b->getChromosome()->name() || $a->start <=> $b->start} map{$project->newTranscript($_)} @list_tr_chrom ;
		foreach my $trans (@list_tr_chrom_trie) {
			push(@transcripts_by_ordered_chrom, $trans) ;	
	}
	} 
	
} 


my @transcripts ; #les noms, pas des objets ? 
if ($order eq "name") {
		@transcripts = @transcripts_by_name;
}
else {@transcripts =@transcripts_by_ordered_chrom ;}

	


	#récupération de la liste des transcrits stockée en base de donnée pour les paramètres déterminés padding et utr pour lesquels un minimum a été trouvé sur un des exons
	#rajouter le paramètre intronic ???
#	my $list_transcripts_min_min = $no->get($project->name,"minimum-".$padding."-".$utr) ;
	#si cette liste n'existe pas, on la génère
#	unless ($list_transcripts_min_min){
#		warn "coucou, pas de liste pre calcule";
#		$list_transcripts_min_min = preload_coverage::computeLowTranscripts($project,\@transcripts,$no,$intronic,$utr,$padding,1);
##		$list_transcripts_min_min = computeLowTranscripts_script($project,\@transcripts,$no,$intronic,$utr,$padding,1);
#		$no->put($project->name,"minimum-".$padding."-".$utr,$list_transcripts_min_min) ;
#		
#	}

 	#liste des min des min des exons pour tous les patients  et tous les transcrits 
	my $list_transcripts_min_min = preload_coverage::computeLowTranscripts($project,\@transcripts,$no,$intronic,$utr,$padding,1);  #boucle sur les transcrits, pas de fork dans la fonction
	$no->put($project->name,"minimum-".$padding."-".$utr,$list_transcripts_min_min) ;
#	warn Dumper $list_transcripts_min_min ;
	
	#patients
	my $lPatients = $project->getPatients() ;

my $list_low_transcripts;
#possible de forker ici sur les transcrits , mais après tests, ça n'accélère pas
#my $pm = new Parallel::ForkManager($fork);

foreach my $t (@transcripts){
	my @z;
#	my $pid = $pm->start and next;
#	warn "start $t->name()";
	foreach my $p (@$lPatients){
		#récupère dans @z la liste des couvertures minimum pour les exons pour chaque transcrit, afin de faire le test du min
		push(@z,$list_transcripts_min_min->{$t->id}->{$p->{id}}) ####problem, on a rien comme valeur =>résolu
	}
	#on prend le min des min
	my $min = min(@z); 
	next if $min >= $limit;
#	warn $t->id." : ".$min;
	#si ce min est < à la limite, on stocke les id des transcrits dans une liste de LowTranscripts
#	$list_transcripts->{$t->id} ++; 
	$list_low_transcripts->{$t->id} ++; 
}
	#liste réduite aux transcripts dont au moins un exon sur au moins un patient présente un min de couv <limite
#	warn Dumper $list_transcripts ;
	warn "nb transcripts avec exons inf cov_limit : ".keys(%$list_low_transcripts) ; 
#écriture des résultats dans des fichiers xls, txt et html
	my $workbook  = Spreadsheet::WriteExcel->new("/data-xfs/dev/cfourrag/GenBo/script/ngs_exome/cecile/Low_exons_".$projectName.".xls");
	my $workbookdigest  = Spreadsheet::WriteExcel->new("/data-xfs/dev/cfourrag/GenBo/script/ngs_exome/cecile/Low_exons_".$projectName."_rendu.xls");
	my $worksheet2 = $workbookdigest->add_worksheet($projectName);
	my $format2 = $workbookdigest->add_format();
    $format2->set_bold();
	my $formatit = $workbookdigest->add_format();
    $formatit->set_italic();

	#faire une double boucle sur cette liste pour extraire les min <limit pour les patients concernés
	my $cpat=0;
	foreach my $patient (sort {$a->{name} cmp $b->{name}} @$lPatients){
#		warn $patient->name();
		
		#fichier texte avec données minimales
		print FILEOUT "\n".$patient->name()."\n" ;
		print FILEOUT2 "<br><b>".$patient->name()."</b><br>";
		#tableur excell avec données complètes
		my $worksheet = $workbook->add_worksheet($patient->name());
		my $format = $workbook->add_format();
    	$format->set_bold();
    	
		

		$worksheet->write( 0, 0, $projectName,  $format);
		$worksheet->write( 1, 0, "Patient :" , $format);
		$worksheet->write( 1, 1,$patient->name(), $format);
		$worksheet2->write( $cpat, 0, "Patient :" , $format2);
		$worksheet2->write( $cpat, 1,$patient->name(), $format2);
		
	# line header
	my $nbline = 4;
	
	my $formatTh= $workbook->add_format(border   => 1);
	$formatTh->set_bold();
	my $formatThe= $workbook->add_format(border   => 1);
	my $nbcol = 0;
	$worksheet->write_string($nbline, $nbcol, "Gene",  $formatTh );
	$nbcol++;
	$worksheet->write_string($nbline, $nbcol, "Transcript",  $formatTh );
	$nbcol++;
	$worksheet->write_string($nbline, $nbcol, "Chromosome",  $formatTh );
	$nbcol++;
	$worksheet->write_string($nbline, $nbcol, "Exon",  $formatTh );
	$nbcol++;
	$worksheet->write_string($nbline, $nbcol, "Start",  $formatTh );
	$nbcol++;
	$worksheet->write_string($nbline, $nbcol, "End",  $formatTh );
	$nbcol++;
	$worksheet->write_string($nbline, $nbcol, "Mean",  $formatTh );
	$nbcol++;
	$worksheet->write_string($nbline, $nbcol, "Min",  $formatTh );
	$nbcol++;
	$nbline++;
	#results
	$nbcol = 0;
	my $nbcol2=0;	

	my $transcripts_cgi = $project->bundle_transcripts() ;

		#essai tri des transcripts par ordre alphabétique des gènes : 
#		my @transcripts_cgi = @{ $projectP->bundle_transcripts() };
#	my @transcripts =	  sort { $a->getGene->external_name cmp $b->getGene->external_name }	  map  { $projectP->newTranscript($_) } @transcripts_cgi;
		my @transcripts_tries = sort { $a->getGene->external_name cmp $b->getGene->external_name }	  map  { $project->newTranscript($_) } keys(%$list_low_transcripts);
#		warn scalar(@transcripts_tries);
#		warn Dumper @transcripts_tries ; #$VAR210 = $VAR1->{'project'}{'objects'}{'transcripts'}{'ENST00000252744_5'};
		
		my $gene_name_precedent;
		my $ENST_precedent;
		my $flag =0 ;
		my $debut = 0 ;
		foreach my $t (@transcripts_tries) {
		my $debug;
		$debug =1 if $t->id eq "ENST00000329235_X";
		
#		warn Dumper $t ;
		my $transcript_chrom = $t->getChromosome()->name;
		my $gene_name =$t->gene_external_name();
		if ($gene_name eq "SRY") {next;}
		#pour n'afficher que le résultat d'un des ENST du gène si il y a plusieurs ENST par gène
		if ($gene_name eq $gene_name_precedent) {if ($flag==0) {print FILEOUT " (".$ENST_precedent.")" ;  $flag = 1 ; next;} else {next ;} }
		else {$flag =0;}
#		warn $transcript_obj->name."-".$gene_name."-chr".$transcript_chrom ;
#		warn $t->name ;
#				my $capture_intspan = $t->getChromosome->getIntSpanCapture();
#		my $exons = $transcript_obj->getExons();
		my $exons = $t->getExons();
		my $c =1;
		my $chaine_num_exons ="exons: ";
		my $hGene_lowExons ;
		my @list_low_exons ;
		foreach my $exon (sort{$a->end*$a->strand <=> $b->end*$b->strand } @$exons){
#			my $s1 = $exon->getGenomicSpan()->intersection($capture_intspan);
			
			my $sstart = $exon->{start};
			my $send = $exon->{end};
			my $intspan_sans_utr = $exon->intspan_no_utr();
#			warn $transcript_obj->name."-exon".$c."   ".$sstart."---".$send;
			
			my ($pos, $intspan_tmp) = $exon->return_start_end_no_utr(padding=>$padding);
#			warn Dumper $pos ;
			unless ($pos){$c +=1 ; next ;} #bizarre 
			my $sstart_sans_utr = $pos->{start};
			my $send_sans_utr = $pos->{end};
			
#			warn "exon ".$c." : ".$sstart."---".$send;
#			warn "exon ".$c." : ".$sstart_sans_utr."---".$send_sans_utr;

			
				my $res2 = $exon->getTranscript->getGene->get_coverage($patient)->coverage($sstart_sans_utr,$send_sans_utr);
				
#				my ($mean, $intspan2, $min) = $exon->compute_mean_min_intstpan ();
				my $mean = $res2->{mean} ;
#				warn $mean ;
				my $min = $res2->{min} ;
#				warn $min ;
				
				
				if ($min <= $cov_limit) {
					warn $res2->{min}." "." $c ".$patient->name if $debug; 
#					warn "exon".$c."-".$min ;
					push(@list_low_exons, $c);
#					warn "liste ".join("," ,@list_low_exons);
					$hGene_lowExons->{$gene_name}=@list_low_exons;
#					warn  Dumper "coucou".$hGene_lowExons ;
					
					my $nbcol = 0;
					$chaine_num_exons .= $c.",";
					$worksheet->write_string($nbline, $nbcol, $gene_name,  $formatThe );
					$nbcol++;
					$worksheet->write_string($nbline, $nbcol, $t->name,  $formatThe );
					$nbcol++;
					$worksheet->write_string($nbline, $nbcol, $transcript_chrom,  $formatThe );
					$nbcol++;
					$worksheet->write_string($nbline, $nbcol, $c,  $formatThe );
					$nbcol++;
					$worksheet->write_string($nbline, $nbcol, $sstart_sans_utr,  $formatThe );
					$nbcol++;
					$worksheet->write_string($nbline, $nbcol, $send_sans_utr,  $formatThe );
					$nbcol++;
					$worksheet->write_string($nbline, $nbcol, $mean,  $formatThe );
					$nbcol++;
					$worksheet->write_string($nbline, $nbcol, $min,  $formatThe );
					$nbcol++;
					$nbline++;
				}
			$c +=1 ;
			
			}
		

if ($hGene_lowExons){ 
#			warn Dumper $hGene_lowExons;
#			warn "liste1 ".join("," ,@list_low_exons);
			warn "coucou ".$patient->name if $debug;
			if (scalar(@list_low_exons) eq 1) {
				if ($debut ==0){ print FILEOUT $gene_name.": exon ".join("," ,@list_low_exons); print FILEOUT2 "<i>".$gene_name."</i>: exon ".join("," ,@list_low_exons)." ; ";}
				else {
				print FILEOUT " ; ".$gene_name.": exon ".join("," ,@list_low_exons);
				print FILEOUT2 "<i>".$gene_name."</i>: exon ".join("," ,@list_low_exons)." ; ";
				$worksheet2->write($cpat+1,$nbcol2,$gene_name, $formatit);
				$worksheet2->write($cpat+1,$nbcol2+1,": exon : ".join("," ,@list_low_exons)." ; ");
				}
			}
			else {if ($debut ==0){ print FILEOUT $gene_name.": exon ".join("," ,@list_low_exons); }
				else {
				print FILEOUT " ; ".$gene_name.": exons ".join("," ,@list_low_exons);
				print FILEOUT2 "<i>".$gene_name."</i>: exons ".join("," ,@list_low_exons)." ; ";
				$worksheet2->write($cpat+1,$nbcol2,$gene_name, $formatit);
				$worksheet2->write($cpat+1,$nbcol2+1,": exons ".join("," ,@list_low_exons)." ; ");
				}
			}
		$debut = 1;	
		$nbcol2++;
		$nbcol2++;
		$gene_name_precedent=$gene_name ;
		$ENST_precedent=$t->name;
		}
		
		}
	$cpat++;
	$cpat++;
	}

	close(FILEOUT);
	close(FILEOUT2);
	$workbook->close();

#retourne cette liste de LowTranscripts (par leur id)
return $list_low_transcripts;

}


getLowExons_by_patient_and_transcript($projectName, $cov_limit, $padding);


