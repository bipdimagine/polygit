#!/usr/bin/perl

use Carp;
use strict;
use JSON;
use Data::Dumper;
use CGI qw/:standard :html3/;
use Set::IntSpan::Fast::XS;
use Set::IntervalTree;
use List::Util qw[min max];
use FindBin qw($Bin);
use Storable qw(store retrieve freeze);

use lib "$Bin/../../../../GenBo/lib/obj-nodb/";


use GBuffer;
use GenBoProject;
use GenBoCache;


############################################################################################################
# Cherche tous les variants structuraux pour chaque patient d'un projet et tous callers confondus.
# Stocke le resultat dans une  table de hash qui est freeze : nom du fihier = project.allSV
# Pour chaque SV la localisation chr_start_end donne accès à la liste des gènes compris dans cet interval
############################################################################################################

my $cgi = new CGI;
my $projectname = $cgi->param('project');
my $type ="BND";
#my $path = "/data-xfs/Manue/Test_SV/".$projectname."/manta/";
#my $path_djv = "/data-xfs/Manue/Test_SV/DejaVu/TransLoc/";



# pour lire les vcfs
my $compteur=0;
	
# pour stocker les translocations 
my $hTransLoc;  
my @listHashRes;

# pour le dejavu
my $hdejavu;


	
###################################
#  lecture des vcf et 
#  gestion du projet via les objets Genbo
###################################

my $caller = "manta";
my $buffer = GBuffer->new();	
my $project = $buffer->newProjectCache( -name => $projectname);

# pour enregister les fichiers allBND
my $path = $project->getSVeqDir();

# pour enregistrer le fichier  dejavu
my $path_djv = $project->DejaVuProjectsSVeq_path;

# boucle sur les patients du projets
my @listPatients = grep {$_->isGenome} @{$project->getPatients()};
foreach my $thePatient (@listPatients)
{
	
	my $patientname = $thePatient->name;
	print $patientname."\n";
	
	#next if ($patientname eq "NPH2363");

	my $file_in = $thePatient->_getCallingSVFileWithMethodName($caller,"variations");
	#my $file_in = $thePatient->getSVFile($caller);
	
	# pour construire le path d'acces au bam pour IGV
	my $patBam = $thePatient->bamUrl;
	my ( $path_Bam, $nothing ) = split( /bwa/, $patBam );
	$path_Bam .= "bwa/";
	
	
	my $bamFile = $path_Bam."/".$patientname.".bam";
	
	
	# pour lire le vcf 
	my $hBND;
	
	# ouverture du fichier manta zippé
	my $fd;
	my $ligne;

	
	
	open($fd," zcat $file_in | ") or die("open: $!");
	
		
	# on skipe le Header et la ligne des champs
	while( defined( $ligne = <$fd> ) &&  $ligne =~ m/^##/ ) {}
	
	# traitement des lignes suivantes pour le fichier en entier
	while( defined( $ligne = <$fd> ) )
	{
			
			my @champs = split(/\t/,$ligne);
			
			
			
			# champ PASS
			next unless ($champs[2] =~ m/MantaBND/);
			next if ( ( $champs[6] ne "PASS")); 	# on ne garde que les lignes ou PASS
			
			
			# champ QUAL
			my $qual = $champs[5];
			
			# position imprecise
			next  if ( $champs[7] =~ m/IMPRECISE/);
		
			# champ chromosome
			my $chrom = $champs[0];
			my $oc = $project->existsChromosome($chrom);
			next unless $oc;
			$chrom = $oc->ucsc_name;
			next  if ($chrom =~ m/GL/); 
			next  if ($chrom =~ m/hs37d5/);
			next  if ($chrom =~ m/M/);
			
			my $chrnum;
			
			if ($chrom =~ m/chr/)
			{
				my @val =split(/r/,$chrom);
				$chrnum=$val[1];
			}
			else
			{
				$chrnum=$chrom;
			}
			
			next  if ($chrnum eq "MT"); 	
			next  if ($chrnum eq "M"); 
			next unless $project->isChromosomeName($chrnum);
			$chrnum = 23 if ($chrnum eq "X");
			$chrnum = 24 if ($chrnum eq "Y");
			
	
			# Position de depart
			my $pos = $champs[1];
			
			# BNDID / MATEID / EVENT
			my $infoBND = $champs[7];
			my @tabinfobnd = split(/;/,$infoBND);
			
			my $MATEID;
			my $EVENT_id;
			my $BNDID = $champs[2];
				
			foreach my $v (@tabinfobnd)
			{
					if ( $v =~ m/MATEID/)
					{
						my ($tag,$val) = split(/=/,$v);
						$MATEID = $val;
					}
			}
			
			# trouver le preffixe commun au deux ID BND et MATE
			my @tab1 = split(/:/,$BNDID);
			my @tab2 = split(/:/,$MATEID);
			my $len = min(scalar(@tab1),scalar(@tab2));
			for( my $i=0; $i< $len; $i++ ) 
			{ 
				$EVENT_id .= $tab1[$i] if ($tab1[$i] eq $tab2[$i]);
			}

			
			
			# REF / ALT
			my $ref = $champs[3];
			my $alt = $champs[4];
			
		
			# pour verifier que la translocation concerne deux chromosome differents
			my $infomateid= $champs[4];
			$infomateid  =~ s/\[/:/g;
			$infomateid =~ s/\]/:/g;
			
			my @tabmateid = split(/:/,$infomateid);
			
			my $chr = $tabmateid[1];
			my $oc = $project->existsChromosome($chr);
			next unless $oc;
			$chr = $oc->ucsc_name;
			
			next  if ($chr =~ m/hs37d5/);
			next  if ($chr =~ m/GL/); 
			
			my $chrnum2;
			
			if ($chr =~ m/chr/)
			{
				my @val =split(/r/,$chr);
				$chrnum2=$val[1];
			}
			else
			{
				$chrnum2=$chr;
			}
			
			$chrnum2 = 23 if ($chrnum2 eq "X");
			$chrnum2 = 24 if ($chrnum2 eq "Y");
			
		
			my $BND_ind_id = $chrnum."_".$pos;
			
			$hBND->{$EVENT_id}->{$BND_ind_id}->{"id"}= $BND_ind_id;
			$hBND->{$EVENT_id}->{$BND_ind_id}->{"chr"}= int($chrnum);
			$hBND->{$EVENT_id}->{$BND_ind_id}->{"pos"}= $pos;
			$hBND->{$EVENT_id}->{$BND_ind_id}->{"refalt"}= $ref."/".$alt;
			$hBND->{$EVENT_id}->{$BND_ind_id}->{"qual"}= $qual;
			
			
			# cytoband autour de la position du BND
			my $deb = $pos-50;
			my $end = $pos+50;
		
			my $hband = $project->getChromosome($chrom)->getCytoband($deb, $end);
	
			my @tb;
			foreach my $b ( keys %{ $hband->{'name'} } ) {
				push( @tb, $b );
			}
			my @band = sort( { $a cmp $b } @tb );

			$hBND->{$EVENT_id}->{$BND_ind_id}->{"cytoband"} = join( ",", @band );
		
			# pour trouver les genes compris dans l'interval
			my $objChr = $project->getChromosome($chrom);
			my $tabGenes = $objChr->getGenesByPosition($deb,$end);
	
			my $genes_liste="";
			my $omim=0;
			my $phenotypes="";
	
			
			foreach my $g (@$tabGenes)
			{
				if (($g->start <= $deb) && ($g->end >= $end))
				{
					$genes_liste .= $g->external_name.":".$g->phenotypes."##";
					$omim = 1 if $g->is_omim_morbid();
				}
			} 
		
			$hBND->{$EVENT_id}->{$BND_ind_id}->{"genes"} = $genes_liste;
			$hBND->{$EVENT_id}->{$BND_ind_id}->{"omim"} = $omim;
			
			
	} #fin de la boucle sur les lignes


	# Pour regrouper les deux parties de la translocation 
	foreach my $eventID ( keys %{ $hBND } ) 
	{
			
		my $nb = keys (%{ $hBND->{$eventID}});
		if ($nb == 2)
		{
			my $chr1 = 0;
			my $chr2 = 0;
			my $pos1 = 0;
			my $pos2 = 0;
			my $cytoband1;
			my $cytoband2;
			my $genes1;
			my $genes2;
			my $omim1;
			my $omim2;
			my $qual1;
			my $qual2;
			my $refalt1;
			my $refalt2;
			
			my $ind_id1;
			my $ind_id2;
					
					
			foreach my $ind_id ( keys %{ $hBND->{$eventID} })
			{
					#warn $chr1." ".$pos1." ".$chr2." ".$pos2;
					
					my ($chr,$pos) = split(/_/,$ind_id);
					
					$chr2 = $chr if ($chr1);
					$pos2 = $pos if ($pos1);
					$chr1 = $chr unless ($chr1);
					$pos1 = $pos unless ($pos1);
					
					#warn $chr1." ".$pos1." ".$chr2." ".$pos2;
			}
			
			
			if (int($chr1) == int($chr2))
			{
				$ind_id1 = $chr1."_".$pos1;
				$ind_id2 = $chr2."_".$pos2;
				
				if (int($pos1) >= int($pos2))
				{
					$ind_id1 = $chr2."_".$pos2;
					$ind_id2 = $chr1."_".$pos1;
				}
			}
			
			if (int($chr1) < int($chr2))
			{
				$ind_id1 = $chr1."_".$pos1;
				$ind_id2 = $chr2."_".$pos2;
			}
			
			if (int($chr1) > int($chr2))
			{
				$ind_id1 = $chr2."_".$pos2;
				$ind_id2 = $chr1."_".$pos1;
			}
				
				
			# recuperer les donnees
			$chr1 = $hBND->{$eventID}->{$ind_id1}->{"chr"};
			$pos1 = $hBND->{$eventID}->{$ind_id1}->{"pos"};
			$cytoband1 = $hBND->{$eventID}->{$ind_id1}->{"cytoband"};
			$genes1 = $hBND->{$eventID}->{$ind_id1}->{"genes"};
			$omim1 = $hBND->{$eventID}->{$ind_id1}->{"omim"};
			$qual1 = $hBND->{$eventID}->{$ind_id1}->{"qual"};
			$refalt1 = $hBND->{$eventID}->{$ind_id1}->{"refalt"};
			
			$chr2 = $hBND->{$eventID}->{$ind_id2}->{"chr"};
			$pos2 = $hBND->{$eventID}->{$ind_id2}->{"pos"};
			$cytoband2 = $hBND->{$eventID}->{$ind_id2}->{"cytoband"};
			$genes2 = $hBND->{$eventID}->{$ind_id2}->{"genes"};
			$omim2 = $hBND->{$eventID}->{$ind_id2}->{"omim"};
			$qual2 = $hBND->{$eventID}->{$ind_id2}->{"qual"};
			$refalt2 = $hBND->{$eventID}->{$ind_id2}->{"refalt"};
		
		
			#creer l'event_id correspondant
			my $event_id = $chr1."_".$pos1."_".$chr2."_".$pos2;
			#warn $event_id if ($chr1== $chr2);
			
			$hTransLoc->{$patientname}->{$event_id}->{"id"}= $event_id;
			$hTransLoc->{$patientname}->{$event_id}->{"TRANSLOC"}= $chr1."to".$chr2;
			$hTransLoc->{$patientname}->{$event_id}->{"IGV"}= $bamFile."et".$chr1."_".$pos1."_".$chr2."_".$pos2;
			$hTransLoc->{$patientname}->{$event_id}->{"QUAL"}= int($qual1);
			
			# info du chromosome de départ
			$hTransLoc->{$patientname}->{$event_id}->{"CHROM1"}= $chr1;
			$hTransLoc->{$patientname}->{$event_id}->{"POS1"}= $pos1;
			$hTransLoc->{$patientname}->{$event_id}->{"CYTOBAND1"}= $cytoband1;
			$hTransLoc->{$patientname}->{$event_id}->{"GENES1"}= $genes1;
			$hTransLoc->{$patientname}->{$event_id}->{"OMIM1"}= $omim1;
			$hTransLoc->{$patientname}->{$event_id}->{"REF/ALT1"}= $refalt1;
			
			# info du chromosome d'arrivée
			$hTransLoc->{$patientname}->{$event_id}->{"CHROM2"}= $chr2;
			$hTransLoc->{$patientname}->{$event_id}->{"POS2"}= $pos2;
			$hTransLoc->{$patientname}->{$event_id}->{"CYTOBAND2"}= $cytoband2;
			$hTransLoc->{$patientname}->{$event_id}->{"GENES2"}= $genes2;
			$hTransLoc->{$patientname}->{$event_id}->{"OMIM2"}= $omim2;
			$hTransLoc->{$patientname}->{$event_id}->{"REF/ALT2"}= $refalt2;	
			
			# pour le calcul du déjavu
			$hdejavu->{$event_id}->{$patientname} = 1;
		}
	}
	
	# on freeze la table hashs correspondant au patient
	my $file_out = $path.$patientname.".allBND";
	store(\ %{ $hTransLoc->{$patientname}}, $file_out) or die "Can't store $file_out for ".$patientname."!\n";
	
} #fin de la boucle patient

# pour le dejavu du projet
my $file_djv = $path_djv."/".$projectname.".SVeqDejavu";
store(\ %{ $hdejavu}, $file_djv) or die "Can't store $file_djv for ".$projectname."!\n";


exit(0);






 	