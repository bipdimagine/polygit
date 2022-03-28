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
use Number::Format qw(:subs);

use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";

use GBuffer;
use GenBoProject;
use GenBoCache;

use layout;
use export_excel;
use export_data;

###########################################################################################################
#  1) recupere les infos concernant tous les SV du projet entre en argument et prealablement freeze : nom du fihier = project.allSV
#  2) regroupe les SV des patients entres en argument et qui ont la même liste de gènes sous un id_global
#  3) Si le project contient des trios on regarde la transmission a partir du SVglobal
#  4) filtre le resultats en fonction des donnéees definies via l'interface
#  5) construit le Json resultant pour affichage
#
############################################################################################################

my $cgi = new CGI;
my $TheProjectName = $cgi->param('projectname');
my $thePatientName     = $cgi->param('filename');
my $minlength   = $cgi->param('minlength');
my $maxlength   = $cgi->param('maxlength');
my $maxfreq     = $cgi->param('maxfreq');
my $select_best = $cgi->param('select_best');
my $chr         = $cgi->param('chrom');
my $cytoband    = $cgi->param('cytoband');
my $listOfGenes = $cgi->param('genes');
my $dejavu      = $cgi->param('dejavu');
my $genotype    = $cgi->param('genotype');
my $transmission    = $cgi->param('transmission');
my $print = $cgi->param('print');
my $fileout = $cgi->param('fileout');

# pour stocker les SV
my $hSVPos;          		# tous par positions chevauchantes

# gestion des callers
my @callers=("wisecondor","canvas","manta");

# pour construire le projet
my $buffer  = GBuffer->new();
my $project = $buffer->newProjectCache( -name => $TheProjectName );

# pour distinguer parmi les autres patients du projet
# ceux qui sont de la même famille
my $thePatient = $project->getPatient($thePatientName);
my $thePatientFamily = $thePatient->getFamily();
my $mothername;
my $fathername;

# pour construire le path d'acces au bam pour IGV
#TODO: a mettre au propre MDR
my $patBam = $thePatient->bamUrl;
my ( $path_Bam, $nothing );
if ($patBam =~ /bwa/) {
	( $path_Bam, $nothing ) = split( /bwa/, $patBam );
	$path_Bam .= "bwa/";
}
else {
	my @lMethods = @{$thePatient->alignmentMethods()};
	my $align_method = $lMethods[0];
	( $path_Bam, $nothing ) = split( /$align_method/, $patBam );
	$path_Bam .= "$align_method/";
}

#my @lMethods = @{$thePatient->alignmentMethods()};
#my $align_method = $lMethods[0];
#my ( $path_Bam, $nothing ) = split( /$align_method/, $patBam );
#$path_Bam .= "$align_method/";


#################################################################
#
# pour le json final recuperation du repertoire de sauvegarde via objet genbo
#
##################################################################

my $varDir = $project->getCNVDir();
my $CNVfile = $varDir.$thePatientName.".allCNV";


my $hPat_CNV = retrieve($CNVfile) or die "Can't retrieve datas from " . $CNVfile . " !\n";
my $hCNV;


my $hGroupedCNV;
my @Filtered_listHashRes;

#####################################################
# pour le dejavu_in_this_project (a part pour etre a jour sans avoir a lancer lesetdejavu)
######################################################

my $dejavuProject_file = $varDir.$TheProjectName."_dejavu.allCNV";
my $hdejavuProject = retrieve($dejavuProject_file) or die "Can't retrieve datas from " . $dejavuProject_file . " !\n";

# pour accelerer acces au dejavu = intervall tree
my $htree_dejavuProject;
foreach my $type (keys %{$hdejavuProject})
{
		foreach my $num (keys %{$hdejavuProject->{$type}})
		{
				
				# 1)  detecter les SV identiques 
				$htree_dejavuProject->{$type}->{$num}= Set::IntervalTree->new;
				
				# remplir les arbres :  regrouper les SV chevauchants
				foreach my $id  (keys %{$hdejavuProject->{$type}->{$num}})
				{
						my ( $t, $c, $d, $f ) = split( /_/, $id );
						$htree_dejavuProject->{$type}->{$num}->insert($id,$d,$f);
				}
		}
}		

#############################
# pour le dejavu global (tous les projets WGS)
#############################
#my $dejavuGlobal_file = "/data-xfs/Manue/Test_SV/DejaVu/newVersion/dejavu_allCNV";
#my $hdejavuGlobal = retrieve($dejavuGlobal_file) or die "Can't retrieve datas from " . $dejavuGlobal_file . " !\n";
# pour accelerer acces au dejavu = intervall tree
#my $htree_dejavuGlobal = $project->dejavuSVIntervalTree();		


#my $lmdb = GenBoNoSqlLmdb->new(dir=>"/data-xfs/Manue/Test_SV/DejaVu/newVersion/",mode=>"r",name=>"dejavu_sv",is_compress=>1);
my $dejavudir = $project->DejaVuCNV_path();
my $lmdb = GenBoNoSqlLmdb->new(dir=>$dejavudir,mode=>"r",name=>"dejavu_sv",is_compress=>1);	

# 1) on regroupe les ids strictement identiques
gather_strict_identicalCNV();

# 2) on regroupe les CNV identiques vus par differents callers 
gatherSV_byPosition();
getInfoCallers();

# 3) on filtre en fonction des infos envoyees par l'interface
filtreCNV();

printJson( \@Filtered_listHashRes, $fileout );
exit(0);


###############################################################################
#
#	methodes
#
#################################################################################


sub printJson {
	my ($listHash, $fileout) = @_;
	my $hash;
	my @t;

	@t = sort {$b->{SCORECNV} <=> $a->{SCORECNV}or $b->{SCORE_GENES} <=> $a->{SCORE_GENES}or $a->{SVTYPE} cmp $b->{SVTYPE}or $a->{POS_INT} <=> $b->{POS_INT}} @$listHash;
			
	my $nb=1;
	foreach my $h (@t)
	{
			$h->{"nb"}=$nb++;
	}

	$hash->{'identifier'} = 'id';
	$hash->{'label'}      = 'id';
	$hash->{'items'}      = \@t;
	
	print $cgi->header('text/json-comment-filtered') unless ($print == 1);
	if ($fileout) {
		open (OUT, ">$fileout");
		print OUT encode_json $hash;
		close(OUT);
	}
	else {
		print encode_json $hash;
	}
	print "\n"  unless ($print == 1);
}

# pour regrouper les ids strictement identiques
sub gather_strict_identicalCNV
{
	foreach my $caller (keys %{$hPat_CNV})
	{
		foreach my $type (keys %{$hPat_CNV->{$caller}})
		{
			foreach my $num (keys %{$hPat_CNV->{$caller}->{$type}})
			{
				foreach my $id  (keys %{$hPat_CNV->{$caller}->{$type}->{$num}})
				{
					$hCNV->{$type}->{$num}->{$id}->{'id'} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'id'};			
					$hCNV->{$type}->{$num}->{$id}->{'TYPE'} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'SVTYPE'};		
					$hCNV->{$type}->{$num}->{$id}->{'CHROM'} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'CHROM'};
					$hCNV->{$type}->{$num}->{$id}->{'START'} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'START'};
					$hCNV->{$type}->{$num}->{$id}->{'END'} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'END'};
					$hCNV->{$type}->{$num}->{$id}->{'LEN'} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'SVLEN'};				
					$hCNV->{$type}->{$num}->{$id}->{'GOLD_G_FREQ'} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'GOLD_G_FREQ'};
					$hCNV->{$type}->{$num}->{$id}->{'GOLD_L_FREQ'} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'GOLD_L_FREQ'};
					$hCNV->{$type}->{$num}->{$id}->{'OMIN_MG'} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'OMIN_MG'};
					$hCNV->{$type}->{$num}->{$id}->{'dbVar_status'} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'dbVar_status'};
					$hCNV->{$type}->{$num}->{$id}->{'RANKAnnot'} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'RANKAnnot'};
					$hCNV->{$type}->{$num}->{$id}->{'DUPSEG'} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'DUPSEG'};
					$hCNV->{$type}->{$num}->{$id}->{'CYTOBAND'} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'CYTOBAND'};
					$hCNV->{$type}->{$num}->{$id}->{'DGV'} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'DGV'};
					
					# pour sauvegarder l'info propre aux differents callers
					$hCNV->{$type}->{$num}->{$id}->{'CALLERS'}->{$caller}++;
					$hCNV->{$type}->{$num}->{$id}->{'ELEMENTARY'}->{$caller} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'ELEMENTARY'};
					$hCNV->{$type}->{$num}->{$id}->{'GT'}->{$caller} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'GT'};
					$hCNV->{$type}->{$num}->{$id}->{'CN'}->{$caller} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'CN'};
					$hCNV->{$type}->{$num}->{$id}->{'RATIO'}->{$caller} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'RATIO'};
					$hCNV->{$type}->{$num}->{$id}->{'QUAL'}->{$caller} = $hPat_CNV->{$caller}->{$type}->{$num}->{$id}->{'QUAL'};
				}
			}
		}
	}
}
	
# pour regrouper les ids recouvrant les mêmes positions

sub gatherSV_byPosition() 
{
	foreach my $type (keys %{$hCNV})
	{
			foreach my $num (keys %{$hCNV->{$type}})
			{
				
				# 1)  detecter les SV identiques 
				my $htree->{$type}->{$num}= Set::IntervalTree->new;
				
				# remplir les arbres :  regrouper les SV chevauchants
				foreach my $id  (keys %{$hCNV->{$type}->{$num}})
				{
						my ( $t, $c, $d, $f ) = split( /_/, $id );
						$htree->{$type}->{$num}->insert($id,$d,$f);
				}
	
				# 2) associer a chaque SV ceux qui le chevauchent 	
				
				my $chr= $project->getChromosome($num);
				
				foreach my $id  (keys %{$hCNV->{$type}->{$num}})
				{	
						my ( $t, $c, $dtheSV, $ftheSV ) = split( /_/, $id );
						my $tab_id = $htree->{$type}->{$num}->fetch($dtheSV,$ftheSV);
						my @tab_id_tmp;
						
						my $lentheSV = abs($ftheSV-$dtheSV);
						
							# pour chacun des SVs chevauchant celui qui nous interresse 
							foreach my $ind_id (@$tab_id)
							{
									# On le garde si il chevauche celui qui nous interresse sur au moins 60% de sa longueur et reciproquement
									my $start1 = $hCNV->{$type}->{$num}->{$ind_id}->{"START"};
									my $end1 = $hCNV->{$type}->{$num}->{$ind_id}->{"END"};
									my $start2 = $hCNV->{$type}->{$num}->{$id}->{"START"};
									my $end2 = $hCNV->{$type}->{$num}->{$id}->{"END"};
									
									#my $same = areSameCNV($start1,$end1,$start2,$end2);
									my $same = $chr->buffer->areSameCNV($start1,$end1,$start2,$end2);
									
									push(@tab_id_tmp,$ind_id) if ($same);
							}
		
							my @tab_id_sorted  = sort ({$a cmp $b} @tab_id_tmp);
							
							my $list_of_id ="";
							foreach my $v (@tab_id_sorted)
							{
								$list_of_id .= $v.",";
							}
							
							$hSVPos->{$list_of_id}++;
							
				} # fin boucle id
			}# fin boucle num
		}# fin boucle type		
		
																				
		#3) creer l'id_global 
		my $nb = 0;
		foreach my $k ( keys( %{$hSVPos} ) ) 
		{
				my $gdeb=0;
				my $gend=0;
				my $genesList;
				my $type;
				my $num;
				my $index;

				my @tab_ind_id = split(/,/,$k);
			
				# creer les ids globaux en privilegiant manta / canvas / wisecondor
				my $canvas =0; 
				my $manta =0;
				
				foreach my $ind_id ( @tab_ind_id ) 
				{	
					my ($t,$n,$d,$e) = split( /_/, $ind_id );
				
					if ( exists($hCNV->{$t}->{$n}->{$ind_id}->{"CALLERS"}->{"manta"}) )
					{
						($type,$num,$gdeb,$gend) = split( /_/, $ind_id );
						$manta++;
					}

					unless ($manta)
					{
						if (exists($hCNV->{$t}->{$n}->{$ind_id}->{"CALLERS"}->{"canvas"}))
						{
							($type,$num,$gdeb,$gend) = split( /_/, $ind_id );
							$canvas++;
						}
						else
						{
							unless ($canvas)
							{
								($type,$num,$gdeb,$gend) = split( /_/, $ind_id );
							}
						}
					}
				}
				
				#my $pos = int($gdeb);
				my $global_id = $type . "_" . $num . "_" . $gdeb . "_" . $gend;
				my $num_int;

				#pour avoir un objet chromosome
				my $chr= $project->getChromosome($num);
								 	
				$num_int = 23 if $num eq "X";			
				$num_int = 24 if $num eq "Y";
				$num_int = int($num) unless (( $num eq "X") || ( $num eq "Y"));
				
				# enregistrer les infos dans la table de hash finale	
				$hGroupedCNV->{$global_id}->{'id'} = $global_id;			
				$hGroupedCNV->{$global_id}->{'TYPE'} = $type;
				$hGroupedCNV->{$global_id}->{'CHROM'} = $num_int;
				$hGroupedCNV->{$global_id}->{'START'} =  format_number($gdeb);
				$hGroupedCNV->{$global_id}->{'END'} =  format_number($gend);
				$hGroupedCNV->{$global_id}->{'LOCUS'} =  format_number($gdeb)."-".format_number($gend);
				$hGroupedCNV->{$global_id}->{'LEN'}  = abs( $gdeb - $gend );
			     
			   
			    # pour le plot
			    #$hGroupedCNV->{$global_id}->{'PLOT'} =  $num_int.":".$gdeb."-".$gend;
			                
				# pour les cytobandes
				my @tb;
				my @band;
				my $hband = $chr->getCytoband( $gdeb, $gend ) if ( $gend > $gdeb );

				foreach my $b ( keys %{ $hband->{'name'} } ) {
					push( @tb, $b );
				}
				@band = sort( { $a cmp $b } @tb );

				$hGroupedCNV->{$global_id}->{'CYTOBAND'} = join( ",", @band );
		
				# pour l'acces a DGV  
				my $gchr = "chr".$num;
				my $url_DGV = getDGV_url( $gchr, $gdeb, $gend );
				$hGroupedCNV->{$global_id}->{'DGV'} = $url_DGV;
		
				# retrouver les informations correspondant aux differents id regroupes 
				# ici il faut recuperer les valeurs max ou les concatenation de liste
					
				my $ggfreq = -1;
				my $glfreq = -1;
				my $mg ="no";
				my $dbVar_status = "";
				my $rankannot = 0;
				my $gt = "";
				my $cn ="";
				my $liste_genes;
				my $ds = -1;
				my $maxdejavu_inproject = -1;
				my $maxdejavu_global = -1;
				my $djprojectmax;
				my $list_other = "-";
				my $list_other_globale = "-";
				my $dj;
				my $sg = -1;
				my $bg;
				
				# pour les infos liees au caller	
				foreach my $caller (@callers)
				{
						$hGroupedCNV->{$global_id}->{$caller."_global"} = 0; 
						$hGroupedCNV->{$global_id}->{$caller."_elementary"}= 0;
				}
				
				# pour les infos liees au caller
				$hGroupedCNV->{$global_id}->{'RATIO'} = 0;
				
				foreach my $ind_id ( @tab_ind_id )
				{	

						$ggfreq = $hCNV->{$type}->{$num}->{$ind_id}->{'GOLD_G_FREQ'} if ($hCNV->{$type}->{$num}->{$ind_id}->{'GOLD_G_FREQ'} > $ggfreq);
						$glfreq = $hCNV->{$type}->{$num}->{$ind_id}->{'GOLD_L_FREQ'} if ($hCNV->{$type}->{$num}->{$ind_id}->{'GOLD_L_FREQ'} > $glfreq);
					
						$mg = "yes" if ($hCNV->{$type}->{$num}->{$ind_id}->{'OMIN_MG'} eq "yes");
						$dbVar_status .= $hCNV->{$type}->{$num}->{$ind_id}->{'dbVar_status'}." " unless ($dbVar_status =~ m/$hCNV->{$type}->{$num}->{$ind_id}->{'dbVar_status'}/);
						$rankannot = $hCNV->{$type}->{$num}->{$ind_id}->{'RANKAnnot'} if ($hCNV->{$type}->{$num}->{$ind_id}->{'RANKAnnot'} > $rankannot);			
						$ds = $hCNV->{$type}->{$num}->{$ind_id}->{'DUPSEG'} if ($hCNV->{$type}->{$num}->{$ind_id}->{'DUPSEG'} > $ds);
					
					
						foreach my $caller (keys %{$hCNV->{$type}->{$num}->{$ind_id}->{'CALLERS'}})
						{
							$hGroupedCNV->{$global_id}->{$caller."_global"} = $ind_id;
							$hGroupedCNV->{$global_id}->{$caller."_elementary"}= $hCNV->{$type}->{$num}->{$ind_id}->{'ELEMENTARY'}->{$caller};		
							$hGroupedCNV->{$global_id}->{'GT'}.= $hCNV->{$type}->{$num}->{$ind_id}->{'GT'}->{$caller}." "  unless ( $hGroupedCNV->{$global_id}->{'GT'} =~ m/$hCNV->{$type}->{$num}->{$ind_id}->{'GT'}->{$caller}/ );
							$hGroupedCNV->{$global_id}->{'CN'}.= $hCNV->{$type}->{$num}->{$ind_id}->{'CN'}->{$caller}." " unless ( $hGroupedCNV->{$global_id}->{'CN'} =~ m/$hCNV->{$type}->{$num}->{$ind_id}->{'CN'}->{$caller}/ );
							$hGroupedCNV->{$global_id}->{'RATIO'} = $hCNV->{$type}->{$num}->{$ind_id}->{'RATIO'}->{$caller} unless ( $hCNV->{$type}->{$num}->{$ind_id}->{'RATIO'}->{$caller} == 0);						
							$hGroupedCNV->{$global_id}->{'QUAL'} .= $caller.":".$hCNV->{$type}->{$num}->{$ind_id}->{'QUAL'}->{$caller}."/";
						}
					
				
						# pour le dejavu in this project
						my ($dj,$dj_list) = getDejavuAndTransmission($type,$num,$ind_id,$global_id);   
					
						 if ($dj > $maxdejavu_inproject)
						{
						 	$maxdejavu_inproject = $dj;
						 	$list_other = $dj_list if($dj);
						}
					 
						 # pour le dejavu global
						 my ($djprojects,$djpatients,$dj_list_globale) = getDejavuGlobal($type,$num,$ind_id,$global_id);   
					
						if ($djpatients > $maxdejavu_global)
						{
					 		$maxdejavu_global = $djpatients;
					 		$djprojectmax = $djprojects;
					 		$list_other_globale = $dj_list_globale;
						}
			} #fin de la boucle sur les ind_id
					
			# pour le score_qual
			my $score_qual = getScoreQual($global_id);
			$hGroupedCNV->{$global_id}->{'QUAL'} .= " ".$score_qual;
			
			$hGroupedCNV->{$global_id}->{'GOLD_G_FREQ'} = $ggfreq;
			$hGroupedCNV->{$global_id}->{'GOLD_L_FREQ'} = $glfreq;
			$hGroupedCNV->{$global_id}->{'OMIN_MG'} = $mg;
			$hGroupedCNV->{$global_id}->{'dbVar_status'} = $dbVar_status;
			$hGroupedCNV->{$global_id}->{'RANKAnnot'} = $rankannot;
			$hGroupedCNV->{$global_id}->{'DUPSEG'} = $ds;
			
			
			$hGroupedCNV->{$global_id}->{'DEJAVU_P'} = $maxdejavu_inproject.";".$list_other;
			$hGroupedCNV->{$global_id}->{'DEJAVU_G'} = $maxdejavu_inproject.";".$list_other.";".$djprojectmax.";".$maxdejavu_global.";".$list_other_globale;
			
			
			#dans les genotypes
			my $liste_finale2=" ";
			my @tabgt = split(/ /,$hGroupedCNV->{$global_id}->{'GT'});
				
			foreach my $gt (@tabgt)
			{
				next  if ($liste_finale2 =~ m/$gt/);
				$liste_finale2 .= $gt." ";
			}
			
			$hGroupedCNV->{$global_id}->{'GT'} = $liste_finale2;
			
								
	} # fin boucle sur k = liste des ind_id regroupes dans un global_id
}


sub getDejavuAndTransmission
{
		my ($type,$num,$ind_id,$global_id) = @_;
	
		# pour le dejavu et la transmission
		$mothername = "none";
		$fathername = "none";
		my $famillymembers;
		my $famillymembersNames;
					
		if ($thePatient->isChild())
		{
			$mothername = $thePatientFamily->mother() if ( $thePatientFamily->mother() );
			$fathername = $thePatientFamily->father() if ( $thePatientFamily->father() );
		}
		
		my $chr= $project->getChromosome($num);
					
		my $transmission=0;
		my $nbdejavuProject = 0;
		my $list_of_other_patient = "";
		
		my ($t,$c,$start1,$end1) = split(/_/,$ind_id);
					
		# on recherche dans le dejavu des  CNV chevauchants
		my $tab_id = $htree_dejavuProject->{$type}->{$num}->fetch($start1,$end1);
		

		# puis pour chacun d'eux on regarde ceux qui sont identiques au sens de areSameCNV
		foreach my $djv_id (@$tab_id)
		{
				my ($t,$c,$start2,$end2) = split(/_/,$djv_id);
				
				my $identity = $chr->buffer->areSameCNV($start1,$end1,$start2,$end2);
				
				if ( $identity )
				{
						$transmission=0;
						
						#dejavu ou transmission ?
						foreach my $pname (keys %{$hdejavuProject->{$type}->{$num}->{$djv_id}})
						{
								next if ($pname eq $thePatientName);
								
								if ($pname eq $mothername)
								{
										$transmission += 2;
								}
								else
								{
										if ($pname eq $fathername)
										{
												$transmission += 1;
										}
										else
										{
											unless($list_of_other_patient =~ m/$pname/ )
											{
												$list_of_other_patient .= $pname.",";
												$nbdejavuProject++;
											}
										}
								}
						}
					
						if ($thePatientFamily->isTrio() && $thePatient->isChild() )
						{
								$hGroupedCNV->{$global_id}->{'TRANSMISSION'} .= "both  " if (($transmission == 3) && ($hGroupedCNV->{$global_id}->{'TRANSMISSION'} !~ m/both/));
								$hGroupedCNV->{$global_id}->{'TRANSMISSION'} .= "mother " if (($transmission == 2) && ($hGroupedCNV->{$global_id}->{'TRANSMISSION'} !~ m/mother/));
								$hGroupedCNV->{$global_id}->{'TRANSMISSION'} .= "father " if (($transmission == 1) && ($hGroupedCNV->{$global_id}->{'TRANSMISSION'} !~ m/father/));
								$hGroupedCNV->{$global_id}->{'TRANSMISSION'} .= "denovo " if (($transmission == 0) && ($hGroupedCNV->{$global_id}->{'TRANSMISSION'} !~ m/denovo/));
						}
						else 
						{
								$hGroupedCNV->{$global_id}->{'TRANSMISSION'} = "-";
						}
				}
		}
		return ($nbdejavuProject,$list_of_other_patient);
}

sub getDejavuGlobal
{
		my ($type,$num,$ind_id,$global_id) = @_;
		
		my $nbproject=0;
		my $nbpatient=0;
		my $list_of_other_patient;
		
		my $hres;
		my $hresnew;
		
		
		my ($t,$c,$start1,$end1) = split(/_/,$ind_id);

		my $chr = $project->getChromosome($c);
			
		# on recherche dans le dejavuGlobal des  CNV chevauchants		
		my $tab_id = $chr->retreiveDejaVuIdSV($start1,$end1,$type);			
		
		# puis pour chacun d'eux on regarde ceux qui sont identiques au sens de areSameCNV
		foreach my $djv_id (@$tab_id)
		{
				my ($t,$c,$start2,$end2) = split(/_/,$djv_id);
				
				#warn $ind_id."       ".$djv_id;
				
				my $identity = $chr->buffer->areSameCNV($start1,$end1,$start2,$end2);
				
				if ( $identity )
				{	
						
						#dejavu new version
						my $hashdv = $lmdb->get($djv_id);
						
						foreach my $project_name (keys %{$hashdv})
						{
							next if ($project_name eq $TheProjectName);
							foreach my $pname (keys %{$hashdv->{$project_name}})
							{
								my $etiq = $project_name.":".$pname;
								$list_of_other_patient .= $project_name.":".$pname."," unless($list_of_other_patient =~ m/$etiq/ );
								$hres->{$project_name}->{$pname}++;
							}
						}
				}
		}
		

		
		my @list = split(/,/,$list_of_other_patient);
		my @sorted_list = sort @list;
		my $theliste;
		
		foreach my $n (@sorted_list)
		{
			$theliste .= $n.","; 
		}
		
		foreach my $proj ( keys %{$hres})
		{
			$nbproject++;
			foreach my $pat (keys %{$hres->{$proj}})
			{
				$nbpatient++;
			}
		}
		
		return ($nbproject,$nbpatient,$theliste);
}




sub getInfoCallers() 
{
	 foreach my $gid  ( keys %{ $hGroupedCNV })
	{
		        my $score = 0;
		
				my $wc = 1;
				my $can = 1;
				my $man = 1;
				
				my ($qual,$scorequal) = split(/ /,$hGroupedCNV->{$gid}->{'QUAL'});
				
				my @tabval = split(/\//,$qual);
				
				my $qwc = 1;
				my $qc = 1;
				my $qm = 1;
				
				my $scoretransmission;
				
				foreach my $val (@tabval)
				{
					my ($k,$v) =split(/:/,$val);
					
					$qwc = 0.75 if ( ($k eq "wisecondor") && ( abs($v) < 40) );
					$qc = 0.75 if ( ($k eq "canvas") && ( abs($v) < 20) );
					$qm = 0.75 if ( ($k eq "manta") && ( abs($v) < 400) );
				}
				
				if ( ( $hGroupedCNV->{$gid}->{"wisecondor_global"} ) && $wc) 
				{ 
						$score += 4*$qwc; $wc=0;
						$hGroupedCNV->{$gid}->{"wisecondor_elementary"} .= ";".$qwc;
				}
				
				if ( ( $hGroupedCNV->{$gid}->{"canvas_global"} ) && $can) 
				{ 
						$score += 2*$qc; $can=0;
						$hGroupedCNV->{$gid}->{"canvas_elementary"} .= ";".$qc;
				}
				
				if ( ( $hGroupedCNV->{$gid}->{"manta_global"} ) && $man) 
				{ 
					$score += 1*$qm; $man=0;
					$hGroupedCNV->{$gid}->{"manta_elementary"} .= ";".$qm;
				}
				
				$hGroupedCNV->{$gid}->{'SCORECALLER'} = $score;
				$hGroupedCNV->{$gid}->{'SCORECNV'} = getScoreCNV($gid);
				
				my $bamFiles = $path_Bam."/".$thePatientName.".bam";
				$bamFiles .= ",".$path_Bam."/".$mothername.".bam" if ($mothername ne "none");
				$bamFiles .= ",".$path_Bam."/".$fathername.".bam" if ($fathername ne "none");
				
				my $bamNames = $thePatientName;
				$bamNames .= ",".$mothername if ($mothername ne "none");
				$bamNames .= ",".$fathername if ($fathername ne "none");
						
				# pour IGV on garde preferentielement les infos de manta puis de canvas et a defaut de wisecondor
				if ( $hGroupedCNV->{$gid}->{"manta_global"} )
				{
						$hGroupedCNV->{$gid}->{"IGV"} = $bamNames.";".$bamFiles.";".$hGroupedCNV->{$gid}->{"manta_global"};
				}
				else 
				{
						if ( $hGroupedCNV->{$gid}->{"canvas_global"} )
		   				{
								$hGroupedCNV->{$gid}->{"IGV"} = $bamNames.";".$bamFiles.";".$hGroupedCNV->{$gid}->{"canvas_global"};
						}
						else 
						{
								$hGroupedCNV->{$gid}->{"IGV"} =  $bamNames.";".$bamFiles.";".$hGroupedCNV->{$gid}->{"wisecondor_global"};
						}
				}
		}
}

sub filtreCNV
{
	
	foreach my $global_id ( keys %{$hGroupedCNV} ) 
	{
		my($type,$num,$gdeb,$gend) = split(/_/,$global_id);
		my $SV_ok = 1;

		# filtre sur le scoreCNV
		# pour les bestones on garde ceux dont le score est > a 2.75 
		next  if ( ( $hGroupedCNV->{$global_id}->{'SCORECNV'} < 2.75 ) && $select_best );

		# filtre sur les chromosomes
		next if ( ( $hGroupedCNV->{$global_id}->{'CHROM'} != $chr ) && ( $chr ne "all" )&& ( $chr ne "noXY" ) );
		next if ( ( $hGroupedCNV->{$global_id}->{'CHROM'} == 23 )&& ( $chr eq "noXY" ) );
		next if ( ( $hGroupedCNV->{$global_id}->{'CHROM'} == 24 )&& ( $chr eq "noXY" ) );

		# filtre sur la taille 
		next if ( $hGroupedCNV->{$global_id}->{"LEN"} < $minlength );
		unless ( $maxlength eq "nomax")
		{
			next if ( $hGroupedCNV->{$global_id}->{"LEN"} > $maxlength );
		}
		
		# filtre sur les frequences
		unless ( $maxfreq eq "nomax" ) {
			next if ( $hGroupedCNV->{$global_id}->{'GOLD_G_FREQ'} > $maxfreq );
			next if ( $hGroupedCNV->{$global_id}->{'GOLD_L_FREQ'} > $maxfreq );
		}

		# filtre sur le genotype
		my $pass = 1;
		my $infoGenotype = $hGroupedCNV->{$global_id}->{"GT"};

		unless ( $genotype eq "both" )
		{
			$pass = 0;

			if (( $genotype eq "ho" ) && (   ( $infoGenotype =~ m/1\/1/ ) || ( $infoGenotype =~ m/1\/2/ ) || ( $infoGenotype =~ m/-/) ) )
			{
				$pass = 1;
			}
			
			if (   ( $genotype eq "he" ) && ( ( $infoGenotype =~ m/0\/1/ ) || ( $infoGenotype =~ m/-/ ) ) )
			{
				$pass = 1;
			}
		}
		next unless ($pass);
		
		# filtre sur la transmission
		$pass = 1;
		my $infoTransmission = $hGroupedCNV->{$global_id}->{"TRANSMISSION"};
	
		unless ( $transmission  eq "all" )
		{
			
			$pass = 0;
			if ( ($transmission eq "mother") &&  ($infoTransmission =~ m/mother/ ) &&  !($infoTransmission =~ m/father/ ) )
			{
				$pass = 1;
			}
			
			if ( ($transmission eq "father") &&  ($infoTransmission =~ m/father/ ) &&  !($infoTransmission =~ m/mother/ ) )
			{
				$pass = 1;
			}
			
			if ( ($transmission eq "both") && ( ($infoTransmission =~ m/both/ ) || (($infoTransmission =~ m/father/ ) &&  ($infoTransmission =~ m/mother/ ) )))
			{
				$pass = 1;
			}
			
			if ( ($transmission eq "denovo") &&  ($infoTransmission =~ m/denovo/ ) &&  !($infoTransmission =~ m/father/ ) &&  !($infoTransmission =~ m/mother/ ) &&  !($infoTransmission =~ m/both/ ) )
			{
				$pass = 1;
			}
		}
		next unless ($pass);

		# filtre sur le dejavu :
		my ($djv_inThisProject,$rien) = split(/;/,$hGroupedCNV->{$global_id}->{"DEJAVU_P"});
		my ($djvp,$list1,$djv_otherProject,$djv_inOtherProject,$list2) = split(/;/,$hGroupedCNV->{$global_id}->{"DEJAVU_G"});
		my $nbotherWG = $djv_inThisProject + $djv_inOtherProject;
		
		unless($dejavu eq "all")
		{next if (  $nbotherWG > $dejavu)};
		
		# filtre sur une liste de cytobandes
		my @tabcytoband = split(/,/,$cytoband);
		$pass = 0;
		foreach my $cb (@tabcytoband) 
		{
			if (   ( $hGroupedCNV->{$global_id}->{"CYTOBAND"} =~ m/$cb/ )|| ( $cb eq "all" ) )
			{
				$pass = 1;
				next;
			}
		}
		next unless ($pass);
		
		$hGroupedCNV->{$global_id}->{"TYPE"} = "TRI" if (($hGroupedCNV->{$global_id}->{"TYPE"} eq "DUP") && ($hGroupedCNV->{$global_id}->{"RATIO"} > 0.7 ) );
		$hGroupedCNV->{$global_id}->{"TYPE"} = "QUA" if (($hGroupedCNV->{$global_id}->{"TYPE"} eq "DUP") && ($hGroupedCNV->{$global_id}->{"RATIO"} > 1 ) );
		
		# au final on calcul pour tous les CNV restant le ratio de couverture ...
		my ($t,$chr,$start,$end) = split(/_/,$global_id);
		my $score_dude = -1;
		#my $score_dude = $thePatient->cnv_value_dude($chr,$start,$end);
		
		$hGroupedCNV->{$global_id}->{"SCORE_DUDE"}=$score_dude;
		
		# au final on cherche les genes et ne calcule le score genes que pour les CNV selectionnés
		#Pour le score_genes
		my $chr = $project->getChromosome($num);
		my $tabGenes = $chr->getGenesByPosition($gdeb,$gend);

		# calculer le score gene max correspondant a la liste de gene associe au variant
		my $max = 0;
		my @names;
		if ( scalar(@$tabGenes) )  # si le variant recouvre des gènes
		{	
				foreach my $g (sort {$b->score <=> $a->score} @$tabGenes){
							$max = $g->score unless $max;
							push (@names,$g->external_name);
				}
		
				$hGroupedCNV->{$global_id}->{'SCORE_GENES'}=$max;
				$hGroupedCNV->{$global_id}->{'BEST_GENE'}= $names[0];
				$hGroupedCNV->{$global_id}->{'GENES'} =join(" ",@names);

			# filtre sur une liste de genes
			my @tabgenes_annot = split(/,/,$hGroupedCNV->{$global_id}->{"GENES"} );
			$pass   = 0;
			my @tabgenes = split(/ /,$listOfGenes);
			foreach my $gene_annot (@tabgenes_annot) {
				foreach my $gene (@tabgenes) {
					if ( ( $gene_annot eq $gene ) || ( $gene eq "all" ) ) {
						$pass = 1;
							last;
					}
				}
				last if ( $pass == 1 );
			}
			next unless ($pass);
		}
		else
		{
				$hGroupedCNV->{$global_id}->{'SCORE_GENES'}="-";
				$hGroupedCNV->{$global_id}->{'BEST_GENE'}= "-";
				$hGroupedCNV->{$global_id}->{'GENES'} ="-";
		}
		push( @Filtered_listHashRes, { %{ $hGroupedCNV->{$global_id} } } );
	}# fin de la boucle sur les global_id

}

sub getDGV_url {
	my ( $chr, $d, $f ) = @_;
	my $url = "http://dgv.tcag.ca/gb2/gbrowse/dgv2_hg19/?name=" . $chr . "%3A". $d . "-"  . $f  . ";search=Search";
	return $url;
}

sub getScoreCNV{
	
	my ($gid ) = @_;
	my $score = $hGroupedCNV->{$gid}->{'SCORECALLER'};
		
	#pour présenter en premier les CNV denovo
	if ( ($hGroupedCNV->{$gid}->{'TRANSMISSION'}  =~ m/denovo/) && !($hGroupedCNV->{$gid}->{'TRANSMISSION'} =~ m/mother/ ) && !($hGroupedCNV->{$gid}->{'TRANSMISSION'}  =~ m/father/ ) && !($hGroupedCNV->{$gid}->{'TRANSMISSION'}  =~ m/both/) )
	{
			$score += 0.6; 
	}
				
	# pour rétrograder le CNV situes au niveau de segmental duplication (chevauchement de plus de 40%)
	if ($hGroupedCNV->{$gid}->{'DUPSEG'} > 40)
	{
			my $sd_score= int($hGroupedCNV->{$gid}->{'DUPSEG'}/40); 
			$score -= $sd_score;		#	(de -1 à -2,5 )
	}
	
	return $score;
}

sub getScoreQual {
	my ($global_id ) = @_;

	my @tval = split(/\//,$hGroupedCNV->{$global_id}->{'QUAL'});
	
	my $scorequal;
	my $nbcaller=0;
	
	
	foreach my $val (@tval)
	{
		
		my ($c,$q) = split(/:/,$val);

		if ($c eq "wisecondor")
		{
			if ($q)
			{
				$nbcaller++;
				if ( abs($q) >= 40)
				{
					$scorequal++;
				}
			}
		}
		
		if ($c eq "canvas")
		{
			if ($q)
			{
				$nbcaller++;
				if  ($q >= 20)
				{
					$scorequal++;
				}
			}
		}
		
		if ($c eq "manta")
		{
			if ($q)
			{
				$nbcaller++;
				if ( $q >= 400 )
				{
					$scorequal++;
				}
			}
		}
	}
	
	$scorequal /= $nbcaller;
	
	return $scorequal; 
}
