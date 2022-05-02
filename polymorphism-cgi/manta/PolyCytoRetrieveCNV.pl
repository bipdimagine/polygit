#!/usr/bin/perl
$|=1;

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
use Digest::MD5::File qw(dir_md5_hex file_md5_hex url_md5_hex file_md5);
use GBuffer;
use GenBoProject;
use GenBoCache;
use Time::HiRes;
use layout;
use export_excel;
use export_data;
use Capture::Tiny ':all';

###########################################################################################################
#  1) recupere les infos concernant tous les SV du projet entre en argument et prealablement freeze : nom du fihier = project.allSV
#  2) regroupe les SV des patients entres en argument et qui ont la même liste de gènes sous un id_global
#  3) Si le project contient des trios on regarde la transmission a partir du SVglobal
#  4) filtre le resultats en fonction des donnéees definies via l'interface
#  5) construit le Json resultant pour affichage
#
############################################################################################################
my $tt = 0;
my $cgi = new CGI;
my $TheProjectName = $cgi->param('project');
my $thePatientName     = $cgi->param('patient');
my $minlength   = $cgi->param('minlength');
my $select_best = $cgi->param('select_best');
my $chrom        = $cgi->param('chrom');
my $listOfGenes = $cgi->param('genes');
my $dejavu      = $cgi->param('dejavu');
my $transmission    = $cgi->param('transmission');
my $omim      = $cgi->param('omim');
my $print = $cgi->param('print');
my $force= $cgi->param( 'force' ) == 1 ;

#$force = 1;

# ceux qui ne sont plus utilisés pour filtrer
my $maxlength   = $cgi->param('maxlength');
my $maxfreq     = $cgi->param('maxfreq');
my $cytoband    = $cgi->param('cytoband');
my $genotype    = $cgi->param('genotype');


# pour construire le projet
my $buffer  = GBuffer->new();
my $project = $buffer->newProjectCache( -name => $TheProjectName );
my $thePatient = $project->getPatient($thePatientName);

my %hkeys = $cgi->Vars;
my @keys;
my $string;


foreach my $k  (sort {$a cmp $b} keys %hkeys){
	next if $k =~ /force/;
	next if $k =~ /user/;
	push(@keys,"$k");
	my $c = $hkeys{$k};
	$c =~ s/\"//g;
	$c =~ s/\+/ /g;
	push(@keys,$c);
}

push(@keys,file_md5_hex($Bin."/PolyCytoRetrieveCNV.pl") );
my $cache_id= join(";",@keys).";"."polycyto_cnv";

my $no_cache;
my $text;
#if ($force == 1) {
#	$no_cache = $thePatient->get_lmdb_cache("c");
#	$no_cache->put("toto","titi");
#	$no_cache->close();
#}
#else {
	$no_cache = $thePatient->get_lmdb_cache("r");
	$text = $no_cache->get_cache($cache_id);
	#warn "----->".$text;
	$no_cache->close();
#}	

$| = 1;
if ($text){
	print $text;
	exit(0);
}
$no_cache->close();
my $stdout = tee_stdout {
print $cgi->header('text/json-comment-filtered');

print "{\"progress\":\".$cache_id";
};

# les seuils pour la comparaison des variants
#my $seuilSameEvent = 60;
#my $seuilDejaVu = 75;
#my $seuilTransmited = 85;

my $seuilSameEvent = 70;
my $seuilDejaVu = 90;
my $seuilTransmited = 90;

# pour stocker les SV
my $hSVPos;          		# tous par positions chevauchantes

# gestion des callers
my @callers=("wisecondor","canvas","manta");


$project->cgi_object(1);

# pour distinguer parmi les autres patients du projet
# ceux qui sont de la même famille

my $thePatientFamily = $thePatient->getFamily();
my $mothername;
my $fathername;

my $trio = 0;
$trio = 1  if $thePatientFamily->isTrio();

# pour la gestion des parents si ils existent
if ($trio && ($thePatient->isChild()) )
{
		$mothername = $thePatientFamily->mother() if ( $thePatientFamily->mother() );
		$fathername = $thePatientFamily->father() if ( $thePatientFamily->father() );
}

# pour construire le path d'acces au bam pour IGV

my $bam_dir  = $thePatient->getProject->getAlignmentUrl($thePatient->alignmentMethod);
	


# pour recuperer les BND Manta associes au CNV
my $pathBND = $project->getSVeqDir();
my $path_djv = $project->DejaVuSVeq_path;

my $PatientBNDfile = $pathBND.$thePatientName.".allBND";
my $MotherBNDfile = $pathBND.$mothername.".allBND" if $mothername;
my $FatherBNDfile = $pathBND.$fathername.".allBND" if $fathername;

my $TranslocDejavufile = $path_djv."SVeqDejavu.all";

my $hTransLoc1;
my $hdejavuBNDProject;
my $htranslocdejavu;
my $bpok = 1;

if (-f $PatientBNDfile)
{
	$hTransLoc1->{"patient"} = retrieve($PatientBNDfile) or die "Can't retrieve datas from " . $PatientBNDfile . " !\n";
}
else
{
	$bpok = 0; 	# les translocs ne sont pas traitées
}

if($bpok)
{
	# pour acceder au  dejavu global des transloc
	$htranslocdejavu = retrieve($TranslocDejavufile) or die "Can't retrieve datas from " . $TranslocDejavufile . " !\n";

	# gestion du dejavu des parents
	if ($trio)
	{
		if (-f $MotherBNDfile){$hTransLoc1->{"mother"} = retrieve($MotherBNDfile) or die "Can't retrieve datas from " . $MotherBNDfile . " !\n";}
		if (-f $FatherBNDfile){$hTransLoc1->{"father"} = retrieve($FatherBNDfile) or die "Can't retrieve datas from " . $FatherBNDfile . " !\n";}
	}
}


foreach my $who ( keys %{$hTransLoc1})
{
	foreach my $event_id (keys %{$hTransLoc1->{$who}})
	{
		my ($chr1, $bp1, $chr2, $bp2) = split(/_/, $event_id);
		$hdejavuBNDProject->{$who}->{$chr1}->{$event_id} = 1;
		$hdejavuBNDProject->{$who}->{$chr2}->{$event_id} = 1; 
	}
}

# pour accelerer acces = intervall tree
my $htree_dejavuBNDProject;
foreach my $who (keys %{$hdejavuBNDProject})
{
		foreach my $num (keys %{$hdejavuBNDProject->{$who}})
		{
				# 1)  detecter les SV identiques 
				$htree_dejavuBNDProject->{$who}->{$num}= Set::IntervalTree->new;
				
				# remplir les arbres :  regrouper les SV chevauchants
				foreach my $id  (keys %{$hdejavuBNDProject->{$who}->{$num}})
				{
						my ( $c1, $bp1, $c2, $bp2 ) = split( /_/, $id );
						
						$htree_dejavuBNDProject->{$who}->{$c1}->insert($id,$bp1-5000,$bp1+5000) if ($c1 == $num);
						$htree_dejavuBNDProject->{$who}->{$c2}->insert($id,$bp2-5000,$bp2+5000) if ($c2 == $num);
				}
		}
}	


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
my $hFilteredCNV;
my @Filtered_listHashRes;

#####################################################
# pour le dejavu_in_this_project (a part pour etre a jour sans avoir a lancer le setdejavu)
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

my $dejavudir = $project->DejaVuCNV_path();
my $lmdb = GenBoNoSqlLmdb->new(dir=>$dejavudir,mode=>"r",name=>"dejavu_sv",is_compress=>1);	

print '.';
# 1) on regroupe les ids strictement identiques
gather_strict_identicalCNV();

# 2) on regroupe les CNV identiques vus par differents callers 
gatherSV_byPosition();

# 3) on filtre en fonction des infos envoyees par l'interface
filtreCNV();

# 4) resultat
my $stdout2 = tee_stdout {
printJson( \@Filtered_listHashRes );
};

$no_cache = $thePatient->get_lmdb_cache("w");
$no_cache->put_cache($cache_id,$stdout.$stdout2,2400);
$no_cache->close();


exit(0);


###############################################################################
#
#	methodes
#
#################################################################################


sub printJson {
	my ($listHash) = @_;
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
	
	
	#print encode_json $hash;
	
	my $json_encode = encode_json $hash;
	print ".\",";
	$json_encode =~ s/{//;
	print $json_encode;
	
	
	print "\n"  unless ($print == 1);
}

# pour regrouper les ids strictement identiques
sub gather_strict_identicalCNV
{
	print "cocuou";
	foreach my $caller (keys %{$hPat_CNV})
	{
		print $caller;
		foreach my $type (keys %{$hPat_CNV->{$caller}})
		{
			print $type;
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
			print '*';
			foreach my $num (keys %{$hCNV->{$type}})
			{
				#pour avoir un objet chromosome
				my $chr= $project->getChromosome($num);
				
				# 1)  detecter les SV identiques 
				my $htree->{$type}->{$num}= Set::IntervalTree->new;
				
				# remplir les arbres :  regrouper les SV chevauchants
				foreach my $id  (keys %{$hCNV->{$type}->{$num}})
				{
						my ( $t, $c, $d, $f ) = split( /_/, $id );
						$htree->{$type}->{$num}->insert($id,$d,$f);
				}
				
				# 2) associer a chaque SV ceux qui le chevauchent 	
				foreach my $id  (keys %{$hCNV->{$type}->{$num}})
				{	
						my ( $t, $c, $dtheSV, $ftheSV ) = split( /_/, $id );
						my $tab_id = $htree->{$type}->{$num}->fetch($dtheSV,$ftheSV);
						my @tab_id_tmp;
						
						my $lentheSV = abs($ftheSV-$dtheSV);
						
							# pour chacun des SVs chevauchant celui qui nous interresse 
							foreach my $ind_id (@$tab_id)
							{
									# On le garde si il chevauche celui qui nous interresse sur au moins 70% de sa longueur et reciproquement
									my $start1 = $hCNV->{$type}->{$num}->{$ind_id}->{"START"};
									my $end1 = $hCNV->{$type}->{$num}->{$ind_id}->{"END"};
									my $start2 = $hCNV->{$type}->{$num}->{$id}->{"START"};
									my $end2 = $hCNV->{$type}->{$num}->{$id}->{"END"};
									
									my $same = 100;
									#TODO: Pq les deux bornes en meme temps ?
									$same = $project->dejavuSV->getIdentityBetweenCNV($start1,$end1,$start2,$end2) unless (($start1==$start2) && ($end1==$end2));
									push(@tab_id_tmp,$ind_id) if ($same >= $seuilSameEvent);
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
		print 'PP';
		foreach my $k ( keys( %{$hSVPos} ) ) 
		{
				$project->print_dot(1);
				my $gdeb;
				my $gend;
				my $genesList;
				my $type;
				my $num;
				my $index;
				
				my @tab_ind_id = split(/,/,$k);
				
				# pour le dejavu global		
				my $flag=0;
				
				# creer les ids globaux en privilegiant manta / canvas / wisecondor
				my $canvas =0; 
				my $manta =0;
				
				foreach my $ind_id ( @tab_ind_id ) 
				{	
					my ($t,$n,$d,$e) = split( /_/, $ind_id );
					$flag = 1 if ( exists($hCNV->{$t}->{$n}->{$ind_id}->{"CALLERS"}->{"wisecondor"}) );
					
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
				
				my $global_id = $type . "_" . $num . "_" . $gdeb . "_" . $gend;
				
				my $num_int;
				$num_int = 23 if $num eq "X";			
				$num_int = 24 if $num eq "Y";
				$num_int = int($num) unless (( $num eq "X") || ( $num eq "Y"));
				
				# filtre sur les chromosomes
				next if ( ( $num_int != $chrom) && ( $chrom ne "all" )&& ( $chrom ne "noXY" ) );
				next if ( ( $num_int == 23 )&& ( $chrom eq "noXY" ) );
				next if ( ( $num_int == 24 )&& ( $chrom  eq "noXY" ) );
		
		
				# filtre sur la taille 
				my $gid_len = abs( $gdeb - $gend );
				next if ( $gid_len < $minlength );
				unless ( $maxlength eq "nomax")
				{
						next if ( $gid_len > $maxlength );
				} 
				
	
				# pour le dejavu global
				my ($dj_otherprojects,$dj_patients_otherproject,$nbw,$nbc,$nbm,$scorecaller,$dj_list_globale) = getDejavuGlobal($type,$num,$global_id,$flag);  

				# pour le dejavu in this project
				my ($dj_patients_inproject,$list_patient_inproject,$theTransmission) = getDejavuAndTransmission($type,$num,$global_id);   
				
	
				# enregistrer les infos dans la table de hash finale	
				$hGroupedCNV->{$global_id}->{'id'} = $global_id;			
				$hGroupedCNV->{$global_id}->{'TYPE'} = $type;
				$hGroupedCNV->{$global_id}->{'CHROM'} = $num_int;
				$hGroupedCNV->{$global_id}->{'START'} =  format_number($gdeb);
				$hGroupedCNV->{$global_id}->{'END'} =  format_number($gend);
				$hGroupedCNV->{$global_id}->{'LOCUS'} =  format_number($gdeb)."-".format_number($gend);
				$hGroupedCNV->{$global_id}->{'LEN'}  = $gid_len;
			    $hGroupedCNV->{$global_id}->{'TRANSMISSION'}  = $theTransmission;
			    $hGroupedCNV->{$global_id}->{'DEJAVU_P'} = $dj_patients_inproject.";".$list_patient_inproject;
				$hGroupedCNV->{$global_id}->{'DEJAVU_G'} = $global_id."+".$dj_patients_inproject.";".$list_patient_inproject.";".$dj_otherprojects.";".$dj_patients_otherproject.";".$dj_list_globale.";".$nbw.";".$nbc.";".$nbm.";".$flag;
			    
			   
			   #pour avoir un objet chromosome
				my $chr= $project->getChromosome($num);
				my $chrLength= $chr->length();
				
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
				
				# pour l'accès à gnomAD
				my $url_gnomAD = getgnomAD_url($num,$gdeb,$gend);
				$hGroupedCNV->{$global_id}->{'gnomAD'} = $url_gnomAD;
				
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

				} #fin de la boucle sur les ind_id


			# filtre sur le dejavu

			unless ($dejavu eq "all")
			{
				unless ( $hGroupedCNV->{$global_id}->{"wisecondor_elementary"} )
				{
						next if (  ($dj_patients_inproject + $dj_patients_otherproject) > $dejavu ) ;
				}	
			}
			
			# pour le score_qual
			my $score_qual = getScoreQual($global_id);
			$hGroupedCNV->{$global_id}->{'QUAL'} .= " ".$score_qual;
			
			$hGroupedCNV->{$global_id}->{'GOLD_G_FREQ'} = $ggfreq;
			$hGroupedCNV->{$global_id}->{'GOLD_L_FREQ'} = $glfreq;
			$hGroupedCNV->{$global_id}->{'OMIN_MG'} = $mg;
			$hGroupedCNV->{$global_id}->{'dbVar_status'} = $dbVar_status;
			$hGroupedCNV->{$global_id}->{'RANKAnnot'} = $rankannot;
			$hGroupedCNV->{$global_id}->{'DUPSEG'} = $ds;
					
			#pour le genotype
			my $liste_finale2=" ";
			my @tabgt = split(/ /,$hGroupedCNV->{$global_id}->{'GT'});
				
			foreach my $gt (@tabgt)
			{
				next  if ($liste_finale2 =~ m/$gt/);
				$liste_finale2 .= $gt." ";
			}
			$hGroupedCNV->{$global_id}->{'GT'} = $liste_finale2;
			
			 # pour le plot
			 my $transmis ="-";
			 my $tr =  $hGroupedCNV->{$global_id}->{'TRANSMISSION'};
			 $transmis = "denovo" if ($tr =~ m/denovo/);
			 $transmis = "mother" if ( (($tr =~ m/mother/) || ($tr =~ m/maybeM/)) && ($tr !~ m/father/) );
			 $transmis = "father" if ( (($tr =~ m/father/) || ($tr =~ m/maybeF/)) && ($tr !~ m/mother/) );
			 $transmis = "both" if ( ($tr =~ m/both/) || ( (($tr =~ m/mother/) || ($tr =~ m/maybM/))  && (($tr =~ m/father/) || ($tr =~ m/maybF/))));
			 
			 $gend -=5000;
			 $hGroupedCNV->{$global_id}->{'PLOT'} =  $num_int.":".$gdeb."-".$gend.";".$transmis; 
			 
			 # pour igv
			my $bamFiles = $bam_dir.$thePatientName.".bam";
			my $bamNames = $thePatientName;
			if ($thePatientFamily->isTrio() )
			{
				my $members = $thePatientFamily->getMembers();

				foreach my $m (@$members)
				{
						my $membername = $m->name();
						$bamFiles .= ",".$bam_dir.$membername.".bam" unless ($membername eq $thePatientName);
						$bamNames .= ",".$membername unless ($membername eq $thePatientName);
				}
			}
			$hGroupedCNV->{$global_id}->{"IGV"} = $bamNames.";".$bamFiles.";".$hGroupedCNV->{$global_id}->{"id"};
			 
			 # pour le score caller 
			 $hGroupedCNV->{$global_id}->{'SCORECALLER'} = getScoreCallers($global_id);
			 
			 # pour rechercher les BND associés au CNV
			 if ($bpok)
			{
					$hGroupedCNV->{$global_id}->{'BPManta'}=";";
					$hGroupedCNV->{$global_id}->{'BPManta_mother'}=";";
					$hGroupedCNV->{$global_id}->{'BPManta_father'}=";";
			
				
					$hGroupedCNV->{$global_id}->{'BPManta'} = getBNDManta($global_id,"patient");
					
					if ($trio)
					{
						if ($hdejavuBNDProject->{'mother'})
						{
							$hGroupedCNV->{$global_id}->{'BPManta_mother'}=";";
							$hGroupedCNV->{$global_id}->{'BPManta_mother'} = getBNDManta($global_id,"mother");
						}
						else
						{
							$hGroupedCNV->{$global_id}->{'BPManta_mother'} = "X";
						}
				
						if ($hdejavuBNDProject->{'father'})
						{
							$hGroupedCNV->{$global_id}->{'BPManta_father'}=";";
							$hGroupedCNV->{$global_id}->{'BPManta_father'} = getBNDManta($global_id,"father");
						}
						else
						{
								$hGroupedCNV->{$global_id}->{'BPManta_father'} = "X";
						}
					}
					else
					{
						$hGroupedCNV->{$global_id}->{'BPManta_mother'} = "X";
						$hGroupedCNV->{$global_id}->{'BPManta_father'} = "X";
					}
			}
			 
			 # pour claculer le score final du CNV
			 $hGroupedCNV->{$global_id}->{'SCORECNV'} = getScoreCNV($global_id);
								
	} # fin boucle sur k = liste des ind_id regroupes dans un global_id

}


sub getDejavuAndTransmission
{
		my ($type,$num,$global_id) = @_;

		my $chr= $project->getChromosome($num);
					
		my $nbdejavuProject = 0;
		my $list_of_other_patient = "";
		
		my ($t,$c,$start1,$end1) = split(/_/,$global_id);
					
		# on recherche dans le dejavu des  CNV chevauchants
		my $tab_id = $htree_dejavuProject->{$type}->{$num}->fetch($start1,$end1);
		

		# puis pour chacun d'eux on regarde ceux qui sont identiques au sens de areSameCNV
		my $theTransmission;
		my $transmission=" ";
		my $maxidM = 0;
		my $maxidF = 0;

		foreach my $djv_id (@$tab_id)
		{
				my ($t,$c,$start2,$end2) = split(/_/,$djv_id);
				my $identity = int($project->dejavuSV->getIdentityBetweenCNV($start1,$end1,$start2,$end2));
				
				if ( $identity >= $seuilTransmited*0.8)
				{
						#dejavu ou transmission ?
						foreach my $pname (keys %{$hdejavuProject->{$type}->{$num}->{$djv_id}})
						{		
								next if ($pname eq $thePatientName);
								
								if ( ($pname eq $mothername) || ($pname eq $fathername))
								{
									if ($pname eq $mothername) 
									{
										 if ($identity >= $seuilTransmited)
										 {
										 	$transmission .= "mother " unless ($transmission =~ m/mother/); 		# identity a plus de 90 = transmission
										 }
										 else
										 {
										 	$transmission .= "maybeM(".$identity."%) " unless ( ($transmission =~ m/mother/) ||  ( ($transmission =~ m/maybeM/) && ($identity>$maxidM)) ); 	# identity entre 70 et 90 = denovo bof
											$maxidM = $identity if ($identity > $maxidM);
										 }
									}
									
									if ($pname eq $fathername)
									{
										if ($identity >= $seuilTransmited)
										{
										 	$transmission .= "father " unless ( $transmission =~ m/father/);
										}
										else
										{
											$transmission .= "maybeF(".$identity."%) " unless ( ($transmission =~ m/father/) ||  ( ($transmission =~ m/maybeF/) && ($identity>$maxidF)) ); 	# identity entre 70 et 90 = denovo bof
											$maxidF = $identity if ($identity > $maxidF);
										}
									}
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
		}
		
		if ($thePatientFamily->isTrio() && $thePatient->isChild() )
		{
				 if ($transmission eq " ")
				 {
				 	$theTransmission = "strict-denovo";
				 }
				 else
				 {
				 	$theTransmission = $transmission;
				 }
		}
		else 
		{
				$theTransmission = "X";
		}
		return ($nbdejavuProject,$list_of_other_patient,$theTransmission);
}

sub getDejavuGlobal
{
	my ($type,$num,$event_id,$flag) = @_;
	

	my $nbproject=0;
	my $nbpatient=0;
	my $nbDJV_Wisecondor=0;
	my $nbDJV_Canvas=0;
	my $nbDJV_Manta=0;
	
	my $dejavuMax = $dejavu;
	$dejavuMax = "all" if $flag;		# cnv vu par wisecondor
	
	my $scorecaller_evt=0;
	my $list_of_other_patient;
	
	my ($t,$c,$start1,$end1) = split(/_/,$event_id);

	my $no = $project->dejavuSV();
	my $hashdv = $no->get_cnv($c,$start1,$end1,$t,$dejavuMax,$seuilDejaVu);
	
	my %dj;

	foreach my $project_name (keys %{$hashdv})
	{
			next if ($project_name eq $TheProjectName);
			$nbproject++;
			foreach my $pname (keys %{$hashdv->{$project_name}})
			{
					my $scorecaller=0;
					my $c;
					$nbpatient++;
					my $res;
					foreach my $caller (keys %{$hashdv->{$project_name}->{$pname}})
					{
							if ($caller eq "wisecondor") {$scorecaller +=4; $c="w";$nbDJV_Wisecondor++;}
							if ($caller eq "canvas") {$scorecaller +=2; $c="c";$nbDJV_Canvas++;}
							if ($caller eq "manta") {$scorecaller +=1; $c="m";$nbDJV_Manta++;}
							
							my ($start,$end,$identity) = split(/_/,$hashdv->{$project_name}->{$pname}->{$caller});
							$res .= int($identity)."%_".$c." ";
					}
					$scorecaller_evt += $scorecaller;
					my $etiq = $project_name.":".$pname.":".$res;
					$dj{$etiq} ++;
				}
		}

		my $theliste;
		$theliste = join(",",sort keys %dj);
		$scorecaller_evt /= $nbpatient if $nbpatient;
		
		return ($nbproject,$nbpatient,$nbDJV_Wisecondor,$nbDJV_Canvas,$nbDJV_Manta,$scorecaller_evt,$theliste);
}



sub getScoreCallers() 
{
				my ($gid) = @_;

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
				
				return $score;
				
				
}


sub getScoreCNV
{
	my ($gid) = @_;
	my $score;
	
		$score = $hGroupedCNV->{$gid}->{'SCORECALLER'};
		
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
	
		# pour remonter  ceux qui presentent un BP dans une region de 1kb autour des positions debut ou fin
		$score += 1 if ($hGroupedCNV->{$gid}->{'BPManta'} ne ";");		
		$score += 0.5 if ( ($hGroupedCNV->{$gid}->{'BPManta_mother'} ne ";") && ($hGroupedCNV->{$gid}->{'BPManta_mother'} ne "X") );	
		$score += 0.5 if ( ($hGroupedCNV->{$gid}->{'BPManta_father'} ne ";")   && ($hGroupedCNV->{$gid}->{'BPManta_father'} ne "X")   );	
	
		return $score;

}


sub filtreCNV
{
	
	foreach my $global_id ( keys %{$hGroupedCNV} ) 
	{
		my($type,$num,$gdeb,$gend) = split(/_/,$global_id);
		my $SV_ok = 1;

		# filtre sur le scoreCNV : pour les bestones on garde ceux dont le score est > a 3 
		next  if ( ( $hGroupedCNV->{$global_id}->{'SCORECNV'} < 3 ) && ($select_best == 0) );
		next  if ( ( $hGroupedCNV->{$global_id}->{'SCORECNV'} < 1 ) && ($select_best == 1) );
		next  if ( ( $hGroupedCNV->{$global_id}->{'SCORECNV'} < 0.5 ) && ($select_best == 2) );	
		
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
		
		# filtre sur les gènes omim
		if ($omim)
		{
			next unless ($hGroupedCNV->{$global_id}->{"OMIN_MG"} eq "yes");
		}
				
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
		
		# pour le plot
		my $ty = $hGroupedCNV->{$global_id}->{'TYPE'};
		$hGroupedCNV->{$global_id}->{"PLOT"}	 = $ty.";".$hGroupedCNV->{$global_id}->{"PLOT"};
		
		# au final on cherche les genes et ne calcule le score genes que pour les CNV selectionnés
		#Pour le score_genes
		my $chr = $project->getChromosome($num);
		my $tabGenes = $chr->getGenesByPosition($gdeb,$gend);

		# calculer le score gene max correspondant a la liste de gene associe au variant
		my $max = 0;
		my $maxname;
		
		my $htabscore;
		my @names;
		if ( scalar(@$tabGenes) )  # si le variant recouvre des gènes
		{	
				foreach my $g (@$tabGenes){
							$htabscore->{$g->external_name} = $g->raw_score;
				}
				
				my $nb = 0;
				foreach my $name (keys %{$htabscore})
				{
							if ($htabscore->{$name} >= $max)
							{
								$max = $htabscore->{$name};
								$maxname = $name;
							}
							$nb++;
							push (@names,$name);
				}
		
				$hGroupedCNV->{$global_id}->{'SCORE_GENES'}=$max;
				$hGroupedCNV->{$global_id}->{'GENES'} =join(" ",@names);
				$hGroupedCNV->{$global_id}->{'GENES'} = $maxname.":".$max.":".$nb." ".$hGroupedCNV->{$global_id}->{'GENES'};
				#locus a la place de genes name pour le chargement de genescout
				$hGroupedCNV->{$global_id}->{'GENES'} .= ';'.$chr->id().':'.$gdeb.'-'.$gend;


			# filtre sur une liste de genes
			my @tabgenes_annot = split(/ /,$hGroupedCNV->{$global_id}->{"GENES"} );
			$pass   = 0;
			my @tabgenes = split(/,/,$listOfGenes);
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
				$hGroupedCNV->{$global_id}->{'SCORE_GENES'}=-2;
				$hGroupedCNV->{$global_id}->{'GENES'} ="-";
		}
		
		
		#################################################################################
		# pour aller rechercher les infos PR /PR dans les vcfs manta si et seulement si evt vu par manta
		
		if ( $hGroupedCNV->{$global_id}->{"manta_global"} ) 
		{
			my $file_in = $thePatient->_getCallingSVFileWithMethodName("manta","variations");
			my $tabix = $buffer->getSoftware("tabix");
			my $cmd1 = "$tabix $file_in chr$num:$gdeb-$gend | grep $gend | cut -f 9,10 ";
			my $res1 = `$cmd1`;
			
 			unless ($res1)
 			{
 				$cmd1 = "$tabix $file_in $num:$gdeb-$gend | grep $gend | cut -f 9,10 ";
				$res1 = `$cmd1`;
 			}
 		
 			chomp($res1);


 			my @champs = split(" ",$res1);
 	
 			my @tabLabel = split(/:/,$champs[0]);
 			my @tabValue = split(/:/,$champs[1]);
 			my $pr;
 			my $sr;
 	
 	
 			for (my $i=0; $i<@tabLabel.length;$i++)
 			{
				$pr = $tabValue[$i] if ($tabLabel[$i] eq "PR"); 
				$sr = $tabValue[$i] if ($tabLabel[$i] eq "SR"); 	
 			}

 			my($ref,$alt)=split(/,/,$pr);
 			$hGroupedCNV->{$global_id}->{'PR'} =$ref."/".$alt;
 	
 			my($ref,$alt)=split(/,/,$sr);
 			$hGroupedCNV->{$global_id}->{'SR'} = $ref."/".$alt;
		} 
		else
		{
			$hGroupedCNV->{$global_id}->{'PR'} = "/";
			$hGroupedCNV->{$global_id}->{'SR'} = "/";
		}
		
		
		# pour la validation
		$hGroupedCNV->{$global_id}->{'VALIDATION'} = $thePatientName.";".$global_id.";".rand(500);
		
		#pour garder ceux qui ont passe les filtres
		$hFilteredCNV->{$global_id}= $hGroupedCNV->{$global_id};
				
		
	}# fin de la boucle sur les global_id
	
	if (scalar(keys(%{ $hFilteredCNV})) == 0) 
	{ 
				my $hash;
				$hash->{'id'} = "-";
				$hash->{'TYPE'} = "-";
				$hash->{'CHROM'} = "-";
				$hash->{'CYTOBAND'} = "-";
				$hash->{'LOCUS'} = "-";
				$hash->{'LEN'} = "-";
				$hash->{'RATIO'} = "-";
				$hash->{'DUPSEG'} = "-";
				$hash->{'SCORECNV'} = "-";
				$hash->{'PLOT'} = "-";
				$hash->{'CN'} = "-";
				$hash->{'IGV'} = "-";
				$hash->{'wisecondor_elementary'} = "Sorry";
				$hash->{'canvas_elementary'} = "No";
				$hash->{'manta_elementary'} = "Results";
				$hash->{'BPManta'} = "-";
				$hash->{'BPManta_mother'} = "-";
				$hash->{'BPManta_father'} = "-";
				$hash->{'TRANSMISSION'} = "-";
				$hash->{'DEJAVU_G'} = "-";
				$hash->{'SCORE_GENES'} = "-";
				$hash->{'RANKAnnot'} = "-";
				$hash->{'OMIM_MG'} = "-";
				$hash->{'GENES'} = "-";
				$hash->{'DGV'} = "-";
				$hash->{'GOLD_G_FREQ'} = "-";
				$hash->{'GOLD_L_FREQ'} = "-";
				$hash->{'dbVar_status'} = "-";
			 	 push(@Filtered_listHashRes, $hash); 
	}
	else
	{
		foreach my $global_id ( keys %{$hFilteredCNV} ) 
		{
				push( @Filtered_listHashRes, { %{ $hFilteredCNV->{$global_id} } } );
		}
	}
}

sub getDGV_url {
	my ( $chr, $d, $f ) = @_;
	my $url = "http://dgv.tcag.ca/gb2/gbrowse/dgv2_hg19/?name=" . $chr . "%3A". $d . "-"  . $f  . ";search=Search";
	return $url;
}

sub getgnomAD_url {
	my ( $chr, $d, $f ) = @_;
	my $url = "https://gnomad.broadinstitute.org/region/".$chr."-".$d."-".$f."?dataset=gnomad_sv_r2_1";
	return $url;
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


sub getBNDManta{
	
	my ($global_id,$who) = @_;
	my $out = ";";
	my $djv=0;

	my($type,$num,$gdeb,$gend) = split(/_/,$global_id);

	$num = "23" if ($num eq "X");
	$num = "24" if ($num eq "Y");
	
	if ( defined ($htree_dejavuBNDProject->{$who}->{$num}) )
	{
		
		my $tab_id1 = $htree_dejavuBNDProject->{$who}->{$num}->fetch($gdeb-5000,$gdeb+5000);
		my $tab_id2 = $htree_dejavuBNDProject->{$who}->{$num}->fetch($gend-5000,$gend+5000);
		
		foreach my $event_id (@$tab_id1)
		{
			$djv = getDejavuEvent($event_id);
			next if ($djv > 20);
			$out .= $event_id.";"  unless ($out =~ m/$event_id/ ) ;
		}
		
		foreach my $event_id (@$tab_id2)
		{
			$djv = getDejavuEvent($event_id);
			next if ($djv > 20);
			$out .= $event_id.";"  unless ($out =~ m/$event_id/ ) ;
		}
	}
		return $out;
}

sub getDejavuEvent {
	
	my ($event_id ) = @_;
	my $hdjv_project;
	my $hdjv_patient_iop;
	
		my ($chr1,$bp1,$chr2,$bp2) = split(/_/,$event_id);
		
		#pour ne garder que ceux avec un dejavu < 20
		foreach my $djv_event_id (keys %{$htranslocdejavu})
		{

			my ($c1,$p1,$c2,$p2) = split(/_/,$djv_event_id);
			
			next if $c1 != $chr1;
			next if $c2 !=  $chr2;
			next if ( ($p1 < $bp1-50) || ($p1 >$bp1+50) );
			next if ( ($p2 < $bp2-50) || ($p2 >$bp2+50) );
		
	
			foreach my $dejavu_project (keys %{$htranslocdejavu->{$djv_event_id}})
			{
				foreach my $dejavu_patient (keys %{$htranslocdejavu->{$djv_event_id}->{$dejavu_project}})
				{
					unless($dejavu_project eq $TheProjectName)
					{
						$hdjv_patient_iop->{$dejavu_project.":".$dejavu_patient}++; 
					}
				}
			}
		}
		
		return scalar(keys %$hdjv_patient_iop);
}



