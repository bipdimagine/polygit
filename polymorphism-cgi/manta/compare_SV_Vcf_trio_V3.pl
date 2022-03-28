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
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";

use GBuffer;
use GenBoProject;
use GenBoCache;


##################################################################################
# 	La version 3 permet de faire l'analyse à partir du vcf d'un des callers de type canvas ou wcondor
# 	puis pour chacun des variants trouvé de rechercher l'existaence d'un point de cassure, 
#  dans une zone de 5kb autour des bornes des CNVs trouvées, dans les resultats de lumpy ou manta
# 	
#  + accès à IGV ?
###################################################################


my $cgi = new CGI;
my $projectname = $cgi->param('projectname');
my $patname = $cgi->param('filename');
my $length = $cgi->param('length');
my $caller = $cgi->param('callers');
my $hide_sd = $cgi->param('hide_sd');
my $transm = $cgi->param('transmission');


my $ref ="NPH_15k";
my $type ="DUP,DEL";

	
my $path = "/data-xfs/Manue/Test_SV/";
my $fichier_cytoband = $path."/cytoband.gff3";
my $fichier_dupseg = $path."/super_Duplication.gff3";

my $path2 = "/data-xfs/Manue/Test_SV/".$projectname."/";
	

	
	
	
	my $fdm;
	my $fdl;
	my $fdc;
	my $fdwc;
	my $fdcb;
	my $fdds;
	my $fdannot;
	
	# pour stocker les coords des SV
	my %hlumpytree;
	my %hmantatree;
	my %hcanvastree;
	my %hwcondortree;
	
	# pour stocker les bandes 
	my %hcytoband;
	
	# pour stocker les duplication segmental
	my %hdupseg;
	
	# pour lire les vcfs
	my $ligne;
	my %hSV;
	my $id=0;
	
	# pour traiter les callers separement
	my %hLumpy;	
	my %hManta;			
	my %hCanvas;			
	my %hWCondor;	
	my %hSVall;
	my %hCanvasCN;
	
	# pour le json final
	my @listHashRes;

	
##########################################
#
#  gestion des trios
#
##########################################
my $buffer = GBuffer->new();	
my $project = $buffer->newProjectCache( -name => $projectname, -typeFilters=>"");
my $listPatients = $project->getPatients();
my $pat;
my $mother;
my $father;

my $mother_status;
my $father_status;

# pour récupérer les objets enfant/mother/father
foreach my $p (@$listPatients)
{
		$pat = $p if ( ($p->isChild() ) &&  ($p->name() eq $patname) );
		$mother = $p if ($p->isMother() );
		$father = $p if ($p->isFather() );
}



# les fichiers vcf et annotation du patient pour les differents callers
#my $fichier_manta = $path2."MANTA/".$patname.".vcf.gz";
#my $fichier_lumpy = $path2."LUMPY/".$patname.".vcf.gz";
#my $fichier_canvas = $path2."CANVAS/".$patname.".vcf.gz";
my $fichier_canvas = $pat->getSVFile("canvas");
my $fichier_manta = $pat->getSVFile("manta");
my $fichier_lumpy = $pat->getSVFile("lumpy");


my $fichier_condor = $path2."WCONDOR/".$ref."/".$patname."_centro_".$ref."_aberrations.bed";

# et les annotations 
my $fichier_Annot = $path2.$caller."/".lc($caller)."_SV_".$patname.".annotated.tsv";
$fichier_Annot = $path2.$caller."/".$ref."/".lc($caller)."_SV_".$patname.".annotated.tsv" if $caller eq "WCONDOR";

# les les fichiers vcf des parents pour les differents callers
	

#my $fichier_manta_mother = $path2."MANTA/".$mother->name.".vcf.gz";
#my $fichier_canvas_mother = $path2."CANVAS/".$mother->name.".vcf.gz";
#my $fichier_manta_father = $path2."MANTA/".$father->name.".vcf.gz";
#my $fichier_canvas_father= $path2."CANVAS/".$father->name.".vcf.gz";
my $fichier_manta_mother = $mother->getSVFile("manta");
my $fichier_manta_father = $father->getSVFile("manta");
my $fichier_canvas_mother = $mother->getSVFile("canvas");
my $fichier_canvas_father = $father->getSVFile("canvas");

my $mother_status = $mother->status();
my $father_status = $father->status();


##########################################
#
#    Lecture de la première ligne du fichier Annot 
#    pour recuperer le format
#
##########################################

	my $hannot;
	open($fdannot,$fichier_Annot) or die("open: $!");
	
	# lecture de la première ligne
	$ligne = <$fdannot>;
	my @champs = split(/\t/,$ligne);
	my $ind=0;
	
	foreach my $c (@champs)
	{
		chomp($c);
		$hannot->{$c}=$ind;   	#on stocke le numero du champs
		$ind++;
	}

	
##########################################
#
#  Lecture du fichier cytoband
#
##########################################

	# ouverture du fichier cytoband 
	open($fdcb,$fichier_cytoband) or die("open: $!");
	
	# lecture ligne par ligne
	while( defined( $ligne = <$fdcb> ) )
	{
	
		next unless ($ligne =~ m/^chr/ );
		my @champs = split(/\t/,$ligne);
		my @info = split(/;/,$champs[8]);
		
		
		my ($n,$band) = split(/=/,$info[0]);
		my ($i,$idband) = split(/=/,$info[1]);
		
		my $chr = $champs[0];
		my $start = $champs[3];
		my $end = $champs[4];
		
		# version avec IntervalTree
		if (!exists $hcytoband{$chr})
		{
		 	my $CytoBandTree = Set::IntervalTree->new;
		 	$hcytoband{$chr}=$CytoBandTree;
		}
		
		$hcytoband{$chr}->insert($band,$start,$end);
		
	}
	
##########################################
#
#  Lecture du fichier dupseg
#
##########################################

	# ouverture du fichier dupseg
	open($fdds,$fichier_dupseg) or die("open: $!");
	
	# lecture ligne par ligne
	while( defined( $ligne = <$fdds> ) )
	{
	
		next unless ($ligne =~ m/^chr/ );
		my @champs = split(/\t/,$ligne);
		my @info = split(/;/,$champs[8]);
		
		
		my ($n,$id) = split(/=/,$info[0]);
		my ($i,$percentid) = split(/=/,$info[1]);
		
		my $chr = $champs[0];
		my $start = $champs[3];
		my $end = $champs[4];
		
		my $idLocation = $id.":".$start."_".$end;
		
		
		# version avec IntervalTree
		if (!exists $hdupseg{$chr})
		{
		 	my $DupSegTree = Set::IntervalTree->new;
		 	$hdupseg{$chr}=$DupSegTree;
		}
		
		$hdupseg{$chr}->insert($idLocation,$start,$end);
		
	}
	
####################################
#
#			LUMPY / CANVAS ou MANTA
#
####################################

if (-e $fichier_condor)
{
# ouverture du fichier wisecondor et parsing
	open($fdwc, $fichier_condor) or die("open: $!");

	# lecture de la première ligne
	$ligne = <$fdwc>;

	
	my $SVtype;
	my $SVend ;
	my $SVscore;
	my $SVlog;
	my $SVchr;
	my $SVdeb;
	
	# traitement des lignes suivantes pour le fichier en entier
	while( defined( $ligne = <$fdwc> ) )
	{
		my @champs = split(/\t/,$ligne);
		
		$SVchr="chr".$champs[0];
		$SVdeb=$champs[1];
		$SVend=$champs[2];
		$SVscore=$champs[3];
		$SVlog=$champs[4];
		$SVtype = $champs[5];
		
		$SVtype = "DEL" if  ($SVtype =~ m/loss/);
		$SVtype = "DUP" if ($SVtype =~ m/gain/);
		
		
#		saveVariant($SVtype,$SVchr,$SVdeb,$SVend,"-",$caller,\%hWCondor) if ($caller eq "WCONDOR");
		
		# version avec IntervalTree
		if ( !exists $hwcondortree{$SVchr}->{$SVtype} ) 
		{
		 	my $newTree = Set::IntervalTree->new;
		 	$hwcondortree{$SVchr}->{$SVtype} = $newTree;
		}
		
		my $tag = "wcondor_".$SVtype."_".$SVchr."_".$SVdeb."_".$SVend;
		$hwcondortree{$SVchr}->{$SVtype}->insert($tag,$SVdeb,$SVend);
		
	}
}

if (-e $fichier_lumpy)
{
	# ouverture du fichier lumpy 
	open($fdl," zcat $fichier_lumpy | ") or die("open: $!");

	# on skipe le Header et la ligne des champs
	while( defined( $ligne = <$fdl> ) &&  $ligne =~ m/^#/ ) {};
	
	my $SVtype;
	my $SVend ;
	my $SVlength;
	my $SVchr;
	my $SVdeb;
	my $SVSU=50;
	my $id=0;
	
	
	# traitement des lignes suivantes pour le fichier en entier
	while( defined( $ligne = <$fdl> ) )
	{
		my @champs = split(/\t/,$ligne);
		my @champsINFO= split(/;/,$champs[7]);
		my @champsPAT= split(/:/,$champs[9]);
		
		$SVchr=$champs[0];
		# on ne tient pas compte des chromosomes particuliers
		#next  if ($SVchr eq "chrX");	
		#next  if ($SVchr eq "chrY"); 
		next  if ($SVchr eq "chrMT"); 
		next  if ($SVchr eq "chrM"); 	
		next  if ($SVchr =~ m/^GL/); 
		next  if ($SVchr =~ m/^hs37d5/);
		
		$SVdeb=$champs[1];

		# recuperation des donnees de INFO : TYPE  END et length
		foreach my $c (@champsINFO)
		{
			my ($key,$val)=split(/=/,$c);
			$SVtype = $val if ($key eq "SVTYPE");
			$SVend = $val if ($key eq "END");
			$SVlength = $val if ($key eq "SVLEN");
			$SVSU = $val if ($key eq "SU");
		}
		
		next unless ($type =~ m/$SVtype/);				# pour l'instant on ne s'interresse que au DUP/DEL

		# on ne garde que les lignes qui de longueur > length
		next if ( abs($SVlength) < $length);		
		
#		my $GT = $champsPAT[0];
#		
#		# si ça passe ....
#		saveVariant($SVtype,$SVchr,$SVdeb,$SVend,$GT,$caller,\%hLumpy) if ($caller eq "LUMPY");
		
		# version avec IntervalTree
		if ( !exists $hlumpytree{$SVchr}->{$SVtype} ) 
		{
		 	my $newTree = Set::IntervalTree->new;
		 	$hlumpytree{$SVchr}->{$SVtype} = $newTree;
		}
		
		my $tag = "lumpy_".$SVtype."_".$SVchr."_".$SVdeb."_".$SVend;
		$hlumpytree{$SVchr}->{$SVtype}->insert($tag,$SVdeb,$SVend);
		
	}
}

if (-e $fichier_canvas)
{
	
	# ouverture du fichier canvas et parsing
	open($fdc," zcat $fichier_canvas | ") or die("open: $!");

	# on skipe le Header et la ligne des champs
	while( defined( $ligne = <$fdc> ) &&  $ligne =~ m/^#/ ) {}

	my $SVtype;
	my $SVend ;
	my $SVlength;
	my $SVchr;
	my $SVdeb;
	my $id=0;
	
	# traitement des lignes suivantes pour le fichier en entier
	while( defined( $ligne = <$fdc> ) )
	{
		
					my @champs = split(/\t/,$ligne);
					my @champsINFO= split(/;/,$champs[7]);
					my @champsPAT= split(/:/,$champs[9]);
					
					$SVchr=$champs[0];
					next  if ($SVchr eq "chrX");	
					next  if ($SVchr eq "chrY"); 
					next  if ($SVchr eq "chrMT"); 	
					next  if ($SVchr =~ m/^GL/); 
					next  if ($SVchr =~ m/^hs37d5/);
					
					
					$SVdeb=$champs[1];
					my $CN = $champs[4];
		 
					my @canvasinfo=split(/:/,$champs[2]);
					$SVtype = $canvasinfo[1];
		
					$SVtype = "DEL" if ($SVtype eq "LOSS");
					$SVtype = "DUP" if ($SVtype eq "GAIN");
					$SVtype = "LOH" if ($SVtype eq "LOH");
		

					next unless ($type =~ m/$SVtype/);				# pour l'instant on ne s'interresse que au DUP/DEL 

					next if ( ( $champs[6] ne "PASS"));				# on ne garde que les lignes ou PASS
		
					# recuperation des donnees de INFO : TYPE  END et length
					foreach my $c (@champsINFO)
					{
							my ($key,$val)=split(/=/,$c);
							$SVend = $val if ($key eq "END");
							$SVlength = $val if ($key eq "CNVLEN");
					}
		
					# on ne garde que les lignes qui de longueur > length
					next if ( abs($SVlength) < $length);		
					
					my $GT = $champsPAT[0];
					
					# si ca passe ....
					# pour conserver le copynumber
					$hCanvasCN{$SVtype.":".$SVchr.":".$SVdeb.":".$SVend} = $CN;
									
					saveVariant($SVtype,$SVchr,$SVdeb,$SVend,$GT,$caller,\%hCanvas)  if ($caller eq "CANVAS");
					
					# version avec IntervalTree
					if ( !exists $hcanvastree{$SVchr}->{$SVtype} ) 
					{
		 				my $newTree = Set::IntervalTree->new;
		 				$hcanvastree{$SVchr}->{$SVtype} = $newTree;
					}
		
					my $tag = "canvas_".$SVtype."_".$SVchr."_".$SVdeb."_".$SVend;
					$hcanvastree{$SVchr}->{$SVtype}->insert($tag,$SVdeb,$SVend);
		}
}
	
if (-e $fichier_manta)
{	
		# ouverture du fichier manta zippé
		open($fdm," zcat $fichier_manta | ") or die("open: $!");
		
		# on skipe le Header et la ligne des champs
		while( defined( $ligne = <$fdm> ) &&  $ligne =~ m/^##/ ) {}

		my $SVtype;
		my $SVend ;
		my $SVlength;
		my $SVchr;
		my $SVdeb;
		my $id=0;
		
		# traitement des lignes suivantes pour le fichier en entier
		while( defined( $ligne = <$fdm> ) )
		{
			
				my @champs = split(/\t/,$ligne);
				my @champsINFO= split(/;/,$champs[7]);
				my @champsPAT = split(/:/,$champs[9]);

				$SVchr=$champs[0];
				
				#next  if ($SVchr eq "chrX");	
				#next  if ($SVchr eq "chrY"); 
				next  if ($SVchr eq "chrMT"); 	
				next  if ($SVchr =~ m/^GL/); 
				next  if ($SVchr =~ m/^hs37d5/);
				
				
				$SVdeb=$champs[1];
				my $qual = $champs[5];
		
				# recuperation des donnees de INFO : TYPE  END et length
				foreach my $c (@champsINFO)
				{
					my ($key,$val)=split(/=/,$c);
					$SVtype = $val if ($key eq "SVTYPE");
					$SVend = $val if ($key eq "END");
					$SVlength = $val if ($key eq "SVLEN");
				}
				
		
				next unless ($qual > 500); 					# on ne garde que les lignes de bonne qualite
				next unless ($type =~ m/$SVtype/);		# pour l'instant on ne s'interresse que au DUP/DEL et INS
				next if ( abs($SVlength) < $length);		# on ne garde que les lignes qui de longueur > length
				next if ( ( $champs[6] ne "PASS")); 		# on ne garde que les lignes ou PASS
	

				my $GT = $champsPAT[0];		

				# si ca passe ...				
				saveVariant($SVtype,$SVchr,$SVdeb,$SVend,$GT,$caller,\%hManta)  if ($caller eq "MANTA");
				
				# version avec IntervalTree
				if ( !exists $hmantatree{$SVchr}->{$SVtype} ) 
				{
		 				my $newTree = Set::IntervalTree->new;
		 				$hmantatree{$SVchr}->{$SVtype} = $newTree;
				}
		
				my $tag = "manta_".$SVtype."_".$SVchr."_".$SVdeb."_".$SVend;
				$hmantatree{$SVchr}->{$SVtype}->insert($tag,$SVdeb,$SVend);
		}
}		

annotVariants($caller,\%hManta)  if ($caller eq "MANTA");
annotVariants($caller,\%hCanvas)  if ($caller eq "CANVAS");

%hSVall = %hCanvas if ($caller eq "CANVAS" );
%hSVall = %hManta if ($caller eq "MANTA" );

 
# ---------------------------------------
# creation des json pour l'interface 
# --------------------------------------

if (scalar(keys(%hSVall)) == 0) 
{ 
		my $hash;
		$hash->{'id'}="-";
		$hash->{'CHROM_POS_END'}="-";
		$hash->{'CALLER'}="-";
		$hash->{'CYTOBAND'}="-";
		$hash->{'SVLEN'}="-";
		$hash->{'SVTYPE'}="-";
		$hash->{'DGV'}="-";
		$hash->{'GT'}="-";
		$hash->{'CN'}="-";
		$hash->{'DUPSEG'}="-";
		$hash->{'GENES'}="-";
		$hash->{'GOLD'}="-";
		$hash->{'GOLD_G_freq'}="-";
		$hash->{'GOLD_L_freq'}="-";
		$hash->{'RANKAnnot'}="-";
		$hash->{'DDD_DEL_FREQ'}="-";
		$hash->{'DDD_DUP_FREQ'}="-";
		$hash->{'OMIN_MG'}="-";
		$hash->{'dbVar_event'}="-";
		$hash->{'dbVar_status'}="-";
		$hash->{'OTHERS'}="-";
		$hash->{'BP_LUMPY'}="-";
		$hash->{'BP_MANTA'}="-";
		$hash->{'TRANSMISSION'}="-";
		$hash->{'ScoreT'} = "-";
		
		push(@listHashRes, $hash);
}
else
{
	foreach my $id (sort(keys(%hSVall))) 
	{ 
			push(@listHashRes, {%{$hSVall{$id}}} ); 
	}	
}

printJson(\@listHashRes);
exit(0);


###############################################################################################		


sub saveVariant()
{
	my ($SVtype,$SVchr,$SVdeb,$SVend,$GT,$Caller,$hresCaller) = @_;	
		
	my $url_DGV = getDGV_url($SVchr,$SVdeb,$SVend);
	my $band = getCytoBand($SVchr,$SVdeb,$SVend);
	my $copyNumber = $hCanvasCN{$SVtype.":".$SVchr.":".$SVdeb.":".$SVend}; 
	
	
	# pour les segmental duplication
	my $dupseg = getDupSeg($SVchr,$SVdeb,$SVend);
	return 0 if ($hide_sd && ($dupseg > 50)); 	# on cache les duplications segmental
	
	my ($ch,$num) =split(/r/,$SVchr);
	
	
	
	$hresCaller->{$id}->{'id'}=int($id);
	$hresCaller->{$id}->{'CHROM'}=$SVchr;
	$hresCaller->{$id}->{'CHROM_POS_END'}=$num."_".$SVdeb."_".$SVend;
	$hresCaller->{$id}->{'CALLER'}=$Caller;
	$hresCaller->{$id}->{'POS'}=$SVdeb;
	$hresCaller->{$id}->{'END'}=$SVend;
	$hresCaller->{$id}->{'CYTOBAND'}=$band;
	$hresCaller->{$id}->{'SVLEN'}=abs($SVdeb-$SVend);
	$hresCaller->{$id}->{'SVTYPE'}=$SVtype;
	$hresCaller->{$id}->{'DGV'}=$url_DGV;
	$hresCaller->{$id}->{'GT'}=$GT;
	$hresCaller->{$id}->{'CN'}=$copyNumber;
	$hresCaller->{$id}->{'DUPSEG'}=$dupseg;
	
		
	$id++;
	return 1;
}

sub annotVariants()
{
	
	my ($Caller,$hresCaller) = @_;	
	
	# boucle sur les variants
	foreach my $id ( keys(%{$hresCaller} ) )
	{
		
		my $SVtype = $hresCaller->{$id}->{'SVTYPE'};
		my $SVchr = $hresCaller->{$id}->{'CHROM'};
		my $SVdeb = $hresCaller->{$id}->{'POS'};
		my $SVend = $hresCaller->{$id}->{'END'};
		
		my $l = abs( $SVdeb-$SVend);
		
		# pour voir si le variant a ete vu avec les autres caller et definir les point de cassure
		my $lumpy_bp = "-";	
		my $manta_bp = "-";	
		$lumpy_bp = getBreakpoint($SVchr,$SVdeb,$SVend,"lumpy",$fichier_lumpy) if (-e $fichier_lumpy);	
		$manta_bp = getBreakpoint($SVchr,$SVdeb,$SVend,"manta",$fichier_manta) if (-e $fichier_manta);		
		my $others = getInfoOtherCaller($SVtype,$SVchr,$SVdeb,$SVend,$Caller);
		
		
		# pour voir si le variant est transmis par le père et/ou la mère
		my %hscoret;
		my $SVGT = $hresCaller->{$id}->{'GT'};
		my $transmission = " ";
		
		$transmission .= getTransmission($SVchr,$SVdeb,$SVend,"mother","manta", $fichier_manta_mother,\%hscoret) if (-e $fichier_manta_mother);
		$transmission .= getTransmission($SVchr,$SVdeb,$SVend,"father","manta", $fichier_manta_father,\%hscoret) if (-e $fichier_manta_father);
		$transmission .= getTransmission($SVchr,$SVdeb,$SVend,"mother","canvas",$fichier_canvas_mother,\%hscoret) if (-e $fichier_canvas_mother);
		$transmission .= getTransmission($SVchr,$SVdeb,$SVend,"father","canvas", $fichier_canvas_father,\%hscoret) if (-e $fichier_canvas_father);
		
		# calcul du score transmission
		my $scoreTransmission = getScoreTransmission(\%hscoret,$SVGT,$mother_status,$father_status);
		#my ($score,$type) = split(/ /,$scoreTransmission);

		# pour n'afficher que ceux repondant a un modele donne
		unless ( $transm =~  m/all/)
		{	
			unless ( $transm =~  m/$scoreTransmission/)
			{
				delete ($hresCaller->{$id});
				next;
			}
		}
		
		
		$hresCaller->{$id}->{'ScoreT'} = $scoreTransmission;
					
		#############################
		#  annotations a partir de AnnotSV
		#############################
		
		#Genes
		my $liste_of_genes = "-";
		
		$liste_of_genes = getAnnots($SVchr,$SVdeb,$SVend,$hannot->{"Gene name"});
		my @tabgenes = split(/\//,$liste_of_genes);
		$liste_of_genes =~ s/\//,/g;
		
		
		#DGV
		my $dgv_gain = getAnnots($SVchr,$SVdeb,$SVend,$hannot->{'DGV_GAIN_IDs'});
		my $dgv_loss = getAnnots($SVchr,$SVdeb,$SVend,$hannot->{'DGV_LOSS_IDs'});
		my $dgv_gain_freq = getAnnots($SVchr,$SVdeb,$SVend,$hannot->{'DGV_GAIN_Frequency'});
		my $dgv_loss_freq = getAnnots($SVchr,$SVdeb,$SVend,$hannot->{'DGV_LOSS_Frequency'});

		my $SV_rank = getAnnots($SVchr,$SVdeb,$SVend,$hannot->{'AnnotSV ranking'});	
		
		#DDD
		my $DDD_DEL_freq =  getAnnots($SVchr,$SVdeb,$SVend,$hannot->{'DDD_DEL_Frequency'});
		my $DDD_DUP_freq =  getAnnots($SVchr,$SVdeb,$SVend,$hannot->{'DDD_DUP_Frequency'});
		
		my $OMIN_MG = getAnnots($SVchr,$SVdeb,$SVend,$hannot->{'morbidGenes'});	
		
		#dbVar
		my $dbVar_event =  getAnnots($SVchr,$SVdeb,$SVend,$hannot->{'dbVar_event'});
		$dbVar_event =~ s/;/ /g;
		my $dbVar_status =  getAnnots($SVchr,$SVdeb,$SVend,$hannot->{'dbVar_status'});
		$dbVar_status =~ s/;/ /g;	
		
		$hresCaller->{$id}->{'GENES'}=$liste_of_genes;
		$hresCaller->{$id}->{'GOLD'}=$dgv_gain." ".$dgv_loss;
		$hresCaller->{$id}->{'GOLD_G_freq'}=$dgv_gain_freq;
		$hresCaller->{$id}->{'GOLD_L_freq'}=$dgv_loss_freq;
		$hresCaller->{$id}->{'RANKAnnot'}=$SV_rank;
		$hresCaller->{$id}->{'DDD_DEL_FREQ'}=$DDD_DEL_freq;
		$hresCaller->{$id}->{'DDD_DUP_FREQ'}=$DDD_DUP_freq;
		$hresCaller->{$id}->{'OMIN_MG'}=$OMIN_MG;
		$hresCaller->{$id}->{'dbVar_event'}=$dbVar_event;
		$hresCaller->{$id}->{'dbVar_status'}=$dbVar_status;
		$hresCaller->{$id}->{'OTHERS'}=$others;
		$hresCaller->{$id}->{'BP_LUMPY'}=$lumpy_bp;
		$hresCaller->{$id}->{'BP_MANTA'}=$manta_bp;
		$hresCaller->{$id}->{'TRANSMISSION'}=$transmission;
		
	}
}


sub printJson {
	my ($listHash) = @_;
	my $hash;
	$hash->{'identifier'} = 'id';
	$hash->{'label'} = 'id';
	$hash->{'items'} = $listHash;
	print $cgi->header('text/json-comment-filtered');
	print encode_json $hash;
	print "\n";
}

sub getDGV_url
{
	my ($chr,$d,$f) = @_;
	my $url = "http://dgv.tcag.ca/gb2/gbrowse/dgv2_hg19/?name=".$chr."%3A".$d."-".$f.";search=Search";
	return $url;
}

# version avec interval_tree
sub getCytoBand
{
	my ($chr,$d,$f) = @_;
	my $cytoband;
	my $result;
	
	$result = $hcytoband{$chr}->fetch($d,$f);
	
	return "noband"	 if (scalar(@$result) == 0);
		
	while (scalar(@$result) > 0)
	{
		my $c = shift(@$result);
		$cytoband .= $c.","; 
	}
	
	return $cytoband;
}

sub getDupSeg
{
	my ($chr,$d,$f) = @_;
	my $dupseg;
	my $result;

	my $cov50 = 0;
	my $cov10 = 0;
	
	# recherche des dupseg chevauchant le fragment
	$result = $hdupseg{$chr}->fetch($d,$f);
	my $nb_dupseg = scalar(@$result);
	
	# si il n y en a pas
	return "-"	 if ( $nb_dupseg == 0);
	
	# sinon on creer un intspan
	my $span = Set::IntSpan::Fast::XS->new();
	
	foreach my $res (@$result)
	{
		 my ($id,$coord) = split(/:/,$res);
		 my ($ds_start, $ds_end) = split(/_/,$coord); 
		 
		 # on ne conserve que la partie comprise entre $d et $f
		 $ds_start = $d if ($ds_start <= $d);
		 $ds_end = $f if ($ds_end >= $f);
		 
	    $span->add_range($ds_start,$ds_end);
	}
	
	# calcul de la couverture globale des dupseg
	my $globale_cov=0;
	my $list_pos = $span->as_string();
	my @tabpos = split(/,/,$list_pos);

	foreach my $coords (@tabpos)
	{
		my ($debut,$fin) =split(/-/,$coords);
		my $taille = ($fin-$debut);
		$globale_cov = $globale_cov + $taille;
	}
	
	return $globale_cov/($f-$d)*100;
}

sub getAnnots
{
	my ($chr,$d,$f,$num) = @_;
	
	my @n = split(/chr/,$chr);
	my $chrnum = $n[1];
	
	my $cmd =  "cat ".$fichier_Annot." | grep ".$chrnum."_".$d."_".$f." | grep full";
 	my $res1 = `$cmd`;
 	chomp($res1);
 	
 	my @champsAnnot = split(/\t/,$res1);

 	return $champsAnnot[$num] if $champsAnnot[$num];
 	return "-";
}

sub getBreakpoint
{
	my ($chr,$deb,$end,$caller,$fichier) = @_;
	
	my $b1 = $deb-5000;
	my $b2 = $deb+5000;
	
	my $cmd =  "tabix ".$fichier." ".$chr.":".$b1."-".$b2;
 	my $res1 = `$cmd`;
 	chomp($res1);
 	 
 	my $out = " ";
 	
 
 	if ($res1)
 	{
		# determiner le nombre de SV present dans l'interval 
		my @champsSV = split(/\n/,$res1);
		my $nbSV = scalar(@champsSV);
		
		$out .=  "";
		
		for(my $i=0 ; $i < $nbSV ; $i++)
		{
			my @champs = split(/\t/,$champsSV[$i]);
			
			my $bpchr=$champs[0];
			my $bppos=$champs[1];
			my $bpalt=$champs[4];
			my $bptype = "?";
			my $bpend = "?";
			my $bplength = "?";
	
			my @champsINFO= split(/;/,$champs[7]);
			
			foreach my $c (@champsINFO)
			{
					my ($key,$val)=split(/=/,$c);
					$bptype = $val if ($key eq "SVTYPE");
					$bpend = $val if ($key eq "END");
					$bplength = $val if ($key eq "SVLEN");
			}
 			$out .= "  ".$bptype.":".$bpchr.":".$bppos."-".$bpend." ".$bpalt."next" if ( abs($bppos - $deb) <= 5000);
 		}
		
 	}
	return $out;
 	
}
	
sub getInfoOtherCaller
{
	my ($type,$chr,$deb,$end,$caller) = @_;
	
	
	my $res1;
	my $res2;
	my $res3;
	
	my $out = " ";
	
	my $set = Set::IntSpan::Fast::XS->new();	
	$set->add_range($deb,$end);
	
	my $setRes = Set::IntSpan::Fast::XS->new();	
	my $setChevauchant = Set::IntSpan::Fast::XS->new();	
	
	
	if ($caller eq "CANVAS")
	{	
 		$res1 = $hlumpytree{$chr}->{$type}->fetch($deb,$end) if exists $hlumpytree{$chr}->{$type};
 		$res2 = $hmantatree{$chr}->{$type}->fetch($deb,$end) if exists $hmantatree{$chr}->{$type};
 		$res3 = $hwcondortree{$chr}->{$type}->fetch($deb,$end) if exists $hwcondortree{$chr}->{$type};
	}
	
	if ($caller eq "LUMPY")
	{
			$res1 = $hcanvastree{$chr}->{$type}->fetch($deb,$end) if exists $hcanvastree{$chr}->{$type};
 			$res2 = $hmantatree{$chr}->{$type}->fetch($deb,$end)  if exists $hmantatree{$chr}->{$type};
 			$res3 = $hwcondortree{$chr}->{$type}->fetch($deb,$end) if exists $hwcondortree{$chr}->{$type};
	}
	
	if ($caller eq "MANTA") 
	{
			$res1 = $hcanvastree{$chr}->{$type}->fetch($deb,$end)  if exists $hcanvastree{$chr}->{$type};
 			$res2 = $hlumpytree{$chr}->{$type}->fetch($deb,$end)	if exists $hlumpytree{$chr}->{$type};
 			$res3 = $hwcondortree{$chr}->{$type}->fetch($deb,$end) if exists $hwcondortree{$chr}->{$type};
	}
	
	if ($caller eq "WCONDOR") 
	{
			$res1 = $hcanvastree{$chr}->{$type}->fetch($deb,$end)  if exists $hcanvastree{$chr}->{$type};
 			$res2 = $hlumpytree{$chr}->{$type}->fetch($deb,$end)	if exists $hlumpytree{$chr}->{$type};
 			$res3 = $hmantatree{$chr}->{$type}->fetch($deb,$end)   if exists $hmantatree{$chr}->{$type};
	}	
	
	
	$out .= join(",",@$res1).","  if ($res1 && scalar(@$res1));
 	$out .= join(",",@$res2).","  if ($res2 && scalar(@$res2));
 	$out .= join(",",@$res3).","  if ($res3 && scalar(@$res3));
 	
 	return "-" unless ($out ne " ");
 	
 	# pour ne garder pour les autres callers que les evenement chevauchant a 70% au moins 
 	# et dont au moins une des bornes est a moins de 5kb des bornes originales 
 	my $lc1 = 0.7 * abs($end-$deb);
	my @res = split(/,/,$out);
	
	my $out_res = "";
	
	foreach my $sv (@res)
	{
		my @ch = split(/_/,$sv);
		my $d = $ch[3];
		my $f = $ch[4];
		
		my $lc2 = 0.7 * abs($d-$f);
		
		next unless ( (abs($d - $deb) <= 5000)  || ( abs($f- $end) <= 5000) );
		
		$setRes = Set::IntSpan::Fast::XS->new();	
		$setRes->add_range($d,$f);
		
		$setChevauchant = $set->intersection($setRes);
		my $l = $setChevauchant->as_string(); 
		
		my ($b1,$b2) = split(/-/,$l);
		$out_res .= $sv."," if ( ( $b2-$b1 >= $lc1) &&  ($b2-$b1 >= $lc2) );
 	}	
 	
 	return $out_res;
}

sub getTransmission
{
	my ($chr,$deb,$end,$parent,$caller,$fichier,$hscore) = @_;
	
	my $b1 = $deb-5000;
	my $b2 = $deb+5000;
	
	my $cmd;
	my $res1;
	
	$cmd =  "tabix ".$fichier." ".$chr.":".$b1."-".$b2;
 	$res1 = `$cmd`;
 	chomp($res1);

 	my $out = " ";
 
 	if ($res1)
 	{
		# determiner le nombre de SV present dans l'interval 
		my @champsSV = split(/\n/,$res1);
		my $nbSV = scalar(@champsSV);
		
		
		for(my $i=0 ; $i < $nbSV ; $i++)
		{
			my @champs = split(/\t/,$champsSV[$i]);
			
			my $bpchr=$champs[0];
			my $bppos=$champs[1];
			my $bpalt=$champs[4];
			my $bptype = "?";
			my $bpend = "?";
			my $bplength = "?";
	
			my @champsINFO= split(/;/,$champs[7]);
			my @champsPAT = split(/:/,$champs[9]);
			
			my $GT= $champsPAT[0];
	
			foreach my $c (@champsINFO)
			{
					my ($key,$val)=split(/=/,$c);
					$bptype = $val if ($key eq "SVTYPE");
					$bpend = $val if ($key eq "END");
					$bplength = $val if ($key eq "SVLEN");
			}
			
			if ( abs($bppos - $deb) <= 5000)
			{
				$out .= "  ".$caller."_".$parent.":".$bptype.":".$bpchr.":".$bppos."-".$bpend." ".$bpalt." ".$GT."next"; 
 				$hscore->{$caller}->{$parent}->{"GT"} = $GT;
			}
		}
 	}
	return $out;
}

sub getScoreTransmission
{
	my ($hscoret,$SVGT,$mother_status,$father_status) = @_;
	my $score=0;
	my $transmission_type = "-";
	my $GT_mother;
	my $GT_father;
	
	my $both = 0;
	
		foreach my $caller (keys(%{$hscoret}))
		{
			$GT_mother = 0;
			$GT_father = 0;
			$GT_mother = $hscoret->{$caller}->{"mother"}->{"GT"} if ($hscoret->{$caller}->{"mother"}->{"GT"});
			$GT_father   = $hscoret->{$caller}->{"father"}->{"GT"}   if ($hscoret->{$caller}->{"father"}->{"GT"});
			

			if (  $GT_mother && $GT_father ) 		# transmission des deux parents
			{
				if ( ($mother_status == 1)  && ( $father_status == 1) ) 		# les deux parents sont sains
				{		
					if ($SVGT eq "0/1")
					{
							if ( ($GT_mother eq "0/0")  &&  ( $GT_father eq "0/0" ) )
							{
									$score = 4; # on a un denovo
									$transmission_type = "denovo";
									$both = 1;
							}
							else
							{
								if ( ($GT_mother eq "0/1") ||  ($GT_father eq "0/1") )
								{
										$score = 0;
										$transmission_type = "both";
										$both = 1;
								}
								else
								{
									if ( ($GT_mother eq "1/1")  &&  ( $GT_father eq "1/1" )  )
									{
										$score = -2; # erreur mendeliene
										$transmission_type = "error";
										$both = 1;
									}
								}
							}
					}
					else		# SVGT = 1/1 ou 1/2
					{	
						if ( ($GT_mother eq "0/1")  &&  ( $GT_father eq "0/1" )  )
						{
							$score = 4; # on a un recessif
							$transmission_type = "recessive";
							$both = 1;
						}
						else
						{
							if ( ($GT_mother eq "0/0") ||  ($GT_father eq "0/0") )
							{
								$score = -2; # erreur mendeliene
								$transmission_type = "error";
								$both = 1;
							}
							else
							{
								if ( ($GT_mother eq "1/1") ||  ($GT_father eq "1/1"))
								{
									$score = 0; 	# ne joue pas sur le phenotype 
									$transmission_type = "both";
									$both = 1;
								}
							}
						}			
					}
			}
			else
			{
				$score += 0  if (( $father_status == 1)   && ( $SVGT ne  $GT_father)); 	    # ne joue pas sur le phenotype 
				$score += 0  if (( $mother_status == 1) && ( $SVGT ne  $GT_mother)); 	    # ne joue pas sur le phenotype 
				$transmission_type = "-";
				$both = 1;
			}
		}
		else #transmission d'un seul parent
		{
			unless ($both)
			{
				if ($GT_mother)  # transmission par la mère uniquement
				{
					if ( $mother_status == 1) # parent  sain
					{		
							$score = 2;
							$transmission_type = "mother"; 
					}
				}
				
				if ($GT_father)  # transmission par le père uniquement
				{
					if ( $father_status == 1) # parent  sain
					{		
						$score = 2;
						$transmission_type = "father";
					}
				}
			}
		}
	}
	#return $score." ".$transmission_type;
	return $transmission_type;
}	


 	