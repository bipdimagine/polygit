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
 
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";

use GBuffer;
use GenBoProject;
use GenBoCache;


############################################################################################################
# Cherche tous les variants structuraux pour chaque patient d'un projet et tous callers confondus.
# Stocke le resultat dans une  table de hash qui est freeze : nom du fihier = project.allSV
# Pour chaque SV la localisation chr_start_end donne accès à la liste des gènes compris dans cet interval
############################################################################################################

my $cgi = new CGI;
my $projectname = $cgi->param('projectname');

warn $projectname;



my $length = 1000;
my $ref ="NPH_15k";

my $type ="DUP,DEL";
#my $type ="DUP,DEL,INS,INV";
	


my $path2 = "/data-xfs/Manue/Test_SV/".$projectname."/";

	
my $fd;
my $fdcb;
my $fdds;
my $fdannot;
	
###############
# tables de hash
###############

my $hannot;
	
# pour lire les vcfs
my $ligne;
my %hSV;
my $compteur=0;
	
# pour stocker les SV
my %hSVall;  # tous par id

# pour sauvegarder les infos patient de chaque SV
my $hpatientInfo;
my %hCanvasCN;


	
	
###################################
#
#  lecture des vcf et 
#  gestion du projet via les objets Genbo
#
###################################


my $buffer = GBuffer->new();	
my $project = $buffer->newProjectCache( -name => $projectname);
my $listPatients = $project->getPatients();

# pour la sauvegarde du fichier de sortie
my $variationsDir = $project->getVariationsDir();
my $file_out = $variationsDir.$projectname.".allSV";
warn $file_out;


my $listCallers = $project->callingSVMethods();
my %hpatientFiles;

# table des noms de fichiers 
foreach my $pat (@$listPatients)
{
	my $patname = $pat->name();
	
	# pour récupérer le bamfile
	my $patientBam = $pat->getBamFileName("bwa");
	
	my @tab = split(/ngs/,$patientBam);
	$patientBam = "/NGS".$tab[1];
	
	$hpatientFiles{$patname}->{'bam'} = $patientBam;
	
	foreach my $caller (@$listCallers)
	{
			$hpatientFiles{$patname}->{$caller}->{'vcf'} = $pat->getSVFile($caller); 		
			$hpatientFiles{$patname}->{$caller}->{'annotSV'} = $pat->getAnnotSVFileName($caller); 	
			$hpatientFiles{$patname}->{$caller}->{'annotSV'} = 0 if ($caller eq "lumpy");
	}
}
	

# pour récupérer l'objet patient
foreach my $pat (@$listPatients)
{
		my $SVtype;
		my $SVend ;
		my $SVlength;
		my $SVchr;
		my $SVdeb;
		my $GT;
		
		my $patname = $pat->name();
		

		#######################################
		# gestion des differents callers SV
		#######################################

		foreach my $caller (@$listCallers)
		{
				my $dir = $project->getVariationsDir($caller);
	
				##########################################
				#
				#    Lecture de la première ligne du fichier Annot 
				#    pour recuperer le format
				#
				##########################################

				my $fichier_Annot = $hpatientFiles{$patname}->{$caller}->{'annotSV'} ;
				
				if ($fichier_Annot)
				{
					open($fdannot,$fichier_Annot) or die("open: $!");
	
					# lecture de la première ligne
					$ligne = <$fdannot>;
					my @champs = split(/\t/,$ligne);
					my $ind=0;
	
					foreach my $c (@champs)
					{
						chomp($c);
						$hannot->{$caller}->{$c}=$ind;   	#on stocke le numero du champs
						$ind++;
					}
				}
				
				my $fichierPatient = $hpatientFiles{$patname}->{$caller}->{'vcf'};
				
				if ($caller  eq "wisecondor")
				{
					
					# ouverture du fichier wisecondor et parsing
					open($fd," zcat $fichierPatient | ") or die("open: $!");
			
					# lecture de la première ligne
					$ligne = <$fd>;

					my $SVscore;
					my $SVlog;
					
					$GT = "-";

	
					# traitement des lignes suivantes pour le fichier en entier
					while( defined( $ligne = <$fd> ) )
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
		
						# on ne garde que les lignes qui de longueur > length
						my $SVlength = $SVend-$SVdeb;
						next if ( abs($SVlength) < $length);		
						
					
						saveVariant($patname,$SVtype,$SVchr,$SVdeb,$SVend,$GT,$caller,\%hSVall);
						annotVariant($patname,$SVtype,$SVchr,$SVdeb,$SVend,$GT,$caller,\%hSVall); 
						
					}
				}
		
		if ($caller eq "lumpy")
		{

				# ouverture du fichier lumpy 
				open($fd," zcat $fichierPatient | ") or die("open: $!");

				# on skipe le Header et la ligne des champs
				while( defined( $ligne = <$fd> ) &&  $ligne =~ m/^#/ ) {};
	
				my $SVSU=50;

	
			# traitement des lignes suivantes pour le fichier en entier
			while( defined( $ligne = <$fd> ) )
			{
				my @champs = split(/\t/,$ligne);
				my @champsINFO= split(/;/,$champs[7]);
				my @champsPAT= split(/:/,$champs[9]);
				
				$SVtype="-";
				$SVend=0;
				$SVlength=0;
				$SVdeb=0;
				$SVchr=$champs[0];
				
				# on ne tient pas compte des chromosomes particuliers
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
		
				$GT = $champsPAT[0];
		
				# si ça passe ....
				#saveVariant($patname,$SVtype,$SVchr,$SVdeb,$SVend,$GT,$caller,\%hSVall);
				#annotVariant($patname,$SVtype,$SVchr,$SVdeb,$SVend,$GT,$caller,\%hSVall); 
			}
	}

	if ($caller eq "canvas")
	{
	
			# ouverture du fichier canvas et parsing
			open($fd," zcat $fichierPatient | ") or die("open: $!");

			# on skipe le Header et la ligne des champs
			while( defined( $ligne = <$fd> ) &&  $ligne =~ m/^#/ ) {}

			# traitement des lignes suivantes pour le fichier en entier
			while( defined( $ligne = <$fd> ) )
			{
		
					my @champs = split(/\t/,$ligne);
					my @champsINFO= split(/;/,$champs[7]);
					my @champsPAT= split(/:/,$champs[9]);
					
					
					$SVchr=$champs[0];

					next  if ($SVchr eq "chrMT"); 	
					next  if ($SVchr =~ m/^GL/); 
					next  if ($SVchr =~ m/^hs37d5/);
					
					
					$SVdeb=$champs[1];
					my $CN = $champs[4];
					$CN =~ s/<//;
					$CN =~ s/>//;
		 
					my @canvasinfo=split(/:/,$champs[2]);
					$SVtype = $canvasinfo[1];
		
					$SVtype = "DEL" if ($SVtype eq "LOSS");
					$SVtype = "DUP" if ($SVtype eq "GAIN");
					$SVtype = "LOH" if ($SVtype eq "LOH");
		

					next unless ($type =~ m/$SVtype/);				# pour l'instant on ne s'interresse que au DUP/DEL 

					next if ( ( $champs[6] ne "PASS"));				# on ne garde que les lignes ou PASS
		
					# recuperation des donnees de INFO : TYPE  END et length
					$SVend=0;
					$SVlength=0;
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
					$hCanvasCN{$patname}->{$SVtype.":".$SVchr.":".$SVdeb.":".$SVend} = $CN;
									

					saveVariant($patname,$SVtype,$SVchr,$SVdeb,$SVend,$GT,$caller,\%hSVall);
					annotVariant($patname,$SVtype,$SVchr,$SVdeb,$SVend,$GT,$caller,\%hSVall); 					
			}
	}
	
	if ($caller eq "manta")
	{	
				# ouverture du fichier manta zippé
				open($fd," zcat $fichierPatient | ") or die("open: $!");
		
				# on skipe le Header et la ligne des champs
				while( defined( $ligne = <$fd> ) &&  $ligne =~ m/^##/ ) {}
	
				# traitement des lignes suivantes pour le fichier en entier
			while( defined( $ligne = <$fd> ) )
			{
			
				my @champs = split(/\t/,$ligne);
				my @champsINFO= split(/;/,$champs[7]);
				my @champsPAT = split(/:/,$champs[9]);

				$SVchr=$champs[0];
				
				next  if ($SVchr eq "chrMT"); 	
				next  if ($SVchr =~ m/^GL/); 
				next  if ($SVchr =~ m/^hs37d5/);
				
				
				$SVdeb=$champs[1];
				my $qual = $champs[5];
		
				# recuperation des donnees de INFO : TYPE  END et length
				$SVend=0;
				$SVlength=0;
				$SVtype="-";
				
				foreach my $c (@champsINFO)
				{
					my ($key,$val)=split(/=/,$c);
					$SVtype = $val if ($key eq "SVTYPE");
					$SVend = $val if ($key eq "END");
					$SVlength = $val if ($key eq "SVLEN");
				}
				
		
				next unless ($qual > 400); 				# on ne garde que les lignes de bonne qualite
				next unless ($type =~ m/$SVtype/);		# pour l'instant on ne s'interresse que au DUP/DEL et INS
				next if ( abs($SVlength) < $length);	# on ne garde que les lignes qui de longueur > length
				next if ( ( $champs[6] ne "PASS")); 	# on ne garde que les lignes ou PASS
	
				$GT = $champsPAT[0];		

				# si ca passe ...				
				saveVariant($patname,$SVtype,$SVchr,$SVdeb,$SVend,$GT,$caller,\%hSVall);
				annotVariant($patname,$SVtype,$SVchr,$SVdeb,$SVend,$GT,$caller,\%hSVall); 			
			}
		}	
		
	} # fin de la boucle sur les callers
	
} #fin de la boucle sur les patients



#####################################
#  Sauvegarde de la table de hash resultante 
#  pour utilisation ulterieure
#####################################

#my $cacheDir = $project->getCacheStructuralVariantsDir();
#my $file_out =$path2.$projectname.".allSV";

store(\ %hSVall, $file_out) or die "Can't store $file_out!\n";
exit(0);


##############################################################################################################################################		


sub saveVariant()
{
	my ($patname,$SVtype,$SVchr,$SVdeb,$SVend,$GT,$Caller,$hresCaller) = @_;	
	
	my $copyNumber = 0;
	$copyNumber = $hCanvasCN{$patname}->{$SVtype.":".$SVchr.":".$SVdeb.":".$SVend}; 
	

	# pour trouver les genes compris dans l'interval
	my $objChr = $project->getChromosome($SVchr);
	my $tabGenes = $objChr->getGenesByPosition($SVdeb,$SVend);
	
	my $genes_liste="";
	
	foreach my $g (@$tabGenes)
	{
		$genes_liste .= $g->external_name.",";
	} 

	my $num;
	my $ch;
	if ($SVchr =~ m/chr/)
	{
			($ch,$num)  = split(/r/,$SVchr);
	}
	else
	{
			$num =$SVchr;
			$SVchr = "chr".$SVchr;
	}		
	my $id = $SVtype."_".$num."_".$SVdeb."_".$SVend;

	
	# ce qui depend du variant
	$hresCaller->{$id}->{'id'}=$id;
	$hresCaller->{$id}->{'SVLEN'}=abs($SVdeb-$SVend);
	$hresCaller->{$id}->{'SVTYPE'}=$SVtype;
	$hresCaller->{$id}->{'GENESLISTE'}="-";	
	$hresCaller->{$id}->{'GENES'}="-";
	$hresCaller->{$id}->{'SCOREMAX_GENES'}=0;
	$hresCaller->{$id}->{'CHROM'}=$SVchr;
	$hresCaller->{$id}->{'LIST_PATIENT'} .= $patname." ";
	
	
	# calculer le score gene max correspondant a la liste de gene associe au variant
	my $chr = $project->getChromosome($SVchr);
	my $genes = $chr->getGenesByPosition($SVdeb,$SVend);
	
	if ( scalar(@$genes) )  # si le variant recouvre des gènes
	{	
		my @names;
		my $max;
		foreach my $g (sort {$b->score <=> $a->score} @$genes){
			$max = $g->score unless $max;
			push (@names,$g->external_name);
		}
		
		$hresCaller->{$id}->{'SCOREMAX_GENES'}=$max;
		$hresCaller->{$id}->{'GENESLISTE'} =join(" ",@names);
		$hresCaller->{$id}->{'GENES'} =join(",",@names);
	}
	
	# ce qui depend du patient	
	my $pat = $project->getPatient($patname);
	
	my $pedinfo;
	$pedinfo = "child" if $pat->isChild();
	$pedinfo ="mother" if $pat->isMother();
	$pedinfo ="father" if $pat->isFather();
	my $status = $pat->status();
	
	$hresCaller->{$id}->{'PATIENT'} .= $id." ".$patname." ".$Caller." ".$GT." ".$pedinfo." ".$status.";";
	$hresCaller->{$id}->{$Caller}->{$patname}->{'GT'} = $GT;
	$hresCaller->{$id}->{$Caller}->{$patname}->{'CN'} = $copyNumber if ($copyNumber);
	
	
}

sub annotVariant()
{
	my ($patname,$SVtype,$SVchrFromVcf,$SVdeb,$SVend,$GT,$Caller,$hresCaller)=  @_;
	
	warn " Annote  $compteur :  $patname $Caller $SVtype $SVchrFromVcf $SVdeb $SVend";
		
	my $l = abs( $SVdeb-$SVend);
		
	my $num;
	my $ch;
	my $SVchr;
	
	
	if ($SVchrFromVcf =~ m/chr/)
	{
			($ch,$num)  = split(/r/,$SVchrFromVcf);
			$SVchr=$SVchrFromVcf;
	}
	else
	{
			$num =$SVchrFromVcf;
			$SVchr = "chr".$SVchrFromVcf;
	}		
	my $id = $SVtype."_".$num."_".$SVdeb."_".$SVend;
		
	# pour voir si le variant a ete vu avec les autres caller et definir les points de cassure
	my $lumpy_bp = "";	
	my $manta_bp = "";	
	$lumpy_bp = getBreakpoint($SVchrFromVcf,$SVdeb,$SVend,"lumpy",$patname,$hpatientFiles{$patname}->{'lumpy'}->{'vcf'}) if ( (-e $hpatientFiles{$patname}->{'lumpy'}->{'vcf'}) );	
	$manta_bp = getBreakpoint($SVchrFromVcf,$SVdeb,$SVend,"manta",$patname,$hpatientFiles{$patname}->{'manta'}->{'vcf'}) if ( (-e $hpatientFiles{$patname}->{'manta'}->{'vcf'}) );	
	
		
		#############################
		#  annotations a partir de AnnotSV
		#############################
		
		#Genes
#		my $liste_of_genes = "-";
#		
#		$liste_of_genes = getAnnots($patname,$Caller,$SVchr,$SVdeb,$SVend,$hannot->{$Caller}->{"Gene name"});
#		my @tabgenes = split(/\//,$liste_of_genes);
#		$liste_of_genes =~ s/\//,/g;
		
		
		#DGV
		my $dgv_gain = getAnnots($patname,$Caller,$SVchr,$SVdeb,$SVend,$hannot->{$Caller}->{'DGV_GAIN_IDs'});
		my $dgv_loss = getAnnots($patname,$Caller,$SVchr,$SVdeb,$SVend,$hannot->{$Caller}->{'DGV_LOSS_IDs'});
		my $dgv_gain_freq = getAnnots($patname,$Caller,$SVchr,$SVdeb,$SVend,$hannot->{$Caller}->{'DGV_GAIN_Frequency'});
		my $dgv_loss_freq = getAnnots($patname,$Caller,$SVchr,$SVdeb,$SVend,$hannot->{$Caller}->{'DGV_LOSS_Frequency'});

		#SV_rank
		my $SV_rank = getAnnots($patname,$Caller,$SVchr,$SVdeb,$SVend,$hannot->{$Caller}->{'AnnotSV ranking'});	
		$hresCaller->{$id}->{'RANKAnnot'}=$SV_rank;
		$hresCaller->{$id}->{'RANKAnnot'} ++ if ($hresCaller->{$id}->{'SCOREMAX_GENES'} > 4);	
		$hresCaller->{$id}->{'RANKAnnot'} .= ",".$hresCaller->{$id}->{'SCOREMAX_GENES'};
		
		#DDD
		my $DDD_DEL_freq =  getAnnots($patname,$Caller,$SVchr,$SVdeb,$SVend,$hannot->{$Caller}->{'DDD_DEL_Frequency'});
		my $DDD_DUP_freq =  getAnnots($patname,$Caller,$SVchr,$SVdeb,$SVend,$hannot->{$Caller}->{'DDD_DUP_Frequency'});
		
		my $OMIN_MG = getAnnots($patname,$Caller,$SVchr,$SVdeb,$SVend,$hannot->{$Caller}->{'morbidGenes'});	
		
		#dbVar
		my $dbVar_event =  getAnnots($patname,$Caller,$SVchr,$SVdeb,$SVend,$hannot->{$Caller}->{'dbVar_event'});
		$dbVar_event =~ s/;/ /g;
		my $dbVar_status =  getAnnots($patname,$Caller,$SVchr,$SVdeb,$SVend,$hannot->{$Caller}->{'dbVar_status'});
		$dbVar_status =~ s/;/ /g;	
		
#		$hresCaller->{$id}->{'GENES'}=$liste_of_genes;
		$hresCaller->{$id}->{'GOLD'}=$dgv_gain." ".$dgv_loss;
		$hresCaller->{$id}->{'GOLD_G_freq'}=$dgv_gain_freq;
		$hresCaller->{$id}->{'GOLD_L_freq'}=$dgv_loss_freq;
		
		$hresCaller->{$id}->{'DDD_DEL_FREQ'}=$DDD_DEL_freq;
		$hresCaller->{$id}->{'DDD_DUP_FREQ'}=$DDD_DUP_freq;
		$hresCaller->{$id}->{'OMIN_MG'}=$OMIN_MG;
		$hresCaller->{$id}->{'dbVar_event'}=$dbVar_event;
		$hresCaller->{$id}->{'dbVar_status'}=$dbVar_status;
		$hresCaller->{$id}->{'BREAKPOINTS'} .= $lumpy_bp.$manta_bp;
		$compteur++;
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



sub getAnnots
{
	my ($patname,$caller,$chr,$d,$f,$num) = @_;
	
	my @n = split(/chr/,$chr);
	my $chrnum = $n[1];
	
	my $fichier_Annot = $hpatientFiles{$patname}->{$caller}->{'annotSV'};
	
	my $cmd =  "cat ".$fichier_Annot." | grep ".$chrnum."_".$d."_".$f." | grep full";
 	my $res1 = `$cmd`;
 	chomp($res1);
 	
 	my @champsAnnot = split(/\t/,$res1);

 	return $champsAnnot[$num] if $champsAnnot[$num];
 	return "-";
}

sub getBreakpoint
{
	my ($chr,$deb,$end,$caller,$patname,$fichier) = @_;
	
	my $b1 = $deb-5000;
	my $b2 = $deb+5000;
	
	my $cmd =  "tabix ".$fichier." ".$chr.":".$b1."-".$b2;
 	my $res1 = `$cmd`;
 	chomp($res1);
 	 
 	my $out = "-";
 	
 
 	if ($res1)
 	{
		# determiner le nombre de SV present dans l'interval 
		my @champsSV = split(/\n/,$res1);
		my $nbSV = scalar(@champsSV);
		
		$out =  "";
		
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
 			$out .= $caller.":".$patname."_".$bptype.":".$bpchr.":".$bppos."-".$bpend."_".$bpalt." " if (( abs($bppos - $deb) <= 5000) && ($bptype eq "BND"));
 		}
		
 	}
	return $out;
 	
}
	


 	