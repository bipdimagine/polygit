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
use Clone qw(clone);
use Parallel::ForkManager;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use Getopt::Long;

use GBuffer;



###################################################################
# Cherche tous les variants structuraux d'un projet pour chaque patient et tous callers confondus.
# Reconstruit les CNV fragmentes (même caller et bornes a moins de 10% de la longeur)
# Construit une table de hash par patient et freeze ces tables  : nom du fihier = patient.allSV
# Pour le dejavu construit egalement une table de hash qui conserve tous les CNV du projet 
# et garde pour chaque CNV l'info du patient et du caller qui l'a détecté.
# Pour chaque CNV la localisation chr_start_end donne accès à la liste des gènes compris dans cet interval
####################################################################

my $fork = 5;
my $cgi = new CGI;


my $limit;
my $projectname;
my $patient_name;
my $fork =1 ;
GetOptions(
	'project=s' => \$projectname,
#	'patient=s' => \$patient_name,
	'fork=s' => \$fork,
);

 #= $cgi->param('projectname');
#my $patient_name = $cgi->param('patient');
#$fork = $cgi->param('fork');

#$fork=1;

# pour récupérer les objets project et patient
my $buffer = GBuffer->new();

my $project = $buffer->newProjectCache( -name => $projectname);

my $listCallers = $project->callingSVMethods();
my $listPatients = $project->getPatients();



# on gardera que les DUP DEL de taille supérieur à 1kb
#my $length = 1000;


# pour le fichier annoSV
my $hannot;
my $fdannot;

# pour lire les vcfs
my $ligne;
my $fd;
my $compteur=0;
	
# pour stocker les SV
my $hPat_elementaryCNV;  # tous les CNV du vcf  (élémentaires = avant regroupement)
my $hPat_CNV; # pour les id_globaux

# pour le dejavu
my $hdejavu;



#########################################
#  Pour les duplications segmentales
#########################################
my $cytodir = $project->dirCytoManue;	
my $fichier_dupseg = $cytodir."/super_Duplication.gff3";
warn $fichier_dupseg;
my $fdds;

#  Lecture du fichier dupseg
	my $hdupseg;
	open($fdds,$fichier_dupseg) or die("open: $!");
	
	# lecture ligne par ligne
	while( defined( $ligne = <$fdds> ) )
	{
		next unless ($ligne =~ m/^chr/ );
		my @champs = split(/\t/,$ligne);
		my @info = split(/;/,$champs[8]);
		
		my ($n,$idd) = split(/=/,$info[0]);
		my ($i,$percentid) = split(/=/,$info[1]);
		
		my $chr = $champs[0];
		my $start = $champs[3];
		my $end = $champs[4];
		
		my $idLocation = $idd.":".$start."_".$end;
		
		# version avec IntervalTree
		if (!exists $hdupseg->{$chr})
		{
		 	my $DupSegTree = Set::IntervalTree->new;
		 	$hdupseg->{$chr}=$DupSegTree;
		}		
		$hdupseg->{$chr}->insert($idLocation,$start,$end);		
	}	
	
	
############
#  lecture des vcf 
############


my $SVtype;
my $SVend ;
my $SVlength;
my $SVchr;
my $SVdeb;


# pour la sauvegarde des fichiers de sortie
#my $variationsDir = "/data-xfs/Manue/Test_SV/".$projectname."/newVersion/";
#my $variationsDir= $project->getVariationsDir();
my $variationsDir= $project->getCNVDir();
my $pm = new Parallel::ForkManager($fork);
$project->getChromosomes();		
my $job_id = time;
my $hjobs ;
$pm->run_on_finish(
    	sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    		my $j = $data->{job};
    		delete $hjobs->{$j};
    		
    }
    );
    $project->buffer->dbh_deconnect();
foreach my $patobj (@$listPatients)
{
	warn "------------";
	warn $patobj->name;
	warn "---------------";
	next if  $patient_name && $patobj->name ne $patient_name ;
	#or $patobj->name ne "dl-2-E-sg-A";
	next unless $patobj->isGenome;
	$job_id ++;
	$hjobs->{$job_id} ++;
	my $pid = $pm->start and next;
	$project->buffer->dbh_reconnect();
	my $patname= $patobj->name();
	my $file_out = $variationsDir.$patname.".allCNV";
	unlink  $file_out if -e $file_out;
	#next if (-e $file_out);
	
	
	foreach my $caller (@$listCallers)
	{
		next if $caller eq "lumpy";
		my $GT = "-";
		my $CN = "-";
		my $WC_ratio=0;
		my $WC_zscore=0;
		my $QUAL=0;
		
		my $dir = $project->getVariationsDir($caller);
	
		#    Lecture de la première ligne du fichier Annot 
		#    pour recuperer le format

		my $fichier_Annot = $patobj->getAnnotSVFileName($caller);
				
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
				
		my $fichierPatient = $patobj->getSVFile($caller);
				
		if ($caller  eq "wisecondor")
		{
			# ouverture du fichier wisecondor et parsing
			open($fd," zcat $fichierPatient | ") or die("open: $!");
			
			# lecture de la première ligne
			$ligne = <$fd>;

			# traitement des lignes suivantes pour le fichier en entier
			while( defined( $ligne = <$fd> ) )
			{
						my @champs = split(/\t/,$ligne);
		
						$SVchr="chr".$champs[0];
						$SVdeb=$champs[1];
						$SVend=$champs[2];
						$WC_ratio=$champs[3];
						$QUAL=$champs[4];
						$SVtype = $champs[5];
		
				  		$SVtype = "DEL" if  ($SVtype =~ m/loss/);
						$SVtype = "DUP" if ($SVtype =~ m/gain/);
		
						# on ne garde que les lignes qui de longueur > length
						#my $SVlength = $SVend-$SVdeb;
						#next if ( abs($SVlength) < $length);		
						
						# si ca passe ...	
						saveVariant($patobj,$caller,$SVtype,$SVchr,$SVdeb,$SVend,$GT,$CN,$WC_ratio,$QUAL);
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
					$CN =~ s/<//g;
					$CN =~ s/>//g;
					
					$QUAL = $champs[5];
		 
					my @canvasinfo=split(/:/,$champs[2]);
					$SVtype = $canvasinfo[1];
		
					$SVtype = "DEL" if ($SVtype eq "LOSS");
					$SVtype = "DUP" if ($SVtype eq "GAIN");
					$SVtype = "LOH" if ($SVtype eq "LOH");
		

					next unless ( ($SVtype eq "DUP") || ($SVtype eq "DEL") );				# pour l'instant on ne s'interresse que au DUP/DEL 
					next if ( ( $champs[6] ne "PASS"));													# on ne garde que les lignes ou PASS
		
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
					#next if ( abs($SVlength) < $length);		
					
					my $GT = $champsPAT[0];
					
					# si ca passe ...	
					saveVariant($patobj,$caller,$SVtype,$SVchr,$SVdeb,$SVend,$GT,$CN,$WC_ratio,$QUAL);
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
				next  if ($SVchr =~ m/^NC_/);
				
				$SVdeb=$champs[1];
				$QUAL = $champs[5];
		
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
				
		
				#next unless ($QUAL > 400); 																# on ne garde que les lignes de bonne qualite
				next unless ( ($SVtype eq "DUP") || ($SVtype eq "DEL") );				# pour l'instant on ne s'interresse que au DUP/DEL 
				#next if ( abs($SVlength) < $length);													# on ne garde que les lignes qui de longueur > length
				next if ( ( $champs[6] ne "PASS")); 													# on ne garde que les lignes ou PASS
	
				$GT = $champsPAT[0];		

				# si ca passe 			
				saveVariant($patobj,$caller,$SVtype,$SVchr,$SVdeb,$SVend,$GT,$CN,$WC_ratio,$QUAL);
				
			} # fin du fichier
		}	# fin du if
		
	} # fin de la boucle sur les callers
	
	#on recopie tels quel sles CNV manta
	$hPat_CNV->{$patname}->{'manta'} = clone($hPat_elementaryCNV->{$patname}->{'manta'});
	
	# gestion des CNV fragmentes par wisecondor et canvas
	gatherCNV_from_samecaller($patname,"wisecondor");
	gatherCNV_from_samecaller($patname,"canvas");
	
	# on recupere la liste des genes et le score max corespondant au CNV
	# et les infos de cytogenetique : duplications segmentaires et cytoband
	setComplementaryInfos($patname);	
	# on freeze la table correspondant à chaque patient individuelement
	store(\ %{$hPat_CNV->{$patname}}, $file_out) or die "Can't store $file_out for ".$patname."!\n";
$pm->finish(0,{job=>$job_id});
} # fin boucle sur les patients
$pm->wait_all_children();
confess() if scalar(keys %{$hjobs});
$project->buffer->dbh_reconnect();

###############
# pour le dejavu
###############

foreach my $patobj (@$listPatients)
{
	next unless $patobj->isGenome;
	my $patname= $patobj->name();
	my $hPat_allCNV;
	
	warn "dejavu ".$patname;
	
	my $CNVfile = $variationsDir.$patname.".allCNV";
	my $hPat_allCNV = retrieve($CNVfile) or die "Can't retrieve datas from " . $CNVfile . " !\n";
	
	foreach my $caller (keys %{$hPat_allCNV})
	{
		foreach my $type (keys %{$hPat_allCNV->{$caller}})
		{
			foreach my $num (keys %{$hPat_allCNV->{$caller}->{$type}})
			{
				foreach my $global_id  (keys %{$hPat_allCNV->{$caller}->{$type}->{$num}})
				{
						$hdejavu->{$type}->{$num}->{$global_id}->{$patname}->{$caller}++;
				}
			}
		}
	}
}
						
						
						
		
# freeze la table du dejavu pour le projet
my $file_dejavu_inthisproject = $variationsDir.$projectname."_dejavu.allCNV";
my $file_dejavu = $project->DejaVuCNVFile();
store(\ %{$hdejavu}, $file_dejavu_inthisproject) or die "Can't store $file_dejavu_inthisproject!\n";
store(\ %{$hdejavu}, $file_dejavu) or die "Can't store $file_dejavu!\n";

exit(0);


##################################
#   methodes
####################################		


sub saveVariant()
{
	my ($patobj,$Caller,$SVtype,$SVchr,$SVdeb,$SVend,$GT,$CN,$WC_ratio,$QUAL) = @_;
	my $patname = $patobj->name();
	
	print $patobj->name." Save : $compteur :  $patname $Caller $SVtype $SVchr $SVdeb $SVend \n" if $compteur %100 ==0 ;
	
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

	# on stocke dans la table $hPat_elementaryCNV les variants élémentaires
	$hPat_elementaryCNV->{$patname}->{$Caller}->{$SVtype}->{$num}->{$id}->{'id'}=$id;
	$hPat_elementaryCNV->{$patname}->{$Caller}->{$SVtype}->{$num}->{$id}->{'SVTYPE'}=$SVtype;
	$hPat_elementaryCNV->{$patname}->{$Caller}->{$SVtype}->{$num}->{$id}->{'CHROM'}=$num;
	$hPat_elementaryCNV->{$patname}->{$Caller}->{$SVtype}->{$num}->{$id}->{'START'}=$SVdeb;
	$hPat_elementaryCNV->{$patname}->{$Caller}->{$SVtype}->{$num}->{$id}->{'END'}=$SVend;
	$hPat_elementaryCNV->{$patname}->{$Caller}->{$SVtype}->{$num}->{$id}->{'SVLEN'}=abs($SVdeb-$SVend);
	$hPat_elementaryCNV->{$patname}->{$Caller}->{$SVtype}->{$num}->{$id}->{'GT'}=$GT;
	$hPat_elementaryCNV->{$patname}->{$Caller}->{$SVtype}->{$num}->{$id}->{'CN'}=$CN;
	$hPat_elementaryCNV->{$patname}->{$Caller}->{$SVtype}->{$num}->{$id}->{'RATIO'}=$WC_ratio;
	$hPat_elementaryCNV->{$patname}->{$Caller}->{$SVtype}->{$num}->{$id}->{'QUAL'}=$QUAL;
	$hPat_elementaryCNV->{$patname}->{$Caller}->{$SVtype}->{$num}->{$id}->{'ELEMENTARY'}= $id;
	
	annotVariant($patobj,$Caller,$SVtype,$SVchr,$SVdeb,$SVend,$GT); 
}

sub annotVariant()
{
	my ($patobj,$Caller,$SVtype,$SVchrFromVcf,$SVdeb,$SVend,$GT) = @_;
	my $patname = $patobj->name();
	my $l = abs( $SVdeb-$SVend);
	my $num;
	my $ch;
	my $SVchr;
	
	#warn " Annot : $compteur :  $patname $Caller $SVtype $SVchrFromVcf $SVdeb $SVend";
	
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
		
	
		#  annotations a partir de AnnotSV

		#DGV
		my $dgv_gain_freq = getAnnots($patobj,$Caller,$SVchr,$SVdeb,$SVend,$hannot->{$Caller}->{'DGV_GAIN_Frequency'});
		my $dgv_loss_freq =  getAnnots($patobj,$Caller,$SVchr,$SVdeb,$SVend,$hannot->{$Caller}->{'DGV_LOSS_Frequency'});
		
	
		$dgv_gain_freq = -1 if ($dgv_gain_freq eq "-");
		$dgv_loss_freq = -1 if ($dgv_loss_freq eq "-");
		
		#SV_rank
		my $SV_rank = getAnnots($patobj,$Caller,$SVchr,$SVdeb,$SVend,$hannot->{$Caller}->{'AnnotSV ranking'});	

		#OMIM
		my $OMIN_MG = getAnnots($patobj,$Caller,$SVchr,$SVdeb,$SVend,$hannot->{$Caller}->{'morbidGenes'});	
		
		#dbVar
		my $dbVar_status =  getAnnots($patobj,$Caller,$SVchr,$SVdeb,$SVend,$hannot->{$Caller}->{'dbVar_status'});
		$dbVar_status =~ s/;/ /g;	
		
		# on stocke
		$hPat_elementaryCNV->{$patname}->{$Caller}->{$SVtype}->{$num}->{$id}->{'GOLD_G_FREQ'}= $dgv_gain_freq;
		$hPat_elementaryCNV->{$patname}->{$Caller}->{$SVtype}->{$num}->{$id}->{'GOLD_L_FREQ'}= $dgv_loss_freq;
		$hPat_elementaryCNV->{$patname}->{$Caller}->{$SVtype}->{$num}->{$id}->{'OMIN_MG'}= $OMIN_MG;
		$hPat_elementaryCNV->{$patname}->{$Caller}->{$SVtype}->{$num}->{$id}->{'dbVar_status'}= $dbVar_status;
		$hPat_elementaryCNV->{$patname}->{$Caller}->{$SVtype}->{$num}->{$id}->{'RANKAnnot'}= $SV_rank;
		$compteur++;
}

sub gatherCNV_from_samecaller()
{
	my ($patname,$caller) = @_;
	
	my $htree;
	my $hintspan;

	# 1)  detecter les SV chevauchants 
	
	# creer les arbres
	foreach my $type (keys %{$hPat_elementaryCNV->{$patname}->{$caller}})
	{
			foreach my $num (keys %{$hPat_elementaryCNV->{$patname}->{$caller}->{$type}})
			{
				foreach my $id  (keys %{$hPat_elementaryCNV->{$patname}->{$caller}->{$type}->{$num}})
				{
						my ($t,$c,$d,$f) = split( /_/, $id );
						$htree->{$t}->{$c}= Set::IntervalTree->new;
						$hintspan->{$t}->{$c}= Set::IntSpan::Fast->new();
				}
			}
	}	

	# remplir les arbres :  regrouper les SV chevauchants
	foreach my $type (keys %{$hPat_elementaryCNV->{$patname}->{$caller}})
	{
			foreach my $num (keys %{$hPat_elementaryCNV->{$patname}->{$caller}->{$type}})
			{
				foreach my $id  (keys %{$hPat_elementaryCNV->{$patname}->{$caller}->{$type}->{$num}})
				{
						my ($t,$c,$d,$f) = split( /_/, $id );
						$htree->{$t}->{$c}->insert($id,$d,$f);
				}
			}
	}
	
	# 2) associer a chaque SV trouve avec le même caller ceux qui lui sont proches	
	foreach my $type (keys %{$hPat_elementaryCNV->{$patname}->{$caller}})
	{
			foreach my $num (keys %{$hPat_elementaryCNV->{$patname}->{$caller}->{$type}})
			{
				foreach my $id  (keys %{$hPat_elementaryCNV->{$patname}->{$caller}->{$type}->{$num}})
				{
		
						#on cherche les regroupements
						my ($t,$c,$dtheSV,$ftheSV) = split( /_/, $id );
		
						my $padding = 0.1 * abs($ftheSV-$dtheSV);
		
						my $tab_id = $htree->{$t}->{$c}->fetch($dtheSV-$padding,$ftheSV+$padding);

						# on regroupe les id dans un intspan
						my $gdeb=0;
						my $gend=0;
			
						foreach my $ind_id ( @$tab_id ) 
						{
							my ($t,$c,$d,$f) = split(/_/, $ind_id);
				
		 				   	$gdeb = $d if ( ($d < $gdeb) || ($gdeb==0));
		 					$gend = $f if ( ($f > $gend) );
 						}
						$hintspan->{$t}->{$c}->add_range($gdeb,$gend);
				}
			}
	}
	gather_id($htree,$hintspan,$patname,$caller);
}

sub gather_id() {
	
		my ($htree,$hintspan,$patname,$caller) = @_;
		my $nb = 0;

		foreach my $type (keys %{$hintspan})
		{
			foreach my $num (keys %{$hintspan->{$type}})
			{
				my $liste_of_bornes = $hintspan->{$type}->{$num}->as_string() if defined( $hintspan->{$type}->{$num});
				$liste_of_bornes .= "," if ($liste_of_bornes !~ m/,/);
				my @tab_bornes = split(/,/,$liste_of_bornes);
				
				foreach my $bornes (@tab_bornes)
				{
					my ($start,$end) = split(/-/,$bornes);
					
					my $tab_ind_id = $htree->{$type}->{$num}->fetch($start,$end);
					
					my $global_id = $type."_".$num."_".$start."_". $end;
		
					# enregistrer les infos dans la table de hash 
					$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'id'} = $global_id;			
					$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'SVTYPE'} = $type;
					$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'CHROM'} = $num;
					$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'START'} = $start;
					$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'END'} = $end;
					$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'SVLEN'}  = abs( $start - $end );
		
					# retrouver les informations correspondant aux differents id regroupes 
					# ici il faut recuperer les valeurs max ou les concatenation de liste
					
					my $ggfreq = -1;
					my $glfreq = -1;
					my $mg ="no";
					my $dbVar_status = "";
					my $rankannot = 0;
					my $gt = "";
					my $elementary_ids = "";
					my $cn;
					my $ratio;
					my $qual;
					my $index;
					
					foreach my $id ( @$tab_ind_id )
					{					
						$ggfreq = $hPat_elementaryCNV->{$patname}->{$caller}->{$type}->{$num}->{$id}->{'GOLD_G_FREQ'} if ($hPat_elementaryCNV->{$patname}->{$caller}->{$type}->{$num}->{$id}->{'GOLD_G_FREQ'} > $ggfreq);
						$glfreq = $hPat_elementaryCNV->{$patname}->{$caller}->{$type}->{$num}->{$id}->{'GOLD_L_FREQ'} if ($hPat_elementaryCNV->{$patname}->{$caller}->{$type}->{$num}->{$id}->{'GOLD_L_FREQ'} > $glfreq);
						$mg = "yes" if ($hPat_elementaryCNV->{$patname}->{$caller}->{$type}->{$num}->{$id}->{'OMIN_MG'} eq "yes");
						$dbVar_status .= $hPat_elementaryCNV->{$patname}->{$caller}->{$type}->{$num}->{$id}->{'dbVar_status'}." ";
						$rankannot = $hPat_elementaryCNV->{$patname}->{$caller}->{$type}->{$num}->{$id}->{'RANKAnnot'} if ($hPat_elementaryCNV->{$patname}->{$caller}->{$type}->{$num}->{$id}->{'RANKAnnot'} > $rankannot);
						$gt .= $hPat_elementaryCNV->{$patname}->{$caller}->{$type}->{$num}->{$id}->{'GT'}." " unless ($gt =~ m/$hPat_elementaryCNV->{$patname}->{$caller}->{$type}->{$num}->{$id}->{'GT'}/) ;
						$cn .= $hPat_elementaryCNV->{$patname}->{$caller}->{$type}->{$num}->{$id}->{'CN'}." " unless ($cn =~ m/$hPat_elementaryCNV->{$patname}->{$caller}->{$type}->{$num}->{$id}->{'CN'}/) ;
						$ratio += $hPat_elementaryCNV->{$patname}->{$caller}->{$type}->{$num}->{$id}->{'RATIO'};
						$qual += $hPat_elementaryCNV->{$patname}->{$caller}->{$type}->{$num}->{$id}->{'QUAL'};
						$elementary_ids .= $id." ";
					}
					
					$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'ELEMENTARY'}= $elementary_ids;
					$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'GOLD_G_FREQ'} = $ggfreq;
					$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'GOLD_L_FREQ'} = $glfreq;
					$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'OMIN_MG'} = $mg;
					$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'dbVar_status'} = $dbVar_status;
					$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'RANKAnnot'} = $rankannot;
					$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'GT'}=$gt;
					$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'CN'}=$cn;
					$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'RATIO'}=$ratio/scalar(@$tab_ind_id);			#on prend la moyenne
					$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'QUAL'}=$qual/scalar(@$tab_ind_id);			#on prend la moyenne
			}
		}
	}
}


sub setComplementaryInfos()
{
	my ($patname) =@_;
	
	foreach my $caller (keys %{$hPat_CNV->{$patname}})
	{
		foreach my $type (keys %{$hPat_CNV->{$patname}->{$caller}})
		{
			foreach my $num (keys %{$hPat_CNV->{$patname}->{$caller}->{$type}})
			{
				foreach my $global_id  (keys %{$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}})
				{
						my ($t,$c,$start,$end) = split(/_/, $global_id);
						
						# pour trouver les genes compris dans l'interval
						die($c." ".$global_id) if ($c ne $num);
						my $chr = "chr".$num;
						
						#TODO: chromosome MT next
						next if ($chr eq 'chrMT');
						my $objChr = $project->getChromosome($chr);
						
						my $tabGenes = $objChr->getGenesByPosition($start,$end);
	
						my $genes_liste="";
	
						foreach my $g (@$tabGenes)
						{
							$genes_liste .= $g->external_name.",";
						} 


						# calculer le score gene max correspondant a la liste de gene associe au variant
						$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'GENES'} = " ";
						$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'SCORE_GENES'} = 0;
						if ( scalar(@$tabGenes) )  # si le variant recouvre des gènes
						{	
							my @names;
							my $max;
							foreach my $g (sort {$b->score <=> $a->score} @$tabGenes){
										$max = $g->score unless $max;
										push (@names,$g->external_name);
							}
		
							$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'SCORE_GENES'}=$max;
							$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'BEST_GENE'}= $names[0];
							$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'GENES'} =join(" ",@names);
						}
						
						# pour detecter la presence de duplication segmentaire
						$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'DUPSEG'}= getDupSeg($chr,$start,$end);
						
						# pour les cytobandes
						my @tb;
						my @band;
						my $hband = $project->getChromosome($num)->getCytoband( $start, $end ) if ( $end> $start);

						foreach my $b ( keys %{ $hband->{'name'} } ) {
								push( @tb, $b );
						}
						@band = sort( { $a cmp $b } @tb );
						$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'CYTOBAND'} = join( ",", @band );
						
						# pour l'acces a DGV
						my $url_DGV = getDGV_url( $chr, $start, $end );
						$hPat_CNV->{$patname}->{$caller}->{$type}->{$num}->{$global_id}->{'DGV'} = $url_DGV;
						
						# pour le dejavu in this project
						#$hdejavu->{$type}->{$num}->{$global_id}->{$patname}->{$caller}++;
				}
			}
		}
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



sub getAnnots
{
	my ($patobj,$caller,$chr,$d,$f,$num) = @_;
	
	my @n = split(/chr/,$chr);
	my $chrnum = $n[1];
	
	
	my $patname=$patobj->name();
	my $fichier_Annot = $patobj->getAnnotSVFileName($caller);
	
	my $cmd =  "cat ".$fichier_Annot." | grep ".$chrnum."_".$d."_".$f." | grep full";
 	my $res1 = `$cmd`;
 	chomp($res1);
 	
 	my @champsAnnot = split(/\t/,$res1);

 	return $champsAnnot[$num] if $champsAnnot[$num];
 	return "-";
}



sub getDGV_url
{
	my ($chr,$d,$f) = @_;
	my $url = "http://dgv.tcag.ca/gb2/gbrowse/dgv2_hg19/?name=".$chr."%3A".$d."-".$f.";search=Search";
	return $url;
}

sub getDupSeg
{
	my ($chr,$d,$f) = @_;
	my $dupseg;
	my $result;

	
	# recherche des dupseg chevauchant le fragment
	# exemple eval pour chrMT qui ne fetch pas	
	eval { $result = $hdupseg->{$chr}->fetch($d,$f); };
	if($@) {
		warn "\n\n\nERROR: can't fetch $chr... NEXT...\n\n\n";
		return '-';
	}
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

 	