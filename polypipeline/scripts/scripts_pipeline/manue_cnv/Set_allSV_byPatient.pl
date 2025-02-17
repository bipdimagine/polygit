#!/usr/bin/perl
use FindBin qw($Bin);
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
use Text::CSV qw( csv );
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use Getopt::Long;
use lib "$Bin";
require  "$Bin/SVParser.pm";
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
	'patient=s' => \$patient_name,
	'fork=s' => \$fork,
);

 #= $cgi->param('projectname');
#my $patient_name = $cgi->param('patient');
#$fork = $cgi->param('fork');

#$fork=1;

# pour récupérer les objets project et patient
my $buffer = GBuffer->new();

my $project = $buffer->newProjectCache( -name => $projectname);


my $listPatients = $project->getPatients();
if ($patient_name) {
	$listPatients = [$project->getPatient($patient_name)];
}



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

my $hPat_CNV; # pour les id_globaux

# pour le dejavu
my $hdejavu;



#########################################
#  Pour les duplications segmentales
#########################################
my $cytodir = $project->get_cytology_directory;	
my $fichier_dupseg = $cytodir."/segmental_duplication.bed";
my $fdds;

#  Lecture du fichier dupseg
	my $hdupseg;
	open($fdds,$fichier_dupseg) or die("open: $!");
	
	# lecture ligne par ligne
	while( defined( $ligne = <$fdds> ) )
	{
		next unless ($ligne =~ m/^chr/ );
		my @champs = split(/\t/,$ligne);
		
		my $chr = $champs[0];
		my $start = $champs[1];
		my $end = $champs[2];
		
		my $idLocation = $chr.":".$start."_".$end;
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
my $variationsDir = $project->getCNVDir();

my $pm = new Parallel::ForkManager($fork);
$project->getChromosomes();		
my $job_id = time;
my $hjobs ;
my $chrs;
foreach my $c (@{$project->getChromosomes}){
	next if $c->name eq "MT";
	$chrs->{$c->name} ++;
	$chrs->{$c->ucsc_name} ++;
}

$pm->run_on_finish(
    	sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    		my $j = $data->{job};
    		delete $hjobs->{$j};
    		
    }
    );
  $project->buffer->hash_genes_omim_morbid();
  $project->setPhenotypes();
    
  $project->disconnect();
  $buffer->dbh_deconnect();
  
  
  foreach my $patobj (@$listPatients)
{
	my $listCallers = $patobj->callingSVMethods();
	 $patobj->isGenome;
}
  $project->disconnect();
  $buffer->dbh_deconnect();

foreach my $patobj (@$listPatients)
{
	warn "------------";
	warn $patobj->name;
	warn "---------------";
	my $listCallers = $patobj->callingSVMethods();
	#or $patobj->name ne "dl-2-E-sg-A";
	next unless $patobj->isGenome;
	$job_id ++;
	$hjobs->{$job_id} ++;
	my $pid = $pm->start and next;
	$listCallers = $patobj->callingSVMethods();
	$project->buffer->dbh_reconnect();
	my $patname= $patobj->name();
	my $file_out = $variationsDir.$patname.".allCNV";
	unlink  $file_out if -e $file_out;
	#next if (-e $file_out);
	my $patient = $patobj;
	my $hPat_elementaryCNV;  # tous les CNV du vcf  (élémentaires = avant regroupement)
	warn Dumper @$listCallers;
	foreach my $caller (@$listCallers)
	{
		next if $caller eq "lumpy";
		warn "+++".$caller;
		
		my $dir = $project->getVariationsDir($caller);
	
		#    Lecture de la première ligne du fichier Annot 
		#    pour recuperer le format

		my $fichier_Annot = $patobj->getAnnotSVFileName($caller);
			warn $fichier_Annot;
		if ($fichier_Annot)
		{
			open($fdannot,$fichier_Annot) or die("*************  open: $!  $fichier_Annot");
	
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
			close($fdannot);
		}
		my $fichierPatient = $patobj->getSVFile($caller);
		my $res;	
		warn $caller;
		if ($caller  eq "wisecondor")
		{
	
				my $hash = SVParser::parse_wisecondor($patient);
				warn Dumper $hash;
				$hPat_CNV->{'wisecondor'} = SVParser::gatherCNV_from_samecaller($patname,$hash);
		}
		elsif ($caller  eq "manta" or $caller  eq "pbsv" or $caller  eq "dragen-sv") {
			$hPat_CNV->{'manta'}  = SVParser::parse_vcf($patient,$caller);
		
		}
		elsif ($caller  eq "hificnv" or $caller  eq "canvas" or $caller  eq "dragen-cnv") {
			my $hash  = SVParser::parse_vcf($patient,$caller);
			$hPat_CNV->{'canvas'} = SVParser::gatherCNV_from_samecaller($patname,$hash);
			warn Dumper $hPat_CNV->{'canvas'};
		}
	}
	# on recupere la liste des genes et le score max corespondant au CNV
	# et les infos de cytogenetique : duplications segmentaires et cytoband
	
	$project->buffer->disconnect();
	#$project->buffer->dbh_reconnect();
	
	
	
	
	warn Dumper keys %$hPat_CNV;
	warn $file_out;
	# on freeze la table correspondant à chaque patient individuelement
	store($hPat_CNV, $file_out) or die "Can't store $file_out for ".$patname."!\n";
	
	
	
#	warn "--- end gather --- ";
##	my $file_out_gather = $variationsDir.$patname.".gatherCNV";
#	unlink  $file_out_gather if -e $file_out_gather;
#	store($gather, $file_out_gather) or die "Can't store $file_out_gather for ".$patname."!\n";
#	warn "end freeze";
	$pm->finish(0,{job=>$job_id});
} # fin boucle sur les patients
$pm->wait_all_children();
warn "END PATIENT !!!";
confess() if scalar(keys %{$hjobs});
$project->buffer->dbh_reconnect();

###############
# pour le dejavu
###############
warn "end fork ";
warn "=================";
warn "=================";
warn "=================";
my $hPat_allCNV;

foreach my $patObj (@{$project->getPatients()})
{
	next unless $patObj->isGenome;
	my $patname= $patObj->name();
	
	
	my $CNVfile = $variationsDir.$patname.".allCNV";
	 $hPat_allCNV->{$patObj->id} = retrieve($CNVfile) or die "Can't retrieve datas from " . $CNVfile . " !\n";
	
	foreach my $caller (keys %{$hPat_allCNV->{$patObj->id}})
	{
		
				foreach my $global_id  (keys %{$hPat_allCNV->{$patObj->id}->{$caller}})
				{
						$hdejavu->{$global_id}->{$patname}->{$caller}++;
						#$csv->print($fh, [$global_id,"$type",$num,$patname,$caller]) or die "Erreur d'écriture dans le fichier CSV : " . $csv->error_diag();
				}
		
	}
}

warn $variationsDir."/dejavu/";	
				
my $nodejavu = GenBoNoSqlDejaVuCNV->new( dir =>$variationsDir."/dejavu/", mode => "c" );	
foreach my $chr (@{$project->getChromosomes}){
	$nodejavu->create_table($chr->name);
}
foreach my $gb (keys %{$hdejavu}){
	$nodejavu->insert_cnv($gb,$hdejavu->{$gb},{});
	
}		


# Fermer le fichier
			
		
# freeze la table du dejavu pour le projet
my $file_dejavu_inthisproject = $variationsDir.$projectname."_dejavu.allCNV";
my $file_dejavu = $project->DejaVuCNVFile();
warn $file_dejavu_inthisproject;
warn $file_dejavu;
store($hdejavu, $file_dejavu_inthisproject) or die "Can't store $file_dejavu_inthisproject!\n";
store($hdejavu, $file_dejavu) or die "Can't store $file_dejavu!\n";




exit(0);


##################################
#   methodes
####################################		



sub setComplementaryInfos
{
	my ($hcnv,$patient)  =@_;
	my $project = $patient->project;
	foreach my $global_id  (keys %{$hcnv}){
		my ($t,$c,$start,$end) = split(/_/, $global_id);
			my $objChr = $project->getChromosome($c);
						
						my $tabGenes = $objChr->getGenesByPosition($start,$end);
						
						my $genes_liste="";
	
						foreach my $g (@$tabGenes)
						{
							warn $objChr->name." $start $end" unless $g;
							warn Dumper $tabGenes unless $g;
							$genes_liste .= $g->external_name.",";
						} 


						# calculer le score gene max correspondant a la liste de gene associe au variant
						$hcnv->{$global_id}->{'GENES'} = " ";
						$hcnv->{$global_id}->{'SCORE_GENES'} = 0;
						if ( scalar(@$tabGenes) )  # si le variant recouvre des gènes
						{	
							my @names;
							my $max;
							foreach my $g (sort {$b->score <=> $a->score} @$tabGenes){
										$max = $g->score unless $max;
										push (@names,$g->external_name);
							}
		
							$hcnv->{$global_id}->{'SCORE_GENES'}=$max;
							$hcnv->{$global_id}->{'BEST_GENE'}= $names[0];
							$hcnv->{$global_id}->{'GENES'} =join(" ",@names);
						}
						
						# pour detecter la presence de duplication segmentaire
						$hcnv->{$global_id}->{'DUPSEG'}= getDupSeg($objChr->ucsc_name,$start,$end);
						
						# pour les cytobandes
						my @tb;
						my @band;
						my $hband = $objChr->getCytoband( $start, $end ) if ( $end> $start);

						foreach my $b ( keys %{ $hband->{'name'} } ) {
								push( @tb, $b );
						}
						@band = sort( { $a cmp $b } @tb );
						$hcnv->{$global_id}->{'CYTOBAND'} = join( ",", @band );
						
						# pour l'acces a DGV
						my $url_DGV = getDGV_url( $objChr->ucsc_name, $start, $end );
						$hcnv->{$global_id}->{'DGV'} = $url_DGV;
						# pour le dejavu in this project
						#$hdejavu->{$global_id}->{$patname}->{$caller}++;
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




sub gather_CNV
{
	
	my ($hcnv,$patient) = @_;
	my $hCNV;
	
	warn "START";
	foreach my $caller (keys %{$hcnv})
	{
		print $caller."++";

				foreach my $id  (keys %{$hPat_CNV->{$caller}})
				{
					$hCNV->{$id}->{'id'} = $hPat_CNV->{$caller}->{$id}->{'id'};			
					$hCNV->{$id}->{'TYPE'} = $hPat_CNV->{$caller}->{$id}->{'SVTYPE'};		
					$hCNV->{$id}->{'CHROM'} = $hPat_CNV->{$caller}->{$id}->{'CHROM'};
					$hCNV->{$id}->{'START'} = $hPat_CNV->{$caller}->{$id}->{'START'};
					$hCNV->{$id}->{'END'} = $hPat_CNV->{$caller}->{$id}->{'END'};
					$hCNV->{$id}->{'LEN'} = $hPat_CNV->{$caller}->{$id}->{'SVLEN'};				
					$hCNV->{$id}->{'GOLD_G_FREQ'} = $hPat_CNV->{$caller}->{$id}->{'GOLD_G_FREQ'};
					$hCNV->{$id}->{'GOLD_L_FREQ'} = $hPat_CNV->{$caller}->{$id}->{'GOLD_L_FREQ'};
					$hCNV->{$id}->{'dbVar_status'} = $hPat_CNV->{$caller}->{$id}->{'dbVar_status'};
					$hCNV->{$id}->{'RANKAnnot'} = $hPat_CNV->{$caller}->{$id}->{'RANKAnnot'};
					$hCNV->{$id}->{'DUPSEG'} = $hPat_CNV->{$caller}->{$id}->{'DUPSEG'};
					$hCNV->{$id}->{'CYTOBAND'} = $hPat_CNV->{$caller}->{$id}->{'CYTOBAND'};
					$hCNV->{$id}->{'DGV'} = $hPat_CNV->{$caller}->{$id}->{'DGV'};
					$hCNV->{$id}->{'KARYOTYPE_ID'}=$hPat_CNV->{$caller}->{$id}->{'KARYOTYPE_ID'};
					# pour sauvegarder l'info propre aux differents callers
					$hCNV->{$id}->{'CALLERS'}->{$caller}++;
					$hCNV->{$id}->{'INFOS'}->{$caller}= $hPat_CNV->{$caller}->{$id}->{'INFOS'};
					$hCNV->{$id}->{'ELEMENTARY'}->{$caller} = $hPat_CNV->{$caller}->{$id}->{'ELEMENTARY'};
					$hCNV->{$id}->{'GT'}->{$caller} = $hPat_CNV->{$caller}->{$id}->{'GT'};
					$hCNV->{$id}->{'CN'}->{$caller} = $hPat_CNV->{$caller}->{$id}->{'CN'};
					$hCNV->{$id}->{'RATIO'}->{$caller} = $hPat_CNV->{$caller}->{$id}->{'RATIO'};
					$hCNV->{$id}->{'QUAL'}->{$caller} = $hPat_CNV->{$caller}->{$id}->{'QUAL'};
					
				}
	}
	return gatherSV_byPosition($hCNV,$patient);
	
}


	
# pour regrouper les ids recouvrant les mêmes positions

#sub getDejavuAndTransmission
#{
#		my ($global_id,$patient) = @_;
#
#		my $chr= $project->getChromosome($num);
#					
#		my $nbdejavuProject = 0;
#		my $list_of_other_patient = "";
#		
#		my ($t,$c,$start1,$end1) = split(/_/,$global_id);
#					
#		# on recherche dans le dejavu des  CNV chevauchants
#		
#		
#		my $tab_id = $htree_dejavuProject->{$type}->{$num}->fetch($start1,$end1);
#		
#
#		# puis pour chacun d'eux on regarde ceux qui sont identiques au sens de areSameCNV
#		my $theTransmission;
#		my $transmission=" ";
#		my $maxidM = 0;
#		my $maxidF = 0;
#
#		foreach my $djv_id (@$tab_id)
#		{
#				my ($t,$c,$start2,$end2) = split(/_/,$djv_id);
#				my $identity = int($project->dejavuCNV->getIdentityBetweenCNV($start1,$end1,$start2,$end2));
#				
#				if ( $identity >= $seuilTransmited*0.8)
#				{
#						#dejavu ou transmission ?
#						foreach my $pname (keys %{$hdejavuProject->{$type}->{$num}->{$djv_id}})
#						{		
#								next if ($pname eq $thePatientName);
#								
#								if ( ($pname eq $mothername) || ($pname eq $fathername))
#								{
#									if ($pname eq $mothername) 
#									{
#										 if ($identity >= $seuilTransmited)
#										 {
#										 	$transmission .= "mother " unless ($transmission =~ m/mother/); 		# identity a plus de 90 = transmission
#										 }
#										 else
#										 {
#										 	$transmission .= "maybeM(".$identity."%) " unless ( ($transmission =~ m/mother/) ||  ( ($transmission =~ m/maybeM/) && ($identity>$maxidM)) ); 	# identity entre 70 et 90 = denovo bof
#											$maxidM = $identity if ($identity > $maxidM);
#										 }
#									}
#									
#									if ($pname eq $fathername)
#									{
#										if ($identity >= $seuilTransmited)
#										{
#										 	$transmission .= "father " unless ( $transmission =~ m/father/);
#										}
#										else
#										{
#											$transmission .= "maybeF(".$identity."%) " unless ( ($transmission =~ m/father/) ||  ( ($transmission =~ m/maybeF/) && ($identity>$maxidF)) ); 	# identity entre 70 et 90 = denovo bof
#											$maxidF = $identity if ($identity > $maxidF);
#										}
#									}
#								}
#								else
#								{
#											unless($list_of_other_patient =~ m/$pname/ )
#											{
#												$list_of_other_patient .= $pname.",";
#												$nbdejavuProject++;
#											}
#								}
#						}
#				}
#		}
#		
#		if ($thePatientFamily->isTrio() && $thePatient->isChild() )
#		{
#				 if ($transmission eq " ")
#				 {
#				 	$theTransmission = "strict-denovo";
#				 }
#				 else
#				 {
#				 	$theTransmission = $transmission;
#				 }
#		}
#		else 
#		{
#				$theTransmission = "X";
#		}
#		return ($nbdejavuProject,$list_of_other_patient,$theTransmission);
#}
#


sub getScoreQual {
	my ($hash ) = @_;
	
	
	my $scorequal;
	my $nbcaller= scalar keys %$hash;
	if (exists $hash->{wisecondor}){
		if ( $hash->{wisecondor} >= 40)
				{
					$scorequal++;
				}
	}
	if (exists $hash->{canvas}){
		if ( $hash->{canvas} >= 20)
				{
					$scorequal++;
				}
	}
		if (exists $hash->{manta}){
		if ( $hash->{manta} >= 400)
				{
					$scorequal++;
				}
	}
	
	
	$scorequal /= $nbcaller;
	return $scorequal; 
}
 	