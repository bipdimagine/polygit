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

my $compteur;
my $cgi = new CGI;
my $projectname = $cgi->param('project');
my $patientname = $cgi->param('patient');
my $chr         = $cgi->param('chrom');
#my $qual = $cgi->param('qual');
my $dejavu      = $cgi->param('dejavu');
my $listOfGenes = $cgi->param('genes');
my $omim      = $cgi->param('omim');
my $transmis = $cgi->param('transmission');
my $scoreEvent = $cgi->param('score_event');
#

my $type;

my $hfilteredEvent;

my @listHashRes;
my $hbreakpoint;
my $hallPatient;


# creer un objet projet
my $buffer = GBuffer->new();	
my $project = $buffer->newProjectCache( -name => $projectname);

# pour distinguer parmi les autres patients du projet
# ceux qui sont de la même famille
my $thePatient = $project->getPatient($patientname);
my $thePatientFamily = $thePatient->getFamily();
my $mothername;
my $fathername;

if ($thePatient->isChild())
{
		$mothername = $thePatientFamily->mother() if ( $thePatientFamily->mother() );
		$fathername = $thePatientFamily->father() if ( $thePatientFamily->father() );
}
	
					
# pour acceder au fichier allBND
my $path = $project->getSVeqDir();
my $BNDfile = $path.$patientname.".allBND";

# pour acceder au  dejavu global
my $path_djv = $project->DejaVuSVeq_path;
my $TranslocDejavufile = $path_djv."SVeqDejavu.all";

###################################
# pour recuperer les evenements DUP DEL associés 
###################################

my $pathCNV = $project->getCNVDir();
my $hCNV;


my $dejavuCNVProject_file = $pathCNV.$projectname."_dejavu.allCNV";
my $hdejavuCNVProject= retrieve($dejavuCNVProject_file) or die "Can't retrieve datas from " . $dejavuCNVProject_file . " !\n";

foreach my $t (keys %{$hdejavuCNVProject})
{
	foreach my $n (keys %{$hdejavuCNVProject->{$t}})
	{
		foreach my $id (keys %{$hdejavuCNVProject->{$t}->{$n}})
		{
			foreach my $name (keys %{$hdejavuCNVProject->{$t}->{$n}->{$id}})
			{
				$hCNV->{$name}->{$id}=1;
			}
		
		}
	}
}

# pour accelerer acces = intervall tree
my $htree_dejavuCNVProject;
foreach my $type (keys %{$hdejavuCNVProject})
{
		foreach my $num (keys %{$hdejavuCNVProject->{$type}})
		{
				
				# 1)  detecter les SV identiques 
				$htree_dejavuCNVProject->{$type}->{$num}= Set::IntervalTree->new;
				
				# remplir les arbres :  regrouper les SV chevauchants
				foreach my $id  (keys %{$hdejavuCNVProject->{$type}->{$num}})
				{
						my ( $t, $c, $d, $f ) = split( /_/, $id );
						$htree_dejavuCNVProject->{$t}->{$c}->insert($id,$d-5000,$d+5000);
						$htree_dejavuCNVProject->{$t}->{$c}->insert($id,$f-5000,$f+5000);
				}
		}
}	
	
## pour acceder au dejavu global  des CNVs
my $dejavuCNVdir = $project->DejaVuCNV_path();
my $lmdbCNV = GenBoNoSqlLmdb->new(dir=>$dejavuCNVdir,mode=>"r",name=>"dejavu_sv",is_compress=>1);	

###########################################################################################

# on recupère les données
my $hTransLoc = retrieve($BNDfile) or die "Can't retrieve datas from " . $BNDfile . " !\n";
my $hdejavu = retrieve($TranslocDejavufile) or die "Can't retrieve datas from " . $TranslocDejavufile . " !\n";
my $hdejavubis;

my $nbPatientTotal = scalar (keys %{$hallPatient});

# pour accelerer le dejavu des BND
foreach my $djv_event_id (keys %{$hdejavu})
{
		my ($c1,$p1,$c2,$p2) = split(/_/,$djv_event_id);
		my $k = $c1."_".$c2;
		$hdejavubis->{$k}->{$djv_event_id}=$hdejavu->{$djv_event_id};
}

foreach my $event_id (keys %$hTransLoc)
{

	# filtre sur les gènes omim
	if ($omim)
	{
		next unless ( ($hTransLoc->{$event_id}->{"OMIM1"}) ||  ($hTransLoc->{$event_id}->{"OMIM2"}));
	}

	
	
	# pour le dejavu
	my $listpat_itp ="";
	my $listpat_iop ="";
	
	# pour les chromosomes
	my ($chr1, $bp1, $chr2, $bp2) = split(/_/, $event_id);
	$type = "transloc" if $chr1 != $chr2;
	$type = "inv" if $chr1 == $chr2;
	

	unless($chr==0) {
		if ( $chr == 25 && ($type eq "transloc"))  #chromosomes acrocentriques = 13,14,15,21,22
		{
			next unless ((($chr1 == 13) || ($chr1 == 14) || ($chr1 == 15) || ($chr1 == 21) || ($chr1 == 22)) && (($chr2 == 13) || ($chr2 == 14) || ($chr2 == 15) || ($chr2 == 21) || ($chr2 == 22)));
		}
		if ( $chr == 25 && ($type eq "inv"))  #chromosomes acrocentriques = 13,14,15,21,22 
		{
			next unless ( ($chr1 == 13) || ($chr1 == 14) || ($chr1 == 15) || ($chr1 == 21) || ($chr1 == 22) );
		}
		
		if ( $chr != 25 )
		{
			next unless (($chr1 == $chr) || ($chr2 == $chr));
		}
	}
	
	my $c1 = $chr1;
	my $c2 = $chr2;
	
	if ($chr1 == 23)
	{
		$c1 ="X";
	}
	if ($chr1 == 24)
	{
		$c1 ="Y";
	}
	if ($chr2 == 23)
	{
		$c2 ="X";
	}
	if ($chr2 == 24)
	{
		$c2 ="Y";
	}
	
	my $cb1;
	my $cb2;
	
	# cas particulier des inversions
	if ($type eq "inv")
	{
		$hTransLoc->{$event_id}->{"LENGTH"} = abs($bp1-$bp2);
	}
	else
	{
		$hTransLoc->{$event_id}->{"LENGTH"} = "-" ;
	}
	
	$hTransLoc->{$event_id}->{"TRANSLOC"} = $type."##".$c1."##".$hTransLoc->{$event_id}->{"CYTOBAND1"}."##".$c2."##".$hTransLoc->{$event_id}->{"CYTOBAND2"};
	
	$hTransLoc->{$event_id}->{"CYTOBAND1"} = $c1.$hTransLoc->{$event_id}->{"CYTOBAND1"};
	$hTransLoc->{$event_id}->{"CYTOBAND2"} = $c2.$hTransLoc->{$event_id}->{"CYTOBAND2"};
	
	
	
	my $hdjv_project;
	my $hdjv_patient_itp;
	my $hdjv_patient_iop;

	#pour la transmission
	my $transmission=0;
	my $infoTransmission="-";
		
	my $identity1;
	my $identity2;
	
	my $k=$chr1."_".$chr2;
	
	foreach my $djv_event_id (keys %{$hdejavubis->{$k}})
	{
		
		my ($c1,$p1,$c2,$p2) = split(/_/,$djv_event_id);

		next if ( ($p1 < $bp1-50) || ($p1 >$bp1+50) );
		next if ( ($p2 < $bp2-50) || ($p2 >$bp2+50) );
				
		$identity1 = $p1-$bp1;
		$identity2 = $p2-$bp2;
		
		$transmission=0;
		
		foreach my $dejavu_project (keys %{$hdejavu->{$djv_event_id}})
		{
			
			$hdjv_project->{$dejavu_project}++  unless ($dejavu_project eq $projectname);
		
			foreach my $dejavu_patient (keys %{$hdejavu->{$djv_event_id}->{$dejavu_project}})
			{
				if ($dejavu_project eq $projectname)
				{

					#dejavu_itp ou transmission
					if ($dejavu_patient eq $mothername)
					{
								$transmission += 2;
					}
					else
					{
							if ( $dejavu_patient eq $fathername)
							{
								$transmission += 1;
							}
							else
							{
									$hdjv_patient_itp->{$dejavu_patient}++ unless($dejavu_patient eq $patientname);
							}
					}
				}
				else
				{
					$hdjv_patient_iop->{$dejavu_project.":".$dejavu_patient.":bp1=".$identity1."_m bp2=".$identity2."_m"}++; 
				}
			}
		}
	
		# etablir la transmission
		if ($thePatientFamily->isTrio() && $thePatient->isChild() )
		{
			$infoTransmission .= "both  " if (($transmission == 3) && ($infoTransmission !~ m/both/));
			$infoTransmission .= "mother " if (($transmission == 2) && ($infoTransmission !~ m/mother/));
			$infoTransmission .= "father " if (($transmission == 1) && ($infoTransmission !~ m/father/));
				
			if ($transmission == 0)
			{
				$infoTransmission .= "denovo " if ($infoTransmission !~ m/denovo/);
			}
		}
		else
		{
			$infoTransmission = "X";
		}
	}

	# filtre sur la transmission
	
		my $pass = 1;
	
		unless ( ($transmis  eq "all")  || ( $infoTransmission eq "X"))
		{
			$pass = 0;
			if ( ($transmis eq "mother") &&  ($infoTransmission =~ m/mother/ ) &&  !($infoTransmission =~ m/father/ ) )
			{
				$pass = 1;
			}
			
			if ( ($transmis eq "father") &&  ($infoTransmission =~ m/father/ ) &&  !($infoTransmission =~ m/mother/ ) )
			{
				$pass = 1;
			}
			
			if ( ($transmis eq "both") && ( ($infoTransmission =~ m/both/ ) || (($infoTransmission =~ m/father/ ) &&  ($infoTransmission =~ m/mother/ ) )))
			{
				$pass = 1;
			}
			
			if ( ($transmis eq "denovo") &&  ($infoTransmission =~ m/denovo/ ) &&  !($infoTransmission =~ m/father/ ) &&  !($infoTransmission =~ m/mother/ ) &&  !($infoTransmission =~ m/both/ ) )
			{
				$pass = 1;
			}
		}
		next unless ($pass);

	$hTransLoc->{$event_id}->{"TRANSMISSION"} = $infoTransmission;	

	my $nbdejavu_project = scalar keys %{$hdjv_project};
	my $nbPatient_itp = scalar keys %{$hdjv_patient_itp };
	my $nbPatient_iop = scalar keys %{$hdjv_patient_iop };
	my $nbdejavuTotal = $nbPatient_itp + $nbPatient_iop;
	
	$hTransLoc->{$event_id}->{"nbdejavu"} = $nbdejavuTotal;
	
	$listpat_itp = join(",",sort keys %{$hdjv_patient_itp});
	$listpat_iop =  join(",",sort keys %{$hdjv_patient_iop});
	$hTransLoc->{$event_id}->{"dejavu"} = $type."(".$event_id.")+".$nbPatient_itp.";".$listpat_itp.",;".$nbdejavu_project.";".$nbPatient_iop.";".$listpat_iop;
	

	# pour afficher tous les membres de la famille dans IGV si trio
	my $bam_dir     = $thePatient->getProject->getAlignmentUrl($thePatient->alignmentMethod);
	my $bamFiles = $bam_dir.$patientname.".bam";
	my $bamNames = $patientname;
	
	if ($thePatientFamily->isTrio() )
	{
		my $members = $thePatientFamily->getMembers();

		foreach my $m (@$members)
		{
			my $membername = $m->name();
			$bamFiles .= ",".$bam_dir.$membername.".bam" unless ($membername eq $patientname);
			$bamNames .= ",".$membername unless ($membername eq $patientname);
		}
	}
	$hTransLoc->{$event_id}->{'IGV'} = $bamFiles.";".$bamNames.";".$chr1."_".$bp1."_".$chr2."_".$bp2;
	
	# filtres
	unless($dejavu eq "all") {next if ( $hTransLoc->{$event_id}->{"nbdejavu"} > $dejavu)};
	
	
	$hTransLoc->{$event_id}->{"SCORE_EVENT"} = getScoreEvent($event_id);
		
	# filtre sur le score
	next if ($hTransLoc->{$event_id}->{"SCORE_EVENT"} < 10) && $scoreEvent==0;
	next if ($hTransLoc->{$event_id}->{"SCORE_EVENT"} < 5) && $scoreEvent== 1;
	
	############################################
	# pour afficher les CNVs lies a l'evenement
	
	$hTransLoc->{$event_id}->{'CNV'}=";";
	
	if ($thePatientFamily->isTrio() )
	{
		$hTransLoc->{$event_id}->{'CNV_mother'}=";";
		$hTransLoc->{$event_id}->{'CNV_father'}=";";
	}
	else
	{
		$hTransLoc->{$event_id}->{'CNV_mother'}="X";
		$hTransLoc->{$event_id}->{'CNV_father'}="X";
	}
	
	getLinkedCNV($event_id);
	
	
	#############################################################
	# pour aller rechercher les infos PR /PR dans les vcfs manta 

	my $file_in = $thePatient->_getCallingSVFileWithMethodName("manta","variations");
	my $tabix = $buffer->getSoftware("tabix");	
	my $cmd1 = "$tabix $file_in $c1:$bp1-$bp1 | grep $c2 | cut -f 9,10 ";
	my $res1 = `$cmd1`;
			
 	unless ($res1)
 	{
 				$cmd1 = "$tabix $file_in chr$c1:$bp1-$bp1 | grep chr$c2 | cut -f 9,10 ";
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
 	$hTransLoc->{$event_id}->{'PR'} =$ref."/".$alt;
 	
 	my($ref,$alt)=split(/,/,$sr);
 	$hTransLoc->{$event_id}->{'SR'} = $ref."/".$alt;
 
 
  	# calculer le score des gènes compris dans l'evenement
	
		my $genesListe;
		my $objChr  = $project->getChromosome($c1);
		my $deb = $bp1-50;
		my $end = $bp1+50; 
		my $genebp1bis = $c1.":".$deb."_".$end."##";
		my $tabgenesbp1 = $objChr->getGenesByPosition($deb,$end);
	
		if ( scalar(@$tabgenesbp1) >=1 )  # si le variant recouvre des gènes
		{	
			foreach my $g (@$tabgenesbp1)
			{
				$genebp1bis .= $g->external_name.";".$g->raw_score."##";
				$genesListe .= $g->external_name;
			}
		}
		else
		{
			$genebp1bis="-";
		}
		$hTransLoc->{$event_id}->{"GENES1"}=$genebp1bis;

		$objChr  = $project->getChromosome($c2) if ($c1 ne $c2);
		$deb = $bp2-50;
		$end = $bp2+50; 
		my $genebp2bis=$c2.":".$deb."_".$end."##";
		my $tabgenesbp2 = $objChr->getGenesByPosition($deb,$end);
	
		if ( scalar(@$tabgenesbp2) >=1 )  # si le variant recouvre des gènes
		{	
			foreach my $g (@$tabgenesbp2)
			{
				$genebp2bis .= $g->external_name.";".$g->raw_score."##";
				$genesListe .= $g->external_name;
			}
		}
		else
		{
			$genebp2bis="-";
		}
		$hTransLoc->{$event_id}->{"GENES2"}=$genebp2bis;
		
		
		# dans le cas de inversions 
		my $geneInv;
		my $scoremax=0;

		
		if ($type eq "inv")
		{
			my $tabgenesinv = $objChr->getGenesByPosition($bp1,$bp2);
			if ( scalar(@$tabgenesinv) >=1 )  # si le variant recouvre des gènes
			{	
				foreach my $g (@$tabgenesinv)
				{
					my $gscore=$g->raw_score; 
					my $gname= $g->external_name;
					
					$genesListe .= $gname;
					
					if($gscore > $scoremax)
					{
						$scoremax  =$gscore;
						$geneInv = $gname.";".$gscore."##".$geneInv;
					}
					else
					{
						$geneInv .= $gname.";".$gscore."##";
					} 
				}
				$geneInv = $c1.":".$bp1."_".$bp2."##".$geneInv;
				$hTransLoc->{$event_id}->{"ALLGENES"}=$geneInv;
			}
			else
			{
				$hTransLoc->{$event_id}->{"ALLGENES"}="-";
			}
		}
		else
		{
			$hTransLoc->{$event_id}->{"ALLGENES"}="-";
		}
		
	
	
	# dernier filtre sur la liste de gènes

	unless ($listOfGenes eq "all")
	{
		$pass   = 0;
		my @tabgenes = split(/,/,$listOfGenes);
		foreach my $gene (@tabgenes) 
		{
			if ( $genesListe =~ m/$gene/ )
			{
				$pass = 1;
			}
			last if ( $pass == 1 );
		}
		next unless ($pass);
	}
	
	
	# pour le json final
	push( @listHashRes, { %{$hTransLoc->{$event_id}} } );
}	



if ( (scalar(@listHashRes) == 0) )
{
	my $hash;
	$hash->{"id"}= "-";
	$hash->{"TRANSLOC"}= "-";
	$hash->{"IGV"}= "-";
	$hash->{"QUAL"}= "-";
	$hash->{"ALLGENES"}= "-";	
	$hash->{"LENGTH"}= "-";
	$hash->{"POS1"}= "-";
	$hash->{"CYTOBAND1"}= "-";
	$hash->{"GENES1"}= "-";
	$hash->{"OMIM1"}= "-";
	$hash->{"REF/ALT1"}= "-";
	$hash->{"TYPE"}= "-";
	$hash->{"POS2"}= "-";
	$hash->{"CYTOBAND2"}= "-";
	$hash->{"GENES2"}= "-";
	$hash->{"OMIM2"}= "-";
	$hash->{"REF/ALT2"}= "-";
	$hash->{"dejavuBP1"} = "-";
	$hash->{"dejavuBP2"} = "-";
	$hash->{"dejavu"} = "-";
	$hash->{'TRANSMISSION'} = "-";
	$hash->{'CNV'}= "-";
	$hash->{'CNV_mother'}= "-";
	$hash->{'CNV_father'}= "-";
	$hash->{"SCORE_EVENT"} = "-";
	$hash->{"PR"} = "-";
	$hash->{"SR"} = "-";
	push( @listHashRes, { %{$hash} } );
}

#close($fd);
printJson( \@listHashRes );

exit(0);


############
#          #
# methodes #
#          #
############


sub printJson {
	my ($listHash) = @_;
	my $hash;
	$hash->{'identifier'} = 'id';
	$hash->{'label'} = 'id';

	my @t = sort {$b->{SCORE_EVENT} <=> $a->{SCORE_EVENT} or $a->{CHROM1} <=> $b->{CHROM1}or $a->{CHROM2} <=> $b->{CHROM2}} @$listHash;
	
	$hash->{'items'} = \@t;
	print $cgi->header('text/json-comment-filtered');
	print encode_json $hash;
	print "\n";
}

sub getScoreEvent()
 {
	my ($event_id) = @_;

	my $score = 0;
	
	# break point dans un gene omim morbid
	$score += 1  if $hTransLoc->{$event_id}->{"OMIM1"};	
	$score += 1  if $hTransLoc->{$event_id}->{"OMIM2"};
	
	# score sur le dejavu
	my $nbdejavu = $hTransLoc->{$event_id}->{"nbdejavu"};
	if ($nbdejavu <= 10)
	{
		$score += 10 - $nbdejavu;
	}
	
	# score sur la QUAL
	$score += int($hTransLoc->{$event_id}->{"QUAL"} /100)/10;
	
	#pour présenter en premier les CNV denovo
	if ( ($hTransLoc->{$event_id}->{'TRANSMISSION'}  =~ m/denovo/) && !($hTransLoc->{$event_id}->{'TRANSMISSION'} =~ m/mother/ ) && !($hTransLoc->{$event_id}->{'TRANSMISSION'}  =~ m/father/ ) && !($hTransLoc->{$event_id}->{'TRANSMISSION'}  =~ m/both/) )
	{
			$score += 1.5; 
	}
	
	return $score;
 }
 
 sub getLinkedCNV{
	
	my ($event_id) = @_;
	my $out = ";";
	my $djv=0;
	
	my ($chr1, $bp1, $chr2, $bp2) = split(/_/, $event_id);
	#warn $event_id;
	
	$chr1 = "X" if ($chr1 == 23);
	$chr1 = "Y" if ($chr1 == 24);
	$chr2 = "X" if ($chr2 == 23);
	$chr2 = "Y" if ($chr2 == 24);
	
	# on recherche dans le dejavu des  CNV chevauchants avec bp1 et bp2
	foreach my $t (keys %$htree_dejavuCNVProject)
	{
		my $htab_id;
		
		if ( defined ($htree_dejavuCNVProject->{$t}->{$chr1}))
		{
			my $tab_id1 = $htree_dejavuCNVProject->{$t}->{$chr1}->fetch($bp1-5000,$bp1+5000);
			foreach my $cnv_id (@$tab_id1)
			{
				$htab_id->{$cnv_id}=1;
			}
		}
		
		if ( defined ($htree_dejavuCNVProject->{$t}->{$chr2}))
		{
			my $tab_id2 = $htree_dejavuCNVProject->{$t}->{$chr2}->fetch($bp2-5000,$bp2+5000);
			foreach my $cnv_id (@$tab_id2)
			{
				$htab_id->{$cnv_id}=1;
			}
		}
		
		foreach my $cnv_id ( keys %{$htab_id})
		{
			$djv = getDejavuCNV($cnv_id);
			next if ($djv > 20);
			
			$hTransLoc->{$event_id}->{'CNV'} .= $cnv_id.";"  if $hCNV->{$patientname}->{$cnv_id} == 1;
			$hTransLoc->{$event_id}->{'CNV_mother'} .= $cnv_id.";"  if $hCNV->{$mothername}->{$cnv_id} == 1;
			$hTransLoc->{$event_id}->{'CNV_father'} .= $cnv_id.";"  if $hCNV->{$fathername}->{$cnv_id} == 1;
			
		}
	}
}

sub getDejavuCNV
{
	my ($event_id) = @_;
	
	my $nbproject=0;
	my $nbpatient=0;
	my $scorecaller_evt=0;
	my $list_of_other_patient;
	
	my ($t,$c,$start1,$end1) = split(/_/,$event_id);

	my $no = $project->dejavuSV();
	my $hashdv = $no->get_cnv($c,$start1,$end1,$t,$dejavu,90);

	foreach my $project_name (keys %{$hashdv})
	{
			$nbproject++;
			foreach my $pname (keys %{$hashdv->{$project_name}})
			{
					$nbpatient++;
			}
	}

	return $nbpatient;
}
