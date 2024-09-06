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
my $patient = $thePatient;
my $thePatientFamily = $thePatient->getFamily();

my $file_in = $thePatient->getSRFile;
my $tabix = Bio::DB::HTS::Tabix->new( filename => $file_in );

my $mothername;
my $fathername;
my $mother_id;
my $father_id;
if ($thePatient->isChild())
{
		
		$mother_id = $thePatientFamily->getMother()->id if ( $thePatientFamily->mother() );
		$father_id = $thePatientFamily->getFather()->id if ( $thePatientFamily->father() );
}
	
					
# pour acceder au fichier allBND

# pour acceder au  dejavu global

###################################
# pour recuperer les evenements DUP DEL associés 
###################################

#my $pathCNV = $project->getCNVDir();
my $hCNV;


#my $dejavuCNVProject_file = $pathCNV.$projectname."_dejavu.allCNV";
#my $hdejavuCNVProject= retrieve($dejavuCNVProject_file) or die "Can't retrieve datas from " . $dejavuCNVProject_file . " !\n";
#
#foreach my $t (keys %{$hdejavuCNVProject})
#{
#	foreach my $n (keys %{$hdejavuCNVProject->{$t}})
#	{
#		foreach my $id (keys %{$hdejavuCNVProject->{$t}->{$n}})
#		{
#			foreach my $name (keys %{$hdejavuCNVProject->{$t}->{$n}->{$id}})
#			{
#				$hCNV->{$name}->{$id}=1;
#			}
#		
#		}
#	}
#}

# pour accelerer acces = intervall tree
#my $htree_dejavuCNVProject;
#foreach my $type (keys %{$hdejavuCNVProject})
#{
#		foreach my $num (keys %{$hdejavuCNVProject->{$type}})
#		{
#				
#				# 1)  detecter les SV identiques 
#				$htree_dejavuCNVProject->{$type}->{$num}= Set::IntervalTree->new;
#				
#				# remplir les arbres :  regrouper les SV chevauchants
#				foreach my $id  (keys %{$hdejavuCNVProject->{$type}->{$num}})
#				{
#						my ( $t, $c, $d, $f ) = split( /_/, $id );
#						$htree_dejavuCNVProject->{$t}->{$c}->insert($id,$d-5000,$d+5000);
#						$htree_dejavuCNVProject->{$t}->{$c}->insert($id,$f-5000,$f+5000);
#				}
#		}
#}	
#	
## pour acceder au dejavu global  des CNVs
#my $dejavuCNVdir = $project->DejaVuCNV_path();

###########################################################################################

# on recupère les données
my $dir = $project->getCacheDir(). "/SV/";
my $nodejavu = GenBoNoSqlDejaVuSV->new( dir => $dir, mode => "r" );
my $hTransLoc = $nodejavu->get_all_sv(10,$patient);


my $hdejavubis;
my $no_sv = $project->dejavuSV();
#my $no1 = $project->dejavuCNV();
my $nbPatientTotal = scalar (keys %{$hallPatient});

# pour accelerer le dejavu des BND


foreach my $sv (@$hTransLoc) {

	# filtre sur les gènes omim
	if ($omim)
	{
		next unless ( ($sv->{"OMIM1"}) ||  ($sv->{"OMIM2"}));
	}
	
	# pour le dejavu
	my $listpat_itp ="";
	my $listpat_iop ="";
	
	# pour les chromosomes
	$type = "transloc";
	$type = "inv" if $sv->{TYPE} eq "INV";
	
#
#	unless($chr==0) {
#		if ( $chr == 25 && ($type eq "transloc"))  #chromosomes acrocentriques = 13,14,15,21,22
#		{
#			next unless ((($chr1 == 13) || ($chr1 == 14) || ($chr1 == 15) || ($chr1 == 21) || ($chr1 == 22)) && (($chr2 == 13) || ($chr2 == 14) || ($chr2 == 15) || ($chr2 == 21) || ($chr2 == 22)));
#		}
#		if ( $chr == 25 && ($type eq "inv"))  #chromosomes acrocentriques = 13,14,15,21,22 
#		{
#			next unless ( ($chr1 == 13) || ($chr1 == 14) || ($chr1 == 15) || ($chr1 == 21) || ($chr1 == 22) );
#		}
#		
#		if ( $chr != 25 )
#		{
#			next unless (($chr1 == $chr) || ($chr2 == $chr));
#		}
#	}
	
	
	
	
	
	
	# cas particulier des inversions

	
	 #20|9852693|9852798|21|20490513|20490613
	#20|54108694|54108694|21|32583601|32583607
	$sv->{"TRANSLOC"} = $sv->{TYPE}."##".$sv->{CHROM1}."##".$sv->{"CYTOBAND1"}."##".$sv->{CHROM2}."##".$sv->{"CYTOBAND2"};
	$sv->{"CYTOBAND1"} = $sv->{CHROM1}.$sv->{"CYTOBAND1"};
	$sv->{"CYTOBAND2"} = $sv->{CHROM2}.$sv->{"CYTOBAND2"};
	
	
	
	my $hdjv_project;
	my $hdjv_patient_itp;
	my $hdjv_patient_iop;

	#pour la transmission
	my $transmission=0;

		
	my $identity1;
	my $identity2;
	
	if (exists $sv->{PATIENTS}->{$mother_id} && exists $sv->{PATIENTS}->{$father_id}){
		$sv->{PATIENTS}->{$patient->id}->{transmission} =  "both";
	}
	elsif (exists $sv->{PATIENTS}->{$mother_id}) {
			$sv->{PATIENTS}->{$patient->id}->{transmission} =  "mother";
	}
	elsif (exists $sv->{PATIENTS}->{$father_id}) {
			$sv->{PATIENTS}->{$patient->id}->{transmission} =  "father";
	}
	elsif ($father_id > 0 && $mother_id > 0 ) {
			$sv->{PATIENTS}->{$patient->id}->{transmission} =  "denovo";
	}
	else {
			$sv->{PATIENTS}->{$patient->id}->{transmission} = "X";
		}
	


	# filtre sur la transmission
	
		my $pass = 1;
	
#		unless ( ($transmis  eq "all")  || ( $sv->{PATIENTS}->{$patient->id}->{transmission} eq "X"))
#		{
#			$pass = 0;
#			if ( ($transmis eq "mother") &&  ($sv->{PATIENTS}->{$patient->id}->{transmission} =~ m/mother/ ) &&  !($sv->{PATIENTS}->{$patient->id}->{transmission} =~ m/father/ ) )
#			{
#				$pass = 1;
#			}
#			
#			if ( ($transmis eq "father") &&  ($sv->{PATIENTS}->{$patient->id}->{transmission} =~ m/father/ ) &&  !($sv->{PATIENTS}->{$patient->id}->{transmission} =~ m/mother/ ) )
#			{
#				$pass = 1;
#			}
#			
#			if ( ($transmis eq "both") && ( ($sv->{PATIENTS}->{$patient->id}->{transmission} =~ m/both/ ) || (($sv->{PATIENTS}->{$patient->id}->{transmission} =~ m/father/ ) &&  ($sv->{PATIENTS}->{$patient->id}->{transmission} =~ m/mother/ ) )))
#			{
#				$pass = 1;
#			}
#			
#			if ( ($transmis eq "denovo") &&  ($sv->{PATIENTS}->{$patient->id}->{transmission} =~ m/denovo/ ) &&  !($sv->{PATIENTS}->{$patient->id}->{transmission} =~ m/father/ ) &&  !($sv->{PATIENTS}->{$patient->id}->{transmission} =~ m/mother/ ) &&  !($sv->{PATIENTS}->{$patient->id}->{transmission} =~ m/both/ ) )
#			{
#				$pass = 1;
#			}
#		}


	my $nbdejavu_project = $sv->{DEJAVU}->{PROJECTS};
	my $nbPatient_itp =  scalar(keys%{$sv->{PATIENTS}}) - 1 ;
	
	my $nbPatient_iop = $sv->{DEJAVU}->{PATIENTS};
	my $nbdejavuTotal = $nbPatient_itp + $nbPatient_iop;
	
	$sv->{"nbdejavu"} = $nbdejavuTotal+"0";
	
	$listpat_itp = join(",",sort map {$sv->{PATIENTS}->{$_}->{PATIENT_NAME}}keys %{$sv->{PATIENTS}});    ### !!!!  join(",",sort keys %{$hdjv_patient_itp});
	$listpat_iop =   join(",",sort @{$sv->{DEJAVU}->{ARRAY}}) if exists $sv->{DEJAVU}->{ARRAY};
	
	$sv->{"DEJAVU"}->{STRING} = $type."(".$sv->{ID}.")+".$nbPatient_itp.";".$listpat_itp.",;".$nbdejavu_project.";".$nbPatient_iop.";".$listpat_iop;
	


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
	$sv->{'IGV'} = $bamFiles.";".$bamNames."et".$sv->{CHROM1}."_".$sv->{POS1}."_".$sv->{CHROM2}."_".$sv->{POS2};
	
	# filtres
	
	$sv->{SCORE_EVENT} = getScoreEvent($sv,$patient);
	die($sv->{SCORE_EVENT}) if 	$sv->{SCORE_EVENT} eq "NaN";
	# filtre sur le score

	
	#next if ($sv->{"SCORE_EVENT"} < 10) && $scoreEvent==0;
	#next if ($sv->{"SCORE_EVENT"} < 5) && $scoreEvent== 1;
	
	############################################
	# pour afficher les CNVs lies a l'evenement
	
	$sv->{'CNV'}=";";
	
	if ($thePatientFamily->isTrio() )
	{
		$sv->{'CNV_mother'}=";";
		$sv->{'CNV_father'}=";";
	}
	else
	{
		$sv->{'CNV_mother'}="X";
		$sv->{'CNV_father'}="X";
	}
	#getLinkedCNV($sv->{ID});
	#############################################################
	# pour aller rechercher les infos PR /PR dans les vcfs manta 

	
	$sv->{'PR'} = join("/",@{$sv->{PATIENTS}->{$patient->id}->{INFOS}->{PR}}) if exists $sv->{PATIENTS}->{$patient->id}->{INFOS}->{PR};
	$sv->{'SR'} = join("/",@{$sv->{PATIENTS}->{$patient->id}->{INFOS}->{SR}}) if exists $sv->{PATIENTS}->{$patient->id}->{INFOS}->{SR};
	$sv->{'PR'} = join("/",@{$sv->{PATIENTS}->{$patient->id}->{INFOS}->{AD}}) if exists $sv->{PATIENTS}->{$patient->id}->{INFOS}->{AD};
 
 
  	# calculer le score des gènes compris dans l'evenement
	
		my $genesListe;
		my $objChr  = $project->getChromosome($sv->{CHROM1});
		my $deb = $sv->{POS1}-50;
		my $end = $sv->{POS1}+50; 
		
		my $genebp1bis = $sv->{CHROM1}.":".$deb."_".$end."##";
		

		
		if ($sv->{GENES1}){
			foreach my $g (keys %{$sv->{GENES1}}){
				$genesListe .= $g;
				$genebp1bis .= $g.";".$sv->{GENES1}->{$g}."##";
			}
		}
		else
		{
			$genebp1bis="-";
		}
		$sv->{"GENES1"}=$genebp1bis;
		
		$deb = $sv->{POS2}-50;
		 $end = $sv->{POS2}+50; 
		my $genebp2bis=$sv->{CHROM2}.":".$deb."_".$end."##";
		if ($sv->{GENES2}){
			foreach my $g (keys %{$sv->{GENES2}}){
				$genesListe .= $g;
				$genebp2bis .= $g.";".$sv->{GENES2}->{$g}."##";
			}
		}
		else
		{
			$genebp2bis="-";
		}
		
		$sv->{"GENES2"}=$genebp2bis;
		
		
		# dans le cas de inversions 
		my $geneInv;
		my $scoremax=0;
		$sv->{"ALLGENES"}="";
		
		get_local_CNV($sv,$sv->{CHROM1},$sv->{POS1},"DEL",$patient->id,$mother_id,$father_id);
		get_local_CNV($sv,$sv->{CHROM1},$sv->{POS1},"DUP",$patient->id,$mother_id,$father_id);
#		if ($type eq "inv")
#		{
#			my $tabgenesinv = $objChr->getGenesByPosition($bp1,$bp2);
#			if ( scalar(@$tabgenesinv) >=1 )  # si le variant recouvre des gènes
#			{	
#				foreach my $g (@$tabgenesinv)
#				{
#					my $gscore= $g->score; 
#					my $gname= $g->external_name;
#					
#					$genesListe .= $gname;
#					
#					if($gscore > $scoremax)
#					{
#						$scoremax  =$gscore;
#						$geneInv = $gname.";".$gscore."##".$geneInv;
#					}
#					else
#					{
#						$geneInv .= $gname.";".$gscore."##";
#					} 
#				}
#				$geneInv = $c1.":".$bp1."_".$bp2."##".$geneInv;
#				$hTransLoc->{$event_id}->{"ALLGENES"}=$geneInv;
#			}
#			else
#			{
#				$hTransLoc->{$event_id}->{"ALLGENES"}="-";
#			}
#		}
#		else
#		{
#			$hTransLoc->{$event_id}->{"ALLGENES"}="-";
#		}
		
	
	
	# dernier filtre sur la liste de gènes
#
#	unless ($listOfGenes eq "all")
#	{
#		$pass   = 0;
#		my @tabgenes = split(/,/,$listOfGenes);
#		foreach my $gene (@tabgenes) 
#		{
#			if ( $genesListe =~ m/$gene/ )
#			{
#				$pass = 1;
#			}
#			last if ( $pass == 1 );
#		}
#		next unless ($pass);
#	}
#	
	# pour le json final
	my $hash;
	$hash->{"id"}= $sv->{ID};
	$hash->{"TRANSLOC"}= $sv->{TRANSLOC};
	$hash->{"IGV"}= $sv->{IGV};
	$hash->{"QUAL"}=$sv->{QUAL};
	$hash->{"ALLGENES"}= "-";	
	$hash->{"LENGTH"}=  $sv->{LENGTH};
	$hash->{"POS1"}= $sv->{POS1};
	$hash->{"CYTOBAND1"}=  $sv->{CYTOBAND1};
	$hash->{"GENES1"}= $sv->{GENES1};
	$hash->{"OMIM1"}= $sv->{OMIM1};
	$hash->{"REF/ALT1"}= "-";
	$hash->{"TYPE"}= $sv->{TYPE};
	$hash->{"POS2"}= $sv->{POS2};
	$hash->{"CYTOBAND2"}=  $sv->{CYTOBAND2};
	$hash->{"GENES2"}= $sv->{GENES2};
	$hash->{"OMIM2"}= $sv->{OMIM2};
	$hash->{"REF/ALT2"}= "-";
	$hash->{"dejavuBP1"} = "-";
	$hash->{"dejavuBP2"} = "-";
	$hash->{"dejavu"} = $sv->{DEJAVU}->{STRING};
	$hash->{'TRANSMISSION'} = $sv->{PATIENTS}->{$patient->id}->{transmission};
	$hash->{'CNV'}= $sv->{CNV};
	$hash->{'CNV_mother'}= $sv->{CNV_mother};
	$hash->{'CNV_father'}= $sv->{CNV_father};
	
	$hash->{"SCORE_EVENT"} = $sv->{SCORE_EVENT};
	$hash->{"PR"} = $sv->{PR};
	$hash->{"SR"} =  $sv->{SR};
	push( @listHashRes, $hash);
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
	my ($sv,$patient) = @_;

	my $score = 0;
	
	# break point dans un gene omim morbid
	$score += 1  if $sv->{"OMIM1"};	
	$score += 1  if $sv->{"OMIM2"};
	# score sur le dejavu
	my $nbdejavu = $sv->{"nbdejavu"};
	
	if ($nbdejavu <= 10)
	{
		$score += 10 - $nbdejavu;
	}
	# score sur la QUAL
	$score += int($sv->{"QUAL"} /100)/10;
	#pour présenter en premier les CNV denovo
	if ( ($sv->{PATIENTS}->{$patient->id}->{'TRANSMISSION'}  =~ m/denovo/) && !($sv->{PATIENTS}->{$patient->id}->{'TRANSMISSION'} =~ m/mother/ ) && !($sv->{PATIENTS}->{$patient->id}->{'TRANSMISSION'}  =~ m/father/ ) && !($sv->{PATIENTS}->{$patient->id}  =~ m/both/) )
	{
			$score += 1.5; 
	}
	return $score;
 }
 
# sub getLinkedCNV{
#	
#	my ($event_id) = @_;
#	my $out = ";";
#	my $djv=0;
#	
#	my ($chr1, $bp1, $chr2, $bp2) = split(/_/, $event_id);
#	
#	$chr1 = "X" if ($chr1 == 23);
#	$chr1 = "Y" if ($chr1 == 24);
#	$chr2 = "X" if ($chr2 == 23);
#	$chr2 = "Y" if ($chr2 == 24);
#	
#	# on recherche dans le dejavu des  CNV chevauchants avec bp1 et bp2
#	foreach my $t (keys %$htree_dejavuCNVProject)
#	{
#		my $htab_id;
#		
#		if ( defined ($htree_dejavuCNVProject->{$t}->{$chr1}))
#		{
#			my $tab_id1 = $htree_dejavuCNVProject->{$t}->{$chr1}->fetch($bp1-5000,$bp1+5000);
#			foreach my $cnv_id (@$tab_id1)
#			{
#				$htab_id->{$cnv_id}=1;
#			}
#		}
#		
#		if ( defined ($htree_dejavuCNVProject->{$t}->{$chr2}))
#		{
#			my $tab_id2 = $htree_dejavuCNVProject->{$t}->{$chr2}->fetch($bp2-5000,$bp2+5000);
#			foreach my $cnv_id (@$tab_id2)
#			{
#				$htab_id->{$cnv_id}=1;
#			}
#		}
#		
#		foreach my $cnv_id ( keys %{$htab_id})
#		{
#			$djv = getDejavuCNV($cnv_id);
#			next if ($djv > 20);
#			
#			$hTransLoc->{$event_id}->{'CNV'} .= $cnv_id.";"  if $hCNV->{$patientname}->{$cnv_id} == 1;
#			$hTransLoc->{$event_id}->{'CNV_mother'} .= $cnv_id.";"  if $hCNV->{$mothername}->{$cnv_id} == 1;
#			$hTransLoc->{$event_id}->{'CNV_father'} .= $cnv_id.";"  if $hCNV->{$fathername}->{$cnv_id} == 1;
#			
#		}
#	}
#}

sub get_local_CNV {
my ($sv,$chr_name,$start,$type,$pid,$mid,$fid) = @_;
my $dir = $project->getCacheDir(). "/CNV/";
my $nodejavu = GenBoNoSqlDejaVuCNV->new( dir => $dir, mode => "r" );

my $chr = $patient->project->getChromosome($chr_name);
my $l = abs ($start- $chr->length);
my $a = $start;
my $b =  $chr->length;
my $del1 = $nodejavu->get_cnv_project($type,$chr_name,$a,$b,30);
$a = 1;
$b =  $start;

my $del2 = $nodejavu->get_cnv_project($type,$chr_name,$a,$b,30);

my @objs;
push(@objs,@$del1);
push(@objs,@$del2);

return  unless @objs;

		foreach my $cnv (@objs)
		{
			next if $cnv->{dejavu}->{nb_patients} > 30;
			$sv->{'CNV'}.=$cnv->{ID}.";" if exists $cnv->{PATIENTS}->{$pid};
			$sv->{'CNV_mother'}.=$cnv->{ID}.";" if exists $cnv->{PATIENTS}->{$mid};
			$sv->{'CNV_father'}.=$cnv->{ID}.";" if exists $cnv->{PATIENTS}->{$fid};
		}
		
return $objs[-1];
}


sub getDejavuCNV
{
	my ($event_id) = @_;
	my $nbproject=0;
	my $nbpatient=0;
	my $scorecaller_evt=0;
	my $list_of_other_patient;
	
	my ($t,$c,$start1,$end1) = split(/_/,$event_id);

	my $no = $project->dejavuCNV();
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
