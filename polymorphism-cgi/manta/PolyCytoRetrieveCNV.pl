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




###################################################################
# Cherche tous les variants structuraux d'un projet pour chaque patient et tous callers confondus.
# Reconstruit les CNV fragmentes (même caller et bornes a moins de 10% de la longeur)
# Construit une table de hash par patient et freeze ces tables  : nom du fihier = patient.allSV
# Pour le dejavu construit egalement une table de hash qui conserve tous les CNV du projet 
# et garde pour chaque CNV l'info du patient et du caller qui l'a détecté.
# Pour chaque CNV la localisation chr_start_end donne accès à la liste des gènes compris dans cet interval
####################################################################

#my $fork = 15;
my $limit;

#my $projectname;
#my $patientname;
#my $print=0;
#
#GetOptions(
#	'project=s' => \$projectname,
#	'patient=s' => \$patientname,
#	'fork=s' => \$fork,
#);


my $cgi = new CGI;
my $projectname = $cgi->param('project');
my $patientname     = $cgi->param('patient');
my $fork = $cgi->param('fork');
my $minlength   = $cgi->param('minlength');
my $select_best = $cgi->param('select_best');
my $chrom        = $cgi->param('chrom');
my $listOfGenes = $cgi->param('genes');
my $dejavu      = $cgi->param('dejavu');
my $transmission    = $cgi->param('transmission');
my $omim      = $cgi->param('omim');
my $print = $cgi->param('print');
my $force= $cgi->param( 'force' ) == 1 ;


# pour récupérer les objets project et patient
my $buffer = GBuffer->new();
my $project = $buffer->newProjectCache( -name => $projectname);

my $dir = $project->getCacheDir(). "/CNV2/";
my $nodejavu = GenBoNoSqlDejaVuCNV->new( dir => $dir, mode => "r" );

#my $patient = $project->getPatients()->[0];
my $patient = $project->getPatient($patientname);
#warn $patient->name;

my $nb;
my $hGroupedCNV;
my @Filtered_listHashRes;
my $global_id;


# pour construire le path d'acces au bam pour IGV

my $bam_dir  = $patient->getProject->getAlignmentUrl($patient->alignmentMethod);
my $bamFiles = $bam_dir.$patientname.".bam";
my $bamNames = $patientname;

my $patientFamily = $patient->getFamily();
my $trio;

if ($patientFamily->isTrio())
{
	my $members = $patientFamily->getMembers();
	$trio = 1;
	foreach my $m (@$members)
	{
		my $membername = $m->name();
		$bamFiles .= ",".$bam_dir.$membername.".bam" unless ($membername eq $patientname);
		$bamNames .= ",".$membername unless ($membername eq $patientname);
	}
}

# pour la gestion des parents si ils existent
my $mothername;
my $fathername;

if ($trio && ($patient->isChild()) )
{
		$mothername = $patientFamily->mother() if ( $patientFamily->mother() );
		$fathername = $patientFamily->father() if ( $patientFamily->father() );
}
			

foreach my $chr (@{$project->getChromosomes}){
	my $x =  $nodejavu->get_all_cnv($chr->name,10,$patient);
	$nb += scalar(@$x);
	foreach my $cnv (@$x){
		my $global_id = $cnv->{"ID"};
		#warn $global_id;	
		
		# cnv version patrick
		#warn Dumper $cnv;
		# pour les colonnes de l'interface
		
		# 1) ce qui ne dépend que du cnv
		
		#les bornes du CNV
		my $gdeb = $cnv->{'START'};
		my $gend = $cnv->{'END'};
		my $chrnum = $chr->name;
 		
		$hGroupedCNV->{$global_id}->{'id'} = $global_id;
		$hGroupedCNV->{$global_id}->{'TYPE'} = $cnv->{'TYPE'};
		$hGroupedCNV->{$global_id}->{'CHROM'} = $cnv->{'CHROM'};
		$hGroupedCNV->{$global_id}->{'CYTOBAND'} = $cnv->{'CYTOBAND'};
		$hGroupedCNV->{$global_id}->{'LOCUS'} =  $gdeb."-".$gend;
		$hGroupedCNV->{$global_id}->{'LEN'}  = $cnv->{'LEN'};
		$hGroupedCNV->{$global_id}->{'DUPSEG'} = $cnv->{'DUPSEG'};
		
		
		#pour les gènes emportés par le CNV
		my @genes = getGenes($chr, $gdeb,$gend);
		$hGroupedCNV->{$global_id}->{'SCORE_GENES'} = $genes[0];;
		$hGroupedCNV->{$global_id}->{'GENES'} = $genes[1];
		
		
		# pour l'accès à gnomAD
		
		$hGroupedCNV->{$global_id}->{'gnomAD'} = "https://gnomad.broadinstitute.org/region/".$chrnum."-".$gdeb."-".$gend."?dataset=gnomad_sv_r2_1";
		
		$hGroupedCNV->{$global_id}->{'OMIM'} = $cnv->{'OMIM'};
		$hGroupedCNV->{$global_id}->{'DGV'} = $cnv->{'DGV'};
		
		$hGroupedCNV->{$global_id}->{'GOLD_G_FREQ'} = $cnv->{'GOLD_G_FREQ'};
		$hGroupedCNV->{$global_id}->{'GOLD_L_FREQ'} = $cnv->{'GOLD_L_FREQ'};
		$hGroupedCNV->{$global_id}->{'dbVar_status'} = $cnv->{'dbVar_status'};
		
		$hGroupedCNV->{$global_id}->{'wisecondor_elementary'} = 0;
		$hGroupedCNV->{$global_id}->{'canvas_elementary'} = 0; 
		$hGroupedCNV->{$global_id}->{'manta_elementary'} = 0;
		
		# pour les infos dépendantes du patient et du caller
		foreach my $id ( keys( %{$cnv->{'PATIENTS'}->{$patient->id}->{'INFOS_ELEMENTARY'}} ) ) 
		{

			$hGroupedCNV->{$global_id}->{'RATIO'} = "-";	   	
			$hGroupedCNV->{$global_id}->{'GT'}=" ";
			$hGroupedCNV->{$global_id}->{'PR'} = "-";
			$hGroupedCNV->{$global_id}->{'SR'} = "-";
			$hGroupedCNV->{$global_id}->{'QUAL'} =  "";
			  
			$hGroupedCNV->{$global_id}->{'SCORECALLER'} = $cnv->{'PATIENTS'}->{$patient->id}->{'SCORE_CALLERS'};
			
			foreach my $caller ( keys( %{$cnv->{'PATIENTS'}->{$patient->id}->{'INFOS_ELEMENTARY'}->{$id}->{"CALLERS"}} ) ) 
			{	
				$hGroupedCNV->{$global_id}->{'GT'} .= $cnv->{'PATIENTS'}->{$patient->id}->{'INFOS_ELEMENTARY'}->{$id}->{"CALLERS"}->{$caller}->{'GT'}." ";
				$hGroupedCNV->{$global_id}->{'QUAL'} .= $caller.":".$cnv->{'PATIENTS'}->{$patient->id}->{'INFOS_ELEMENTARY'}->{$id}->{"CALLERS"}->{$caller}->{'QUALITY'}."/";
				
				if ($caller eq "manta")
				{
					
					if (exists $cnv->{'PATIENTS'}->{$patient->id}->{'INFOS_ELEMENTARY'}->{$id}->{'CALLERS'}->{$caller}->{'INFOS'}->{'PR'})
					{
						my $tabPR = $cnv->{'PATIENTS'}->{$patient->id}->{'INFOS_ELEMENTARY'}->{$id}->{'CALLERS'}->{$caller}->{'INFOS'}->{'PR'};
						$hGroupedCNV->{$global_id}->{'PR'} = @$tabPR[0]."/".@$tabPR[1];
					}
					
					if (exists $cnv->{'PATIENTS'}->{$patient->id}->{'INFOS_ELEMENTARY'}->{$id}->{'CALLERS'}->{$caller}->{'INFOS'}->{'SR'})
					{
						my $tabSR = $cnv->{'PATIENTS'}->{$patient->id}->{'INFOS_ELEMENTARY'}->{$id}->{'CALLERS'}->{$caller}->{'INFOS'}->{'SR'};
						$hGroupedCNV->{$global_id}->{'SR'} = @$tabSR[0]."/".@$tabSR[1];
					}
					
					$hGroupedCNV->{$global_id}->{'manta_elementary'} = $caller.";".$cnv->{'PATIENTS'}->{$patient->id}->{'INFOS_ELEMENTARY'}->{$id}->{"CALLERS"}->{$caller}->{'SCORE'};
				}
				if ($caller eq "wisecondor")
				{
					$hGroupedCNV->{$global_id}->{'RATIO'} = $cnv->{'PATIENTS'}->{$patient->id}->{'INFOS_ELEMENTARY'}->{$id}->{'CALLERS'}->{$caller}->{'RATIO'};
					$hGroupedCNV->{$global_id}->{'wisecondor_elementary'} = $caller.";".$cnv->{'PATIENTS'}->{$patient->id}->{'INFOS_ELEMENTARY'}->{$id}->{"CALLERS"}->{$caller}->{'SCORE'};		
				}
				if ($caller eq "canvas")
				{
					$hGroupedCNV->{$global_id}->{'canvas_elementary'} = $caller.";".$cnv->{'PATIENTS'}->{$patient->id}->{'INFOS_ELEMENTARY'}->{$id}->{"CALLERS"}->{$caller}->{'SCORE'};
				}
			}
		}
				
		# pour la transmission
		$hGroupedCNV->{$global_id}->{'TRANSMISSION'}  = "-";
		if ($trio)
		{
			$hGroupedCNV->{$global_id}->{'TRANSMISSION'}  = getTransmission($cnv,$mothername,$fathername); 	
		}
		
		# pour le dejavu global
		my $nb_projects = $cnv->{'dejavu'}->{'nb_projects'};
		my $nb_patients = $cnv->{'dejavu'}->{'nb_patients'};
		my $list_patients = $cnv->{'dejavu'}->{'string'};
		my $nbcoverage =  $cnv->{'dejavu'}->{'caller_coverage'};
		my $nbdepth =  $cnv->{'dejavu'}->{'caller_depth'};
		my $nbsr =  $cnv->{'dejavu'}->{'caller_sr'};
		
		my $flag = 0;
		$flag = 1 if ($nbcoverage > 0);
		
		
		$hGroupedCNV->{$global_id}->{'DEJAVU_G'} = $global_id."+".getDejaVu_project($cnv).";".$nb_projects.";".$nb_patients.";".$list_patients.";".$nbcoverage.";".$nbdepth.";".$nbsr.";".$flag;
		$hGroupedCNV->{$global_id}->{'PLOT'} =  $cnv->{'TYPE'}.";".$cnv->{'CHROM'}.":".$cnv->{'START'}."-".$cnv->{'END'}.";".$hGroupedCNV->{$global_id}->{'TRANSMISSION'};	 
		
		
		$hGroupedCNV->{$global_id}->{"IGV"} = $bamNames.";".$bamFiles.";".$hGroupedCNV->{$global_id}->{"id"};
		$hGroupedCNV->{$global_id}->{'BPManta'}="-";	      # manquant
		$hGroupedCNV->{$global_id}->{'BPManta_mother'}="-";	   # manquant
   		$hGroupedCNV->{$global_id}->{'BPManta_father'}="-";	   # manquant
		
		# scoreCNV 
		$hGroupedCNV->{$global_id}->{'SCORECNV'}= getScoreCNV($hGroupedCNV->{$global_id}->{'SCORECALLER'},$hGroupedCNV->{$global_id}->{'TRANSMISSION'},$hGroupedCNV->{$global_id}->{'DUPSEG'});
				
		# utile mais non affiché	
		$hGroupedCNV->{$global_id}->{'START'} = $cnv->{'START'};
		$hGroupedCNV->{$global_id}->{'END'} =  $cnv->{'END'};
		$hGroupedCNV->{$global_id}->{'DEJAVU_P'} = getDejaVu_project($cnv); 	
		#warn Dumper $hGroupedCNV->{$global_id};
		#die;
		
	}
}

# pour l'affichage
if (scalar(keys(%{ $hGroupedCNV})) == 0) 
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
				$hash->{'wisecondor_elementary'} = "-";
				$hash->{'canvas_elementary'} = "-";
				$hash->{'manta_elementary'} = "-";
				$hash->{'BPManta'} = "-";
				$hash->{'BPManta_mother'} = "-";
				$hash->{'BPManta_father'} = "-";
				$hash->{'TRANSMISSION'} = "-";
				$hash->{'DEJAVU_G'} = "-";
				$hash->{'SCORE_GENES'} = "-";
				$hash->{'RANKAnnot'} = "-";
				$hash->{'OMIM'} = "-";
				$hash->{'GENES'} = "-";
				$hash->{'DGV'} = "-";
				$hash->{'GOLD_G_FREQ'} = "-";
				$hash->{'GOLD_L_FREQ'} = "-";
				$hash->{'dbVar_status'} = "-";
			 	 push(@Filtered_listHashRes, $hash); 
}
else
{
	foreach my $global_id ( keys %{$hGroupedCNV} ) 
	{
			push( @Filtered_listHashRes, { %{ $hGroupedCNV->{$global_id} } } );
	}
}

printJson( \@Filtered_listHashRes );
#warn Dumper $hGroupedCNV;			
exit(0);


################################################@

sub getDejaVu_project
{
	my ($cnvHash) = @_;
	my $nbdjv_inproject = 0;
	my $listpatients_djv_inproject = "";
	foreach my $patient_id ( keys( %{$cnvHash->{'PATIENTS'}} ) ) 
	{
		 $nbdjv_inproject++;
		 $listpatients_djv_inproject .= $cnvHash->{'PATIENTS'}->{$patient_id}->{'NAME'}.",";
	}
	
	my $result = $nbdjv_inproject.";".$listpatients_djv_inproject;
	return $result;
}

sub getTransmission
{
	my ($cnvHash,$mothername,$fathername) = @_;
	my $transmission = " ";
	my @tab = ("denovo","father","mother","both");
	
	
	foreach my $patient_id ( keys( %{$cnvHash->{'PATIENTS'}} ) ) 
	{
		 if ($cnvHash->{'PATIENTS'}->{$patient_id}->{'NAME'} eq $mothername)
		 {
		 	$transmission .= "mother"." ";
		 }
		 if ($cnvHash->{'PATIENTS'}->{$patient_id}->{'NAME'} eq $fathername)
		 {
		 		$transmission .= "father"." ";
		 }
	}
	return $transmission;
}

sub getGenes
{
	my ($chr,$gdeb,$gend) = @_;
	my $tabGenes = $chr->getGenesByPosition($gdeb,$gend);
	
	my @tab;
	
	# au final on cherche les genes et ne calcule le score genes que pour les CNV selectionnés
		my $genes_liste = "";
		my @genes_names;
		my $gname;
		my $gscore = 0;
		my $scoremax = 0;
	
		if ( scalar(@$tabGenes) >=1 )  # si le variant recouvre des gènes
		{
			foreach my $g (@$tabGenes)
			{
				$gscore=$g->score; 
				$gname= $g->external_name;
							 
				if($gscore > $scoremax)
				{
					$scoremax = $gscore;
					$genes_liste = $gname.";".$gscore."##".$genes_liste; # placé en tête de liste
				}
				else
				{
					$genes_liste .= $gname.";".$gscore."##";
				} 
			}
			$tab[0] = $scoremax;
			$tab[1] = $chr->name.":".$gdeb."_".$gend."##".$genes_liste;
		}
		else
		{
			$tab[0] = 0;
			$tab[1] = "-";
		}
		return @tab;
}

sub getScoreCNV
{
	my ($score_caller,$transmission,$dupseg) = @_;
	my $score;
	
		$score = $score_caller;
		
		#pour présenter en premier les CNV denovo
		if ( ($transmission  =~ m/denovo/) && !($transmission =~ m/mother/ ) && !($transmission =~ m/father/ ) && !($transmission =~ m/both/) )
		{
			$score += 0.6; 
		}
		
		# pour rétrograder le CNV situes au niveau de segmental duplication (chevauchement de plus de 40%)
		if ($dupseg > 40)
		{
			my $sd_score= int($dupseg/40); 
			$score -= $sd_score;		#	(de -1 à -2,5 )
		}
	
		# pour remonter  ceux qui presentent un BP dans une region de 1kb autour des positions debut ou fin
		#$score += 1 if ($hGroupedCNV->{$gid}->{'BPManta'} ne ";");		
		#$score += 0.5 if ( ($hGroupedCNV->{$gid}->{'BPManta_mother'} ne ";") && ($hGroupedCNV->{$gid}->{'BPManta_mother'} ne "X") );	
		#$score += 0.5 if ( ($hGroupedCNV->{$gid}->{'BPManta_father'} ne ";")   && ($hGroupedCNV->{$gid}->{'BPManta_father'} ne "X")   );	
	
		return $score;

}

sub printJson {
	my ($listHash) = @_;
	my $hash;
	my @t;

	#@t = sort { $a->{SVTYPE} cmp $b->{SVTYPE} } @$listHash;
	@t = sort {$b->{SCORECNV} <=> $a->{SCORECNV}or $b->{SCORE_GENES} <=> $a->{SCORE_GENES}or $a->{SVTYPE} cmp $b->{SVTYPE}} @$listHash;
			
	my $nb=1;
	foreach my $h (@t)
	{
			$h->{"nb"}=$nb++;
	}

	$hash->{'identifier'} = 'id';
	#$hash->{'label'}      = 'id';
	$hash->{'items'}      = \@t;
	
	
	print $cgi->header('text/json-comment-filtered');
	print encode_json $hash;
	print "\n";
	
	
}





