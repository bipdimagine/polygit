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
use lib "$Bin/../../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";

use GBuffer;
use GenBoProject;
use GenBoCache;

use layout;
use export_excel; 
use export_data;

#------------------------------------------------------------
# prepare les données pour le plot des variants He=0 Ho=1
#------------------------------------------------------------

# Instanciation d'un objet buffer pour gérer les connexions aux bases de données  
my $buffer = GBuffer->new();		

# recupere les options  
my $cgi = new CGI;
my $projectname = $cgi->param('project');
my $patient_name = $cgi->param('patient');
my $chr_name = $cgi->param('chr');
my $seuil = $cgi->param('minlength');

die("\n\nNo -project option... Die...\n\n") unless ($projectname);

$chr_name = "X" if ($chr_name == 23);
$chr_name = "Y" if ($chr_name == 24);

#  instanciation d'un objet project_cache et des objets patient et chromosome
my $project = $buffer->newProjectCache( -name => $projectname, -typeFilters=>"");
my $patient = $project->getPatient("$patient_name");
my $chr = $project->getChromosome($chr_name);


# pour savoir si on a ou pas les parents
my $patientFamily = $patient->getFamily();
my $trio = 0;
$trio=1 if ( ($patient->isChild()) &&  ($patientFamily->mother()) );
$trio=2 if ( ($patient->isChild()) &&  ($patientFamily->father()) );
$trio=3 if ( ($patient->isChild()) &&  ($patientFamily->mother()) && ($patientFamily->father()) );

# pour le json final
my $hres;
my @listHashRes;



##############################
# 	Regions Ho du patient   #
##############################
warn "variant Ho du patient debut";

my $vect= $patient->getVectorOrigin($chr);		# tous les variants du patient
$vect &= $chr->getVectorSubstitutions();		# qui sont des substitutions
	
my $list = to_array($vect,$chr->name);	
$project->setListVariants($list);


my $nb=0;
my $RegionHo="0 ";
my $RegionHoEnCours = " ";
my $BornesRegionHoStart=" ";
my $BornesRegionHoStop=" ";
my $RegionHoStart="0";
my $RegionHoStop ="0";


my $HeCourant;
my $HoCourant;

warn "boucle sur les variants du patient";
my $i=0;
while(my $v = $project->nextVariant)
{	
	warn $i++;
	if ($v->isHeterozygote($patient)) 
	{
		$HeCourant= $v->start;
		
 		if ($RegionHoEnCours eq " ") # progression dans une region He
		{
			next; 
		}
		else	# fin d'une région Ho
		{
			
			if ($nb >= $seuil)
			{
				$RegionHoStop = $HeCourant;				# he qui stop la region Ho = debut d'une region he
				$RegionHo .= $RegionHoEnCours.$HeCourant." ";		# on stocke les variants Ho si nb > seuil	$RegionHoEnCours = " ";
				$BornesRegionHoStart .= $RegionHoStart." ";
				$BornesRegionHoStop  .= $RegionHoStop." ";
				
				#$RegionHoStart = "0";
				#$RegionHoStop = "0";
			}
			$RegionHoEnCours = " ";
			$nb=0;
		}
	}
	else
	{
		$HoCourant=$v->start;
		if ($RegionHoEnCours eq " ")  #debut d'un region Ho
		{
			$RegionHoStart = $HoCourant;	
			$nb++;
		}
		else	# progression dans une region Ho
		{
			$nb++;
		}
		$RegionHoEnCours .= $v->start." ";
	}
}
	
$hres->{'PLOT'} = $trio;
$hres->{'PatPOSHo'} = $RegionHo;
$hres->{'PatPOSHoStart'} = $BornesRegionHoStart;
$hres->{'PatPOSHoStop'} = $BornesRegionHoStop;

if( ($trio==1) || ($trio==3) ) # on a la mère
{
	my $mother = $patient->getFamily()->getMother();

	##############################
	# 	Regions Ho de la mere  #
	##############################
	warn "variant Ho de la mere";
	
	$vect= $mother->getVectorOrigin($chr);		# tous les variants de la mere
	$vect &= $chr->getVectorSubstitutions();	# qui sont des substitutions
	
	my $list_mom = to_array($vect,$chr->name);	
	$project->setListVariants($list_mom);
	
	
	 $nb=0;
	 my $RegionHo_mom=" ";
	 my $RegionHoEnCours_mom = " ";
	 my $BornesRegionHoStart_mom=" ";
	 my $BornesRegionHoStop_mom=" ";
	 my $RegionHoStart_mom="0";
	 my $RegionHoStop_mom ="0";


	 my $HeCourant_mom;
	 my $HoCourant_mom;
	 
	warn "boucle sur les variants du mother";
	my $i=0;
	while(my $v = $project->nextVariant)
	{	
		warn $i++;
		if ($v->isHeterozygote($mother)) 
		{
			$HeCourant_mom= $v->start;
		
 			if ($RegionHoEnCours_mom eq " ") # progression dans une region He
			{
				next; 
			}
			else	# fin d'une région Ho
			{
			
				if ($nb >= $seuil)
				{
					$RegionHoStop_mom = $HeCourant_mom;					# he qui stop la region Ho = debut d'une region he
					$RegionHo_mom .= $RegionHoEnCours_mom." ";		#  on stocke les variants Ho si nb > seuil	$RegionHoEnCours = " ";
					$BornesRegionHoStart_mom .= $RegionHoStart_mom." ";
					$BornesRegionHoStop_mom  .= $RegionHoStop_mom." ";
				
					$RegionHoStart_mom = "0";
					$RegionHoStop_mom = "0";
				}
				$RegionHoEnCours_mom = " ";
				$nb=0;
			}
		}
		else
		{
			$HoCourant_mom=$v->start;
			if ($RegionHoEnCours_mom eq " ")  #debut d'un region Ho
			{
				$RegionHoStart_mom = $HoCourant_mom;	
				$nb++;
			}
			else	# progression dans une region Ho
			{
				$nb++;
			}
			$RegionHoEnCours_mom .= $v->start." ";
		}
	}
	
	$hres->{'MotherPOSHo'} = $RegionHo_mom;
	$hres->{'MotherPOSHoStart'} = $BornesRegionHoStart_mom;
	$hres->{'MotherPOSHoStop'} = $BornesRegionHoStop_mom;
}

if( ($trio==2) || ($trio==3) ) # on a le père
{	
	my $father = $patient->getFamily()->getFather();

	##############################
	# 	Regions Ho du pere       #
	##############################
	warn "variant Ho du pere";
	
	$vect= $father->getVectorOrigin($chr);		# tous les variants de la mere
	$vect &= $chr->getVectorSubstitutions();	# qui sont des substitutions
	
	my $list_dad = to_array($vect,$chr->name);	
	$project->setListVariants($list_dad);
	
	$nb=0;
	 my $RegionHo_dad=" ";
	 my $RegionHoEnCours_dad = " ";
	 my $BornesRegionHoStart_dad=" ";
	 my $BornesRegionHoStop_dad=" ";
	 my $RegionHoStart_dad="0";
	 my $RegionHoStop_dad ="0";


 	my $HeCourant_dad;
	my $HoCourant_dad;

	warn "boucle sur les variants du father";
	my $i=0;
	while(my $v = $project->nextVariant)
	{	
		warn $i++;
		if ($v->isHeterozygote($father)) 
		{
			$HeCourant_dad= $v->start;
		
 			if ($RegionHoEnCours_dad eq " ") # progression dans une region He
			{
				next; 
			}
			else	# fin d'une région Ho
			{
			
				if ($nb >= $seuil)
				{
					$RegionHoStop_dad = $HeCourant_dad;					# he qui stop la region Ho = debut d'une region he
					$RegionHo_dad .= $RegionHoEnCours_dad." ";			#  on stocke les variants Ho si nb > seuil	$RegionHoEnCours = " ";
					$BornesRegionHoStart_dad .= $RegionHoStart_dad." ";
					$BornesRegionHoStop_dad  .= $RegionHoStop_dad." ";
				
					$RegionHoStart_dad = "0";
					$RegionHoStop_dad = "0";
				}
				$RegionHoEnCours_dad = " ";
				$nb=0;
			}
		}
		else
		{
			$HoCourant_dad=$v->start;
			if ($RegionHoEnCours_dad eq " ")  #debut d'un region Ho
			{
				$RegionHoStart_dad = $HoCourant_dad;	
				$nb++;
			}
			else	# progression dans une region Ho
			{
				$nb++;
			}
			$RegionHoEnCours_dad .= $v->start." ";
		}
	}
	
	$hres->{'FatherPOSHo'} = $RegionHo_dad;
	$hres->{'FatherPOSHoStart'} = $BornesRegionHoStart_dad;
	$hres->{'FatherPOSHoStop'} = $BornesRegionHoStop_dad;
}
	
push( @listHashRes, { %{ $hres } } );
	
printJson( \@listHashRes );
exit(0);





################
#  methods
################

sub printJson {
	my ($listHash) = @_;
	my $hash;
	my @t;
	
	$hash->{'identifier'} = 'PLOT';
	$hash->{'items'} = \@$listHash;
	
	print $cgi->header('text/json-comment-filtered');
	print encode_json $hash;
	print "\n";
}


sub to_array {
	my ( $v, $name ) = @_;
	my $set  = Set::IntSpan::Fast::XS->new( $v->to_Enum );
	my $iter = $set->iterate_runs();
	my @t;
	my $x=0;
	while ( my ( $from, $to ) = $iter->() ) {
		for my $member ( $from .. $to ) {
			$x++;
			next if (($x%50 != 0) && $project->isGenome);
			if ($name) {
				push( @t, $name . "!" . $member );
			}
			else {
				push( @t, $member );
			}
		}
	}
	return \@t;
}


