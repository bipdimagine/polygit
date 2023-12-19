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

#----------------------------------------------------------------------------------------------------
# 
#-----------------------------------------------------------------------------------------------------

# Instanciation d'un objet buffer pour gérer les connexions aux bases de données  
my $buffer = GBuffer->new();		

# recupere les options  
my $cgi = new CGI;
my $projectname = $cgi->param('project');
my $patient_name = $cgi->param('patient');
my $chr_name = $cgi->param('chr');


die("\n\nNo -project option... Die...\n\n") unless ($projectname);

$chr_name = "X" if ($chr_name == 23);
$chr_name = "Y" if ($chr_name == 24);



#  instanciation d'un objet project_cache et des objets patient et chromosome
my $project = $buffer->newProjectCache( -name => $projectname, -typeFilters=>"");
my $patient = $project->getPatient($patient_name);
my $chr = $project->getChromosome($chr_name);


# pour savoir si on a ou pas les parents
my $patientFamily = $patient->getFamily();
my $trio =0;
$trio=1 if ( ($patient->isChild()) &&  ($patientFamily->mother()) && ($patientFamily->father()) );


# pour le json final
my $hres;
my @listHashRes;


$hres->{'PLOT'} = $trio;


if($trio)
{
	# acces aux fichiers des parents
	my $mother = $patient->getFamily()->getMother();
	my $father = $patient->getFamily()->getFather();
	my $mothername = $mother->name();
	my $fathername = $father->name();
	
	
	#	########################################
	#											#
	#     calcul de la fréquence allélique		#
	#											#
	#	########################################

	my $v1;
	my $v2;


	##########	Version on garde tout ce qui est unique à un des deux parents et retrouvé chez l'enfant #################

	$v1 = $mother->getVectorOrigin($chr);			# tous les variants de la mère
	$v1 &= $chr->getVectorSubstitutions();			# qui sont des substitutions
 	$v1 -= $father->getVectorOrigin($chr);			# qui ne sont pas chez le père
 	
 	# dont le ratio est supérieur à 20% chez la mère
	my $vector_ratio_name_mother = $mother->name . "_ratio_20";
	my $vquality_mother = $chr->getVectorScore($vector_ratio_name_mother);
 	$v1 &= $vquality_mother;
 	
	$v1 &= $patient->getVectorOriginHe($chr);			# qui sont chez l'enfant

	##########################################################################################################

	$v2 = $father->getVectorOrigin($chr);			# tous les variants du père
	$v2 &= $chr->getVectorSubstitutions();			# qui sont des substitutions
 	$v2 -= $mother->getVectorOrigin($chr);			# qui ne sont pas chez la mère
 	
 	# dont le ratio est supérieur à 20% chez le père
	my $vector_ratio_name_father = $father->name . "_ratio_20";
	my $vquality_father = $chr->getVectorScore($vector_ratio_name_father);
 	$v2 &= $vquality_father;
	$v2 &= $patient->getVectorOriginHe($chr);			# qui sont chez l'enfant


	my $v3 = $v1+$v2;
	my $list3 = to_array($v3,$chr->name);	
	$project->setListVariants($list3);

	my @array;
	my $x;

	 #boucle sur les variants 
 
	while(my $v = $project->nextVariant){
		$x++;		
		my $p1 = $v->getPourcentAllele($mother);
		my $p2 = $v->getPourcentAllele($father);
		my $value = $v->getPourcentAllele($patient);
		my $vdepth = $v->getDepth($patient);   # à comparer avec la profondeur moyenne du patient 
		
	
		next if $value eq "-";
		#next if $vdepth < 10;
		
		if ($p1 ne "-"){
			my $res = $v->start.",".$value.",null";
			push(@array,$res);
		}

		if ($p2 ne "-"){
			my $res = $v->start.",null,".$value;
			push(@array,$res);
		}
	}
	

	my $lpos="0 ";
	my $lmother="0 ";
	my $lfather="0 ";
	
	##############################
	#
	# calcul de la fraction foetal
	#
	##############################
	my @foetal_frac;
	foreach my $values (@array)
	{
		my ($p,$m,$f)=split(/,/,$values);
		
		$lpos .= $p." ";
		$lmother .= $m." ";
		$lfather .= $f." "; 
		
		push(@foetal_frac,$f) unless $f eq "null";
	}
	

	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@foetal_frac);
	my $ff_mean = $stat->mean;
	$ff_mean = sprintf ("%.2f", $ff_mean);

	my $sv = $stat->standard_deviation();
	$sv = sprintf ("%0.2f", $sv);

	
	$hres->{'POSball'} = $lpos;
	$hres->{'MOTHER'} = $lmother;
	$hres->{'FATHER'} = $lfather;
	$hres->{'mean'} = $ff_mean. " +/- ". $sv;

	

	
	####################  balance allelique calculée chez le pere ##########################
	my $vf;
	
	$vf = $father->getVectorOrigin($chr);
	$vf &= $chr->getVectorSubstitutions();
	
	#on retire de l'affichage les variants dont le ratio est inferieur à 20
	$vf &= $vquality_father;
	

	my $list = to_array($vf,$chr->name);
	my @array;
		
	my $lpos_dad="";
	my $ltrans_dad="";
	my $ldad="";
	
	if( scalar(@$list) > 0)
	{
		$project->setListVariants($list);
		my $x;

		# boucle sur les variants 
		while(my $v = $project->nextVariant){
			$x++;
			my $value = $v->getPourcentAllele($father);
			my $transmis = 0;

			next unless $v->getPourcentAllele($mother) eq "-";		# on passe les variants existants chez l'autre parent
			
			$transmis = 1 if( ($v->getPourcentAllele($patient) ne "-")  && ($v->getPourcentAllele($mother) eq "-"));
			next if $value eq "-";
	
			my $res = $v->start.",".$value.",".$transmis;
			push(@array,$res);
		}
	}
	else
	{
		my $res = "0,0,0";
		push(@array,$res);
	}
	
	foreach my $values (@array)
	{
		my ($pos,$vfather,$transmis)=split(/,/,$values);
		$lpos_dad .= $pos." ";
		$ldad .= $vfather." ";
		$ltrans_dad .= $transmis." ";
	}
	
	$hres->{'POSball_dad'} = $lpos_dad;
	$hres->{'BAll_dad'} = $ldad;
	$hres->{'transmission_dad'} = $ltrans_dad;
	
	
	###################  balance allelique calculée chez la mère ###########################
	my $vm;
	
	$vm = $mother->getVectorOrigin($chr);
	$vm &= $chr->getVectorSubstitutions();
	
	#on retire de l'affichage les variants dont le ratio est inferieur à 20
	$vm &= $vquality_mother;

	my $list = to_array($vm,$chr->name);
	$project->setListVariants($list);
	my @array;
	my $x;
	
	my $lpos_mom="";
	my $ltrans_mom;
	my $lmom="";

	if( scalar(@$list) > 0)
	{
		$project->setListVariants($list);
		my $x;

		# boucle sur les variants 
		while(my $v = $project->nextVariant){
			$x++;
			my $value = $v->getPourcentAllele($mother);
			my $transmis = 0;
		
			next unless $v->getPourcentAllele($father) eq "-";	# on passe les variants existants chez l'autre parent
			
			$transmis=1 if( ($v->getPourcentAllele($patient) ne "-") && ($v->getPourcentAllele($father) eq "-"));
			next if $value eq "-";
	
			my $res = $v->start.",".$value.",".$transmis;
			push(@array,$res);
		}
	}
	else
	{
		my $res = "0,0,0";
		push(@array,$res);
	}
	
	foreach my $values (@array)
	{
		my ($pos,$vmother,$transmis)=split(/,/,$values);
		$lpos_mom .= $pos." ";
		$lmom .= $vmother." ";
		$ltrans_mom .= $transmis." ";
	}
	
	$hres->{'POSball_mom'} = $lpos_mom;
	$hres->{'BAll_mom'} = $lmom;
	$hres->{'transmission_mom'} = $ltrans_mom;
	
}
else
{
	confess("Not trio");
	exit(0);
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
	
	#@t = sort { $a->{POS} <=> $b->{POS}} @$listHash;
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
			next if (($x%20 != 0) && $project->isGenome);
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