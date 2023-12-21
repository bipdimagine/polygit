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



die("\n\nNo -project option... Die...\n\n") unless ($projectname);





#  instanciation d'un objet project_cache et des objets patient et chromosome
my $project = $buffer->newProjectCache( -name => $projectname, -typeFilters=>"");
my $patient = $project->getPatient("$patient_name");
#my $chr = $project->getChromosome($chr_name);


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
	
	my $lpos=" ";
	my $lmother=" ";
	my $lfather=" ";
	
	my $lpos_dad="";
	my $ltrans_dad="";
	my $ldad="";
	
	my $lpos_mom="";
	my $ltrans_mom;
	my $lmom="";
	
	my @array_plasma;
	my @array_father;
	my @array_mother;
	
	my @foetal_frac;
	
	my $posP = 0 ;
	my $posF = 0 ;
	my $posM = 0 ;
	
	# boucle sur les chromosomes
	foreach my $chr ( @{$project->getChromosomes()} )
	{
	
	my $chrname = $chr->name();
	warn $chrname;
	
	next if $chrname eq "MT";
	
	$chrname = "X" if ($chrname == 23);
	$chrname = "Y" if ($chrname == 24);

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



	 #boucle sur les variants 
 
	while(my $v = $project->nextVariant){
		$posP += 100;		
		my $p1 = $v->getPourcentAllele($mother);
		my $p2 = $v->getPourcentAllele($father);
		my $value = $v->getPourcentAllele($patient);
		my $vdepth = $v->getDepth($patient);   
		
	
		next if $value eq "-";
		next if $vdepth < 20;
		
		#my $pos = $v->start + $x;
		
		if ($p1 ne "-"){
			my $res = $posP.",".$value.",null";
			push(@array_plasma,$res);
		}

		if ($p2 ne "-"){
			my $res = $posP.",null,".$value;
			push(@array_plasma,$res);
		}
	}
	

	
	##############################
	#
	# calcul de la fraction foetal
	#
	##############################
	
	foreach my $values (@array_plasma)
	{
		my ($p,$m,$f)=split(/,/,$values);
	
		$lpos .= $p." ";
		$lmother .= $m." ";
		$lfather .= $f." "; 
		
		push(@foetal_frac,$f) unless $f eq "null";
	}
	

	
#	####################  balance allelique calculée chez le pere ##########################
#	my $vf;
#	
#	$vf = $father->getVectorOrigin($chr);
#	$vf &= $chr->getVectorSubstitutions();
#	
#	#on retire de l'affichage les variants dont le ratio est inferieur à 20
#	$vf &= $vquality_father;
#	
#
#	my $list = to_array($vf,$chr->name);
#
#	if( scalar(@$list) > 0)
#	{
#		$project->setListVariants($list);
#
#
#		# boucle sur les variants 
#		while(my $v = $project->nextVariant){
#			$posF += 100;
#			my $value = $v->getPourcentAllele($father);
#			my $transmis = 0;
#
#			next unless $v->getPourcentAllele($mother) eq "-";		# on passe les variants existants chez l'autre parent
#			
#			$transmis = 1 if( ($v->getPourcentAllele($patient) ne "-")  && ($v->getPourcentAllele($mother) eq "-"));
#			next if $value eq "-";
#	
#			my $res = $posF.",".$value.",".$transmis;
#			push(@array_father,$res);
#		}
#	}
#	else
#	{
#		my $res = "0,0,0";
#		push(@array_father,$res);
#	}
#	
#	foreach my $values (@array_father)
#	{
#		my ($p,$vfather,$transmis)=split(/,/,$values);
#		$lpos_dad .= $p." ";
#		$ldad .= $vfather." ";
#		$ltrans_dad .= $transmis." ";
#	}
#	
#	
#	
#	###################  balance allelique calculée chez la mère ###########################
#	my $vm;
#	
#	$vm = $mother->getVectorOrigin($chr);
#	$vm &= $chr->getVectorSubstitutions();
#	
#	#on retire de l'affichage les variants dont le ratio est inferieur à 20
#	$vm &= $vquality_mother;
#
#	my $list = to_array($vm,$chr->name);
#	$project->setListVariants($list);
#
#	
#	if( scalar(@$list) > 0)
#	{
#		$project->setListVariants($list);
#
#		# boucle sur les variants 
#		while(my $v = $project->nextVariant){
#			$posM += 100;
#			my $value = $v->getPourcentAllele($mother);
#			my $transmis = 0;
#		
#			next unless $v->getPourcentAllele($father) eq "-";	# on passe les variants existants chez l'autre parent
#			
#			$transmis=1 if( ($v->getPourcentAllele($patient) ne "-") && ($v->getPourcentAllele($father) eq "-"));
#			next if $value eq "-";
#	
#			my $res = $posM.",".$value.",".$transmis;
#			push(@array_mother,$res);
#		}
#	}
#	else
#	{
#		my $res = "0,0,0";
#		push(@array_mother,$res);
#	}
#	
#	foreach my $values (@array_mother)
#	{
#		my ($p,$vmother,$transmis)=split(/,/,$values);
#		$lpos_mom .= $p." ";
#		$lmom .= $vmother." ";
#		$ltrans_mom .= $transmis." ";
#	}
#	
} # fin de la boucle sur les chromosomes
	

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

#	$hres->{'POSball_dad'} = $lpos_dad;
#	$hres->{'BAll_dad'} = $ldad;
#	$hres->{'transmission_dad'} = $ltrans_dad;
#
#	
#	$hres->{'POSball_mom'} = $lpos_mom;
#	$hres->{'BAll_mom'} = $lmom;
#	$hres->{'transmission_mom'} = $ltrans_mom;
	
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