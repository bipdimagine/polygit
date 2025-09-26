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
my $plasma_name = $cgi->param('plasma');
my $printres = $cgi->param('printres');

# pour suprimer le bruit de fond chez le père  
my $minRatio = "_ratio_20";
my $rm_noise="yes";

die("\n\nNo -project option... Die...\n\n") unless ($projectname);


#  instanciation d'un objet project_cache et des objets plasma et bt
my $project = $buffer->newProjectCache( -name => $projectname, -typeFilters=>"");
my $plasma = $project->getPatient($plasma_name);

# fichiers de resultats
my $pathDPNI = "/data-isilon/download/manue/DPNI/";
my $fd;
my $file = $pathDPNI.$plasma_name.".txt";

open($fd,'>>',$file);
print($fd $projectname."\t".$plasma_name."\n");
print($fd "ID\t Depth\t VAF\t Freq\t DJV\t Score(Depth*FF)\n");				# on sauvegarde les caractéristiques du variant

#die;

# pour savoir si on a ou pas les parents
my $plasmaFamily = $plasma->getFamily();
my $trio =0;
$trio=1 if ( ($plasma->isChild()) &&  ($plasmaFamily->mother()) && ($plasmaFamily->father()) );

# pour le json final
my $hres;
my @listHashRes;

# pour le fichier de resultat
my @array_variant_plasma;

$hres->{'PLOT'} = $trio;
	
if($trio)
{
	# acces aux fichiers des parents
	my $mother = $plasma->getFamily()->getMother();
	my $father = $plasma->getFamily()->getFather();
	my $mothername = $mother->name();
	my $fathername = $father->name();
	
	# pour compter les variants chez le père avant et après élimination du bruit de fond
	my $nbvarFather1;
	my $nbvarFather2;
	
	
	# Profondeur moyenne du père et du plasma
	my $meanFatherCov=$father->coverage()->{"mean"};
	my $meanPlasmaCov=$plasma->coverage()->{"mean"};


	my $vmother;
	my $vfather;
	my $vplasmaFromFather;
	
	my $lpos=" ";
	my $lmother=" ";
	my $lfather=" ";
	

	my @foetal_frac;	
	my $posP = 0 ;
	my $posVariants=0;

	# boucle sur les chromosomes
	foreach my $chr ( @{$project->getChromosomes()} )
	{
		
		my $chrname = $chr->name();
		warn $chrname;
		
		next if $chrname eq "MT";
		
		$chrname = "X" if ($chrname == 23);
		$chrname = "Y" if ($chrname == 24);
		
		
		#$posP = $posP+25000-$posVariants;
		#$posVariants=0; #remise à zero

		### selection des variants uniques à un des deux parents

		# variant de la mère absent chez le père
		$vmother = $mother->getVectorOrigin($chr);								# tous les variants de la mère
		$vmother &= $chr->getVectorSubstitutions();								# qui sont des substitutions
 		$vmother -= $father->getVectorOrigin($chr);								# qui ne sont pas chez le père
 	

		# variants du père absent chez la mère
		$vfather = $father->getVectorOrigin($chr);								# tous les variants du père
		$vfather &= $chr->getVectorSubstitutions();								# qui sont des substitutions
		$nbvarFather1 += count_values($vfather,$chr->name);
 		$vfather -= $mother->getVectorOrigin($chr);								# qui ne sont pas chez la mère
 		
 		# élimination du bruit de fond chez le père
 		#===========================================
		my $vector_ratio_name_father = $father->name.$minRatio;
		unless ($rm_noise eq "no")
		{
			my $vquality_father = $chr->getVectorScore($vector_ratio_name_father);
			$vfather &= $vquality_father;
		}		

		$nbvarFather2 += count_values($vfather,$chr->name);
		
		
		
		######################################
		#
		# (1) Calcul de la fraction foetal
		#
		######################################
	
		$vplasmaFromFather = $vfather & $plasma->getVectorOriginHe($chr);					# variants du père > minRatio, absent chez la mère retrouvés dans le plasma
	
		my $v3 = $vmother+$vplasmaFromFather;												# liste = somme des variants de la mère non presents chez le père  + variants du plasma transmis par le père
		my $list3 = to_array($v3,$chr->name);	
		$project->setListVariants($list3);

		my @array_plasma;
		
		#boucle sur les variants 
		while(my $v = $project->nextVariant)
		{
			next if $v->getDP($mother) < 10;
			my $p1 = $v->getPourcentAllele($mother);							# valeur ou - si le variant n'est pas présent chez la mère
			my $p2 = $v->getPourcentAllele($father);							# valeur ou - si le variant retrouvé dans le plasma ne vient pas du père 
			my $value = $v->getPourcentAllele($plasma);							# la fréquence allélique du variant 
			my $vdepth = $v->getDepth($plasma);   
		
	
			next if $value eq "-";
			
			
			$posP += 100;
			#$posVariants += 100;
			
		
			if ($p1 ne "-"){													# le variant est présent chez la mère
				my $res = $posP.",".$value.",null";
				push(@array_plasma,$res);
			}

			if ($p2 ne "-"){													# le variant n'est pas présent chez la mère et vient du père 
				my $res = $posP.",null,".$value;
				push(@array_plasma,$res);
				
				# on sauvegarde le variant
				push(@array_variant_plasma,$v);				
			}
		}
	
	
		foreach my $values (@array_plasma)
		{
			my ($p,$m,$f)=split(/,/,$values);
	
			$lpos .= $p." ";
			$lmother .= $m." ";
			$lfather .= $f." "; 
		
			push(@foetal_frac,$f) if $f ne "null";							# dans foetal_frac la balance allélique des variants du plasma provenant du père
		}
		


	} # fin de la boucle sur les chromosomes


	# calcul de la fraction foetal
	my $size = $#foetal_frac+1;
	warn "size=".$size;
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@foetal_frac);
	my $ff_mean = $stat->mean;
	$ff_mean = sprintf ("%.2f", $ff_mean);
	
	warn "fraction_foetal_moyenne = ".$ff_mean;
	warn "From ".$size." variants";
	warn "min = ".min(@foetal_frac);
	warn "max = ".max(@foetal_frac);
	
	# et de la deviation standard
	my $sv = $stat->standard_deviation();
	$sv = sprintf ("%.2f", $sv);

	
	# on stocke les valeurs pour le plot
	$hres->{'POSball'} = $lpos;
	$hres->{'MOTHER'} = $lmother;
	$hres->{'FATHER'} = $lfather;
	$hres->{'mean'} = $ff_mean. " +/- ". $sv;
	$hres->{'score'} = ($meanPlasmaCov*$ff_mean)/100;
	
	# Pour le plot
	push( @listHashRes, { %{ $hres } } );
	printJson( \@listHashRes );
	
	# Sauvegarde les caractéristiques du variant
	if ($printres)
	{
		foreach my $v (@array_variant_plasma)
		{
			my $vid = $v->id;
			my $vdepth = $v->getDepth($plasma);
			my $score = ($vdepth*$ff_mean)/100;
			my $vfreq = $v->frequency*100;
			my $djv = "-";
			$djv = $v->nb_deja_vu_samples;
			$vfreq = sprintf ("%.3f", $vfreq);
			warn $vid;
			print($fd $vid."\t".$vdepth."\t".$v->getPourcentAllele($father).";".$v->getPourcentAllele($plasma)."\t".$vfreq."%\t".$djv."\t".$score."\n");
			
		}
	}
	
}
else
{
	confess("missing parents");
	exit(0);
}




exit(0);





################
#  methods
################

sub printJson {
	my ($listHash) = @_;
	my $hash;

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

sub count_values {
	my ( $v, $name ) = @_;
	my $set  = Set::IntSpan::Fast::XS->new( $v->to_Enum );
	my $iter = $set->iterate_runs();
	my $x=0;
	while ( my ( $from, $to ) = $iter->() ) {
		for my $member ( $from .. $to ) {
			$x++;
		}
	}
	return $x;
}
