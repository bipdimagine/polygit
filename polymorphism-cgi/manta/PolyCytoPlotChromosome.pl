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

#----------------------------------------------------------------------------------------------------------------------------------------------
#  recupère les valeurs du patient.bins.bed correspondant au CNV pour faire le plot des ratios
#-----------------------------------------------------------------------------------------------------------------------------------------------

# Instanciation d'un objet buffer pour gérer les connexions aux bases de données  
my $buffer = GBuffer->new();		

# recupere les options  (ici -project echantillon_condition1 echantillon_condition2) 
my $cgi = new CGI;
my $projectname = $cgi->param('project');
my $patient_name = $cgi->param('patient');
my $chr_name = $cgi->param('chr');


die("\n\nNo -project option... Die...\n\n") unless ($projectname);

$chr_name = "X" if ($chr_name == 23);
$chr_name = "Y" if ($chr_name == 24);

#  instanciation d'un objet project_cache et des objets patient et chromosome
my $project = $buffer->newProjectCache( -name => $projectname, -typeFilters=>"");
my $patient = $project->getPatient("$patient_name");
my $chr = $project->getChromosome($chr_name);

my $deb = 1;
my $fin= $chr->length();

# pour savoir si on a ou pas les parents
my $patientFamily = $patient->getFamily();
my $trio =0;
$trio=1 if ( ($patient->isChild()) &&  ($patientFamily->mother()) && ($patientFamily->father()) );

# pour construire le path d'acces au bins.bed de wisecondor 
my $dir = $project->getVariationsDir("wisecondor");

my $filein_pat_bin = $dir."/".$patient_name."_bins.bed.gz";
my $filein_pat_abs  = $dir."/".$patient_name."_aberrations.bed.gz";
my $filein_mom_bin;
my $filein_dad_bin;
my $filein_mom_abs;
my $filein_dad_abs;
my $hres;
my @listHashRes;

# faire un tabix pour recuperer dans le fichier les lignes correspondant au chromosome
my $tabix=$buffer->software("tabix");
my $cmd1 = $tabix." ".$filein_pat_bin. " ".$chr_name.":".$deb."-".$fin." | awk '{print \$2,\$5,\$6}' ";

my $res1 = `$cmd1`;
chomp($res1);

my @tabValPat = split(/\n/,$res1);
my $lpos="";
my $lratio="";
my $lzscore="";
my $i=0;

foreach my $ligne (@tabValPat)
{
	my ($x,$y,$z)=split(/ /,$ligne);
	if (($z<5) && ($z>-5))
	{
			$i++;
			next if ($i%20 != 0);
	}
	$lpos .= $x." ";
	$lratio .= $y." ";
	$lzscore .= $z." ";
}

# pour regarder si le chromosome porte un evenement wisecondor
my $cmd2 = "zmore ".$filein_pat_abs." | awk '{print \$1,\$6}' ";
my $res2 = `$cmd2`;
chomp($res2);
my $blue = 0;
my $red = 0;

my @tabValPat2 = split(/\n/,$res2);
foreach my $ligne (@tabValPat2)
{
	my ($chr,$type)=split(/ /,$ligne);
	next if $chr ne $chr_name;
	
	$blue=1 if $type eq "gain";
	$red=2 if $type eq "loss";
}

$hres->{'TYPE'} = int($red)+int($blue);
$hres->{'PLOT'} = $trio;
$hres->{'POSratio'} = $lpos;
$hres->{'RATIO'} = $lratio;
$hres->{'Zscore'} = $lzscore;

if($trio)
{
	# acces aux fichiers bed des parents si ils existent
	my $mother = $patient->getFamily()->getMother();
	my $father = $patient->getFamily()->getFather();
	my $mothername = $mother->name();
	my $fathername = $father->name();
	
	$filein_mom_bin = $dir."/".$mothername."_bins.bed.gz";
	$filein_dad_bin = $dir."/".$fathername."_bins.bed.gz";
	$filein_mom_abs = $dir."/".$mothername."_aberrations.bed.gz";
	$filein_dad_abs = $dir."/".$fathername."_aberrations.bed.gz";

	# ratio WC pour la mere 
	my $cmd = $tabix." ".$filein_mom_bin. " ".$chr_name.":".$deb."-".$fin." | awk '{print \$2,\$5,\$6}' ";
	my $res1 = `$cmd`;
	chomp($res1);

	my @tabValMom = split("\n",$res1);
	my $lposmum="";
	my $lratiomum="";
	my $lzscoremum="";
	$i=0;

	foreach my $ligne (@tabValMom)
	{
		my ($x,$y,$z)=split(/ /,$ligne);
		if (($z<5) && ($z>-5))
		{
			$i++;
			next if ($i%20 != 0);
		}
		$lposmum .= $x." ";
		$lratiomum .= $y." ";
		$lzscoremum .= $z." ";
	}

	$hres->{'POSratio_mom'} = $lposmum;
	$hres->{'RATIO_mom'} = $lratiomum;
	$hres->{'Zscore_mom'} = $lzscoremum;
	
	# ratio WC pour le pere
	my $cmd2 = $tabix." ".$filein_dad_bin. " ".$chr_name.":".$deb."-".$fin." | awk '{print \$2,\$5,\$6}' ";
	my $res2 = `$cmd2`;
	chomp($res2);

	my @tabValDad = split("\n",$res2);
	my $lposdad="";
	my $lratiodad="";
	my $lzscoredad="";
	$i=0;

	foreach my $ligne (@tabValDad)
	{
		my ($x,$y,$z)=split(/ /,$ligne);
		if (($z<5) && ($z>-5))
		{
			$i++;
			next if ($i%20 != 0);
		}
		$lposdad .= $x." ";
		$lratiodad .= $y." ";
		$lzscoredad .= $z." ";
	}

	$hres->{'POSratio_dad'} = $lposdad;
	$hres->{'RATIO_dad'} = $lratiodad;
	$hres->{'Zscore_dad'} = $lzscoredad;
	
	
	#######################
	#     calcul dela fréquence allélique
	#######################
	my $vector;
	my $nb = $chr->getVariantsVector->Size();
	#my $no = $chr->cache_lmdb_variations();

	my $vector_pos = $chr->getVectorByPosition($deb,$fin);

	my $v1;
	my $v2;
	my $v1bis;
	my $v2bis;

	#######  version patrick avec les variants Ho   ########

	# dup
	$v1 = $mother->getVectorOriginHo($chr);
	$v1 &= $chr->getVectorSubstitutions();
 	$v1 -= $father->getVectorOrigin($chr);
	$v1 &= $patient->getVectorOriginHe($chr);

	$v2 = $father->getVectorOriginHo($chr);
	$v2 &= $chr->getVectorSubstitutions();
 	$v2 -= $mother->getVectorOrigin($chr);
	$v2 &= $patient->getVectorOriginHe($chr);
	
	#del
	$v1bis = $mother->getVectorOriginHe($chr);
	$v1bis &= $chr->getVectorSubstitutions();
 	$v1bis -= $father->getVectorOrigin($chr);
	$v1bis &= $patient->getVectorOrigin($chr);

	$v2bis = $father->getVectorOriginHe($chr);
	$v2bis &= $chr->getVectorSubstitutions();
 	$v2bis -= $mother->getVectorOrigin($chr);
	$v2bis &= $patient->getVectorOrigin($chr);
	
	
	#my $v3 = $v1+$v2;
	my $v3 = $v1+$v2+$v1bis+$v2bis;
		
	my $list3 = to_array($v3,$chr->name);	
	$project->setListVariants($list3);

	my @array;
	my $x;

	# boucle sur les variants 
	while(my $v = $project->nextVariant){
		$x++;		
		my $p1 = $v->getPourcentAllele($mother);
		my $p2 = $v->getPourcentAllele($father);
		my $value = $v->getPourcentAllele($patient);
		next if $value eq "-";
		
		if ($p1 ne "-"){
			my $res = $v->start.",".$value.",null";
			push(@array,$res);
		}
		$p1= "" if $p1 eq "-";

		if ($p2 ne "-"){
			my $res = $v->start.",null,".$value;
			push(@array,$res);
		}
		next if  $p2 eq "-";
		$p2= "" if $p2 eq "-";
	}
	
	my $lpos="";
	my $lmother="";
	my $lfather="";

	foreach my $values (@array)
	{
		my ($p,$m,$f)=split(/,/,$values);
		
		$lpos .= $p." ";
		$lmother .= $m." ";
		$lfather .= $f." ";
		
	}

	$hres->{'POSball'} = $lpos;
	$hres->{'MOTHERball'} = $lmother;
	$hres->{'FATHERball'} = $lfather;
	
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
			next if ($x%5 != 0);
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

