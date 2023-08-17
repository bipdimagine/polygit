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
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";

use GBuffer;
use GenBoProject;
use GenBoCache;
use Digest::MD5::File qw(dir_md5_hex file_md5_hex url_md5_hex file_md5);
use layout;
use export_excel; 
use export_data;
use Capture::Tiny ':all';
#----------------------------------------------------------------------------------------------------------------------------------------------
#  recupère les valeurs du patient.bins.bed correspondant au CNV pour faire le plot des ratios
#-----------------------------------------------------------------------------------------------------------------------------------------------

# Instanciation d'un objet buffer pour gérer les connexions aux bases de données  
my $buffer = GBuffer->new();		

# recupere les options  (ici -project echantillon_condition1 echantillon_condition2) 
my $cgi = new CGI;
my $project_name = $cgi->param('project');
my $patient_name = $cgi->param('patient');
my $chr_name = $cgi->param('chr');
my $deb_cnv = $cgi->param('debcnv');
my $fin_cnv = $cgi->param('fincnv');
my $type =  $cgi->param('type');

die("\n\nNo -project option... Die...\n\n") unless ($project_name);

$chr_name = "X" if ($chr_name == 23);
$chr_name = "Y" if ($chr_name == 24);

#  instanciation d'un objet project_cache et des objets patient et chromosome
my $project = $buffer->newProjectCache( -name => $project_name, -typeFilters=>"");
my $patient = $project->getPatient("$patient_name");
my $chr = $project->getChromosome($chr_name);


my %hkeys = $cgi->Vars;
my @keys;
my $string;


foreach my $k  (sort {$a cmp $b} keys %hkeys){
	next if $k =~ /force/;
	next if $k =~ /user/;
	push(@keys,"$k");
	my $c = $hkeys{$k};
	$c =~ s/\"//g;
	$c =~ s/\+/ /g;
	push(@keys,$c);
}
push(@keys,file_md5_hex($Bin."/PolyCytoPlotCNV.pl") );
my $cache_id= "$chr_name;$project_name;$patient_name;$type"."polycyto_plot_cnv".file_md5_hex($Bin."/PolyCytoPlotCNV.pl");
#warn $cache_id;
my $no_cache;
my $text;
#
#$no_cache = $patient->get_lmdb_cache("r");
#warn $no_cache->filename;
$text = $no_cache->get_cache($cache_id);
$no_cache->close();


$| = 1;
if ($text){
	warn " CACHE";
	print $text;
	exit(0);
}



my $deb = 1;
my $fin= $chr->length();

# pour savoir si on a ou pas les parents
my $patientFamily = $patient->getFamily();

# trio = 3 si on a les deux parents 2 si on a que la mere ou 1  si on a que le pere
my $trio =0;		
$trio = 1 if ( $patient->isChild() &&  $patientFamily->father() );
$trio +=2 if ( $patient->isChild()  && $patientFamily->mother() );


# pour construire le path d'acces au bins.bed de wisecondor 
my $dir = $project->getVariationsDir("wisecondor");
my $filein_pat = $dir."/".$patient_name."_bins.bed.gz";
my $filein_mom;
my $filein_dad;
my $hres;
my @listHashRes;

# faire un tabix pour recuperer dans le fichier les lignes correspondant au chromosome
my $tabix=$buffer->software("tabix");
my $cmd = $tabix." ".$filein_pat. " ".$chr_name.":".$deb."-".$fin." | awk '{print \$2,\$5}' ";
my $res1 = `$cmd`;
chomp($res1);

my @tabValPat = split(/\n/,$res1);
my $lpos="";
my $lratio="";
my $i=0;


foreach my $ligne (@tabValPat)
{
	my ($x,$y)=split(/ /,$ligne);
	if (($y<0.2) && ($y>-0.2))
	{
		$i++;
		next if ($i%10 != 0);
	}
	$lpos .= $x." ";
	$lratio .= $y." ";
}


$hres->{'PLOT'} = $trio;
$hres->{'POSratio'} = $lpos;
$hres->{'RATIO'} = $lratio;

my $mother;
my $father;

if (($trio == 2) || ($trio==3))	# only mom or both
{
	# acces aux fichiers bed des parents si ils existent
	$mother = $patient->getFamily()->getMother();
	my $mothername = $mother->name() if defined($mother);
	
	$filein_mom = $dir."/".$mothername."_bins.bed.gz";
	
	# ratio WC pour la mere 
	my $cmd = $tabix." ".$filein_mom. " ".$chr_name.":".$deb."-".$fin." | awk '{print \$2,\$5}' ";
	my $res1 = `$cmd`;
	chomp($res1);

	my @tabValMom = split(/\n/,$res1);
	my $lpos="";
	my $lratio="";
	my $i = 0;
	
	foreach my $ligne (@tabValMom)
	{
		my ($x,$y)=split(/ /,$ligne);
		if (($y<0.2) && ($y>-0.2))
		{
			$i++;
			next if ($i%10 != 0);
		}
		$lpos .= $x." ";
		$lratio .= $y." ";
	}

	$hres->{'POSratio_mom'} = $lpos;
	$hres->{'RATIO_mom'} = $lratio;
}

if (($trio ==1)	|| ($trio==3)	) # only dad or both
{	
	# acces aux fichiers bed des parents si ils existent
	$father = $patient->getFamily()->getFather();
	my $fathername = $father->name() if defined($father);
	
	$filein_dad = $dir."/".$fathername."_bins.bed.gz";
	
	# ratio WC pour le pere
	my $cmd = $tabix." ".$filein_dad. " ".$chr_name.":".$deb."-".$fin." | awk '{print \$2,\$5}' ";
	my $res1 = `$cmd`;
	chomp($res1);

	my @tabValDad = split(/\n/,$res1);
	my $lpos="";
	my $lratio="";
	my $i=0;
	
	foreach my $ligne (@tabValDad)
	{
		my ($x,$y)=split(/ /,$ligne);
		if (($y<0.2) && ($y>-0.2))
		{
			$i++;
			next if ($i%10 != 0);
		}
		$lpos .= $x." ";
		$lratio .= $y." ";
	}

	$hres->{'POSratio_dad'} = $lpos;
	$hres->{'RATIO_dad'} = $lratio;
}	

	
if ($trio == 3)	#both parents
{
	
	#-----------------------------------------------------------
	#  et calcul de la balance allelique si trio
	#------------------------------------------------------------
	
	my $chr_fasta_name = $chr->fasta_name;
	my $res;
	
	my $v1;
	my $v2;
	my $v1bis;
	my $v2bis;
	
	# (1) selection des snps informatifs a partir des variants Ho : version patrick
	
	if ($type == 1)	# DEL on regarde tous les variants de l'enfant / au niveau de la deletion tous les variants He du parent non délété deviennent Ho
	{
		$v1 = $mother->getVectorOriginHe($chr);
		$v1 &= $chr->getVectorSubstitutions();
 		$v1 -= $father->getVectorOrigin($chr);
		$v1 &= $patient->getVectorOrigin($chr);

		$v2 = $father->getVectorOriginHe($chr);
		$v2 &= $chr->getVectorSubstitutions();
 		$v2 -= $mother->getVectorOrigin($chr);
		$v2 &= $patient->getVectorOrigin($chr);
	}

	if ($type >= 2 ) # DUP les snps informatifs sont ceux Ho chez un des parents non présents chez l'autre et He chez l'enfant 
	{
		$v1 = $mother->getVectorOriginHo($chr);
		$v1 &= $chr->getVectorSubstitutions();
 		$v1 -= $father->getVectorOrigin($chr);
		$v1 &= $patient->getVectorOriginHe($chr);

		$v2 = $father->getVectorOriginHo($chr);
		$v2 &= $chr->getVectorSubstitutions();
 		$v2 -= $mother->getVectorOrigin($chr);
		$v2 &= $patient->getVectorOriginHe($chr);
	}
	
	
	my $v3 = $v1+$v2;	
	my $list3 = to_array($v3,$chr->name);	
	$project->setListVariants($list3);

	my @array;
	my $x;

	# boucle sur les variants 
	my $hvar;
	#my $t =time;
	while(my $v = $project->nextVariant) {
		$hvar->{$v->id} ++;
		$x++;			
		my $p1 = $v->getPourcentAllele($mother);			# 100 si mother Ho value si He  - sinon
		my $p2 = $v->getPourcentAllele($father);			# 100 si father Ho value si He - sinon
		my $value = $v->getPourcentAllele($patient);		# freaquence allelique de l'enfant
		
		next if $value eq "-";
		
		if ($p1 ne "-"){									# la mere donne allele
			my $res = $v->start.",".$value.",null";			# frequence allelique de l'enfant attribue a la mere 
			push(@array,$res);
		}
		$p1= "" if $p1 eq "-";

		if ($p2 ne "-"){									# le pere donne allele
			my $res = $v->start.",null,".$value;			# frequence allelique de l'enfant attribue au pere
			push(@array,$res);
		}
		next if  $p2 eq "-";
		$p2= "" if $p2 eq "-";
	}
	#warn abs(time-$t);
	#die($x);
	
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
	$hres->{'MOTHER'} = $lmother;
	$hres->{'FATHER'} = $lfather;
}
else
{
	my $chr_fasta_name = $chr->fasta_name;
	my $res;

	my $v1;
	
	$v1 = $patient->getVectorOrigin($chr);
	$v1 &= $chr->getVectorSubstitutions();

	my $list = to_array($v1,$chr->name);
	$project->setListVariants($list);
	my @array;
	my $x;
	
	my $lpos="";
	my $lpatient="";

	# boucle sur les variants 
	while(my $v = $project->nextVariant){
		$x++;
		my $value = $v->getPourcentAllele($patient);
		next if $value eq "-";
	
		my $res = $v->start.",".$value;
			push(@array,$res);
		}
	
	foreach my $values (@array)
	{
		my ($pos,$vpatient)=split(/,/,$values);
		$lpos .= $pos." ";
		$lpatient .= $vpatient." ";
	}
	
	$hres->{'POSball'} = $lpos;
	$hres->{'PATIENT'} = $lpatient;
}



push( @listHashRes, { %{ $hres } } );


my $stdout2 = tee_stdout {
printJson( \@listHashRes );
};
 
$no_cache = $patient->get_lmdb_cache("w");
$no_cache->put_cache_text($cache_id,$stdout2,2400);
$no_cache->close();
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
	my ( $v, $name ,$step) = @_;
	my $set  = Set::IntSpan::Fast::XS->new( $v->to_Enum );
	my $iter = $set->iterate_runs();
	my @t;
	my $x=0;
	
	while ( my ( $from, $to ) = $iter->() ) {
		for my $member ( $from .. $to ) {
			$x++;
			next if ($x%10 != 0);
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


