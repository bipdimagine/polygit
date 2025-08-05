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
#  recupère les valeurs du patient.bins.bed correspondant au CNV pour faire le plot des ratios
#-----------------------------------------------------------------------------------------------------

# Instanciation d'un objet buffer pour gérer les connexions aux bases de données
my $buffer = GBuffer->new();

# recupere les options
my $cgi          = new CGI;
my $projectname  = $cgi->param('project');
my $patient_name = $cgi->param('patient');
my $chr_name     = $cgi->param('chr');

die("\n\nNo -project option... Die...\n\n") unless ($projectname);

$chr_name = "X" if ( $chr_name == 23 );
$chr_name = "Y" if ( $chr_name == 24 );

#  instanciation d'un objet project_cache et des objets patient et chromosome
my $project =
  $buffer->newProjectCache( -name => $projectname, -typeFilters => "" );
my $patient = $project->getPatient("$patient_name");
my $chr     = $project->getChromosome($chr_name);

# pour savoir si on a ou pas les parents
my $patientFamily = $patient->getFamily();
my $trio          = 0;
$trio = 1
  if ( ( $patient->isChild() )
	&& ( $patientFamily->mother() )
	&& ( $patientFamily->father() ) );

# pour construire le path d'acces au bins.bed de wisecondor
my $dir = $project->getVariationsDir("wisecondor");

my $filein_pat_bin = $dir . "/" . $patient_name . "_bins.bed.gz";
my $filein_pat_abs = $dir . "/" . $patient_name . "_aberrations.bed.gz";

# pour le json final
my $hres;
my @listHashRes;

# faire un tabix pour recuperer dans le fichier les lignes correspondant au chromosome
my $tabix = $buffer->software("tabix");
my $cmd1 =
	$tabix . " "
  . $filein_pat_bin . " "
  . $chr_name
  . " | awk '{print \$2,\$5,\$6}' ";

my $res1 = `$cmd1`;
chomp($res1);

my @tabValPat = split( /\n/, $res1 );
my $lpos      = "";
my $lratio    = "";
my $lzscore   = "";
my $i         = 0;

foreach my $ligne (@tabValPat) {
	my ( $x, $y, $z ) = split( / /, $ligne );
	if ( ( $z < 5 ) && ( $z > -5 ) ) {
		$i++;
		next if ( $i % 20 != 0 );
	}
	$lpos    .= $x . " ";
	$lratio  .= $y . " ";
	$lzscore .= $z . " ";
}

$hres->{'PLOT'}     = $trio;
$hres->{'POSratio'} = $lpos;
$hres->{'RATIO'}    = $lratio;
$hres->{'Zscore'}   = $lzscore;

if ($trio) {

	# acces aux fichiers bed des parents si ils existent
	my $mother     = $patient->getFamily()->getMother();
	my $father     = $patient->getFamily()->getFather();
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

	$v1 = $mother->getVectorOrigin($chr);    # tous les variants de la mère
	$v1 &= $chr->getVectorSubstitutions();   # qui sont des substitutions
	$v1 -= $father->getVectorOrigin($chr);   # qui ne sont pas chez le père
	$v1 &= $patient->getVectorOrigin($chr);  # qui sont chez l'enfant

	$v2 = $father->getVectorOrigin($chr);    # tous les variants du père
	$v2 &= $chr->getVectorSubstitutions();   # qui sont des substitutions
	$v2 -= $mother->getVectorOrigin($chr);   # qui ne sont pas chez la mère
	$v2 &= $patient->getVectorOrigin($chr);  # qui sont chez l'enfant

	my $v3 = $v1 + $v2;

	my $list3 = to_array( $v3, $chr->name );
	$project->setListVariants($list3);

	my @array;
	my $x;

	#boucle sur les variants

	while ( my $v = $project->nextVariant ) {
		$x++;
		my $p1    = $v->getPourcentAllele($mother);
		my $p2    = $v->getPourcentAllele($father);
		my $value = $v->getPourcentAllele($patient);
		
		next if $value == 0;

		if ( $p1 > 0 ) {
			my $res = $v->start . "," . $value . ",null";
			push( @array, $res );
		}
		$p1 = "" if $p1  == 0;

		if ( $p2> 0 ) {
			my $res = $v->start . ",null," . $value;
			push( @array, $res );
		}
		next if $p2 == "0";
	}
	my $lpos    = "";
	my $lmother = "";
	my $lfather = "";

	foreach my $values (@array) {
		my ( $p, $m, $f ) = split( /,/, $values );

		$lpos    .= $p . " ";
		$lmother .= $m . " ";
		$lfather .= $f . " ";

	}

	$hres->{'POSball'} = $lpos;
	$hres->{'MOTHER'}  = $lmother;
	$hres->{'FATHER'}  = $lfather;

	####################  balance allelique calculée chez le pere ##########################
	my $vf;

	$vf = $father->getVectorOrigin($chr);
	$vf &= $chr->getVectorSubstitutions();

	my $list = to_array( $vf, $chr->name );
	$project->setListVariants($list);
	my @array;
	my $x;

	my $lpos_dad = "";
	my $ltrans_dad;
	my $ldad = "";

	# boucle sur les variants
	while ( my $v = $project->nextVariant ) {
		$x++;
		my $value    = $v->getPourcentAllele($father);
		my $transmis = 0;
		next unless $v->getPourcentAllele($mother) == 0;
		$transmis = 1
		  if ( ( $v->getPourcentAllele($patient) >  0)
			&& ( $v->getPourcentAllele($mother)== 0 ) );
		next if $value == 0;

		my $res = $v->start . "," . $value . "," . $transmis;
		push( @array, $res );
	}

	foreach my $values (@array) {
		my ( $pos, $vfather, $transmis ) = split( /,/, $values );
		$lpos_dad   .= $pos . " ";
		$ldad       .= $vfather . " ";
		$ltrans_dad .= $transmis . " ";
	}

	$hres->{'POSball_dad'}      = $lpos_dad;
	$hres->{'BAll_dad'}         = $ldad;
	$hres->{'transmission_dad'} = $ltrans_dad;

	###################  balance allelique calculée chez la mère ###########################
	my $vm;

	$vm = $mother->getVectorOrigin($chr);
	$vm &= $chr->getVectorSubstitutions();

	my $list = to_array( $vm, $chr->name );
	$project->setListVariants($list);
	my @array;
	my $x;

	my $lpos_mom = "";
	my $ltrans_mom;
	my $lmom = "";

	# boucle sur les variants
	while ( my $v = $project->nextVariant ) {
		$x++;
		my $value    = $v->getPourcentAllele($mother);
		my $transmis = 0;

		next unless $v->getPourcentAllele($father)  == 0;
		$transmis = 1
		  if ( ( $v->getPourcentAllele($patient) > 0 )
			&& ( $v->getPourcentAllele($father) == 0 ) );
		next if $value  == 0;

		my $res = $v->start . "," . $value . "," . $transmis;
		push( @array, $res );
	}

	foreach my $values (@array) {
		my ( $pos, $vmother, $transmis ) = split( /,/, $values );
		$lpos_mom   .= $pos . " ";
		$lmom       .= $vmother . " ";
		$ltrans_mom .= $transmis . " ";
	}

	$hres->{'POSball_mom'}      = $lpos_mom;
	$hres->{'BAll_mom'}         = $lmom;
	$hres->{'transmission_mom'} = $ltrans_mom;

}
else {
	my $chr_fasta_name = $chr->fasta_name;
	my $res;

	my $v1;

	$v1 = $patient->getVectorOrigin($chr);
	$v1 &= $chr->getVectorSubstitutions();

	my $list = to_array( $v1, $chr->name );
	$project->setListVariants($list);
	my @array;
	my $x;

	my $lpos     = "";
	my $lpatient = "";

	# boucle sur les variants
	while ( my $v = $project->nextVariant ) {
		$x++;
		my $value = $v->getPourcentAllele($patient);
		next if $value == 0;

		my $res = $v->start . "," . $value;
		push( @array, $res );
	}

	foreach my $values (@array) {
		my ( $pos, $vpatient ) = split( /,/, $values );
		$lpos     .= $pos . " ";
		$lpatient .= $vpatient . " ";
	}

	$hres->{'POSball'} = $lpos;
	$hres->{'PATIENT'} = $lpatient;
}

push( @listHashRes, { %{$hres} } );

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
	$hash->{'items'}      = \@$listHash;

	print $cgi->header('text/json-comment-filtered');
	print encode_json $hash;
	print "\n";
}

sub to_array {
	my ( $v, $name ) = @_;
	my $set  = Set::IntSpan::Fast::XS->new( $v->to_Enum );
	my $iter = $set->iterate_runs();
	my @t;
	my $x = 0;
	while ( my ( $from, $to ) = $iter->() ) {
		for my $member ( $from .. $to ) {
			$x++;
			next if ( ( $x % 20 != 0 ) && $project->isGenome );
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

