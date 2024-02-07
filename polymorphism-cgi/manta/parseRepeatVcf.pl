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
#use Number::Format qw(:subs);

use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";

use GBuffer;
use GenBoProject;
use GenBoCache;



########################################################################
# Parse les fichiers produits par l'option ExpansionHunter de dragen
########################################################################

my $cgi = new CGI;
my $projectname = $cgi->param('project');
my $patientname = $cgi->param('patient');

my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => $projectname);
my $dir = $project->getSVDir()."TRIPLETS/";

#print $cgi->header('text/json-comment-filtered');


my @listHashRes;
	
###################################
#  lecture des vcf et 
#  gestion du projet via les objets Genbo
###################################

my $caller = "expansionhunter";



#	my $file_in = $thePatient->_getCallingSVFileWithMethodName($caller,"variations");
	my $file_in = $dir.$patientname.".repeats.vcf.gz";
	
	# pour lire le vcf 
	my $hRepeats;
	
	# ouverture du fichier zippé
	my $fd;
	my $ligne;
	my $chrom;
	my $pos;
	my $id;
	my $ref;
	my $alt;
	my $qual;
	my $filter;
	my $info;
	my $format;
	my $values;
	
	my $genotype;

	
	open($fd," zcat $file_in | ") or die("open: $!");
	my  @champsDescription;
	 
	while( defined( $ligne = <$fd> )) 
	{
		chomp($ligne);
		
		# on skipe le Header
		next if ($ligne =~ m/^#/ );
		
		# traitement des lignes suivantes pour le fichier en entier
		# pour supprimer le retour à la ligne
		#$ligne =~ s/\n//;
			
			my @champs = split(/\t/,$ligne);
			$chrom = $champs[0];
			$pos = $champs[1];
			$id =$champs[2];
			$ref =$champs[3];
			$alt =$champs[4];
			$qual =$champs[5];
			$filter =$champs[6];
			$info =$champs[7];
			$format =$champs[8];
			$values =$champs[9];
			
			# aller chercher les valeurs du champs info
			my @champsInfo = split(/;/,$info);
			
			my ($a,$end) = split(/=/,$champsInfo[0]);
			my ($b,$nbRepeatRef) = split(/=/,$champsInfo[1]);
			my ($c,$lengthRepeatRef) = split(/=/,$champsInfo[2]);
			my ($d,$motif) = split(/=/,$champsInfo[3]);
			my ($e,$varId) = split(/=/,$champsInfo[4]);
			my ($f,$repId) = split(/=/,$champsInfo[5]);
			
			
			# aller chercher le genotype et autres informations du dernier champs
			my @champsFormat = split(/:/,$format);
			my @champsValues = split(/:/,$values);
			
			my $ind=0;
			my $val;
			foreach my $f (@champsFormat)
			{ 
				chomp($f);
				$val = $champsValues[$ind];
				$ind++;
				$hRepeats->{$repId}->{$f} = $val if ( ($f eq "SO") ||  ($f eq "GT"));

			}
			
			# aller chercher le nombre de repetitions
			my @champsAlt = split(/,/,$alt);			
			
			
			#remplir la table de hash
			$hRepeats->{$repId}->{"id"}= $varId;
			$hRepeats->{$repId}->{"gene"}= $repId;
			$hRepeats->{$repId}->{"chrom"}= $chrom;
			$hRepeats->{$repId}->{"start"}= $pos;
			$hRepeats->{$repId}->{"end"}= $end;
			$hRepeats->{$repId}->{"nbRepeatRef"}= $nbRepeatRef;
			$hRepeats->{$repId}->{"lengthRepeatRef"}= $lengthRepeatRef;
			$hRepeats->{$repId}->{"motif"}= $motif;
			$hRepeats->{$repId}->{"nbRepeatPat"}= $alt;
			

			# pour le json final
			push( @listHashRes, { %{$hRepeats->{$repId}} } );
	} #fin de la boucle sur les lignes
	


if ( (scalar(@listHashRes) == 0) )
{
	my $hash;
	$hash->{"id"}= "-";
	$hash->{"gene"}= "-";
	$hash->{"chrom"}= "-";
	$hash->{"start"}= "-";
	$hash->{"end"}= "-";
	$hash->{"nbRepeatRef"}= " Ask bioinformatics platform";	
	$hash->{"lengthRepeatRef"}= "-";
	$hash->{"motif"}= " Expansion Hunter did'd run on this project";
	$hash->{"nbRepeatPat"}= "-";
	$hash->{"SO"}= "-";
	$hash->{"GT"}= "-";
	push( @listHashRes, { %{$hash} } );
}



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
	my @t = sort {$a->{gene} cmp $b->{gene}} @$listHash;

	$hash->{'identifier'} = 'id';
	$hash->{'label'}      = 'id';
	$hash->{'items'}      = \@t;
	
	print $cgi->header('text/json-comment-filtered');
	print encode_json $hash;
	print "\n";
}


