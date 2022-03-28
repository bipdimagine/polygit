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
use Number::Format qw(:subs);

use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";

use GBuffer;
use GenBoProject;
use GenBoCache;

use layout;
use export_excel;
use export_data;

###########################################################################################################
#  1) recupere les infos concernant tous les SV du projet entre en argument et prealablement freeze : nom du fihier = project.allSV
#  2) regroupe les SV des patients entres en argument et qui ont la même liste de gènes sous un id_global
#  3) Si le project contient des trios on regarde la transmission a partir du SVglobal
#  4) filtre le resultats en fonction des donnéees definies via l'interface
#  5) construit le Json resultant pour affichage
#
#  Deux cas possibles pour l'affichage des resultats selon si le projet contient ou non  des trios
#  Si on a des trios on regarde en plus la transmission des variants
############################################################################################################

my $cgi = new CGI;
my $projectname = $cgi->param('projectname');
my $thePatientName     = $cgi->param('filename');
my $minlength   = $cgi->param('minlength');
my $maxlength   = $cgi->param('maxlength');
my $maxfreq     = $cgi->param('maxfreq');
my $genotype    = $cgi->param('genotype');
my $select_best = $cgi->param('select_best');
my $hide_sd     = $cgi->param('hide_sd');
my $transm =	$cgi->param('transmission');    # 0 si le projet ne contient pas de trio
my $chr         = $cgi->param('chrom');
my $cytoband    = $cgi->param('cytoband');
my $listOfGenes = $cgi->param('genes');
my $caller      = $cgi->param('caller');
my $dejavu      = $cgi->param('dejavu');

# pour stocker les SV
my $hSVPos;          		# tous par positions chevauchantes
my $hSVCore; 				# les noyaux, chevauchant à au moins 80%
my $hSVallGrouped;    # id commun pour ceux touchant les même gènes

# pour stocker les bandes
#my %hcytoband;

# pour stocker les duplication segmental
my %hdupseg;


# pour construire le projet
my $buffer  = GBuffer->new();
my $project = $buffer->newProjectCache( -name => $projectname );

# pour distinguer parmi les autres patients du projet
# ceux qui sont de la même famille
my $thePatient = $project->getPatient($thePatientName);
my $thePatientFamily = $thePatient->getFamily();

my $patients_all = $project->getPatients();
my @p_family;
my @p_other;
my $thePatient_motherName="-";
my $thePatient_fatherName="-";
my $nb_parents = 0;

foreach my $p (@$patients_all) {
	my $pfam = $p->getFamily();
	
	if ($pfam eq $thePatientFamily)
	{
		push(@p_family,$p);
		if ($p->isMother() && ($p->name() ne $thePatientName))
		{
			$thePatient_motherName = $p->name();
			$nb_parents++;
		}
		if ($p->isFather() && ($p->name() ne $thePatientName))
		{
			$thePatient_fatherName = $p->name();
			$nb_parents++;
		}
	}
	else
	{
		push( @p_other, $p );
	}
}


# pour construire le path d'acces au bam pour IGV
#TODO: a mettre au propre MDR
my $patBam = $thePatient->bamUrl;
my ( $path_Bam, $nothing );
if ($patBam =~ /bwa/) {
	( $path_Bam, $nothing ) = split( /bwa/, $patBam );
	$path_Bam .= "bwa/";
}
else {
	my @lMethods = @{$thePatient->alignmentMethods()};
	my $align_method = $lMethods[0];
	( $path_Bam, $nothing ) = split( /$align_method/, $patBam );
	$path_Bam .= "$align_method/";
}

my $cytodir        = $project->dirCytoManue;
my $fichier_dupseg = $cytodir . "/super_Duplication.gff3";
my $ligne;
my $fdds;

# pour stocker les coords des SV
my $hSVcoord;
my $htree;

##########################################
#
#  Lecture du fichier dupseg
#
##########################################

# ouverture du fichier dupseg
open( $fdds, $fichier_dupseg ) or die("open: $!");

# lecture ligne par ligne
while ( defined( $ligne = <$fdds> ) ) {

	next unless ( $ligne =~ m/^chr/ );
	my @champs = split( /\t/, $ligne );
	my @info   = split( /;/,  $champs[8] );

	my ( $n, $idd )       = split( /=/, $info[0] );
	my ( $i, $percentid ) = split( /=/, $info[1] );

	my $chr   = $champs[0];
	my $start = $champs[3];
	my $end   = $champs[4];

	my $idLocation = $idd . ":" . $start . "_" . $end;

	# version avec IntervalTree
	if ( !exists $hdupseg{$chr} ) {
		my $DupSegTree = Set::IntervalTree->new;
		$hdupseg{$chr} = $DupSegTree;
	}

	$hdupseg{$chr}->insert( $idLocation, $start, $end );

}

#################################################################
#
# pour le json final recuperation du repertoire de sauvegarde via objet genbo
#
##################################################################

my $varDir = $project->getVariationsDir();
my $SVFile = $varDir . "/" . $projectname . ".allSV";


my $hSVAll = retrieve($SVFile) or die "Can't retrieve datas from " . $SVFile . " !\n";

foreach my $id (keys %{$hSVAll})
{
	my @infospatients  = split(/;/,$hSVAll->{$id}->{'PATIENT'});
	foreach my $infopat (@infospatients) 
	{
		$hSVAll->{$id}->{'PATINFO'}->{$infopat}++;
	}
}

my $hSVAllThePatient;

my @Filtered_listHashRes;
my @tabgenes = split( ',', $listOfGenes );
my $resultats_filtres;

my @tabcallers=("wisecondor","canvas");

# traiter le cas particulier des SV fragmentés (wc et canvas)
foreach my $p (@$patients_all) 
{
	my %hSVThePatient_tmp;
	selectSV_from_ThePatient($p,\%hSVThePatient_tmp);
	gatherSV_from_samecaller(\%hSVThePatient_tmp,"wisecondor");
	gatherSV_from_samecaller(\%hSVThePatient_tmp,"canvas");
}

gatherSV_byPosition();
makeGlobal_id_byPosition();
getInfoPatient();

#getGlobalModelTransmission() if $transm;

my $n;
foreach my $global_id ( keys %{$hSVallGrouped} ) {
	my $SV_ok = 1;


	# filtre sur le score caller
	next  if ( ( $hSVallGrouped->{$global_id}->{'SCORECALLER'} < 3 ) && $select_best );

	# filtre sur les chromosomes
	next if ( ( $hSVallGrouped->{$global_id}->{'CHROM'} != $chr ) && ( $chr ne "all" )&& ( $chr ne "noXY" ) );
	next if ( ( $hSVallGrouped->{$global_id}->{'CHROM'} == 23 )&& ( $chr eq "noXY" ) );
	next if ( ( $hSVallGrouped->{$global_id}->{'CHROM'} == 24 )&& ( $chr eq "noXY" ) );
	
	# filtre sur la taille et sur les dupseg
	next if ( $hSVallGrouped->{$global_id}->{"SVLEN"} < $minlength );
	unless ( $maxlength eq "nomax" ) { next if ( $hSVallGrouped->{$global_id}->{"SVLEN"} > $maxlength ); }
	next if ( ( $hSVallGrouped->{$global_id}->{"DUPSEG"} > 40 ) && ( $hide_sd == 1 ) );
	
	# filtre sur les frequences
	my $goldLfreq= -1;
	my $goldGfreq = -1;
	
	my @Gtab = split(/ /,$hSVallGrouped->{$global_id}->{'GOLD_G_freq'});
	my @Ltab = split(/ /,$hSVallGrouped->{$global_id}->{'GOLD_L_freq'});
	
	foreach my $f (@Gtab)
	{
		$goldGfreq = $f if ($f>$goldGfreq && $f ne "-");
	}
	
	foreach my $f (@Ltab)
	{
		$goldLfreq = $f if ($f>$goldLfreq && $f ne "-");
	}
	
	$goldGfreq = "-" if $goldGfreq <0 ;
	$goldLfreq = "-" if $goldLfreq <0 ;
		
	$hSVallGrouped->{$global_id}->{'GOLD_G_freq'} = $goldGfreq;
	$hSVallGrouped->{$global_id}->{'GOLD_L_freq'} = $goldLfreq;
	
	unless ( $maxfreq eq "nomax" ) {
		next if ( $hSVallGrouped->{$global_id}->{'GOLD_G_freq'} > $maxfreq );
		next if ( $hSVallGrouped->{$global_id}->{'GOLD_L_freq'} > $maxfreq );
	}

	# filtre sur une liste de genes
	my @tabgenes_annot = split( ',', $hSVallGrouped->{$global_id}->{"GENES"} );
	my $pass   = 0;

	foreach my $gene_annot (@tabgenes_annot) {
		foreach my $gene (@tabgenes) {
			if ( ( $gene_annot eq $gene ) || ( $gene eq "all" ) ) {
				$pass = 1;
				last;
			}
		}
		last if ( $pass == 1 );
	}
	next unless ($pass);

	# filtre sur le dejavu in this project
	$pass = 0;
	my $countother = 0;
	if ( $hSVallGrouped->{$global_id}->{"PATIENT"} =~ m/$thePatientName/ ) 
	{				
				$pass = 1;
				my $theothers;
				foreach my $pother (@p_other) 
				{
						my $name = $pother->name();
						
						if ( $hSVallGrouped->{$global_id}->{"PATIENT"} =~m/$name/ )
						{
							$countother++;
							$theothers .= $name.",";
						}
				}
				$hSVallGrouped->{$global_id}->{"OTHER"} = $countother.";".$theothers;
				$pass = 0 if ( $countother > $dejavu );
	}
	next unless($pass);

	# filtre sur le genotype
	$pass = 1;
	my $infoGenotype = $hSVallGrouped->{$global_id}->{"GENOTYPE"};

	unless ( $genotype eq "both" ) {
		$pass = 0;

		if (
			( $genotype eq "ho" )
			&& (   ( $infoGenotype eq "1/1" )
				|| ( $infoGenotype eq "1/2" )
				|| ( $infoGenotype eq "1" )
				|| ( $infoGenotype eq "-" ) )
		  )
		{
			$pass = 1;
		}
		if (   ( $genotype eq "he" )&& ( ( $infoGenotype eq "0/1" ) || ( $infoGenotype eq "-" ) ) )
		{
			$pass = 1;
		}

	}
	next unless ($pass);

	# filtre sur une liste de cytobandes
	my @tabcytoband = split( ',', $cytoband );
	my $pass = 0;
	foreach my $cb (@tabcytoband) {
		if (   ( $hSVallGrouped->{$global_id}->{"CYTOBAND"} =~ m/$cb/ )|| ( $cb eq "all" ) )
		{
			$pass = 1;
			next;
		}
	}
	next unless ($pass);

	# filtre sur les modeles de transmission
	if ($transm) {
		my @tabtransm = split( ',', $transm );

		my $pass = 0;
		foreach my $tm (@tabtransm) {
			if (   ( $hSVallGrouped->{$global_id}->{'MODEL_GLOBAL'} =~ m/$tm/ )
				|| ( $tm eq "all" ) )
			{
				$pass = 1;
				next;
			}
		}
		next unless ($pass);
	}
	
	#on suprime ce qui ne sert plus
	delete($hSVallGrouped->{$global_id}->{'PATIENT'});
	delete($hSVallGrouped->{$global_id}->{'PATINFO'});
	push( @Filtered_listHashRes, { %{ $hSVallGrouped->{$global_id} } } );
}

printJson( \@Filtered_listHashRes );
exit(0);

###############################################################################
#
#	methodes
#
#################################################################################

sub printRes {
	my $fres;
	my $fname = $thePatientName . "_" . $minlength . "_" . $maxlength . "_" . $genotype;
	open( $fres, '>', $fname );

	my $titre = "patient=". $thePatientName . " minlength=". $minlength. " maxlength=". $maxlength. " genotype=". $genotype . "\n";
	my $col = "Event Locus Genotype Length Genes";
	printf( $fres $titre . "\n" );
	printf( $fres $col . "\n" );

	my ($listres) = @_;
	my @tabres = split( /next/, $listres );
	foreach my $SV (@tabres) {
		my ( $rien, $gsv ) = split( /et/, $SV );
		my @tabisv = split( /;/, $gsv );

		my $nbchamps = scalar(@tabisv);
		for ( my $i = 0 ; $i < $nbchamps - 1 ; $i++ ) {
			my @tabinfo = split( / /, $tabisv[$i] );
			my ( $type, $chr, $deb, $fin ) = split( /_/, $tabinfo[0] );
			my ( $a, $b, $genes ) = split( / /, $tabisv[-1] );
			my $length = $fin - $deb;
			my $svres =$type . " chr" . $chr . ":" . $deb . "-" . $fin . " " . $tabinfo[3] . " " . $length . " " . $genes . "\n";
			printf( $fres $svres );
		}
	}
}

sub printJson {
	my ($listHash) = @_;
	my $hash;
	my @t;

#	@t = sort {$b->{SCORECALLER} <=> $a->{SCORECALLER}or $a->{SVTYPE} cmp $b->{SVTYPE}} @$listHash if !$select_best;
#	@t = sort {$a->{CHROM} <=> $b->{CHROM}or $a->{SVTYPE} cmp $b->{SVTYPE}or $a->{POS_INT} <=> $b->{POS_INT}} @$listHash if $select_best;
	
	@t = sort {$b->{SCORECALLER} <=> $a->{SCORECALLER}or $a->{CHROM} <=> $b->{CHROM}or $a->{SVTYPE} cmp $b->{SVTYPE}or $a->{POS_INT} <=> $b->{POS_INT}} @$listHash;
	$hash->{'identifier'} = 'id';
	$hash->{'label'}      = 'id';
	$hash->{'items'}      = \@t;
	print $cgi->header('text/json-comment-filtered');
	print encode_json $hash;
	print "\n";
}

sub selectSV_from_ThePatient()
{
	my ($pat,$hpat) = @_;
	my $patname = $pat->name();

	foreach my $id ( keys %{$hSVAll} ) 
	{
			next unless  ($hSVAll->{$id}->{"PATIENT"} =~ m/$patname/);
			$hpat->{$id}= $hSVAll->{$id};
	}
}

sub gatherSV_from_samecaller()
{
	my ($hSVThePatient_tmp,$thecaller) = @_;
	
	my $htree;
	my $hintspan;

	# 1)  detecter les SV chevauchants 
	
	# creer les arbres
	foreach my $id ( keys %{$hSVThePatient_tmp} ) 
	{
		next unless ($hSVThePatient_tmp->{$id}->{"PATIENT"} =~ m/$thecaller/ );
		my ($t,$c,$d,$f) = split( /_/, $id );
		$htree->{$t}->{$c}= Set::IntervalTree->new;
		$hintspan->{$t}->{$c}= Set::IntSpan::Fast->new();
	}

	# remplir les arbres :  regrouper les SV chevauchants
	foreach my $id ( keys %{$hSVThePatient_tmp} ) 
	{
		next unless ($hSVThePatient_tmp->{$id}->{"PATIENT"} =~ m/$thecaller/ );
		my ($t,$c,$d,$f) = split( /_/, $id );
		
		$htree->{$t}->{$c}->insert($id,$d,$f);
	}
	
	# 2) associer a chaque SV trouve avec le même caller ceux qui lui sont proches	
	foreach my $id ( keys %{ $hSVThePatient_tmp} )
	{	
		if ( $hSVThePatient_tmp->{$id}->{"PATIENT"} =~ m/$thecaller/ )
		{
			#on cherche les regroupements
			my ($t,$c,$dtheSV,$ftheSV) = split( /_/, $id );
		
			my $padding = 0.1 * abs($ftheSV-$dtheSV);
		
			my $tab_id = $htree->{$t}->{$c}->fetch($dtheSV-$padding,$ftheSV+$padding);

			# on regroupe les id dans un intspan
			my $gdeb=0;
			my $gend=0;
			
			foreach my $ind_id ( @$tab_id ) 
			{
				my ($t,$c,$d,$f) = split(/_/, $ind_id);
				
		    	$gdeb = $d if ( ($d < $gdeb) || ($gdeb==0));
		 		$gend = $f if ( ($f > $gend) );
 			}
			$hintspan->{$t}->{$c}->add_range($gdeb,$gend);
		}
		
		if ( ($hSVThePatient_tmp->{$id}->{"PATIENT"} !~ m/wisecondor/) && ($hSVThePatient_tmp->{$id}->{"PATIENT"} !~ m/canvas/) )
		{
			# onrecopie tel quel ceux qui ne sont vu ni par canvas ni par wisecondor
			$hSVAllThePatient->{$id} =  $hSVThePatient_tmp->{$id};
		}
	}
	gather_id($hSVThePatient_tmp,$hSVAllThePatient,$htree,$hintspan,$thecaller);
}

sub gather_id() {
	
		my ($hSV_in,$hSV_out,$htree,$hintspan,$thecaller) = @_;
		my $nb = 0;
		
		
		
		foreach my $type (keys %{$hintspan})
		{
			foreach my $chr (keys %{$hintspan->{$type}})
			{
				my $liste_of_bornes = $hintspan->{$type}->{$chr}->as_string() if defined( $hintspan->{$type}->{$chr});
				$liste_of_bornes .= "," if ($liste_of_bornes !~ m/,/);
				warn $liste_of_bornes if ( ($type eq "DUP") && ($chr eq "8") && $thecaller eq "canvas");
				my @tab_bornes = split(/,/,$liste_of_bornes);
				
				foreach my $bornes (@tab_bornes)
				{
					
					my ($deb,$end) = split(/-/,$bornes);
					
					my $tab_ind_id = $htree->{$type}->{$chr}->fetch($deb,$end);
					
					my $global_id = $type."_".$chr."_".$deb."_". $end;
		
					# enregistrer les infos dans une nouvelle table de hash
					$hSV_out->{$global_id}->{'id'}     = $global_id;
					$hSV_out->{$global_id}->{'SVLEN'}  = abs( $deb - $end );
					$hSV_out->{$global_id}->{'SVTYPE'} = $type;
		
		
					# retrouver les informations correspondant aux differents id regroupes
					$hSV_out->{$global_id}->{'RANKAnnot'} = 0;
					$hSV_out->{$global_id}->{'ScoreGene'} = 0;
		
					foreach my $id ( @$tab_ind_id )
					{
							# ce qui depend de la liste de gene et est donc identique pour les differents id
							$hSV_out->{$global_id}->{'CHROM'} = $hSV_in->{$id}->{'CHROM'};
							$hSV_out->{$global_id}->{$thecaller} = $hSV_in->{$id}->{$thecaller};
			
							$hSV_out->{$global_id}->{'GOLD'} .= " ".$hSV_in->{$id}->{'GOLD'} unless ( $hSV_out->{$global_id}->{'GOLD'} =~ m/$hSV_in->{$id}->{'GOLD'}/ );
							$hSV_out->{$global_id}->{'GOLD_G_freq'} .= " ".$hSV_in->{$id}->{'GOLD_G_freq'}, unless ( $hSV_out->{$global_id}->{'GOLD_G_freq'} =~ m/$hSV_in->{$id}->{'GOLD_G_freq'}/ );
							$hSV_out->{$global_id}->{'GOLD_L_freq'} .= " ".$hSV_in->{$id}->{'GOLD_L_freq'} unless ( $hSV_out->{$global_id}->{'GOLD_L_freq'} =~ m/$hSV_in->{$id}->{'GOLD_L_freq'}/ ); 
							$hSV_out->{$global_id}->{'OMIN_MG'} .= " ".$hSV_in->{$id}->{'OMIN_MG'} unless ( $hSV_out->{$global_id}->{'OMIN_MG'} =~ m/$hSV_in->{$id}->{'OMIN_MG'}/ );
							$hSV_out->{$global_id}->{'dbVar_event'} .= $hSV_in->{$id}->{'dbVar_event'} unless ( $hSV_out->{$global_id}->{'dbVar_event'} =~ m/$hSV_in->{$id}->{'dbVar_event'}/ );
							$hSV_out->{$global_id}->{'dbVar_status'} .= $hSV_in->{$id}->{'dbVar_status'} unless ( $hSV_out->{$global_id}->{'dbVar_status'} =~ m/$hSV_in->{$id}->{'dbVar_status'}/ );
							$hSV_out->{$global_id}->{'GENES'} .= $hSV_in->{$id}->{'GENES'}.",";
			
							if ( $hSV_in->{$id}->{'BREAKPOINT'} ne  "-" )
							{
									$hSV_out->{$global_id}->{'BREAKPOINTS'} .= $hSV_in->{$id}->{'BREAKPOINTS'} 
									unless ( $hSV_out->{$global_id}->{'BREAKPOINT'} =~ m/$hSV_in->{$id}->{'BREAKPOINT'}/ )
							}

							# score annotation lié aux gènes couverts par le variant
							my $scores = $hSV_in->{$id}->{'RANKAnnot'};
							my ( $rankAnnot, $scoreGene ) = split( /,/, $scores );
			
							$hSV_out->{$global_id}->{'RANKAnnot'} = max($rankAnnot , $hSV_out->{$global_id}->{'RANKAnnot'});
							$hSV_out->{$global_id}->{'ScoreMaxGene'} = max ($rankAnnot , $hSV_out->{$global_id}->{'ScoreMaxGene'});

							# ce qui depend du patient et/ou des id
	
							$hSV_out->{$global_id}->{'PATIENT'} .= $hSV_in->{$id}->{'PATIENT'}; 
							foreach my $k (keys %{$hSV_in->{$id}->{'PATINFO'}})
							{
									$hSV_out->{$global_id}->{'PATINFO'}->{$k}++; 
							}
			
					}
	
					#pour suprimer les doublons
					my $theliste;
		
					my @tabgenes = split(/,/,$hSV_out->{$global_id}->{'GENES'});
		
					foreach my $gene (@tabgenes)
					{
							next  if ($theliste =~ m/$gene/);
							$theliste .= $gene.",";
					}
					$hSV_out->{$global_id}->{'GENES'}=$theliste;
				}
			}
		}
}
	
# pour regrouper les ids recouvrant les mêmes positions
sub gatherSV_byPosition() {

	#my $idtree;
	my $htree;
	
	
	# 1)  detecter les SV chevauchants 
		
	# creer les arbres
	foreach my $id ( keys %{$hSVAllThePatient} ) 
	{
		my $SVtype = $hSVAllThePatient->{$id}->{"SVTYPE"};
		my $chrom  = $hSVAllThePatient->{$id}->{"CHROM"};
		$htree->{$SVtype}->{$chrom}= Set::IntervalTree->new;
	}

	# remplir les arbres :  regrouper les SV chevauchants
	foreach my $id ( keys %{$hSVAllThePatient} ) 
	{
		my $SVtype = $hSVAllThePatient->{$id}->{"SVTYPE"};
		my $chrom  = $hSVAllThePatient->{$id}->{"CHROM"};
		
		my ( $t, $c, $d, $f ) = split( /_/, $id );
		$htree->{$SVtype}->{$chrom}->insert($id,$d,$f);
	}
	
	# 2) associer a chaque SV ceux qui le chevauchent 	
	foreach my $id ( keys %{$hSVAllThePatient} )
	{	
		my $SVtype = $hSVAllThePatient->{$id}->{"SVTYPE"};
		my $chrom  = $hSVAllThePatient->{$id}->{"CHROM"};
		
		my ( $t, $c, $dtheSV, $ftheSV ) = split( /_/, $id );
		
		my $tab_id = $htree->{$SVtype}->{$chrom}->fetch($dtheSV,$ftheSV);
		my @tab_id_tmp;
		
		my $lentheSV = abs($ftheSV-$dtheSV);
	
		# le nombre de gene portes par le SV qui nous interresse
		my @SVtabgenes = split(/,/,$hSVAllThePatient->{$id}->{"GENES"});
		my $SVnbGenes = scalar(@SVtabgenes);
		
		# pour chacun des SVs chevauchant celui qui nous interresse 
		foreach my $ind_id (@$tab_id)
		{

					# On le garde si il chevauche celui qui nous interresse sur au moins 60% de sa longueur et reciproquement
					my ( $t, $c, $dother, $fother) = split( /_/, $ind_id );
				
					my $overlap = min($ftheSV,$fother) - max($dtheSV,$dother);
					next if ( $overlap < (0.6 * $lentheSV));
				
					# et reciproquement
					my $lenother = abs($fother-$dother);
					next if ( $overlap < (0.6 * $lenother) );
            	
					# et si il partage au moins 70% de gènes couverts par le SV
		
					my $nbGenes;
					my @tabgenes = split(/,/,$hSVAllThePatient->{$ind_id}->{"GENES"});
					foreach my $g (@tabgenes)
					{
						$nbGenes++ if ($hSVAllThePatient->{$id}->{"GENES"} =~ m/$g/);
					}
				
					my $p1 = 0;
					$p1 = int(($nbGenes/$SVnbGenes) *  100) unless ($SVnbGenes ==0);
				
					my $p2 = 0;
					$p2 = int(($SVnbGenes/$nbGenes) *  100) unless ($nbGenes ==0);
								
					push(@tab_id_tmp,$ind_id) if (( $p1 > 70)  ||  ($p2 > 70)) ;
					push(@tab_id_tmp,$ind_id);
		}
		
		my @tab_id_sorted  = sort ({$a cmp $b} @tab_id_tmp);
		
		my $list_of_id ="";
		foreach my $v (@tab_id_sorted)
		{
			$list_of_id .= $v.",";
		}
		
		$hSVPos->{$list_of_id}++;
	}
}

sub makeGlobal_id_byPosition() {
	
	my $nb = 0;

	foreach my $k ( keys( %{$hSVPos} ) ) {
		
		my $gdeb=0;
		my $gend=0;
		my $gnum;
		my $gtype;
		
		my $genesList;
		my $index;

		
		my @tab_ind_id = split( /,/, $k);
	
		
		# creer l'id global version union
		foreach my $ind_id ( @tab_ind_id ) 
		{
			my ($t,$c,$d,$f) = split( /_/, $ind_id );
			$gnum =$c;
			$gtype =$t;
		    
		    $gdeb =$d if ($gdeb==0);
		    $gdeb = $d if ($d < $gdeb);
		 	$gend = $f if ( ($f > $gend) );
 		}
	
		my $pos = int($gdeb);
		my $global_id = $gtype . "_" . $gnum . "_" . $gdeb . "_" . $gend;
		
		
		# enregistrer les infos dans une nouvelle table de hash
		$hSVallGrouped->{$global_id}->{'id'}     = $global_id;
		$hSVallGrouped->{$global_id}->{'SVLEN'}  = abs( $gdeb - $gend );
		$hSVallGrouped->{$global_id}->{'SVTYPE'} = $gtype;
		$hSVallGrouped->{$global_id}->{'POSITIONS'} = format_number($gdeb) . "-" . format_number($gend);

		#pour pouvoir trier sur lea position en numerique
		$hSVallGrouped->{$global_id}->{'POS_INT'} = $pos;

		# enregistrer le nom du chromosome sous forme numerique pour le tri
		my $gchr    = "chr" . $gnum;
		my $gnumInt = $gnum;
		$gnumInt = int(23) if ( $gnum eq "X" );
		$gnumInt = int(24) if ( $gnum eq "Y" );
		$hSVallGrouped->{$global_id}->{'CHROM'} = int($gnumInt);

		# pour l'acces a DGV
		my $url_DGV = getDGV_url( $gchr, $gdeb, $gend );
		$hSVallGrouped->{$global_id}->{'DGV'} = $url_DGV;

		# pour les segmental duplication
		$gchr .= "chr" . $gchr unless ( $gchr =~ m/chr/ );
		my $dupseg = "-";
		$dupseg = getDupSeg( $gchr, $gdeb, $gend ) if ( $gend > $gdeb );
		$hSVallGrouped->{$global_id}->{'DUPSEG'} = $dupseg;
		
		# pour les cytobandes
		my @tb;
		my @band;
		my $hband = $project->getChromosome($gnum)->getCytoband( $gdeb, $gend ) if ( $gend > $gdeb );

		foreach my $b ( keys %{ $hband->{'name'} } ) {
			push( @tb, $b );
		}
		@band = sort( { $a cmp $b } @tb );

		$hSVallGrouped->{$global_id}->{'CYTOBAND'} = join( ",", @band );

		# retrouver les informations correspondant aux differents id regroupes
		$hSVallGrouped->{$global_id}->{'PATIENT'} .= $path_Bam . 'et';
		
		$hSVallGrouped->{$global_id}->{'RANKAnnot'} = 0;
		$hSVallGrouped->{$global_id}->{'ScoreGene'} = 0;
		
		foreach my $id ( @tab_ind_id )
		{
			# ce qui depend de la liste de gene et est donc identique pour les differents id
			$hSVallGrouped->{$global_id}->{'GOLD'} .= $hSVAllThePatient->{$id}->{'GOLD'} unless ( $hSVallGrouped->{$global_id}->{'GOLD'} =~ m/$hSVAllThePatient->{$id}->{'GOLD'}/ );
			$hSVallGrouped->{$global_id}->{'GOLD_G_freq'} .= $hSVAllThePatient->{$id}->{'GOLD_G_freq'} unless ( $hSVallGrouped->{$global_id}->{'GOLD_G_freq'} =~ m/$hSVAllThePatient->{$id}->{'GOLD_G_freq'}/ );
			$hSVallGrouped->{$global_id}->{'GOLD_L_freq'} .= $hSVAllThePatient->{$id}->{'GOLD_L_freq'} unless ( $hSVallGrouped->{$global_id}->{'GOLD_L_freq'} =~ m/$hSVAllThePatient->{$id}->{'GOLD_L_freq'}/ ); 
			$hSVallGrouped->{$global_id}->{'OMIN_MG'} .= $hSVAllThePatient->{$id}->{'OMIN_MG'} unless ( $hSVallGrouped->{$global_id}->{'OMIN_MG'} =~ m/$hSVAllThePatient->{$id}->{'OMIN_MG'}/ );
			$hSVallGrouped->{$global_id}->{'dbVar_event'} .= $hSVAllThePatient->{$id}->{'dbVar_event'} unless ( $hSVallGrouped->{$global_id}->{'dbVar_event'} =~ m/$hSVAllThePatient->{$id}->{'dbVar_event'}/ );
			$hSVallGrouped->{$global_id}->{'dbVar_status'} .= $hSVAllThePatient->{$id}->{'dbVar_status'} unless ( $hSVallGrouped->{$global_id}->{'dbVar_status'} =~ m/$hSVAllThePatient->{$id}->{'dbVar_status'}/ );
			$hSVallGrouped->{$global_id}->{'GENES'} .= $hSVAllThePatient->{$id}->{'GENES'}.",";
			
			if ( $hSVAllThePatient->{$id}->{'BREAKPOINT'} ne  "-" )
			{
				$hSVallGrouped->{$global_id}->{'BREAKPOINTS'} .= $hSVAllThePatient->{$id}->{'BREAKPOINTS'} 
				unless ( $hSVallGrouped->{$global_id}->{'BREAKPOINT'} =~ m/$hSVAllThePatient->{$id}->{'BREAKPOINT'}/ )
			}

			# score annotation lié aux gènes couverts par le variant
			my $scores = $hSVAllThePatient->{$id}->{'RANKAnnot'};
			my ( $rankAnnot, $scoreGene ) = split( /,/, $scores );
			
			$hSVallGrouped->{$global_id}->{'RANKAnnot'} = max($rankAnnot , $hSVallGrouped->{$global_id}->{'RANKAnnot'});
			$hSVallGrouped->{$global_id}->{'ScoreGene'} = max ($scoreGene , $hSVallGrouped->{$global_id}->{'ScoreGene'});

			# ce qui depend du patient et/ou des id
			$hSVallGrouped->{$global_id}->{'PATIENT'} .= $hSVAllThePatient->{$id}->{'PATIENT'};
			foreach my $k (keys %{$hSVAllThePatient->{$id}->{'PATINFO'}})
			{
				$hSVallGrouped->{$global_id}->{'PATINFO'}->{$k}++; 
			}
		}
		
		#pour suprimer les doublons
		my $theliste;

		my @tabgenes = split(/,/,$hSVallGrouped->{$global_id}->{'GENES'});
	
		foreach my $gene (@tabgenes)
		{
			next  if ($theliste =~ m/$gene/);
			$theliste .= $gene.",";
		}
		$hSVallGrouped->{$global_id}->{'GENES'}=$theliste;
	}
	
}

sub getInfoPatient() {
	
	my $hT; # pour la transmission
	
	foreach my $gid ( keys %{$hSVallGrouped} ) 
	{
		my $score = 0;
		
		my @tabPatInfo =  keys( %{$hSVallGrouped->{$gid}->{'PATINFO'}});
		my $info_globalvariant  = join(";",@tabPatInfo);
		
		next unless ($info_globalvariant =~  m/$thePatientName/ );
	
		$hSVallGrouped->{$gid}->{"wisecondor"} = "";
		$hSVallGrouped->{$gid}->{"canvas"} = "";
		$hSVallGrouped->{$gid}->{"manta"} = "";
		$hSVallGrouped->{$gid}->{'GENOTYPE'}   = "-";
		
	
		my $wc = 1;
		my $can = 1;
		my $man = 1;
		
		foreach my $patInfo (@tabPatInfo) 
		{
			my ( $pid, $patname, $caller, $genotype, $ped, $statut ) = split( / /, $patInfo );
			
			# les infos du patient passees en argument pour chaque caller
			if ( $patname eq $thePatientName ) {
				# score associe au variant global
				if ( ( $caller eq "wisecondor") && $wc) { $score +=5; $wc=0;}
				if ( ($caller eq "canvas") && $can) { $score +=2; $can=0;}
				if ( ($caller eq  "manta") && $man) { $score +=1; $man=0;}

				$hSVallGrouped->{$gid}->{$caller}  .= $path_Bam . "et" . $patInfo."ou";
				$hT->{$gid}->{$caller}->{'thechild'}=1;
				$hSVallGrouped->{$gid}->{'GENOTYPE'} = $genotype unless ( $caller eq "wisecondor" );
			}
			
			# pour la transmission si thePatient est un enfant
			if ($thePatient->isChild())
			{
				if ( $patname eq $thePatient_motherName ) 
				{
						$hT->{$gid}->{$caller}->{'mother'}=1;
				}
				else
				{
						if ( $patname eq $thePatient_fatherName )
						{
							$hT->{$gid}->{$caller}->{'father'}=1;
						}
				}
			}
		}
		$hSVallGrouped->{$gid}->{'SCORECALLER'} = $score;
	
		# pour la transmission si $thePatientName est un enfant
		if ( $thePatient->isChild())
		{
			foreach my $caller (keys %{$hT->{$gid}})
			{
				$hT->{$gid}->{$caller}->{"TRANSMISSION"}= "-";
				if ($nb_parents > 0)
				{
					$hT->{$gid}->{$caller}->{"TRANSMISSION"} = "mother" if ($hT->{$gid}->{$caller}->{"thechild"}) && ($hT->{$gid}->{$caller}->{"mother"}) && !($hT->{$gid}->{$caller}->{"father"});	
					$hT->{$gid}->{$caller}->{"TRANSMISSION"} = "father"   if ($hT->{$gid}->{$caller}->{"thechild"}) && ($hT->{$gid}->{$caller}->{"father"}) && !($hT->{$gid}->{$caller}->{"mother"});	
					$hT->{$gid}->{$caller}->{"TRANSMISSION"} = "both"     if ($hT->{$gid}->{$caller}->{"thechild"}) && ($hT->{$gid}->{$caller}->{"mother"}) && ($hT->{$gid}->{$caller}->{"father"});	
					$hT->{$gid}->{$caller}->{"TRANSMISSION"} = "Denovo"     if ($hT->{$gid}->{$caller}->{"thechild"}) && !($hT->{$gid}->{$caller}->{"mother"}) && !($hT->{$gid}->{$caller}->{"father"});	
				}
				$hSVallGrouped->{$gid}->{"TRANSMISSION"} .= $hT->{$gid}->{$caller}->{"TRANSMISSION"};
			}
		}
		
		# pour IGV on garde preferentielement les infos de manta puis de canvas et a defaut de wisecondor
		if ( $hSVallGrouped->{$gid}->{"manta"} ) {
			$hSVallGrouped->{$gid}->{"IGV"} = $hSVallGrouped->{$gid}->{"manta"};
		}
		else 
		{
			if ( $hSVallGrouped->{$gid}->{"canvas"} )
		   {
				$hSVallGrouped->{$gid}->{"IGV"} = $hSVallGrouped->{$gid}->{"canvas"};
			}
			else 
			{
				$hSVallGrouped->{$gid}->{"IGV"} =  $hSVallGrouped->{$gid}->{"wisecondor"};
			}
		}
	}
}

sub getDGV_url {
	my ( $chr, $d, $f ) = @_;
	my $url = "http://dgv.tcag.ca/gb2/gbrowse/dgv2_hg19/?name=" . $chr . "%3A". $d . "-"  . $f  . ";search=Search";
	return $url;
}

sub getDupSeg {
	my ( $chr, $d, $f ) = @_;
	my $dupseg;
	my $result;

	my $cov50 = 0;
	my $cov10 = 0;

	# recherche des dupseg chevauchant le fragment
	$result = $hdupseg{$chr}->fetch( $d, $f );
	my $nb_dupseg = scalar(@$result);

	# si il n y en a pas
	return "-" if ( $nb_dupseg == 0 );

	# sinon on creer un intspan
	my $span = Set::IntSpan::Fast::XS->new();

	foreach my $res (@$result) {
		my ( $id,       $coord )  = split( /:/, $res );
		my ( $ds_start, $ds_end ) = split( /_/, $coord );

		# on ne conserve que la partie comprise entre $d et $f
		$ds_start = $d if ( $ds_start <= $d );
		$ds_end   = $f if ( $ds_end >= $f );

		$span->add_range( $ds_start, $ds_end );
	}

	# calcul de la couverture globale des dupseg
	my $globale_cov = 0;
	my $list_pos    = $span->as_string();
	my @tabpos      = split( /,/, $list_pos );

	foreach my $coords (@tabpos) {
		my ( $debut, $fin ) = split( /-/, $coords );
		my $taille = ( $fin - $debut );
		$globale_cov = $globale_cov + $taille;
	}

	return $globale_cov / ( $f - $d ) * 100;
}

#sub getGlobalModelTransmission( ) {
#
#	my $htransmission;
#	my $infos;
#
#	foreach my $gid ( keys %{$hSVallGrouped} ) {
#		$infos = $hSVallGrouped->{$gid}->{'PATIENT'};
#		my ( $nothing, $info_globalvariant ) = split( /et/, $infos );
#		my @tab_info_individual_variant = split( /;/, $info_globalvariant );
#
#		my $child         = 0;
#		my $healthy_child = 0;
#
#		foreach my $infosPatient (@tab_info_individual_variant) 
#		{
#
#			my ( $id, $patname, $c, $GT ) = split( / /, $infosPatient );
#
#			#next if (($GT eq "-") || ($GT eq " "));	# on passe quand on a pas d'information sur le genotype
#			my $pat    = $project->getPatient($patname);
#			my $family = $pat->getFamily()->name();
#
#			my $mother        = 0;
#			my $father        = 0;
#			my $mother_status = 0;
#			my $father_status = 0;
#			my $mothername    = 0;
#			my $fathername    = 0;
#			my $GT_mother;
#			my $GT_father;
#
#			# valeur par defaut pour les parents
#			$htransmission->{$gid}->{$c}->{$family}->{"mother"}->{"GT"} = "0/0" if ( !exists $htransmission->{$gid}->{$c}->{$family}->{"mother"}->{"GT"} );
#			$htransmission->{$gid}->{$c}->{$family}->{"father"}->{"GT"} = "0/0"  if ( !exists $htransmission->{$gid}->{$c}->{$family}->{"father"}->{"GT"} );
#
#			if ( $pat->isChild() && ( $pat->status() == 2 ) ) {
#				$htransmission->{$gid}->{$c}->{$family}->{$patname}->{"GT"} = $GT;
#				$child = 1;
#			}
#
#			if ( $pat->isMother() ) {
#				$mothername    = $patname;
#				$mother_status = $pat->status();
#				$htransmission->{$gid}->{$c}->{$family}->{"mother"}->{"GT"} = $hSVAll->{$id}->{$c}->{$mothername}->{'GT'} if ( exists $hSVAll->{$id}->{$c}->{$mothername}->{'GT'} );
#			}
#
#			if ( $pat->isFather() ) {
#				$fathername    = $patname;
#				$father_status = $pat->status();
#				$htransmission->{$gid}->{$c}->{$family}->{"father"}->{"GT"} = $hSVAll->{$id}->{$c}->{$fathername}->{'GT'} if ( exists $hSVAll->{$id}->{$c}->{$fathername}->{'GT'} );
#			}
#
#			if ( $pat->isChild() && ( $pat->status() == 1 ) ) {
#				$htransmission->{$gid}->{$c}->{$family}->{"healthy_child"} ->{$GT} = 1;   # permet de savoir si il existe un enfant sain 0/1 ou 1/1
#				$healthy_child = 1;
#			}
#		}
#
#		my $transmission_type = "_";
#		foreach my $infosPatient (@tab_info_individual_variant) {
#
#			my ( $id, $p, $c, $GT ) = split( / /, $infosPatient );
#
#			#next if (($GT eq "-") || ($GT eq " "));	# on passe quand on a pas d'information sur le genotype
#
#			my $pat    = $project->getPatient($p);
#			my $family = $pat->getFamily()->name();
#
#			# genotype des enfants atteints et transmission
#
#			if ( $pat->isChild() && ( $pat->status == 2 ) ) {
#
#				my $GT_casIndex = $htransmission->{$gid}->{$c}->{$family}->{$p}->{"GT"};
#				my $GT_mother = $htransmission->{$gid}->{$c}->{$family}->{"mother"}->{"GT"};
#				my $GT_father = $htransmission->{$gid}->{$c}->{$family}->{"father"}->{"GT"};
#
#				if ( $GT_casIndex eq "0/1" ) {
#					if ( ( $GT_mother eq "0/0" ) && ( $GT_father eq "0/0" ) ) {
#						$transmission_type = $p . " : denovo;";
#					}
#					if ( ( $GT_mother eq "0/1" ) && ( $GT_father eq "0/1" ) ) {
#						$transmission_type = $p . " : both";
#					}
#					if ( ( $GT_mother eq "0/1" ) && ( $GT_father eq "1/1" ) ) {
#						$transmission_type = $p . " : both";
#					}
#					if ( ( $GT_mother eq "1/1" ) && ( $GT_father eq "0/1" ) ) {
#						$transmission_type = $p . " : both";
#					}
#					if ( ( $GT_mother eq "0/0" ) && ( $GT_father eq "0/1" ) ) {
#						$transmission_type = $p . " : father";
#					}
#					if ( ( $GT_mother eq "0/1" ) && ( $GT_father eq "0/0" ) ) {
#						$transmission_type = $p . " : mother";
#					}
#					if ( ( $GT_mother eq "0/0" ) && ( $GT_father eq "1/1" ) ) {
#						$transmission_type = $p . " : father";
#					}
#					if ( ( $GT_mother eq "1/1" ) && ( $GT_father eq "0/0" ) ) {
#						$transmission_type = $p . " : mother";
#					}
#					if ( ( $GT_mother eq "1/1" ) && ( $GT_father eq "1/1" ) ) {
#						$transmission_type = $p . " : error";
#					}
#				}
#				
#				if ( ( $GT_casIndex eq "1/1" ) || ( $GT_casIndex eq "1/2" ) ) 
#				{
#					if (   ( $GT_mother eq "0/1" ) && ( $GT_father eq "0/1" ) )    # parents sains He
#					{
#						$transmission_type = $p. " : recessive"
#						  if (!exists($htransmission->{$gid}->{$c}->{$family}->{"healthy_child"}->{"1/1"}));    # pas de frere ou soeur sains Ho
#						$transmission_type = $p. " : unrelevant" if (exists(
#								$htransmission->{$gid}->{$c}->{$family}->{"healthy_child"}->{"1/1"}));    # au moins un frere ou une soeur saine Ho
#					}
#					if ( ( $GT_mother eq "0/1" ) && ( $GT_father eq "0/0" ) ) {
#						$transmission_type = $p . " : mother?";
#					}
#					if ( ( $GT_mother eq "0/0" ) && ( $GT_father eq "0/1" ) ) {
#						$transmission_type = $p . " : father?";
#					}
#					if ( ( $GT_mother eq "1/1" ) || ( $GT_father eq "1/1" ) ) {
#						$transmission_type =
#						  $p . " : unrelevant";   # ne joue pas sur le phenotype
#					}
#					if ( ( $GT_mother eq "0/0" ) && ( $GT_father eq "0/0" ) ) {
#						$transmission_type = $p . " : denovo";
#					}
#				}
#				if ( $GT_casIndex eq "-" ) # presence du variant pour wisecondor
#				{
#
#					if ( ( $GT_mother eq "-" ) && ( $GT_father eq "-" ) ) {
#						$transmission_type = $p . " : both ";
#					}
#					if ( ( $GT_mother eq "0/0" ) && ( $GT_father eq "-" ) ) {
#						$transmission_type = $p . " : father";
#					}
#					if ( ( $GT_mother eq "-" ) && ( $GT_father eq "0/0" ) ) {
#						$transmission_type = $p . " : mother";
#					}
#					if ( ( $GT_mother eq "0/0" ) && ( $GT_father eq "0/0" ) ) {
#						$transmission_type = $p . " : denovo";
#					}
#				}
#				$hSVallGrouped->{$gid}->{'MODEL_GLOBAL'} .=$c . " : " . $transmission_type . ";";
#			}    # boucle sur les enfants atteints
#		}    # boucle sur les callers
#	}    # boucle sur les gid
#}

