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
my $caller_type_flag = {
	"caller_sr" => 1,
	"caller_depth" => 2,
	"caller_coverage" => 4,
};
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
my $dbh = DBI->connect("dbi:ODBC:Driver=DuckDB;Database=:memory:", "", "", { RaiseError => 1 , AutoCommit => 1});
my $parquet_file = $project->getCacheCNV()."/".$project->name.".".$project->id.".parquet";
my $dir = $project->getCacheCNV(). "/rocks/";
my $rocks = GenBoNoSqlRocks->new(dir=>"$dir",mode=>"r",name=>"cnv");



#my $patient = $project->getPatients()->[0];
my $patient = $project->getPatient($patientname);
my $patient_id = $patient->id;
my $colpatient = "p".$patient->id;
my $query = qq{CREATE TABLE cnvs  AS
                           SELECT * 
                           FROM '$parquet_file'
                           WHERE start > -1  and patient = $patient_id  order by start;
	};
$dbh->do($query);
my $nb;
my $hGroupedCNV;
my @Filtered_listHashRes;
my $global_id;


# pour construire le path d'acces au bam pour IGV

my $bamFiles = $patient->bamUrl();
my $bamNames = $patientname;

my $trio;
my $mother;
my $father;
if ($patient->isChild){
	 $mother = $patient->getFamily->getMother();
	 $father = $patient->getFamily->getFather();
	
}




my $sql = qq{select id from cnvs where nb_dejavu_patients <= 20 and len > 1000;};			

my $sth = $dbh->prepare($sql);
$sth->execute();
my $nb = 0;
while (my $row = $sth->fetchrow_hashref) {

	 my $cnv = $rocks->get($row->{id});

		my $global_id = $cnv->{id};
		# pour les colonnes de l'interface
		
		# 1) ce qui ne dépend que du cnv
		
		#les bornes du CNV
		my $gdeb = $cnv->{start};
		my $gend = $cnv->{end};
		my $chr = $project->getChromosome($cnv->{chromosome});
		
		my $chrnum = $chr->name;
 		
 		my  $t = abs ($cnv->{'end'}-$cnv->{'start'});
			my $taille = 1;
	
	 		if ($t > 1000000 ) {
	 				$taille = int( ($t / 1_000_000) + 0.5 ) ."Mb";
	 		}
	 		else
	 		{
				 if ($t > 1000 ) {
	 					$taille = int( ($t / 1_000) + 0.5 ) ."Kb"; 
	 			}
	 			else
	 			{
	 					$taille = $t . "pb"
	 			}
	 		}	
 		
		$hGroupedCNV->{$global_id}->{'id'} = $global_id;
		$hGroupedCNV->{$global_id}->{'TYPE'} = $cnv->{type};
		
		
		$hGroupedCNV->{$global_id}->{'CHROM'} = $cnv->{type}.";".$chr->name;
		
		
		$hGroupedCNV->{$global_id}->{'LEN'}  = $cnv->{type}.";".(abs($gend -$gdeb) +1);
		$hGroupedCNV->{$global_id}->{'DUPSEG'} = $cnv->{dupseg};
		
		
		#pour les gènes emportés par le CNV
		$hGroupedCNV->{$global_id}->{'SCORE_GENES'} = $cnv->{genes}->[0]->{score};
		my $genes_liste;
			foreach my $g (@{$cnv->{genes}})
			{
				my $gname= $g->{name};
				unless ($genes_liste)
				{
					$genes_liste = $g->{name}.";".$g->{score}."##".$genes_liste; # placé en tête de liste
				}
				else
				{
					$genes_liste .= $g->{name}.";".$g->{score}."##";
				} 
			}
			$hGroupedCNV->{$global_id}->{'GENES'} = listGenes($cnv->{genes},$cnv->{id});
			#$cnv->{chromosome}.":".$cnv->{start}."_".$cnv->{end}."##".$genes_liste;
		
		# pour l'accès à gnomAD
		
		$hGroupedCNV->{$global_id}->{'gnomAD'} = "https://gnomad.broadinstitute.org/region/".$chr->name."-".$gdeb."-".$gend."?dataset=gnomad_sv_r2_1" ;
		
		$hGroupedCNV->{$global_id}->{'gnomAD'} = "https://gnomad.broadinstitute.org/region/".$chr->name."-".$gdeb."-".$gend."?dataset=gnomad_cnv_r4" if ($project->current_genome_version eq "HG38");
		$hGroupedCNV->{$global_id}->{'OMIM'} = $cnv->{'OMIM'};
		$hGroupedCNV->{$global_id}->{'DGV'} = "https://dgv.tcag.ca/gb2/gbrowse/dgv2_hg19/?name=".$chr->ucsc_name.":".$gdeb."-".$gend.';search=Search';
		$hGroupedCNV->{$global_id}->{'DGV'} = "https://dgv.tcag.ca/gb2/gbrowse/dgv2_hg38/?name=".$chr->ucsc_name.":".$gdeb."-".$gend.';search=Search' if ($project->current_genome_version eq "HG38");
		
		$hGroupedCNV->{$global_id}->{'GOLD_G_FREQ'} = $cnv->{'GOLD_G_FREQ'};
		$hGroupedCNV->{$global_id}->{'GOLD_L_FREQ'} = $cnv->{'GOLD_L_FREQ'};
		$hGroupedCNV->{$global_id}->{'dbVar_status'} = $cnv->{'dbVar_status'};
		
		$hGroupedCNV->{$global_id}->{'wisecondor_elementary'} = 0;
		$hGroupedCNV->{$global_id}->{'canvas_elementary'} = 0; 
		$hGroupedCNV->{$global_id}->{'manta_elementary'} = 0;
		# pour les infos dépendantes du patient et du caller
		
		$hGroupedCNV->{$global_id}->{'RATIO'} = "-";	   	
		$hGroupedCNV->{$global_id}->{'RATIO'} = "1";	
	
		$hGroupedCNV->{$global_id}->{'PR'} = "-";
		$hGroupedCNV->{$global_id}->{'SR'} = "-";
		$hGroupedCNV->{$global_id}->{'SCORECALLER'} = 0;
		$hGroupedCNV->{$global_id}->{'SCORECALLER'} = $cnv->{score_caller};
		$hGroupedCNV->{$global_id}->{'SCORECNV'} = 1;#$cnv->{score_caller};
		$hGroupedCNV->{$global_id}->{'GT'}=" ";
		
		
		###################################
		# CALLER COVERAGE WISECONDOR STYLE
		###################################
		
		$hGroupedCNV->{$global_id}->{'caller_coverage'} = "-";
		
		if (test_type($cnv,"coverage")){
			$hGroupedCNV->{$global_id}->{'QUAL'} =  "";
			$hGroupedCNV->{$global_id}->{'QUAL'} =   $cnv->{coverage_zscore};
			my ($ratio) = map {$_->{coverage_ratio}} grep{test_type($_,"coverage")} @{$cnv->{cnv_origin}};
			
			 $ratio = sprintf("%.1f", $ratio);
			 
			my $z_qual =max(map {$_->{coverage_zscore}} grep{test_type($_,"coverage")}@{$cnv->{cnv_origin}});
			
			 $z_qual = sprintf("%.1f", $z_qual);
			$hGroupedCNV->{$global_id}->{'caller_coverage'} = "wisecondor;$ratio;$z_qual;".$cnv->{type};
		
		}
		
		
		###################################
		# CALLER DEPTH CANVAS STYLE
		###################################
		$hGroupedCNV->{$global_id}->{'caller_depth'} = "-";
		
		$hGroupedCNV->{$global_id}->{'caller_depth'} = "-";
		if (test_type($cnv,"depth")){
			$hGroupedCNV->{$global_id}->{'QUAL'} =  "";
			$hGroupedCNV->{$global_id}->{'GT'}= $cnv->{gt};
			$hGroupedCNV->{$global_id}->{'QUAL'} =   $cnv->{sr_qual};
			my $cn = max(map {($_->{depth_CN})} grep{test_type($_,"depth")} @{$cnv->{cnv_origin}});
			$cn = log($cn/2) / log(2) if $cn > 0;
			 $cn = sprintf("%.1f", $cn);
			my $depth_qual =max(map {$_->{depth_qual}} @{$cnv->{cnv_origin}});
			$hGroupedCNV->{$global_id}->{'caller_depth'} ="canvas;$cn;$depth_qual;".$cnv->{type};
			
		
		}
		
		
		###################################
		# CALLER SR MANTA STYLE
		###################################
		$hGroupedCNV->{$global_id}->{'caller_sr'} = "-";
		if (test_type($cnv,"sr")){
			$hGroupedCNV->{$global_id}->{'QUAL'} =  "";
			$hGroupedCNV->{$global_id}->{'GT'}= $cnv->{gt};
			$hGroupedCNV->{$global_id}->{'QUAL'} =   $cnv->{sr_qual};
			my $sr1 = 0;
			my @sr_cnv =  grep{test_type($_,"sr")} @{$cnv->{cnv_origin}};
			$sr1 = max(map {$_->{sr1}}  @sr_cnv);
			my $sr2 = max(map {$_->{sr2}}@sr_cnv) ;
			$hGroupedCNV->{$global_id}->{'SR'}= $sr1."/".$sr2;
		
			my $pr1 = max(map {$_->{pr1}}@sr_cnv);
			my $pr2 = max(map {$_->{pr2}} @sr_cnv);
			$hGroupedCNV->{$global_id}->{'PR'}= $pr1."/".$pr2;
			my $sr_qual =max(map {$_->{sr_qual}} @{$cnv->{cnv_origin}});
			my $toto =  $cnv->{'score'}->{'score_caller_sr'}."";
			my $score =  $hGroupedCNV->{$global_id}->{caller_sr} = "manta;$sr1/$sr2;$pr1/$pr2;$sr_qual;".$cnv->{type} ;
		}
	
		
		
		
		
		
		
		$hGroupedCNV->{$global_id}->{'wisecondor_elementary'} = "-";
		#$cnv_origin->{elementary_caller_coverage};
		$hGroupedCNV->{$global_id}->{'manta_elementary'} = "-";#$cnv_origin->{elementary_caller_sr};
		$hGroupedCNV->{$global_id}->{'canvas_elementary'} = "-";#$cnv_origin->{elementary_caller_depth};
		
				
		# pour la transmission
		$hGroupedCNV->{$global_id}->{'TRANSMISSION'}  = "-";
		$hGroupedCNV->{$global_id}->{'TRANSMISSION'}  = getTransmission($cnv,$mother,$father); 	
		
		# pour le dejavu global
		my $nb_projects = $cnv->{'dejavu'}->{'nb_projects'} +0;
		my $nb_patients = $cnv->{'dejavu'}->{'nb_patients'}+0;
		
		my $list_patients = $cnv->{'dejavu'}->{'string'};
		my $nbcoverage = 0;
		 $nbcoverage =  $cnv->{'dejavu'}->{'caller_coverage'} if $cnv->{'dejavu'}->{'caller_coverage'};
		my $nbdepth = 0;
		 $nbdepth =  $cnv->{'dejavu'}->{'caller_depth'}  if $cnv->{'dejavu'}->{'caller_depth'};
		 my $nbsr = 0;
		$nbsr =  $cnv->{'dejavu'}->{'caller_sr'}+0 if  $cnv->{'dejavu'}->{'caller_sr'};
		
		my $flag = 0;
		$flag = 1 if ($nbcoverage > 0);
		my $nb_patient = scalar(@{$project->getPatients}) -1;
		my $itp = scalar (keys %{$cnv->{patients}}) - 1;
		$hGroupedCNV->{$global_id}->{'DEJAVU_G'} = $global_id.";".$nb_projects.";".$nb_patients.";".$nbcoverage.";".$nbdepth.";".$nbsr.";"."$itp/$nb_patient".";".$project->name();
		$hGroupedCNV->{$global_id}->{'PLOT'} =  $cnv->{'type'}.";".$cnv->{'chromosome'}.":".$cnv->{'start'}."-".$cnv->{'end'}.";".$hGroupedCNV->{$global_id}->{'TRANSMISSION'};	 
		$hGroupedCNV->{$global_id}->{"viewer"};
		
		##################
		## IGV 
		###################
		my $viewer;
		
		$viewer->{dgv}->{url} = "https://dgv.tcag.ca/gb2/gbrowse/dgv2_hg19/?name=".$chr->ucsc_name.":".$gdeb."-".$gend.';search=Search';
		$viewer->{dgv}->{url} = "https://dgv.tcag.ca/gb2/gbrowse/dgv2_hg38/?name=".$chr->ucsc_name.":".$gdeb."-".$gend.';search=Search' if ($project->current_genome_version eq "HG38");
		
		$viewer->{gnomad}->{url} = "https://gnomad.broadinstitute.org/region/".$chr->name."-".$gdeb."-".$gend."?dataset=gnomad_sv_r2_1" ;
		$viewer->{gnomad}->{url} = "https://gnomad.broadinstitute.org/region/".$chr->name."-".$gdeb."-".$gend."?dataset=gnomad_cnv_r4" if ($project->current_genome_version eq "HG38");
	
		#PLOT
		#var out   = '<button class="btn btn-classic btn-m" style="border: 1px solid black;padding:3px" id="' + boutonID + '" onclick="launch_plot(\''+ transmission + '\',' + type + ',' + chr + ',' + debcnv + ',' + fincnv + ' )" style="font-size:16px;">&#128200;</button>';
 		$hGroupedCNV->{$global_id}->{'PLOT'} =  $cnv->{'type'}.";".$cnv->{'chromosome'}.":".$cnv->{'start'}."-".$cnv->{'end'}.";".$hGroupedCNV->{$global_id}->{'TRANSMISSION'};	 
		$viewer->{plot}->{type} = 1;
		$viewer->{plot}->{chromosome} = $cnv->{'chromosome'} ;
		$viewer->{plot}->{start} = $cnv->{start} ;
		$viewer->{plot}->{end} = $cnv->{end} ;
		$viewer->{plot}->{transmission} = $hGroupedCNV->{$global_id}->{'TRANSMISSION'};
		
		$viewer->{igv}->{bam_names} = $bamNames;
		$viewer->{igv}->{bam_files} = $bamFiles;
		$viewer->{igv}->{id} = $hGroupedCNV->{$global_id}->{"id"};
		$viewer->{igv}->{locus_start} = $chr->ucsc_name  . ':' . ($cnv->{'start'} -500). '-' . ($cnv->{'start'}+500);
		$viewer->{igv}->{locus_end} = $chr->ucsc_name   . ':' . ($cnv->{'end'} -500). '-' . ($cnv->{'end'}+500);
		
		
	 		
		$viewer->{igv}->{l} = $taille ;
		$hGroupedCNV->{$global_id}->{"VIEWERS"} = encode_json $viewer;
			##################
		## Locus 
		###################
		
		my @cyto = split(",",$cnv->{cytoband});
		my $locus;
		$locus->{start} = $cnv->{start};
		$locus->{end} = $cnv->{end};
		$locus->{type} = $cnv->{type};
		$locus->{len} = $taille;
		$locus->{cytoband} = $cyto[0];
		$locus->{cytoband} = $cyto[0]."..".$cyto[-1] if (@cyto>1);
		$locus->{viewers} = $viewer;
		$hGroupedCNV->{$global_id}->{'LOCUS'} =  encode_json $locus;
		
		
		
		
		$hGroupedCNV->{$global_id}->{"IGV"} = $bamNames.";".$bamFiles.";".$hGroupedCNV->{$global_id}->{"id"};
		$hGroupedCNV->{$global_id}->{'BPManta'}="-";	      # manquant
		$hGroupedCNV->{$global_id}->{'BPManta_mother'}="-";	   # manquant
   		$hGroupedCNV->{$global_id}->{'BPManta_father'}="-";	   # manquant
		 
		# scoreCNV 
		$hGroupedCNV->{$global_id}->{'SCORECNV'}= getScoreCNV($cnv->{genes},$hGroupedCNV->{$global_id}->{'SCORECALLER'},$hGroupedCNV->{$global_id}->{'TRANSMISSION'},$hGroupedCNV->{$global_id}->{'DUPSEG'});
				
		# utile mais non affiché	
		$hGroupedCNV->{$global_id}->{'START'} = $cnv->{'START'};
		$hGroupedCNV->{$global_id}->{'END'} =  $cnv->{'END'};
		$hGroupedCNV->{$global_id}->{'DEJAVU_P'} = getDejaVu_project($cnv); 	
		#die;
		
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
exit(0);


################################################@

 sub listGenes {
 	my ($agenes,$id,$debug) = @_;
 	my $debug ;
 	my @merged;

	my $res;
	$res->{genes} = [];
	$res->{nb_genes} = scalar(@$agenes);
	my $n =0;
	foreach my $g (@$agenes){
		delete $g->{phenotypes};
		push(@{$res->{genes}},$g);
		last if $n >= 4;
		$n ++;
	}
	$res->{project} = $project->name;
	$res->{id} = $id;
	$res->{patient} = $patient->name;
	$res->{nbh} = 0;
	$res->{nbm} = 0;
	$res->{nbl} = 0;
	
 	return encode_json $res;
 }

sub getDejaVu_project
{
	my ($cnvHash) = @_;
	return "";
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
	my ($cnv,$mother,$father) = @_;
	

	return "?" unless $mother && $father;
		
	my $trm;
	$trm = "?";
	if ($mother){
		$trm = "-";
		$trm = "+"  if (exists $cnv->{patients}->{$mother->id});
	}
	my $trf;
	$trf = "?";
	if ($father){
		$trf = "-";
		$trf = "+"  if (exists $cnv->{patients}->{$father->id});
	}
	
	if ($trf eq "?" && $trm eq "?") {
		return "?";
	}
	if ($trf eq "?" && $trm eq "m") {
		return "m?";
	}
	if ($trf eq "+" && $trm eq "?") {
		return "f?";
	}
	if ($trf eq "+" && $trm eq "+") {
		return "fm";
	}
	if ($trf eq "+" && $trm eq "-") {
		return "f";
	}
	if ($trf eq "-" && $trm eq "+") {
		return "m";
	}
	if ($trf eq "-" && $trm eq "-") {
		return "d";
	}
	return "?";
}



sub getScoreCNV
{
	my ($listgenes,$score_caller,$transmission,$dupseg) = @_;
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
		warn $score;
		$score -- unless @$listgenes;
		
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

sub test_type {
	my($hash,$flag) = @_;
	$flag = "caller_".$flag;
	confess() unless $caller_type_flag->{$flag};
	return $hash->{caller_type_flag} & $caller_type_flag->{$flag};
}



