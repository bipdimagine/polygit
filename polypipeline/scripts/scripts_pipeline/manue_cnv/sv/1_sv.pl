#!/usr/bin/perl

use Carp;
use strict;
use Data::Dumper;
use CGI qw/:standard :html3/;
use Set::IntSpan::Fast::XS;
use Set::IntervalTree;
use List::Util qw[min max  sum];
use FindBin qw($Bin);
use Storable qw(store retrieve freeze);
use Storable 'dclone';
use lib "$Bin/../../../../../GenBo/lib/obj-nodb/";
use Hash::Merge qw/ merge /;
use MCE::Loop;
use Getopt::Long;
use GBuffer;
use GenBoProject;
use Text::CSV;
#require  "$Bin/../../../../GenBo/lib/obj-nodb/GBuffer.pm";
#require "$Bin/../../../../GenBo/lib/obj-nodb/GenBoProject.pm";
require  "$Bin/../SVParser.pm";
use lib "$Bin/../../dejavu/utility/";
use liftOverRegions;

############################################################################################################
# Cherche tous les variants structuraux pour chaque patient d'un projet et tous callers confondus.
# Stocke le resultat dans une  table de hash qui est freeze : nom du fihier = project.allSV
# Pour chaque SV la localisation chr_start_end donne accès à la liste des gènes compris dans cet interval
############################################################################################################

my $projectname;
my $fork;
GetOptions(
	'project=s' => \$projectname,
#	'patient=s' => \$patient_name,
	'fork=s' => \$fork,
);

my $type ="BND";
#my $path = "/data-xfs/Manue/Test_SV/".$projectname."/manta/";
#my $path_djv = "/data-xfs/Manue/Test_SV/DejaVu/TransLoc/";



# pour lire les vcfs
my $compteur=0;
	
# pour stocker les translocations 
my @listHashRes;

# pour le dejavu
my $hdejavu;


	
###################################
#  lecture des vcf et 
#  gestion du projet via les objets Genbo
###################################

my $caller = "manta";
my $buffer = GBuffer->new();	

my $project = $buffer->newProject( -name => $projectname);


my $dir_tmp = $buffer->config_path("root","tmp");
my $filename = "$dir_tmp/".$project->name.".csv";

	my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
	
	my $fh;		
	open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
	$csv->print($fh, ["project","type","chr1_38","position1_38","chr2_38","position2_38","chr1_19","position1_19","chr2_19","position2_19","caller","patient","sr1","sr2","pr1","pr2","gq","qual"]); 
	$csv->print($fh, [0,"NONE","N",-1,"N",-1,"N",-1,"N",-1,0,0,0,0,0,0,0,0]); ; 
	


# pour enregister les fichiers allBND
my $path = $project->getSVeqDir();

# pour enregistrer le fichier  dejavu
my $path_djv = $project->DejaVuProjectsSVeq_path;

# boucle sur les patients du projets
my @listPatients = grep {$_->isGenome} @{$project->getPatients()};

my @all_sv;
my $all_transloc;
foreach my $thePatient (@listPatients)
{
	my $patientname = $thePatient->name;
	my $htransloc ={};
	foreach my $caller (@{$thePatient->callingSRMethods}){
		my ($htransloc1) = SVParser::parse_vcf_bnd($thePatient,$caller);
		$htransloc = merge $htransloc, $htransloc1;
	
	}
	push(@all_sv,values %$htransloc);
	foreach my $event_id (keys %{$htransloc}){
			$hdejavu->{$event_id}->{$patientname} = 1;
	}
	warn scalar (keys %$hdejavu)." ".$thePatient->name;
	my $file_out = $path.$patientname.".allBND.store";
	liftover_sv($thePatient,$htransloc);
}
close($fh);


my $dbh = DBI->connect("dbi:ODBC:Driver=DuckDB;Database=:memory:", "", "", { RaiseError => 1 , AutoCommit => 1});
my $dir_dejavu =  $buffer->config_path("dejavu","sv");
die("not exists dejavu dir : -".$dir_dejavu) unless -e $dir_dejavu;
my $parquet_file = $dir_dejavu.$project->name.".parquet";
my $query = "
	COPY (
        SELECT * from read_csv_auto(['$filename']) order by type
    )
    TO '$parquet_file' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE,ROW_GROUP_SIZE 1000000);";
    warn $query;
	$dbh->do($query);
	$dbh->disconnect;
# pour le dejavu du projet
exit(0);



sub liftover_sv {
	my ($patient,$hdejavu) = @_;
	my $global;
	my $variationsDir = $project->getSVeqDir();();


		foreach my $id (keys %$hdejavu){
			my ($c1,$pos1,$c2,$pos2) = split("_",$id);
			unless (exists $global->{$id}){
				$c1 = return_chr_name($c1);
				next unless $project->isChromosomeName($c1);
				$c2 = return_chr_name($c2);
				next unless $project->isChromosomeName($c2);
			$global->{$id}->{chromosome1} =$project->getChromosome($c1)->ucsc_name;
			$global->{$id}->{chromosome2} = $project->getChromosome($c2)->ucsc_name;;
			$global->{$id}->{position1} = $pos1;
			$global->{$id}->{position2} = $pos2;
			$global->{$id}->{type} = $hdejavu->{$id}->{TYPE};
			
			}
			$global->{$id}->{patient} = $patient->id;
			if (exists $hdejavu->{$id}->{INFOS}){
			$global->{$id}->{sr1} = $hdejavu->{$id}->{INFOS}->{SR}->[0]+0;
			$global->{$id}->{sr2} = $hdejavu->{$id}->{INFOS}->{SR}->[1]+0;
			$global->{$id}->{pr1} = $hdejavu->{$id}->{INFOS}->{PR}->[0]+0;
			$global->{$id}->{pr2} = $hdejavu->{$id}->{INFOS}->{PR}->[1]+0;
			$global->{$id}->{gq} = $hdejavu->{$id}->{INFOS}->{GQ}->[0]+0;
			}
			else {
			$global->{$id}->{sr1} = -1;
			$global->{$id}->{sr2} = -1;
			$global->{$id}->{pr1} = -1;
			$global->{$id}->{pr2} = -1;
			$global->{$id}->{gq} = -1;
			}
			warn  $hdejavu->{$id}->{QUAL};
			$global->{$id}->{qual} = $hdejavu->{$id}->{QUAL};
			warn $global->{$id}->{qual};
		}
	my $lift = liftOverRegions->new(project=>$project,version=>$project->lift_genome_version);
	foreach my $id (keys %$global){
		my $h = $global->{$id};
		$h->{id} = $id;
		$lift->add_region_bnd($h);
	}
	my $lifto =  $lift->liftOver_bnd($project->name);
	 save_csv_translocation($project,$lifto,$dir_tmp);
	
		
}

sub return_chr_name {
	my ($n) = @_;
	return "X" if $n eq "23";
	return "Y" if $n eq "24";
	return "MT" if $n eq "25";
	return $n;
}

sub save_csv_translocation {
	my ($project,$snps,$dir_tmp) = @_;
	
		my $cnv_callers = {
    "wisecondor"          => 1 << 0,  # 2^0 = 1
    "canvas"        => 1 << 1,  # 2^1 = 2
    "manta"        => 1 << 2,  # 2^2 = 4
};
	
	
	
	my $mt ;
	$mt = 1 if $project->getChromosome("MT")->length() > 16571 ;
	#or $project->current_genome_version eq "HG38";	
	
	foreach my $vhh (values %$snps){
			my $chr1 = $project->getChromosome($vhh->{chromosome1})->name;
			my $chr2 = $project->getChromosome($vhh->{chromosome2})->name;
			my $position1 =$vhh->{position1};
			my $position2 =$vhh->{position2};
			my $chrlift1 = "Z";
			my $chrlift2 = "Z";
			my $position1lift = -1;
			my $position2lift = -1;
			if ($vhh->{NB_LIFT} == 2 ){
				my @alift = sort{ $project->getChromosome($vhh->{chromosome1})->karyotypeId <=>  $project->getChromosome($vhh->{chromosome1})->karyotypeId } @{$vhh->{LIFT}};
				if ( $project->isChromosomeName($alift[0]->{chromosome}) && $project->isChromosomeName($alift[1]->{chromosome})){
					$chrlift1 = $project->getChromosome($alift[0]->{chromosome})->name;;
					$chrlift2 = $project->getChromosome($alift[1]->{chromosome})->name;;
					$position1lift = $alift[0]->{position};
					$position2lift = $alift[1]->{position};
				}
		}
		
		
		if ($chr1 eq "MT")  {
			if ($mt == 1){
					$position1 = $position1lift if $position1lift > -1;
			}
			else {
				$position1lift = $position1;
			}
		}
		
		if ($chr2 eq "MT" ) {
			if ($mt == 1){
					$position2 = $position2lift if $position2lift > -1;
			}
			else {
				$position2lift = $position2;
			}
		}
		
		
		my $type =  $vhh->{type};
		
		my $patient = $vhh->{patient};
		my $sr1 = $vhh->{sr1};
		my $sr2 = $vhh->{sr2};
		my $pr1 = $vhh->{pr1};
		my $pr2 = $vhh->{pr2};
		my $gq = $vhh->{gq};
		my $qual = $vhh->{qual};
		if ($project->current_genome_version eq "HG19") {
			$csv->print($fh, [$project->id,$type,$chrlift1,$position1lift,$chrlift2,$position2lift,$chr1,$position1,$chr2,$position2,$cnv_callers->{manta},$patient,$sr1,$sr2,$pr1,$pr2,$gq,$qual]); 
		}
		else {
			$csv->print($fh, [$project->id,$type,$chr1,$position1,$chr2,$position2,$chrlift1,$position1lift,$chrlift2,$position2lift,$cnv_callers->{manta},$patient,$sr1,$sr2,$pr1,$pr2,$gq,$qual]); 
		}
		
	}
	return "\'".$filename."\'";
	
}



sub save {
	my ($merged,$project) = @_;
	confess();
		my $dir = $project->getCacheDir(). "/SV/";
				
		system ("mkdir -p $dir") unless -e $dir;
		my $nodejavu = GenBoNoSqlDejaVuSV->new( dir => $dir, mode => "c" );
		foreach my $sv (@$merged){
			warn Dumper $sv;
			 $nodejavu->insert_local_sv($sv,scalar( keys %{$sv->{'PATIENTS'}}));
		}
		warn $dir;
				$nodejavu->close();
	
}
sub dejavu {
	my ($svs ,$project,$limit) = @_;
	my $no_sv = $project->dejavuSV();
	my $dir_HG38 = $buffer->config_path("dejavu")."/HG38/SVeq/";
	
	foreach my $sv (@$svs){
		my $dvlite = $no_sv->get_sv_dejavu($sv,$limit,$project->name);
		$sv->{DEJAVU}->{PATIENTS} =  $dvlite->{dv_patients};
		$sv->{DEJAVU}->{PROJECTS} =  $dvlite->{dv_projects};
		$sv->{DEJAVU}->{DETAIL} =  dclone $dvlite->{infos};
		my $identity1 = delete  $sv->{DEJAVU}->{DETAIL}->{identity1};
		my $identity2 =  delete  $sv->{DEJAVU}->{DETAIL}->{identity2};
		foreach my $p (keys %{$sv->{DEJAVU}->{DETAIL}}){
			next unless $sv->{DEJAVU}->{DETAIL}->{$p};
			foreach my $dejavu_patient (@{$sv->{DEJAVU}->{DETAIL}->{$p}}){
				push(@{$sv->{DEJAVU}->{ARRAY}},$p.":".$dejavu_patient.":bp1=".$identity1."_m bp2=".$identity2."_m");
			}
			#push(@{$sv->{DEJAVU}->{ARRAY}},map{$p.":".$_}@{$sv->{DEJAVU}->{DETAIL}->{$p}});
		}
		
	}
	
}


sub getGenesList {
	my ($chr,$start,$end) = @_;
			my $tabGenes = $chr->getGenesByPosition($start,$end);
	
			my $genes_liste="";
			my $omim=0;
			my $phenotypes="";
			my $genes = {};
			
			foreach my $g (@$tabGenes)
			{
					$genes->{$g->external_name} = $g->score;
					$genes_liste .= $g->external_name.":".$g->phenotypes."##";
					$omim = 1 if $g->is_omim_morbid();
			}
		return 	($genes,$omim);

}

sub getCytoband {
	my ($chr,$start,$end) = @_;
		my $hband = $chr->getCytoband($start, $end);
	
			my @tb;
			foreach my $b ( keys %{ $hband->{'name'} } ) {
				push( @tb, $b );
			}
			my @band = sort( { $a cmp $b } @tb );
			return  join( ",", @band );
} 

sub genesInfos {
	my ($svs ,$project,$limit) = @_;
	foreach my $sv (@$svs){
			my $chr1 = $project->getChromosome($sv->{CHROM1});
			my ($hgenes,$omim1) = getGenesList($chr1,$sv->{POS1_MIN}-$limit,$sv->{POS1_MAX}+$limit);
			my $BND1_ind_id = $sv->{CHROM1}."_".$sv->{POS1};
			$sv->{CYTOBAND1} =  getCytoband($chr1,$sv->{POS1_MIN}-$limit,$sv->{POS1_MAX}+$limit);
			$sv->{GENES1} = $hgenes;
			$sv->{OMIM1} = $omim1;
				my $chr2 = $project->getChromosome($sv->{CHROM2});
			
			my ($hgenes2,$omim2) = getGenesList($chr2,$sv->{POS2_MIN}-$limit,$sv->{POS2_MAX}+$limit);
			my $BND2_ind_id = $sv->{CHROM2}."_".$sv->{POS2};
			$sv->{CYTOBAND2} =  getCytoband($chr2,$sv->{POS2_MIN}-$limit,$sv->{POS2_MAX}+$limit);
			$sv->{GENES2} = $hgenes2;
			$sv->{OMIM2} = $omim2;
	}
			
}



	