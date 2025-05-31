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
use GenBoDuckDejaVuSv;

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
# pour enregister les fichiers allBND

# pour enregistrer le fichier  dejavu
# boucle sur les patients du projets
my @all_sv;
my @listPatients = grep {$_->isGenome} @{$project->getPatients()};
my $duck = GenBoDuckDejaVuSv->new( project => $project );
foreach my $patient(@listPatients){
	my $psv = $duck->get_sv_project($patient);
	 annot_bnd($psv,$project,$patient);
	
	push(@all_sv,@$psv);
}
save_parquet_rocksdb($project,\@all_sv);
#save($project,\@all_sv);

sub return_chr_name {
	my ($n) = @_;
	return "X" if $n eq "23";
	return "Y" if $n eq "24";
	return "MT" if $n eq "25";
	return $n;
}


sub annot_bnd {
	my ($aSV,$project,$patient) = @_;
	my $distance = 100;
	for my $sv (@$aSV) {
		$sv->{patient} = $patient->id;
		$sv->{id} = $patient->id."!".$sv->{type}."_".$sv->{chrom1}."_".$sv->{pos1}."_".$sv->{chrom2}."_".$sv->{pos2};
		genesInfos($sv,$project,$distance);
		
		dejavu($sv,$project,$distance);
		getScoreEvent($sv,$patient);
	}
	return $aSV
}


sub dejavu {
	my ($sv ,$project,$limit) = @_;
	
		my $dvlite = $duck->get_dejavu($sv->{chrom1},$sv->{pos1},$sv->{chrom2},$sv->{pos2},$limit);
		# if $dvlite->{nb_patients} <10;
		$sv->{dejavu}->{nb_patients} =  $dvlite->{nb_patients};
		$sv->{dejavu}->{nb_projects} =  $dvlite->{nb_projects};
		$sv->{dejavu}->{this_project} =  $dvlite->{this_project};
		
}


sub getGenesList {
	my ($chr,$start,$end,$bp) = @_;
			my $tabGenes = $chr->getGenesByPosition($start,$end);
	
			my $genes_liste="";
			my $omim=0;
			my $phenotypes="";
			my $agenes = [];
			my $max =-1;
			foreach my $g (sort {$b->score <=> $a->score} @$tabGenes){
									my $h;
									$h->{name} = $g->external_name;
									$h->{id} = $g->id;
									$h->{score} = $g->score;
									$h->{bp} = $bp;
									$h->{omim} ++ if $g->is_omim_morbid();
									$h->{omim_id} = $g->omim->{omim_id} if $g->is_omim_morbid();
									$max = $g->score unless $max;
									$omim++;
									$h->{phenotypes} = $g->array_phenotypes;
									my $study = $g->project->getStudy();
									foreach my $p (@{$g->getPanels()}){
									if ($study eq "glucogen") {
										$h->{score} = 8 if ($p->name =~ /glucogen/i);
									}
									}
									push(@{$agenes},$h);	
			}
		my @t= sort {$b->{score} <=> $a->{score}} @{$agenes};
			
		return 	(\@t,$omim);

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
	my ($sv ,$project,$limit) = @_;
			my $chr1 = $project->getChromosome($sv->{chrom1});
			
			if ($sv->{type} eq "TL"){
			my ($hgenes,$omim1) = getGenesList($chr1,$sv->{pos1}-$limit,$sv->{pos1}+$limit,1);
			
			$sv->{cytoband1} =  getCytoband($chr1,$sv->{pos1}-$limit,$sv->{pos1}+$limit);
			$sv->{omim1} = $omim1;
			my $chr2 = $project->getChromosome($sv->{chrom2});
			
			my ($hgenes2,$omim2) = getGenesList($chr2,$sv->{pos2}-$limit,$sv->{pos2}+$limit,2);

			my @toto = sort {$b->{score} <=> $a->{score}}(@$hgenes,@$hgenes2);
			$sv->{genes} = \@toto;
			my $BND2_ind_id = $sv->{chrom2}."_".$sv->{pos2};
			$sv->{cytoband2} =  getCytoband($chr2,$sv->{pos2}-$limit,$sv->{pos2}+$limit);
			$sv->{omim2} = $omim2;
			$sv->{omim} = $omim2 + $omim1;
			}
			elsif ($sv->{type} eq "INV") {
				my ($hgenes,$omim1) = getGenesList($chr1,$sv->{pos1},$sv->{pos2},0);
				$sv->{genes} = $hgenes;
				$sv->{omim} = $omim1;
				$sv->{cytoband1} =  getCytoband($chr1,$sv->{pos1},$sv->{pos1}+2);
				$sv->{cytoband2} =  getCytoband($chr1,$sv->{pos2},$sv->{pos2}+2);
			}
			else {
				die();
			}
			$sv->{max_score_gene}  = max( map{$_->{score}} @{$sv->{genes}});
			$sv->{max_score_gene} =0 unless $sv->{max_score_gene};
			
}



sub save_parquet_rocksdb {
	my ($project,$svs) = @_;
	my $dir_tmp = "/data-beegfs/tmp/";
	my $filename = "$dir_tmp".$project->name.".csv";
	my $fh;		
	open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
	my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
	my $header = ["id","type","chr1","pos1","chr2","pos2","patient","nb_dejavu_patients","nb_dejavu_projects","genes"];
	my $tab_pid;
	my $col_empty =["Z","Z","Z",0,"z",0,0,0,0,"Z"];
	
	
	my $cnv_patients; 
	my $dir = $project->getCacheSV(). "/rocks/";
	my $rocks = GenBoNoSqlRocks->new(dir=>"$dir",mode=>"c",name=>"sv");
	
	$csv->print($fh, $header); 
	$csv->print($fh, $col_empty); 
	
	foreach my $sv (@$svs) {
		my $line;
		$rocks->put($sv->{id},$sv);
		push(@$line,$sv->{id});
		push(@$line,$sv->{type});  
		push(@$line,$sv->{chrom1});
		push(@$line,$sv->{pos1});
		push(@$line,$sv->{chrom2});
		push(@$line,$sv->{pos2});
		push(@$line,$sv->{patient});
		push(@$line,$sv->{dejavu}->{nb_patients});
		push(@$line,$sv->{dejavu}->{nb_projects});
		my $ts;
		foreach my $gene (@{$sv->{genes}}){
			$ts->{$gene->{name}} ++;
		}
		
		my $st_genes = join(",", keys %{$ts});
		$st_genes= "" unless $st_genes;
		push(@$line,$st_genes);
		
		$csv->print($fh, $line); 
	}
	$rocks->close;
	close($fh);
	my $dbh = DBI->connect("dbi:ODBC:Driver=DuckDB;Database=:memory:", "", "", { RaiseError => 1 , AutoCommit => 1});
	my $parquet_file = $project->getCacheSV()."/".$project->name.".".$project->id.".parquet";
	my $query = "
	COPY (
        SELECT * from read_csv_auto(['$filename']) order by type,chr1,pos1
    )
    TO '$parquet_file' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE);";
    warn $query;
	$dbh->do($query);
	$dbh->disconnect;
	warn $parquet_file;
	warn $project->getCacheSV();
	
}

sub getScoreEvent
 {
	my ($sv,$patient) = @_;
	
	my $score = 0;
	
	# break point dans un gene omim morbid
	$score += 1  if $sv->{"max_score_gene"} >= 5 ;	
	$score += 1  if $sv->{"max_score_gene"} >= 4 ;
	$score += 0.5  if $sv->{"max_score_gene"} >= 3 ;	
	# score sur le dejavu
	my $nbdejavu = $sv->{dejavu}->{nb_patients};#$sv->{"nbdejavu"};
	
	if ($nbdejavu <= 10)
	{
		$score += 10 - $nbdejavu;
	}
	# score sur la QUAL
	$score += int($sv->{qual} /100)/10;
	my $dsr = ($sv->{sr2}/($sv->{sr2}+$sv->{sr1}+1));
	my $dpr = ($sv->{pr2}/($sv->{pr2}+$sv->{pr1}+1));
	$score += $dsr;
	$score += $dpr;
	$score = sprintf("%.1f", $score);
	#pour présenter en premier les CNV denovo
#	if ( ($sv->{PATIENTS}->{$patient->id}->{'TRANSMISSION'}  =~ m/denovo/) && !($sv->{PATIENTS}->{$patient->id}->{'TRANSMISSION'} =~ m/mother/ ) && !($sv->{PATIENTS}->{$patient->id}->{'TRANSMISSION'}  =~ m/father/ ) && !($sv->{PATIENTS}->{$patient->id}  =~ m/both/) )
#	{
#			$score += 1.5; 
#	}
$sv->{score} = $score;
	return $score;
 }

#sub gatherSV_byPosition_byPatient
#{
#	my ($aSV,$project,$patient) = @_;
#
#
#	my $all;		
#	my $nbc = 0 ;
#	my $distance = 10000;
#	my @clusters;
#warn "start cluster";
#for my $sv (@$aSV) {
#	my $find;
#	warn Dumper $sv;
#    foreach my $cluster (@clusters) {
#    	next unless $cluster->{chrom1} ne $sv->{CHROM1};
#    	next unless $cluster->{chrom2} ne $sv->{CHROM2};
#        foreach my $member (@{$cluster->{sv}}) {
#            if (
#                $sv->{CHROM1}  eq $member->{CHROM1}  &&
#                $sv->{CHROM2} eq $member->{CHROM2} &&
#                abs($sv->{POS1}  - $member->{POS1})  <= $distance &&
#                abs($sv->{POS2} - $member->{POS2}) <= $distance
#            ) {
#   				push(@{$cluster->{sv}},$sv);
#   				die();
#                $find =1;
#                last ;
#            }
#           
#        }
#         last if $find;
#    }
#    next  if $find;
#    # Si aucun cluster trouvé ? nouveau
#    my $hsv;
#    $hsv->{sv} = [$sv];
#   $hsv->{chrom1} = $sv->{CHROM1};
#   $hsv->{chrom2} = $sv->{CHROM2};
#    $hsv->{chrom2} = $sv->{CHROM2};
#   
#    push @clusters, $hsv ;
#}
#warn "end cluster";
#my $merged_sv;
# foreach my $cluster (@clusters) {
# 	my $SV;
# 	my @p1 = map {$_->{POS1}} @{$cluster->{sv}};
# 	$SV->{pos1} = int(sum(@p1)/scalar(  @{$cluster->{sv}}));
# 	my @p2 = map {$_->{POS2}}  @{$cluster->{sv}};
# 	$SV->{pos2} = int(sum(@p2)/scalar(  @{$cluster->{sv}}));
# 	$SV->{sr1} =  max(map {$_->{sr1}}  @{$cluster->{sv}});
# 	$SV->{sr2} =   max(map {$_->{sr2}}  @{$cluster->{sv}});
# 	$SV->{pr1} = max(map {$_->{pr1}}  @{$cluster->{sv}});
# 	$SV->{pr2} = max(map {$_->{pr2}}  @{$cluster->{sv}});
# 	my $tsv = $cluster->{sv}->[0];
#	$SV->{id} = $patient->id."!".$SV->{chr2}."_".$SV->{pos2}."_".$SV->{chr2}."_".$SV->{pos1};
#		
# 	push(@{$merged_sv},$SV);
# }
# 	warn $project->name;
#	genesInfos($merged_sv,$project,$distance);
#	dejavu($merged_sv,$project,$distance);
#	#save($all_merged_cnv,$project);
#		
#	return 	$merged_sv;			
#}
	