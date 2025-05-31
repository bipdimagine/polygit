#!/usr/bin/perl
use FindBin qw($Bin);
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
use Clone qw(clone);
use Parallel::ForkManager;
use lib "$Bin/../../../../../GenBo/lib/obj-nodb/";
use Hash::Merge qw/ merge /;
use MCE::Loop;
use Getopt::Long;
use Text::CSV;
require  "$Bin/../SVParser.pm";
use lib "$Bin/../../dejavu/utility/";
require  "$Bin/../../../../../GenBo/lib/obj-nodb/GBuffer.pm";
require  "$Bin/../../../../../GenBo/lib/obj-nodb/GenBoDuckDejaVuCnv.pm";
#use GBuffer;



my $caller_type_flag = {
	"caller_sr" => 1,
	"caller_depth" => 2,
	"caller_coverage" => 4,
};
my $type_by_caller = {
	"wisecondor" => 4,
	"manta" => 1,
	"pbsv" => 1,
	"dragen-sv" =>1,
	"hificnv" =>2,
	"canvas" =>2,
	"dragen-cnv" =>2,
	"cnvnator" =>2,
};

my $cnv_callers = {
    "wisecondor"          => 1 << 0,  # 2^0 = 1
    "canvas"        => 1 << 1,  # 2^1 = 2
    "manta"        => 1 << 2,  # 2^2 = 4
    "pbsv"        => 1 << 3,  
    "dragen-sv"        => 1 << 4,  
    "hificnv"        => 1 << 5,  
    "dragen-cnv"        => 1 << 6, 
     "cnvnator"        => 1 << 7, 
};


###################################################################
# Cherche tous les variants structuraux d'un projet pour chaque patient et tous callers confondus.
# Reconstruit les CNV fragmentes (même caller et bornes a moins de 10% de la longeur)
# Construit une table de hash par patient et freeze ces tables  : nom du fihier = patient.allSV
# Pour le dejavu construit egalement une table de hash qui conserve tous les CNV du projet 
# et garde pour chaque CNV l'info du patient et du caller qui l'a détecté.
# Pour chaque CNV la localisation chr_start_end donne accès à la liste des gènes compris dans cet interval
####################################################################

my $fork = 5;
my $dir_tmp = "/data-beegfs/tmp/";

my $limit;
my $projectname;


GetOptions(
	'project=s' => \$projectname,
);

my $col = ["project","patient","chr38","start38","end38","chr19","start19","end19","type","callers","caller_type_flag",'sr1','sr2','sr_qual','pr1','pr2',"depth_qual",'depth_CN','coverage_zscore','coverage_ratio','elementary'];
 #= $cgi->param('projectname');
#my $patient_name = $cgi->param('patient');
#$fork = $cgi->param('fork');

#$fork=1;

# pour récupérer les objets project et patient
my $buffer = GBuffer->new();

my $project = $buffer->newProjectCache( -name => $projectname);

my $duck = GenBoDuckDejaVuCNV->new( project => $project );

my $raw_cnv; 
	my $hp;
	foreach my $p (@{$project->getPatients}){
		my $hash = $duck->get_cnvs_by_project($p);
		
		push(@$raw_cnv,@$hash);
	}
	
	my $cytodir = $project->get_cytology_directory;	
	my $fichier_dupseg = $cytodir."/segmental_duplication.bed";
	my $fdds;
	my $hdupseg;
	open($fdds,$fichier_dupseg) or die("open: $!");
	
	# lecture ligne par ligne
	while( defined( my $ligne = <$fdds> ) )
	{
		next unless ($ligne =~ m/^chr/ );
		my @champs = split(/\t/,$ligne);
		
		my $chr1 = $champs[0];
		next unless $project->isChromosomeName($chr1);
		my $chr = $project->getChromosome($chr1)->name;
		my $start = $champs[1];
		my $end = $champs[2];
		
		my $idLocation = $chr.":".$start."_".$end;
		# version avec IntervalTree
		if (!exists $hdupseg->{$chr})
		{
		 	my $DupSegTree = Set::IntervalTree->new;
		 	$hdupseg->{$chr}=$DupSegTree;
		}		
		$hdupseg->{$chr}->insert($idLocation,$start,$end);		
	}
	
	my $gather = gatherSV_by_Interval($raw_cnv);
	my $file2 = "/data-beegfs/tmp/toto.store";
	store($gather, $file2) or die "Can't store $file2 for ".$file2."!\n";
	save_parquet_rocksdb($gather);
	exit(0);
#my $all_cnv;
#
#	my $hash_patient = read_duckdb(); 
#	my $gather = gatherSV_by_Interval($hash_patient,$p);
#
#
#close($fh);

#



sub gatherSV_by_Interval  
{
	my ($aCNV) = @_;
	 my $hGroupedCNV;
	 
	 my $seuilSameEvent = 0.9;
	 
	 my $hSVPos;

	my $all;		
				#pour avoir un objet chromosome
			
				
				# 1)  detecter les SV identiques 
				my $htree;
				# remplir les arbres :  regrouper les SV chevauchants
				my $interval;
	
				#my $nodejavu = GenBoNoSqlDejaVuCNV->new( dir => $dir, mode => "c" );
				my $hCNV;
				foreach my $cnv  (@$aCNV)
				{
						$hCNV->{$cnv->{uid}} = $cnv;
						push(@{$interval->{$cnv->{type}}->{$cnv->{chromosome}}},[$cnv->{start},$cnv->{end},$cnv->{uid}]);
				
				}
				
				foreach my $type (keys %$interval){
					foreach my $chr_name (keys %{$interval->{$type}}){
						my $chr = $project->getChromosome($chr_name);
						my $merged = merge_intervals($interval->{$type}->{$chr_name},0.75);
						my $complete_merge = merge_hash($type,$chr,$merged,$hCNV);
						
						foreach my $cnv (@$complete_merge){
							push(@{$all},$cnv);
						
					}
						
					}
					
					
					
				}
				annnot_sv($all);
				#$nodejavu->close();
				return $all;
}

warn "ok";
sub merge_hash {
	my ($type,$chr,$merged,$hcnv) = @_;
	 my $total;
	 warn "merge";
	foreach my  $interval (@$merged){
			my @cnvs;
			foreach my $id (split(";",$interval->[2])) {
					push(@cnvs,$hcnv->{$id});
			}
			
			my $real_start;
			my $real_end;
			my @mids = sort {$a->{length} <=> $b->{length} } grep {test_type($_->{caller_type_flag},"caller_sr")} @cnvs;
			my @cids = sort {$a->{start} <=> $b->{start} } grep {test_type($_->{caller_type_flag},"caller_depth")} @cnvs;
			
			 if (@mids) {
			 	$real_start = $mids[-1] ->{start};
			 	$real_end = $mids[-1] ->{end};
			 }
			 elsif (@cids) {
			 	$real_start = $cids[0]->{start};
			 	$real_end = $cids[-1]->{end};
			 	
			 }
			 else {
			 	$real_start = $interval->[0];
			 	$real_end = $interval->[1];
			 }
			
				my $global_id = $type . "_" . $chr->name . "_" . $interval->[0] . "_" . $interval->[1];
			 	my $hfinal;
			 	
			 	
				$hfinal->{id} = $global_id;			
				$hfinal->{type} = $type;
				$hfinal->{karyotypeId} = $chr->karyotypeId;
				$hfinal->{chromosome} = $chr->name;
				$hfinal->{start} =  $real_start;
				$hfinal->{end} =  $real_end;
#				$hfinal->{'LOCUS'} =  format_number($interval->[0])."-".format_number($interval->[1]);
				$hfinal->{len}  = abs($real_end- $real_start)+1;
				
				
				foreach my $cnv (@cnvs){
					my $patient_id = $cnv->{patient};
					$hfinal->{list_patients}->{$patient_id} ++;
					push(@{$hfinal->{patients}->{$patient_id}->{cnv_origin}},$cnv);
					$hfinal->{patients}->{$patient_id}->{callers} =  $hfinal->{patients}->{$patient_id}->{callers} | $cnv->{callers};
					$hfinal->{patients}->{$patient_id}->{caller_type_flag} =  $hfinal->{patients}->{$patient_id}->{caller_type_flag} | $cnv->{caller_type_flag};
				}
				
				#getScoreCallers($hfinal);
				warn "-----";
				dejavu($hfinal);
				genesInfos($hfinal);
				warn "+++++";
				getDupSeg($hfinal);
				
				push(@$total,$hfinal);
	}
	
	warn "ANNOT";
return $total;	
}

exit(0);

sub dejavu {
	my ($cnv) = @_;
	my $nbproject=0;
	my $nbpatient=0;
	my $nbDJV_Wisecondor=0;
	my $nbDJV_Canvas=0;
	my $nbDJV_Manta=0;
	my $nb =0;
	
	my $scorecaller_evt=0;
	my $list_of_other_patient;
	$cnv->{dejavu} = $duck->dejavu($cnv->{type},$cnv->{chromosome},$cnv->{start},$cnv->{end},90);
	
}

 sub genesInfos
{
	my ($cnv)  =@_;
		
			my $objChr = $project->getChromosome($cnv->{chromosome});
						
						my $tabGenes = $objChr->getGenesByPosition($cnv->{start},$cnv->{end});
						
						my $omim = 0;
						
						


						# calculer le score gene max correspondant a la liste de gene associe au variant
						$cnv->{genes} = [];
						$cnv->{omim} = 0;
					
							my @names;
							my $max;
							my $genes;
							foreach my $g (sort {$b->score <=> $a->score} @$tabGenes){
									my $h;
									$h->{name} = $g->external_name;
									$h->{id} = $g->name;
									$h->{score} = $g->score;
									$h->{omim} ++ if $g->is_omim_morbid();;
									$max = $g->score unless $max;
									$cnv->{omim} ++ if $g->is_omim_morbid();
									push(@{$cnv->{genes}},$h);	
							}
						# pour les cytobandes
						my @tb;
						my @band;
						
						my $hband = $objChr->getCytoband($cnv->{start},$cnv->{end});

						foreach my $b ( keys %{ $hband->{'name'} } ) {
								push( @tb, $b );
						}
						@band = sort( { $a cmp $b } @tb );
						$cnv->{cytoband} = join( ",", @band );
						
						# pour l'acces a DGV
						my $url_DGV = getDGV_url( $objChr->ucsc_name, $cnv->{start},$cnv->{end},$project);
						$cnv->{dgv} = $url_DGV;
	
	
}




sub getDGV_url
{
	my ($chr,$d,$f) = @_;
	my $genome = "hg19";
	$genome = "hg38" if $project->annotation_genome_version =~/38/; 
	my $url = "http://dgv.tcag.ca/gb2/gbrowse/dgv2_$genome/?name=".$chr."%3A".$d."-".$f.";search=Search";
	return $url;
}

	
	
sub test_type {
	my($value,$flag) = @_;
	confess() unless $caller_type_flag->{$flag};
	return $value & $caller_type_flag->{$flag};
}



#sub getScoreCallers {
#		my ($gcnv) = @_;
#		my $score = 0;
#		my $sc_total = 0;
#			foreach my $patient (keys %{$gcnv->{patients}}){
#				my $iscore =0;
#				foreach my $icnv (@{$gcnv->{patients}->{$patient}->{cnv_origin}}){
#					
#				if (test_type($icnv->{caller_type_flag},"caller_coverage")){
#					$gcnv->{patients}->{$patient}->{score_callers}->{nb} ++ unless exists 	$gcnv->{patients}->{$patient}->{score_callers}->{caller_coverage};
#					next if $gcnv->{patients}->{$patient}->{score_callers}->{caller_coverage} == 1;
#					$gcnv->{patients}->{$patient}->{score_callers}->{caller_coverage} = 0.75;
#					$gcnv->{patients}->{$patient}->{score_callers}->{caller_coverage} += 0.25 if abs($gcnv->{patients}->{$patient}->{coverage_zscore}) > 40;
#				}
#				if (test_type($icnv->{caller_type_flag},"caller_depth")){
#					$gcnv->{patients}->{$patient}->{score_callers}->{nb} ++  unless exists 	$gcnv->{patients}->{$patient}->{score_callers}->{caller_depth};
#					next if $gcnv->{patients}->{$patient}->{score_callers}->{caller_depth} == 1;
#					$gcnv->{patients}->{$patient}->{score_callers}->{caller_depth} = 0.75;
#					$gcnv->{patients}->{$patient}->{score_callers}->{caller_depth} += 0.25 if abs($gcnv->{patients}->{$patient}->{depth_qual}) > 30;
#				}
#				if (test_type($icnv->{caller_type_flag},"caller_sr")){
#					$gcnv->{patients}->{$patient}->{score_callers}->{nb} ++unless exists 	$gcnv->{patients}->{$patient}->{score_callers}->{caller_sr};
#					next if $gcnv->{patients}->{$patient}->{score_callers}->{caller_sr} == 1;
#					$gcnv->{patients}->{$patient}->{score_callers}->{caller_sr} = 0.75;
#					my $limit = 400;
#					if ($icnv->{callers} & $cnv_callers->{pbsv}){
#						$limit = 30;
#					}
#					$gcnv->{patients}->{$patient}->{score_callers}->{caller_sr} += 0.25 if abs($gcnv->{patients}->{$patient}->{sr_qual}) > $limit;
#				}
#				}
#				#$gcnv->{patients}->{$patient}->{score_callers}->{total} = $iscore;
#			}
#}






sub getIdentityBetweenCNV {
	my (  $start1, $end1, $start2, $end2) = @_;
	
	#retourne le recouvrement en % de la longueur du plus long des deux evenements 
	
	my $overlap = min( $end1, $end2 ) - max( $start1, $start2 );
	confess if abs( $start1 - $end1 ) ==0;
	my $overlap1 = $overlap / abs( $start1 - $end1 );
	my $overlap2 = $overlap / abs( $start2 - $end2 );
	
	return min($overlap1*100,$overlap2*100);
}



sub merge_intervals {
    my ($intervals,$limit) = @_;
    return () unless @$intervals;

    # Trier les intervalles par leur point de début
   my @sorted_intervals = sort { $a->[0] <=> $b->[0] } @$intervals;

    my @merged;
    my $current = shift @sorted_intervals;

    foreach my $interval (@sorted_intervals) {
        # Vérifier le chevauchement de 70% ou plus
        my $xc = getIdentityBetweenCNV($current->[0], $current->[1], $interval->[0], $interval->[1]);
        
        if ($xc >= 70 ) {
            # Fusionner les intervalles
            $current->[1] = $interval->[1] if $interval->[1] > $current->[1];
            $current->[2] .= ";".$interval->[2];
        } else {
            push @merged, $current;
            $current = $interval;
        }
    }

    push @merged, $current;
    return \@merged;
}

	
	



sub getDupSeg
{
	my ($cnv) = @_;
	my $chr = $cnv->{chromosome};
	my $d = $cnv->{start};
	my $f = $cnv->{end};
	my $result;
	
	# recherche des dupseg chevauchant le fragment
	# exemple eval pour chrMT qui ne fetch pas	
	 $result = $hdupseg->{$chr}->fetch($d,$f); 
#	if($@) {
#		warn "\n\n\nERROR: can't fetch $chr... NEXT...\n\n\n";
#		$cnv->{'DUPSEG'} = "-";
#		die();
#		return '-';
#	}
	my $nb_dupseg = scalar(@$result);
	
	# si il n y en a pas
	return "-"	 if ( $nb_dupseg == 0);
	
	# sinon on creer un intspan
	my $span = Set::IntSpan::Fast::XS->new();
	
	foreach my $res (@$result)
	{
		 my ($id,$coord) = split(/:/,$res);
		 my ($ds_start, $ds_end) = split(/_/,$coord); 
		 
		 # on ne conserve que la partie comprise entre $d et $f
		 $ds_start = $d if ($ds_start <= $d);
		 $ds_end = $f if ($ds_end >= $f);
		 
	    $span->add_range($ds_start,$ds_end);
	}
	
	# calcul de la couverture globale des dupseg
	my $globale_cov=0;
	my $list_pos = $span->as_string();
	my @tabpos = split(/,/,$list_pos);

	foreach my $coords (@tabpos)
	{
		my ($debut,$fin) =split(/-/,$coords);
		my $taille = ($fin-$debut);
		$globale_cov = $globale_cov + $taille;
	}
	
	$cnv->{'DUPSEG'} = $globale_cov/($f-$d)*100;
}


 sub annnot_sv {
	my ($merged) = @_;
	
	my $dir_tmp= $project->getCallingPipelineDir("annot_sv1");
	my $file_bed = $dir_tmp."/".$project->name.".bed";
	my $load = qq{module load tcltk/8.6.9};
	my $annotsv = qq{/software/distrib/AnnotSV_2.0/bin/AnnotSV};
	unlink $file_bed if -e $file_bed;
	my $file_tmp = $dir_tmp."/".$project->name.".annotsv.out.tsv";
	unlink $file_tmp if -e $file_tmp;
	my $hh;
	open(BED,">",$file_bed) or die();
	foreach my  $interval (@$merged){
		my $atype = "duplication";
		$atype = "deletion" if $interval->{type} =~  /DEL/;
		my $c = $interval->{chromosome};
		
		$hh->{$c."_".$interval->{start}."_".$interval->{end}."_".$atype} = $interval;
		print BED $c."\t".$interval->{start}."\t".$interval->{end}."\t"."$atype"."\n";
	}
	close BED;
	my $opt = "-svtBEDcol 4";
	my $opt2 = "-genomeBuild GRCh37";
	$opt2 = "-genomeBuild GRCh38" if $project->genome_version_generic() =~/HG38/;
	my $cmd = qq{ export ANNOTSV=/software/distrib/AnnotSV_2.0 && $load && $annotsv $opt2 -SVinputFile $file_bed -outputFile $file_tmp $opt 2>/dev/null };
	system($cmd);
	open(ANNOT,$file_tmp) or die("*************  open: $!  $file_tmp");
	
			# lecture de la première ligne
		my $ligne = <ANNOT>;
		chomp($ligne);
		my @champs = split(/\t/,$ligne);
		my $hannotsv;
		while (my $line = <ANNOT>){
			chomp($line);
			my @data =  split(/\t/,$line);
			my $hash;
			for  (my $i=0;$i<@champs;$i++ ){
				chomp($data[$i]);
				$hannotsv->{$champs[$i]} = $data[$i];
			}
		my $sid = $hannotsv->{'AnnotSV ID'};
		my $cnv = $hh->{$sid};
		die() unless exists $hh->{$sid};
		
		 my $dgv_gain_freq = $hannotsv->{'DGV_GAIN_Frequency'};
		my $dgv_loss_freq =  $hannotsv->{'DGV_LOSS_Frequency'};
		
	
		$dgv_gain_freq = -1 if ($dgv_gain_freq eq "-");
		$dgv_loss_freq = -1 if ($dgv_loss_freq eq "-");
		
		#SV_rank
		my $SV_rank = $hannotsv->{'AnnotSV ranking'};	

		#OMIM
		my $OMIN_MG = $hannotsv->{'morbidGenes'};	
		
		#dbVar
		my $dbVar_status =  $hannotsv->{'dbVar_status'};
		$cnv->{'GOLD_G_FREQ'}= $dgv_gain_freq;
		$cnv->{'GOLD_L_FREQ'}= $dgv_loss_freq;
		$cnv->{'OMIN_MG'}= $OMIN_MG;
		$cnv->{'dbVar_status'}= $dbVar_status;
		$cnv->{'RANKAnnot'}= $SV_rank;
		}
		close(ANNOT);
	unlink $file_bed;
	unlink $file_tmp;
}
sub save_parquet_rocksdb {
	my ($cnvs) = @_;
	my $dir_tmp = "/data-beegfs/tmp/";
	my $filename = "$dir_tmp".$project->name.".csv";
	my $fh;		
	open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
	my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
	my $col = ["id","type","chr","start","end","len", "nb_dejavu_patients","nb_dejavu_projects","genes"];
	my $tab_pid;
	my $col_empty =["Z","Z","Z",0,0,0,0,0,"Z"];
	
	
	my $cnv_patients; 
	my $dir = $project->getCacheCNV(). "/rocks/";
	 my $rocks = GenBoNoSqlRocks->new(dir=>"$dir",mode=>"c",name=>"cnv");


	my @sorted_ids =  map {$_->id} sort {$a->id <=> $b->id} @{$project->getPatients};
	
	
	foreach my $p (@sorted_ids){
		push(@$tab_pid,"p".$p);
		push(@$col_empty,0);
	}
	push(@$col,@$tab_pid);
	$csv->print($fh, $col); 
	$csv->print($fh, $col_empty); 
	
	foreach my $cnv (@$cnvs) {
		my $line;
		$rocks->put($cnv->{id},$cnv);
		push(@$line,$cnv->{id});
		push(@$line,$cnv->{type});  
		push(@$line,$cnv->{chromosome});
		push(@$line,$cnv->{start});
		push(@$line,$cnv->{end});
		push(@$line,abs($cnv->{end}-$cnv->{start})+1);
		push(@$line,$cnv->{dejavu}->{nb_patients});
		push(@$line,$cnv->{dejavu}->{nb_projects});
		my $st_genes = join(",",map {$_->{name}} @{$cnv->{genes}});
		$st_genes= "" unless $st_genes;
		push(@$line,$st_genes);
		foreach my $pid (@sorted_ids) {
			if (exists $cnv->{patients}->{$pid}){
				push(@$line,1);
			}
			else {
				push(@$line,0);
			}
		}
		$csv->print($fh, $line); 
	}
	$rocks->close;
	close($fh);
	my $dbh = DBI->connect("dbi:ODBC:Driver=DuckDB;Database=:memory:", "", "", { RaiseError => 1 , AutoCommit => 1});
	my $parquet_file = $project->getCacheCNV()."/".$project->name.".".$project->id.".parquet";
	my $query = "
	COPY (
        SELECT * from read_csv_auto(['$filename']) order by chr,type,start
    )
    TO '$parquet_file' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE);";
    warn $query;
	$dbh->do($query);
	$dbh->disconnect;
	warn $filename;
	warn $parquet_file;
	warn $project->getCacheDir();
	
}
