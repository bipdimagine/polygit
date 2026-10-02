#!/usr/bin/env perl
use FindBin qw($Bin);
use Carp;
use strict;
use JSON;
use Data::Dumper;
use CGI qw/:standard :html3/;
use Set::IntSpan::Fast;
use Set::IntervalTree;
use List::Util qw[min max];
use FindBin qw($Bin);
use Storable qw(store retrieve freeze);
use Clone qw(clone);
use Parallel::ForkManager;
use strict;
use Text::CSV qw( csv );
use lib "$Bin/../../../../../GenBo/lib/obj-nodb/";
use Getopt::Long;
use lib "$Bin";
use Set::IntSpan::Fast;
require  "$Bin/../SVParser.pm";
require  "$Bin/../parser/parse_pbsv.pm";
require  "$Bin/../parser/parse_hificnv.pm";
require  "$Bin/../parser/parse_sniffles2.pm";
require  "$Bin/../parser/parse_wisecondor.pm";
use lib "$Bin/../../dejavu/utility/";
use liftOverRegions;
use GBuffer;
use Text::CSV;
use Bio::DB::HTS::Tabix;
use MCE::Loop;
#
#je vais filtrer les cnvs en utilsant wisecondor le ratio 
#

my $fork = 5;
my $cgi = new CGI;


my $limit;
my $project_name;
my $patient_name;
my $fork =1 ;

GetOptions(
	'project=s' => \$project_name,
	'patient=s' => \$patient_name,
	'fork=s' => \$fork,
);
my $caller_type_flag = {
	"caller_sr" => 1,
	"caller_depth" => 2,
	"caller_coverage" => 4,
};
my $buffer = GBuffer->new();	
my $project = $buffer->newProject( -name => $project_name);
my $dir = $project->getCacheCNV(). "/rocks/";
my $rocks = GenBoNoSqlRocks->new(dir=>"$dir",mode=>"r",name=>"cnv");
my $parquet_file_quality = $project->getCacheCNV()."/".$project->name.".".$project->id.".cnv_quality.parquet";
my $parquet_file = $project->getCacheCNV()."/".$project->name.".".$project->id.".parquet";

my $dir_tmp = $buffer->config_path("tmp");
my $filename = "$dir_tmp".$project->name.time.".quality.csv";

	#my $col = ["id","patient_id","zscore","ratio","na","nr","pna","len","nb_caller","caller_type","dejavu"];
	my $col = ["id","type","chr","start","end","patient","len","nb_dejavu_patients","nb_dejavu_projects","genes","mz","mr","na","nr","pa","nb_caller","caller_sr","caller_depth","caller_coverage","blacklist","perc_freq_allelic"];


	my $tab_pid;
	
	my @tab;
	my $error =0;
	MCE::Loop->init(
   max_workers => '8', chunk_size => '1',
    gather => sub {
        my ($mce, $data) = @_;
        $error ++ if (exists $data->{error});
        
       
       foreach my $row  (@{$data->{cnvs}}) {
       	 my @t;
		foreach my $k (@$col){
			push(@t,$row->{$k});
			
		}
		   push(@tab,\@t);
       }
     
    }
);
$project->preload_patients();
$buffer->disconnect();
mce_loop {
   my ($mce, $chunk_ref, $chunk_id) = @_;
   	warn "start ".$chunk_id;
   if (ref($chunk_ref) ne "ARRAY") {
   	confess();
   }
   else {
   foreach my $patient (@$chunk_ref){
   	#next if $patient->name ne "CHAL";
   		my $hash;
   		$hash->{error};
   		eval{
   			
		 	$hash->{cnvs} = run($patient);
		 	delete $hash->{error};
   		};
		MCE->gather($chunk_id,$hash);
	}
   }
   	
}   @{$project->getPatients};
#grep{$_->name eq "DYSPH47_EA1"}
MCE::Loop->finish;
	
die() if $error;
my $fh;
open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
	my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
	$csv->print($fh, $col); 
	
	foreach my $line (@tab){
			$csv->print($fh, $line);
	}


$fh->close();
	 my $query = "
	COPY (
        SELECT * from read_csv_auto(['$filename']) order by patient,type,chr,type,start
    )
    TO '$parquet_file_quality'  (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE);";
    print "\n# $query\n";
     
     system("duckdb :memory: -c \"$query \"");	 
     
	
	
sub run {
		my ($patient) = @_;
		my @tres = ();
		my $pid= $patient->id;
		my $sql =qq{select * from '$parquet_file' where id <> 'Z' and patient=$pid ; };
		my $cmd = qq{duckdb -json -c "$sql"};
		my $res =`$cmd`;
		my $array_ref = [];
		$array_ref  = decode_json $res if $res;
		my $bl_name = "spectre.blacklist.bed.gz";
		$bl_name = "encode.blacklist.bed.gz";
		my $tabix_blacklist = Bio::DB::HTS::Tabix->new( filename => "/data-pure/public-data/repository/HG38/blacklist/$bl_name"); 
		my $file = $patient->rawDataWiseCondor();
		my $tabix = Bio::DB::HTS::Tabix->new( filename => $file );
		my $nb = 0;
			my $methods = $patient->getCallingCNVMethodsType("caller_depth");
		my $nb_patients = scalar(@{$project->getPatients});
			
	  	my $bitmask  = 0;
	  	foreach my $m1 (@$methods){
	  			$bitmask |= $buffer->bitmask_caller_sv->{lc($m1)};
	  	}
		foreach my $row (@$array_ref) {
				$nb++;
				warn $patient->name.":".$nb."/".scalar(@$array_ref) if $nb%300 ==0;
	 			my $cnv = $rocks->get($row->{id});
	 			
	 			my $value = 0;
	 			foreach my $c (@{$cnv->{cnv_origin}}){
	 				#die() if $c->{callers} == 512;
	 				$value = $value | $c->{callers};
	 			}
#	 			warn $value if $value ne 1024 and $value ne 32;
	 		#	if (($value & 512) && ($value & 32)) {

    		#		print "Spectre + hificnv\n";
    			#	die();

			#	}
	  	my @t = keys %{$cnv->{score}};
	  	my $start = $cnv->{start};
	  	my $end = $cnv->{end};
	  	if( abs($start-$end) < 1_000_000 && scalar(@t) == 1 &&  exists $cnv->{score}->{score_caller_depth}  ){
	  	#	warn "+++++++++++++++++++++++++++" if  (($value & 512) && ($value & 32));
	  		next unless ($value & $bitmask) == $bitmask;
	  		
	  		
	  	#	warn "------------------------------------";
	  	}
	  	#	next;
	  	#pour region blacklist 
	  	
		
		my $region = $cnv->{chromosome}.":$start-$end";
		my $res_blacklist = $tabix_blacklist->query("chr".$region);
		my $cnv_span = Set::IntSpan::Fast->new("$start-$end");
		my $cnv_len = $end - $start + 1;
		my $bl_span = Set::IntSpan::Fast->new();
		if ($res_blacklist){
			while (my $line = $res_blacklist->next) {
    		my ($b_chr, $b_start, $b_end) = split /\t/, $line;
    			$bl_span->add_range($b_start, $b_end);
			}
		}
		my $intersect = $cnv_span->intersection($bl_span);
		my $overlap_bp = scalar($intersect->as_array);   # taille totale des bases recouvertes
		my $overlap_pct = $overlap_bp / $cnv_len * 100;
		$row->{blacklist} = int($overlap_pct);
			#end region blacklist 
			# frequence allelic 
		my $perc = 100;	 
		if ( $patient->isRevio &&  scalar(@t) == 1){
		
			my $chrom =$project->getChromosome($cnv->{chromosome});
			 my $mean1 = $patient->meanDepth($chrom->fasta_name,$start-500,$start-50);
    		 my $mean2 = $patient->meanDepth($chrom->fasta_name,$end+50,$end+500);
    		 my $m1 = ($mean1+$mean2) /2;
    		 $m1 += 0.001;
    		 my $m = $patient->meanDepth($chrom->fasta_name,$start,$end);
    		 my $p = (abs($m1 - $m) / $m1) * 100;
    		 $row->{blacklist} = 80 if $p < 25;
    		 $row->{blacklist} = 70 if $p < 10;
    		# warn $mean1." ".$mean2." ".$patient->meanDepth($chrom->fasta_name,$start,$end);
			#	warn Dumper $v;
			
		}
		#filer in this project 
		if ($nb_patients> 7){
			my $nb_itp = scalar(values %{$cnv->{patients}});
			if ($nb_itp >= 0.5 * $nb_patients){
				$row->{blacklist} = 70;
			}
			if ($nb_itp >= 0.8 * $nb_patients){
				$row->{blacklist} = 80;
			}
		}
		
		if ($row->{blacklist} < 80 ){
		
		if ($cnv->{type} eq 'DEL' && scalar(@t) == 1  ){
			$perc = analyze_deletion($patient,$cnv);
			
			 #warn $perc."%" if $value == 544;
			if ($perc < 70) {
					$row->{blacklist} = 80;
			}
			if ($perc < 50) {
					$row->{blacklist} = 70;
			}
			
		}
		if ($cnv->{type} eq 'DUP' && scalar(@t) == 1  ){
			$perc = analyze_duplication($patient,$cnv);
			# warn $perc."%" if $value == 544;
			if ($perc < 70) {
					$row->{blacklist} = 80;
			}
			if ($perc < 60) {
					$row->{blacklist} = 70;
			}
		}
		$row->{perc_freq_allelic} = $perc;
			#end frequence allelic  
		if (abs($start-$end) < 1000000 && $patient->isRevio &&  scalar(@t) == 1  ){
			my ($v,$break,$delta)  = analyze_cnv_noise($patient,$cnv);
		
			$row->{blacklist} = 70 if $v > 30 or $break > 20;
			$row->{blacklist} = 80 if( $v > 40 or $break > 30) or $delta > 50 ;
			
#			if ($row->{blacklist} == 0){
#			my $chrom =$project->getChromosome($cnv->{chromosome});
#			 my $mean1 = $patient->meanDepth($chrom->fasta_name,$start-5000,$start-200);
#    		 my $mean2 = $patient->meanDepth($chrom->fasta_name,$end+200,$end+5000);
#    		 my $m1 = ($mean1+$mean2) /2;
#    		 my $m = $patient->meanDepth($chrom->fasta_name,$start,$end);
#    		 warn $m1." ".$m;
#    		 my $p = (abs($m1 - $m) / $m1) * 100;
#    		 warn int($p);
#    		# warn $mean1." ".$mean2." ".$patient->meanDepth($chrom->fasta_name,$start,$end);
#			}
			#	warn Dumper $v;
			
		}
	 	#my $m = $patient->meanDepth();
		}
		$row->{mz} = 0;
		$row->{mr} = 0;
		$row->{na} = 0;
		$row->{nr} = 0;
		$row->{pa} = 0;
		$row->{nb_caller} = scalar(@t);
		$row->{caller_sr} = 0;
		$row->{caller_depth} = 0;
		$row->{caller_coverage} = 0;
		$row->{caller_sr} = 1 if exists $cnv->{score}->{score_caller_sr};
		$row->{caller_depth}= 1 if exists $cnv->{score}->{score_caller_depth};
		$row->{caller_coverage} = 1 if exists $cnv->{score}->{score_caller_coverage};
		push(@tres,$row);
		} #end for each row
		warn "@@@@@@@@@@@@@@@@@@@@@@@@@@@ ".scalar(@tres);
		return \@tres;
}
	
	

sub test_type {
	my($hash,$flag) = @_;
	$flag = "caller_".$flag;
	confess() unless $caller_type_flag->{$flag};
	return $hash->{caller_type_flag} & $caller_type_flag->{$flag};
}

sub analyze_deletion {
	   my ($patient, $cnv) = @_;
	   my $start = $cnv->{start};
	my $end = $cnv->{end};
	my $chr_name =$project->getChromosome($cnv->{chromosome})->ucsc_name;;
	
	my $parquet = $project->parquet_cache_variants();
	my $table_patient = "patient_".$patient->id."_type";
	my $cp2 = "patient_".$patient->id."_type";
	#CAST(patient_55661_transmission AS INTEGER)
	my $asql_patient = [];
	my $asql_patient_only_transmission = [];
	my $t = time;
	#my $sql = qq{ SELECT  variant_start,$table_patient  FROM '$parquet' a where  ${table_patient} <> 0 and variant_type=1 and variant_chromosome='${chr_name}' and variant_start>${start} and variant_end < ${end} order by variant_start;};
	my $sql= qq{SELECT

      count_if(${table_patient} = 1) AS nb_type_1,

        count_if(${table_patient} = 2) AS nb_type_2
FROM '$parquet'

WHERE variant_type = 1
  AND ${table_patient} <> 0
  AND variant_chromosome = '${chr_name}'

  AND variant_start > ${start}

  AND variant_end < ${end};
	};
	
	#warn $sql if $start == 29691640;
	#21:29691640-30238026 
	my $cmd = qq{duckdb -column -noheader -c "$sql"};
	
	my $st = `$cmd`;
	
	chomp($st);
	my $ho = 0;
	my $he = 0;
	
	 ($he,$ho) = split(" ",$st);
if ($chr_name eq "chrX" && $patient->isMale){
	my $debug;
	$debug = 1 if $start == 151548000 && $patient->name eq "DYSPH47_EA1";#-156040000 
	warn "x male $start ".$chr_name ;#if $debug;
	my $intspan = $project->getChromosome($cnv->{chromosome})->intspan_pseudo_autosomal;
	my $intspan2 = Set::IntSpan::Fast->new("$start-$end");
	my $i = $intspan->intersection($intspan2);
	warn $i->as_string if $debug;
	warn $intspan->as_string() if $debug;
	warn $intspan2->as_string() if $debug;
	if ($i->is_empty){
		warn "********************************* coucou ".$ho if $debug;
		return 0 if $ho+$he >= 10;
	}
}	

return 100 if ($ho+$he) < 10;
my $perc = int ( ($ho/($he+$ho)) *100);

#warn "\t\t".$perc;
return $perc;	
}
sub analyze_duplication {
	   my ($patient, $cnv) = @_;
	   my $start = $cnv->{start};
	my $end = $cnv->{end};
	my $chr_name =$project->getChromosome($cnv->{chromosome})->ucsc_name;;

	my $parquet = $project->parquet_cache_variants();
	my $table_patient2 = "patient_".$patient->id."_ratio";
	my $table_patient = "patient_".$patient->id."_type";
	my $cp2 = "patient_".$patient->id."_type";
	#CAST(patient_55661_transmission AS INTEGER)
	my $asql_patient = [];
	my $asql_patient_only_transmission = [];
	my $t = time;
	my $sql = qq{ SELECT  variant_start,$table_patient2  FROM '$parquet' a where  ${table_patient} == 1 and variant_type=1 and variant_chromosome='${chr_name}' and variant_start>${start} and variant_end < ${end} order by variant_start;};
	my $sql = qq{SELECT

     count_if(${table_patient2} >= 60 OR ${table_patient2} <= 35) AS nb_type_1,
      count_if(${table_patient2} < 60 AND ${table_patient2} > 35) AS nb_type_2
	FROM '$parquet'

	WHERE variant_type = 1
  	AND ${table_patient} = 1
  		AND variant_chromosome = '${chr_name}'

  	AND variant_start > ${start}

  	AND variant_end < ${end};
	};
	
	#warn $sql if $chr_name eq "chr1";
	
	my $cmd = qq{duckdb -column -noheader -c "$sql"};
	my $st = `$cmd`;
	chomp($st);
	my ($ok,$nok) = split(" ",$st);
	$ok += 0;
	$nok += 0;
#	open my $fh, '-|', $cmd or die "Cannot execute $cmd: $!";
#	my @snp;
#	my $ok = 0;
#	my $nok = 0 ;
#while (my $line = <$fh>) {
#	 	my ($pos, $value) = split(" ",$line);
#	 	next if ($value > 95);
#	 	if ($value > 60 or $value< 40){
#	 		$ok ++ ;
#	 		next;
#	 	}
#	 		$nok ++;
#	}
#close($fh);	
my $nb = $ok+$nok;
#warn $nb;
return 100 if $nb < 20;
my $perc = int ( ($ok/($ok+$nok)) *100);
warn $perc." ".$patient->name." ".$chr_name.":".$start."-".$end ." DUP" if $chr_name eq "chr1";
return $perc;	
#warn $ok." ".$nok;
	
}
sub analyze_cnv_noise {
    my ($patient, $cnv) = @_;
    
    my $bin_size = 500;
    my $chrom = $project->getChromosome($cnv->{chromosome})->fasta_name;
    my $pos_start = $cnv->{start};
    my $pos_end = $cnv->{end};
	my $bw = $patient->getBigWigFile();
	
    # 1. Calculer combien de fenêtres (bins) tiennent dans ce CNV
    my $region_length = $pos_end - $pos_start;
    my $num_bins = int($region_length / $bin_size);
    
    # Sécurité : si le CNV est plus petit qu'un bin, on prend au moins 1 valeur
    $num_bins = 1 if $num_bins < 1;

    # 2. Lancer la commande UCSC bigWigSummary
    # Rediriger STDERR vers /dev/null pour cacher les warnings de l'UCSC
    my @values;
    for (my $i =$pos_start;$i <= $pos_end-$bin_size;$i+=$bin_size){
    $patient->depth($chrom,$i,$i+$bin_size);
    	push(@values,$patient->meanDepth($chrom,$i,$i+($bin_size*2)));
    }
    
    my $mean1 = $patient->meanDepth($chrom,$pos_start-5000,$pos_start-200);
     my $mean2 = $patient->meanDepth($chrom,$pos_end+200,$pos_end+5000);
    #my $cmd = "/software/bin/bigWigSummary $bw $chrom $pos_start $pos_end $num_bins 2>/dev/null";
    #my $output = `$cmd`;
    #chomp $output;

    # Si la commande échoue ou ne renvoie rien
    #if (!$output) {
     #   return (0, 0, 0, "ERREUR_LECTURE (ou région sans aucune donnée)");
    #}

    # 3. Récupérer les valeurs (elles sont séparées par des tabulations)
    #my @values = split(/\t/, $output);

    # 4. Nettoyer les données : bigWigSummary renvoie "n/a" s'il y a un trou absolu (0 reads)
   # for my $val (@values) {
    #    if ($val eq 'n/a') {
     #       $val = 0;
      #  }
    #}

    # 5. Calculs statistiques
    my $n = scalar(@values);
    my $sum = 0;
    $sum += $_ for @values;
    my $mean = $sum / $n;

    # Si la moyenne est à 0 (aucun read du tout sur tout le CNV)
    if ($mean == 0) {
        return (0, 0, 0, "REGION_Morte (0X)");
    }

    # Calcul de l'écart-type
    my $sq_diff_sum = 0;
    $sq_diff_sum += ($_ - $mean)**2 for @values;
    
    # Division par (n-1) pour la variance de l'échantillon (sécurité pour division par 0)
    my $variance = ($n > 1) ? $sq_diff_sum / ($n - 1) : 0;
    my $stddev = sqrt($variance);

    # Calcul du Coefficient de Variation (en pourcentage)
    my $cv = ($stddev / $mean) * 100;
   my $noise = 0;

for (my $i = 1; $i < @values; $i++) {

    $noise += abs($values[$i] - $values[$i-1]);

}

$noise /= (@values-1);
my $threshold = 0.20;   # 20 %

my $count = 0;

for (my $i=1;$i<@values;$i++){

    my $local = ($values[$i]+$values[$i-1])/2;

    next if $local == 0;

    my $delta = abs($values[$i]-$values[$i-1])/$local;

    $count++ if $delta > $threshold;

}

my $fraction = $count/(@values-1);

	my $r = plateau_breaks(\@values,3,0.20);


my $nbbp = $r->{nb_breaks};
	my $pb = ($r->{nb_breaks}/(@values))*100;
	$mean1+=0.0001;
	$mean2+=0.0001;
#	warn $mean1." ".$mean2;
	my $outside_delta = abs($mean1-$mean2)/(($mean1+$mean2)/2);
#    warn $cv." - noide: $noise - $pb - ".$mean." ".$mean1." ".$mean2." : ".($outside_delta*100);# if $pos_start == 55675001;
  #  warn Dumper @values if $pb> 20;
    
  #  die()   if $noise > 0.3;
   # die();
	return ($mean,$cv,$pb,$outside_delta);
	
}


sub plateau_breaks {

    my ($values,$window,$threshold) = @_;

     $window    //= 4;      # taille des fenêtres

     $threshold //= 0.25;   # 25%

    my @v = @$values;

    my $nb_breaks = 0;

    my @breaks;

    for (my $i=$window;$i<@v-$window;$i++) {

        my ($left,$right) = (0,0);

        for(my $j=1;$j<=$window;$j++) {

            $left  += $v[$i-$j];

            $right += $v[$i+$j];

        }

        $left  /= $window;

        $right /= $window;

        next if $left==0 && $right==0;

        my $ref = ($left+$right)/2;

        next if $ref==0;

        my $delta = abs($left-$right)/$ref;
		
        if($delta > $threshold){

            push @breaks,{

                pos   => $i,

                left  => $left,

                right => $right,

                delta => $delta,

            };

            $nb_breaks++;

        }

    }

    return {

        nb_breaks => $nb_breaks,

        breaks    => \@breaks,

    };

}

