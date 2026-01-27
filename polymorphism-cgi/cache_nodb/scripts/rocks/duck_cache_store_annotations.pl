#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use Parallel::ForkManager;
use Storable qw(store retrieve freeze thaw);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Cwd 'abs_path';
use Digest::MD5 qw(md5_hex);
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use lib "$RealBin/../";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime); 
use JSON::XS;
use Compress::Snappy;
use Getopt::Long;
use Carp; 
use GBuffer;
use Cache_Commons;
use Sys::Hostname;
use GenBoNoSqlRocksVector;
use Text::CSV;
use Devel::Size qw(size total_size);
use MCE::Loop;
use Storable qw(dclone);  # pour cloner le project par process si nécessaire
use Devel::Size qw(total_size);
use Storable qw(dclone);
 my $host = hostname();
use Scalar::Util qw(looks_like_number);
use lib "$RealBin/../../../GenBo/lib/obj-nodb/polyviewer/";
use PolyviewerVariant;
use MCE::Shared;

my $htr = {
	 
                      'mane' => -99,
                      'end' => '-',
                      'spliceAI_cat' => '-',
                      'appris' => '-',
                      'codons_AA' => '-',
                      'name' => '-',
                      'spliceAI' => -99,
                      'dbscsnv' => -99,
                      'impact_score_text' => '-',
                      'consequence' => '-',
                      'polyphen' => -99,
                      'start' => '-',
                      'codons' => '-',
                      'alphamissense' => -99,
                      'sift' => -99,
                      'prot' => '-',
                      'nomenclature' => 'c.-1--07dupT',
                      'impact_score' => -99,
                      'cadd' => -99,
                      'nm' => "-",
                      'exon' => '-',
                      'revel' => -99,
                      'ccds' => "-",
                      'main' => -99
}; 


warn "*_*_*_*_*_".$host."*_*_*_*_*_";

my $fork = 1;
my ($project_name, $chr_name, $annot_version);
my $ok_file;
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'annot_version=s'    => \$annot_version,
	'chr=s'        => \$chr_name,
	'file=s' => \$ok_file,
);

 if ($ok_file && -e $ok_file) {
 	system("rm $ok_file");
 }

my $parquets =[];



unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }


my $buffer  = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProjectCache( -name => $project_name );
my $dir_tmp_cvs =  $project->getCallingPipelineDir($project->name.".parquet.".time);

##################
#construct_sql 
##################
$project->preload_patients();
$project->validations();
my @pids =  map {$_->id} @{$project->getPatients};

my %cast;


my @column_variant = ("start","end","chromosome","index","vector_id","rocksdb_id","type","gnomad_id","name","hgmd_id","clinvar_id","hgmd_class","clinvar_class","isDM","promoterAI","isClinvarPathogenic","keepPathogenic");
my $position_keepPathogenic;
for (my $i=0;$i<@column_variant;$i++){
	$position_keepPathogenic = $i if $column_variant[$i] eq "keepPathogenic";
}

die() unless $position_keepPathogenic;
my @column_frequences_gnomad = ("gnomad_ac","gnomad_an","gnomad_min_pop_name","gnomad_min_pop_freq","gnomad_max_pop_name","gnomad_max_pop_freq","gnomad_ho","getGnomadAC_Male");
my @column_frequences_dejavu = ("other_projects","other_patients","other_patients_ho","similar_projects","similar_patients","similar_patients_ho","in_this_run_patients");

my @integer = ("variant_start","variant_end","variant_gnomad_ac","variant_gnomad_an","variant_gnomad_ho","variant_getGnomadAC_Male","variant_other_projects","variant_other_patients","variant_other_patients_ho","variant_similar_projects","variant_similar_patients","variant_similar_patients_ho","variant_in_this_run_patients");
my @varchar = ("variant_chromosome","variant_index","variant_vector_id","variant_rocksdb_id","variant_type","variant_gnomad_id","variant_name","variant_hgmd_id","variant_clinvar_id","variant_hgmd_class","variant_clinvar_class","variant_gnomad_max_pop_name","variant_gnomad_min_pop_name"); 
my @boolean = ("variant_isDM","variant_isClinvarPathogenic","variant_keepPathogenic");
my @float = ("variant_gnomad_min_pop_freq","variant_gnomad_max_pop_freq","variant_promoterAI");

my @column_variant_names  = map{"variant_".$_} (@column_variant,@column_frequences_gnomad,@column_frequences_dejavu);
my @column_patient;
my @column_gene = ("gene_name","gene_mask");
push(@varchar,"gene_name");
push (@integer,"gene_mask");
foreach my $c ( sort {$a <=> $b} @pids){
	push(@column_patient,"patient_".$c."_ref");
	push(@integer,"patient_".$c."_ref");
	push(@column_patient,"patient_".$c."_alt");
	push(@integer,"patient_".$c."_alt");
	push(@integer,"patient_".$c."_ratio");
	push(@column_patient,"patient_".$c."_ratio");
	
	push(@column_patient,"patient_".$c."_type");
	push(@varchar,"patient_".$c."_type");
	push(@column_patient,"patient_".$c."_transmission");
	push(@integer,"patient_".$c."_transmission");
	#push(@float,"score_".$c);
	#push(@column_gene,"score_".$c);
}



my @batches;
my $aselect;

foreach my $c (@integer){
	  push(@$aselect, "CAST($c AS INTEGER) AS $c");
}

foreach my $c (@varchar){
	  push(@$aselect, "CAST($c AS VARCHAR) AS $c");
}

foreach my $c (@boolean){
	  push(@$aselect, "CAST($c AS BOOLEAN) AS $c");
}
foreach my $c (@float){
	  push(@$aselect, "CAST($c AS FLOAT) AS $c");
}
my $select = join(",",@$aselect);
#	
my $tume = time;
my $diro = $project->rocks_directory();
my $nb = 0;
foreach my $chr (@{$project->getChromosomes}){
	#next if $chr->name ne "21";
	my $no = $chr->get_rocks_variations("r");
   	my $nproc = 10;	
	my $N = $no->size;
 	$no->close;
 
 	my $chunk_size = int($N / $nproc);
	my $remainder  = $N % $nproc;
	my $index = 0;
	for my $i (0 .. $nproc-1) {
    	my $extra = ($i < $remainder) ? 1 : 0;
    	my $size  = $chunk_size + $extra;
    	my $filename = $dir_tmp_cvs."/".$chr->name.".$index.".$project->name.".csv";
    	##$dir_tmp_cvs/".$project->name."_".$chr->name.".csv";
    	$nb ++;
    	push @batches, {start=>$index,end=>$index+$size,filename=>$filename,chromosome=>$chr,nb=>$nb};
    	$index += $size;
	}
}
$project->disconnect();
delete $project->{rocks};
my $shared_hash = MCE::Shared->hash();
MCE::Loop::init {
    chunk_size => 1,  # chaque worker recevra un seul élément
    max_workers => 20, 
     gather => sub {
        my ($mce) = @_;
        delete $shared_hash->{$mce};
     }
};
mce_loop {
    my ($mce, $chra,$chunk_id) = @_;
      $shared_hash->{$chunk_id} ++;
    my $region = $chra->[0];
    warn "start ".$region->{nb}."/".$nb ." ".$region->{chromosome}->ucsc_name;
     warn "end ".$region->{nb}."/".$nb ." ".$region->{chromosome}->ucsc_name;
   	compute($region,$region->{chromosome});
  	MCE->gather($chunk_id);
  
} @batches;
MCE::Loop->finish;
confess()  if %$shared_hash;
my $filename = $dir_tmp_cvs."/*.csv";# join(",",map{"\'".$_->{filename}."\'"} @batches);

my $parquet_file_final = $project->parquet_cache_variants; 
unlink $parquet_file_final if -e $parquet_file_final;
my $query = qq{
	COPY (
        SELECT $select from read_csv_auto([\'$filename\']) order by variant_gnomad_ac
    )
    TO '$parquet_file_final' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE,ROW_GROUP_SIZE 1000000);};
 my $cmd = "duckdb  -c \"$query\"";
 
warn $cmd;

system($cmd) == 0 or die "Erreur DuckDB : $?";  
die() unless -e $parquet_file_final;
warn "ok";
exit(0);


sub compute {
	my ($region,$chr) =@_;
    my $no = $chr->get_rocks_variations("r");
    my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
	my $filename = $region->{filename};
#$dir_tmp_cvs."/".$chr->name.".".$project->name.".csv";#$dir_tmp_cvs/".$project->name."_".$chr->name.".csv";
	my $fh;		
	open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
# my $diro = $project->rocks_directory();
	
		my $t   = time;
		my $res;
	#	my $final_polyviewer_all = GenBoNoSqlRocks->new(dir=>$diro,mode=>"r",name=>"polyviewer_objects",cache=>1);
	my @col = (@column_variant_names,@column_patient,@column_gene);
	$csv->print($fh, \@col); 
    my $nb = 0;
    for (my $i = $region->{start}; $i < $region->{end}; $i++) {
        my $variation = $no->get_index($i);
        $variation->{buffer}  = $buffer;
        $variation->{project} = $project;
        my $duck_variants =[];
     #my @column_variant = ("start","end","chromosome","index","vector_id","rocksdb_id","type","gnomad_id","name","hgmd_id","clinvar_id","hgmd_class","clinvar_class","isDM","isClinvarPathogenic","keepPathogenic");
        my $vh  ={};
        $vh->{start}           = $variation->start;
        $vh->{end}           = $variation->end;
      	$vh->{chromosome}           = $chr->ucsc_name;	
      	$vh->{index} = $chr->name."!".$i;
      	my $index = $chr->name."!".$i;
      	$vh->{vector_id} = $i;
        $vh->{rocksdb_id}  = $variation->genomic_rocksdb_id;
         
        if    ($variation->isVariation())         { $vh->{type} = 1; }
        elsif ($variation->isInsertion())         { $vh->{type} = 2; }
        elsif ($variation->isDeletion())          { $vh->{type} = 3; }
        elsif ($variation->isLargeDeletion())     { $vh->{type} = 4; }
        elsif ($variation->isLargeDuplication())  { $vh->{type} = 5; }
        elsif ($variation->isLargeInsertion())    { $vh->{type} = 6; }
        elsif ($variation->isJunction())          { $vh->{type} = 7; }
        else  { confess("Unknown variation type"); }
         
        $vh->{gnomad_id}           = $variation->name;
        $vh->{name}           = $variation->name;
        $vh->{hgmd_id}  = return_char($variation->hgmd_id);
        $vh->{clinvar_id}  = return_integer($variation->clinvar_id);
        $vh->{clinvar_class} = return_char($variation->clinvar_class);
        $vh->{clinvar_class}   ="-" unless  $variation->clinvar_class;
        $vh->{hgmd_class} = return_char($variation->hgmd_class);
        $vh->{isDM}             = $variation->isDM;
      	$vh->{isDM} =0 unless $vh->{isDM};
        if ($variation->score_clinvar == 5){
        	 $vh->{isClinvarPathogenic} = 1;
        }
        else {
        	 $vh->{isClinvarPathogenic} = 0;
        }
        $vh->{keepPathogenic} = 0;
         if ($vh->{isDM} ==1  or $vh->{isClinvarPathogenic} == 1 ){
         	$vh->{keepPathogenic} = 1 if $variation->getGnomadHO < 100;
         	
         }
         
         $vh->{promoterAI} = -99;
         my $h_promoterAI = $variation->promoterAI();
         if ($h_promoterAI) {
         	my @list_tr = keys %$h_promoterAI;
         	$vh->{promoterAI} = $h_promoterAI->{$list_tr[0]}->{score};
         }
         
         foreach my $k (@column_variant) {
        	die($k) unless exists $vh->{$k};
        	 $vh->{$k} = -1 unless defined $vh->{$k};
        	push(@$duck_variants,$vh->{$k});
        	delete $vh->{$k};
        }
        warn Dumper $vh if keys %$vh;
        die() if keys %$vh;
        #les frequences  GNOMAD
        #@column_frequences_gnomad = ("gnomad_ac","gnomad_an","gnomad_min_pop_name","gnomad_min_pop_freq","gnomad_max_pop_name","gnomad_max_pop_freq","gnomad_ho","getGnomadAC_Male");
        
         #my $vp = $final_polyviewer_all->get($index,1);
         #bless $vp , 'PolyviewerVariant';
        # dejavu_other_projects
        $vh->{gnomad_ac}           = $variation->getGnomadAC;
        $vh->{gnomad_an}        = $variation->getGnomadAN;
       	$vh->{gnomad_min_pop_name}        = $variation->min_pop_name;
       	$vh->{gnomad_min_pop_name}        =  "-" unless 	$vh->{gnomad_min_pop_name} ;
        $vh->{gnomad_min_pop_freq}        = $variation->min_pop_freq;
        $vh->{gnomad_min_pop_freq}        =  -1 if $vh->{gnomad_min_pop_freq} eq "-";
       	$vh->{gnomad_max_pop_name}        = $variation->max_pop_name;
       	$vh->{gnomad_max_pop_name}        =  "-" unless 	$vh->{gnomad_max_pop_name} ;
        $vh->{gnomad_max_pop_freq}        = $variation->max_pop_freq;
         $vh->{gnomad_max_pop_freq}        =  -1 if $vh->{gnomad_max_pop_freq} eq "-";
       	$vh->{gnomad_ho}        = 			$variation->getGnomadHO;
        $vh->{getGnomadAC_Male}        = $variation->getGnomadAC_Male;
          foreach my $k (@column_frequences_gnomad) {
        	die($k) unless exists $vh->{$k};
        	 $vh->{$k} = -1 unless defined $vh->{$k};
        	push(@$duck_variants,$vh->{$k});
        	delete $vh->{$k};
        }
         warn Dumper $vh if keys %$vh;
        die() if keys %$vh;
        #les frequences  DEJAVU
        #my @column_frequences_dejavu = ("other_projects","other_patients","other_patients_ho","similar_projects","similar_patients","similar_patients_ho","in_this_run_patients");
         $vh->{other_projects}   = $variation->other_projects;
#         warn $vp->dejavu_other_projects ." ".$vp->dejavu_other_patients_ho if $vp->dejavu_other_projects > 20000;
         die() if $variation->other_projects > 20000;
        $vh->{other_patients}        = $variation->other_patients;
       	$vh->{other_patients_ho}        = $variation->other_patients_ho;
        $vh->{similar_projects}        = $variation->similar_projects;
       	$vh->{similar_patients}        = $variation->similar_patients;
        $vh->{similar_patients_ho}        = $variation->similar_patients_ho;
       	$vh->{in_this_run_patients}   =  $variation->in_this_run_patients;
        foreach my $k (@column_frequences_dejavu) {
        	die("die :".$k) unless exists $vh->{$k};
        	 $vh->{$k} = -1 unless defined $vh->{$k};
        	push(@$duck_variants,$vh->{$k});
        	delete $vh->{$k};
        }
        confess(Dumper $vh) if keys %$vh;
 		
       my $hpatients;
       my $duck_patients;
        
       
       my $hscore;
        foreach my $pid ( sort {$a <=> $b} @pids){
        	unless (exists $variation->{sequencing_infos}->{$pid}){
        		push(@$duck_patients,0);#ref
        		push(@$duck_patients,-1);#alt
        		push(@$duck_patients,-1);#ratio
        		push(@$duck_patients,0);#type
        		push(@$duck_patients,0);#transmission
        	}
        	else {
        		push(@$duck_patients, $variation->{sequencing_infos}->{$pid}->{max}->[0]);#ref
        		push(@$duck_patients, $variation->{sequencing_infos}->{$pid}->{max}->[1]);#alt
        		my $s =$variation->{sequencing_infos}->{$pid}->{max}->[0]+$variation->{sequencing_infos}->{$pid}->{max}->[1];
        		$s = 1 if $s ==0;
        		my $r = $variation->{sequencing_infos}->{$pid}->{max}->[1]/($s);
        		$r = int($r*100);
        		push(@$duck_patients,$r);#ratio
        		my $v =2 ;
        		$v =1 if $variation->{sequencing_infos}->{$pid}->{max}->[2] eq "he";
        		push(@$duck_patients,$v);#type
        		#warn return_transmissions($project->{objects}->{patients}->{ $pid },$variation);
        		push(@$duck_patients,return_transmissions($project->{objects}->{patients}->{ $pid },$variation));#transmission
        		}
        	}
        	if (@{$variation->{genes_object}}){
				foreach my $gid (@{$variation->{genes_object}}){
					
					my $g = $project->newGene($gid);
					    if ($variation->score_validations($g) > 0){
							#petite astuce pour changer le keep pathogenic dans le tableau 					    	
					    	$duck_variants->[$position_keepPathogenic]  = 1;
					    	
         					#$vh->{keepPathogenic} = 1 ;
		         	
		        		}
						my $duck_gene;
						push(@$duck_gene,$gid);
						push(@$duck_gene,$variation->{annotation}->{$gid}->{mask});
						#foreach my $pid ( sort {$a <=> $b} @pids){
			 			#	my $sc = $variation->scaledScoreVariant($g,$project->{objects}->{patients}->{ $pid });
			 			#	push(@$duck_gene,$sc);
			 			#}
						my $score;
			 			my @values = (@$duck_variants,@$duck_patients,@$duck_gene);
			 			confess() unless scalar(@col) == scalar(@values);
			 			$csv->print($fh, \@values);
			 			if ($index eq "1!127"){
			 				#warn Dumper @values;
			 			}
			 			
			 		
				}
        	}
        	else{
        		my $duck_gene;
        		push(@$duck_gene,"-");
				push(@$duck_gene,2);
				#foreach my $pid ( sort {$a <=> $b} @pids){
			 	#			my $sc = 0;
			 	#			push(@$duck_gene,$sc);
			 	#		}
				my @values = (@$duck_variants,@$duck_patients,@$duck_gene);
			 	confess() unless scalar(@col) == scalar(@values);
			 	$csv->print($fh, \@values);
        	}
	
 	delete $variation->{project};
 	delete $variation->{buffer};
			
    }
	$fh->close;
	$no->close();
	delete $project->{rocks};
   # warn "FIN ".$chr->name." ".$region->{start}."-".$region->{end};

	
} 


sub return_transmissions {
	my ($patient,$variant) = @_;
	my $h_models_ids = {
    solo          => 1 << 0,  # 2^0 = 1
    father        => 1 << 1,  # 2^1 = 2
    mother        => 1 << 2,  # 2^2 = 4
    both          => 1 << 3,  # 2^3 = 8
    is_parent     => 1 << 4,  # 2^4 = 16
    recessif      => 1 << 5,  # 2^5 = 32
    dominant      => 1 << 6,  # 2^6 = 64
    denovo        => 1 << 7,  # 2^7 = 128
    strict_denovo => 1 << 8,  # 2^8 = 256
    error         => 1 << 9,  # 2^9 = 512
    mosaic         => 1 << 10,  # 2^10 = 1024
    uniparental     => 1 << 11,  # 2^11 = 2048
    
};
	
	my $fam = $patient->getFamily;
	if ($fam->isSolo){
		return $h_models_ids->{solo};
	}
	if ($patient->isParent) {
		return $h_models_ids->{is_parent};
	}
	my $model = 0;
	$model |= $h_models_ids->{dominant} if $variant->isDominantTransmission($fam,$patient);

	if ($variant->isBothTransmission($fam,$patient)){
		$model |= $h_models_ids->{both};
		return $model;
	}
	$model |= $h_models_ids->{father} if  $variant->isFatherTransmission($fam,$patient);
	$model |= $h_models_ids->{mother} if $variant->isMotherTransmission($fam,$patient);
	
	$model |= $h_models_ids->{recessif} if $variant->isRecessiveTransmission($fam,$patient);
	$model |= $h_models_ids->{denovo} if $variant->isDenovoTransmission($fam,$patient);
	$model |= $h_models_ids->{strict_denovo} if $variant->isStrictDenovoTransmission($fam,$patient);
	$model |= $h_models_ids->{mosaic} if $variant->isMosaicTransmission($fam,$patient);
	$model |= $h_models_ids->{uniparental} if $variant->isUniparentalDisomyTransmission($fam,$patient);
	die($variant->name." ".$variant->id." ".$patient->name) if $model == 0;	
	return $model;
	
	
	
}

sub return_char {
	my ($value) =@_;
	return "-" unless $value;
	if (looks_like_number($value)){
		return "-";
	}
	return $value;
}
sub return_integer {
	my ($value) =@_;
	return -99 unless $value;
	if (looks_like_number($value)){
		return $value;
	}
	return -99;
}

		
