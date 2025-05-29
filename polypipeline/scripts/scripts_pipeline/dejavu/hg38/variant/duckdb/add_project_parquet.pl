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
use lib "$RealBin/../../../../../../../GenBo/lib/obj-nodb/";
use lib "$RealBin/../../../../../../../GenBo/lib/obj-nodb/packages";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use DBI;
use Compress::Snappy;
use Storable qw/thaw freeze/;
use GenBoNoSqlRocksGenome;
use GenBoNoSqlRocksAnnotation;
use File::Slurp qw(write_file);
use List::Util qw/shuffle/;
use lib "$RealBin/../../../utility";
use liftOverRegions;
use chunks;
use Sereal qw(sereal_encode_with_object sereal_decode_with_object write_sereal_file read_sereal_file);
use GenBoNoSqlDuckDBDejaVu;
use Text::CSV;
use MIME::Base64;
use Storable qw(dclone);
use MCE::Loop;
use MCE::Flow;
use dejavu_duckdb;
my $fork = 1;
my $project_name;
my $force;

GetOptions(
	'fork=s' => \$fork,
	'project=s' => \$project_name,
	'force=s' => \$force,
);


my $f2;


my $dir_tmp = "/tmp/pipeline";
system("mkdir /tmp/pipeline") unless -e "/tmp/pipeline";
my $buffer = new GBuffer;

my $project = $buffer->newProjectCache( -name => "$project_name");


my $lift = liftOverRegions->new(project=>$project,version=>$project->lift_genome_version); 

my $data_lift;	
my @chromosomes = shuffle (@{$project->getChromosomes});
#my $dir_parquet = "/data-beegfs/parquet_HG19_HG38/";
#my $dir_parquet = "/data-beegfs/projects.parquet/";
my $dir_parquet = $buffer->dejavu_parquet_dir();

my $parquet_file = $dir_parquet."/".$project->name.".".$project->id.".parquet";

my $can_dejavu = 1;
foreach my $pname (@{$buffer->getQuery->listProjectsWithoutDejaVu()}) {
	next if $pname ne $project_name;
	$can_dejavu = undef;
	last;
}

if (not $can_dejavu) { $parquet_file .= '.no_dejavu'; }

warn $parquet_file;
#exit(0) if -e $parquet_file and not $force;
$project->getPatients;
$project->preload_patients();
MCE::Loop->init(
   max_workers => $fork, chunk_size => 'auto'
);
$project->lift_genome_version();
$project->disconnect;
my @results = mce_loop {
   my ($mce, $tregions, $chunk_id) = @_;
   
   my $local_project;
   my $files;
   foreach my $chr (@$tregions){
   	my $t = time;
   	my $lift = liftOverRegions->new(project=>$project,version=>$project->lift_genome_version);
   	my $file;
   	if($project->isRocks){ 
     $file = save_and_lift_rocks($chr,$lift);
   	}
    else {
    	$file = save_and_lift_lmdb($chr,$lift);
    }
    push(@$files,$file);
   }
	 MCE->gather($files);
} @chromosomes;	

my @files;
foreach my $t (@results) {
	push(@files,@$t);
}
@files = sort {
    ($b =~ /_X\.csv$/) <=> ($a =~ /_X\.csv$/)  # "_X" en premier
    || $a cmp $b                                # Trie alphabétique normal après
} @files;
my $filename = join(",",grep {$_}@files);

exit(0) if $filename eq "";

my $dbh = DBI->connect("dbi:ODBC:Driver=DuckDB;Database=:memory:", "", "", { RaiseError => 1 , AutoCommit => 1});
my $query = "
	COPY (
        SELECT * from read_csv_auto([$filename]) order by chr38,pos38,chr19,pos19,allele
    )
    TO '$parquet_file' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE,ROW_GROUP_SIZE 1000000);";
    warn $query;
$dbh->do($query);
$filename =~s/,/ /g;
system("rm $filename");
exit(0);



sub save_and_lift_lmdb {
	my($chr,$lift) =@_;
	warn $chr->name." --- ";
	my $dir1 = $buffer->deja_vu_public_dir("HG19","projects_lite");
	my $f1 = "$dir1".$project->name.".lite";
	my $no = GenBoNoSql->new( dir => $dir1 , mode => "r" );
	my $h_hg19 = $no->get($project_name, $chr->name);
	my $patients = $no->get($project_name, "patients");
	$no->close;
	
	
	my $hpat ={};
	my $hpat2 ={};
	my $snps;
	my $lt = time;
	foreach my $name (keys %$patients){
		$hpat->{$patients->{$name}} = $project->getPatient($name)->id;
		$hpat2->{$project->getPatient($name)->id} = $patients->{$name};
	}
	my $t =time;
	my $snps;
	my $variations;
	my $hregions;
	my $cpt =0;
	
	my $lmdb_dir = $project->lmdb_cache_variations_dir();
	$lmdb_dir =~ s/rocks\///;
	warn $lmdb_dir;
	my $lmdb = GenBoNoSqlLmdb->new(
		dir         => $lmdb_dir,
		mode        => 'r',
		is_index    => 1,
		name        => $chr->name,
		is_compress => 1,
		vmtouch     => 1,
	);
	
	my $time = time;
	
	unless ($h_hg19 ) {
		my $file = dejavu_duckdb::save_csv($chr,undef,$dir_tmp);
		$no->close();
		return $file;
	}
	unless (keys %$h_hg19) {
		my $file = dejavu_duckdb::save_csv($chr,undef,$dir_tmp);
		$no->close();
		return $file;
	}
	
	my $vectors = dejavu_duckdb::get_hash_model_variant($chr);
	foreach my $vid (keys %{$h_hg19}){
		$cpt ++ ;

		my $vh;
		my ($c,$p,$a,$b) = split("_",$vid);
		my $vhh;
		$vhh->{chromosome} = $chr->ucsc_name;
		$vhh->{start} = $p;
		$vhh->{end} = $p+length($b);
		my $obj = $lmdb->get($vid);
		next unless $obj;
		my $annex = $obj->{annex};
		my $index_lmdb = $lmdb->get($vid)->{index_lmdb};
		my @ho;
		my @he;
		my @ho_ratio;
		my @he_ratio;
		my $nbhe = 0;
		my $nbho = 0;
		my $max_dp = 0;
		my $max_ratio = 0;
		my $bit_models = 0;
		foreach my $pid (keys %$annex){
			my @values;
			my $dp = $annex->{$pid}->{dp};
			my $alt = $annex->{$pid}->{nb_all_mut};
			$dp =1 if $dp == 0;
			my $ratio = int ($alt/$dp)*100;
			
			if ($dp =~ /,/){
			 my ($a,$b) = split(",",$dp);
			 $dp = $a;
			 }
			$max_dp = $dp if $dp > $max_dp;
			$max_ratio = int($ratio) if int($ratio) > $max_ratio;
			my $model= dejavu_duckdb::find_variant_model($vectors, $index_lmdb, $pid);
			$bit_models = $bit_models | $model;
			
			push(@values,$dp);
			push(@values,$alt);
			push(@values,$model);
			if ($annex->{$pid}->{he} == 1){
				push(@he,$pid);		
				push(@he_ratio,@values);
				$nbhe ++;
			}
			else {
				push(@ho,$pid);		
				push(@ho_ratio,@values);
				$nbho ++;
			}
		}
		
		
		my  $rocksdb_id = chunks::return_rocks_id_from_genbo_id($vid);
	

		my ($s,$all) = split("!",$rocksdb_id);
		$s = $s +0;
		my $vhh;
		$vhh->{chromosome} = $chr->ucsc_name;
		$vhh->{allele} = $all;
		$vhh->{start} = int($s);
		$vhh->{end} = int($s)+length($b);
		$vhh->{rocksid} = $rocksdb_id;
		$lift->add_region($vhh);
		$vhh->{max_ratio} = $max_ratio;
		$vhh->{max_dp} = $max_dp;
		$vhh->{transmissions} = $bit_models;
		$vhh->{patients} =  pack("w*",@he,@ho);;
		$vhh->{he} = scalar(@he);
		$vhh->{ho} = scalar(@ho);
		$vhh->{value} = pack("w*",@he_ratio,@ho_ratio);;
		
		$snps->{$rocksdb_id} = dclone($vhh);
		$cpt++;
		if ($cpt %100000 ==0){
			warn $chr->name."::"."::".$cpt. "time ".	abs(time - $time);
			$time = time;
		}
	}
	my $hash = $lift->liftOver_regions($project->name.".".$chr->name) ;
	my $tab;
	foreach my $a (values %$hash){
		foreach my $vl (@$a){
			my $rid = $vl->{rocksid};
			$snps->{$rid}->{LIFT} =  delete $vl->{LIFT};
		}
	}
	my $file = dejavu_duckdb::save_csv($chr,$snps,$dir_tmp);
	$no->close();
	return $file;
}



sub save_and_lift_rocks {
	my($chr,$lift) =@_;
	
	my $diro = $project->rocks_directory();
	my $final_polyviewer_all = GenBoNoSqlRocks->new(dir=>$diro,mode=>"r",name=>"polyviewer_objects");
	my $vectors = dejavu_duckdb::get_hash_model_variant($chr);
	
	my $snps;
	my $cpt;
	my $time = time;
	for (my $i =0;$i<$chr->size_vector();$i++) {
	my $v = $final_polyviewer_all->get( $chr->name."!".$i);
		my @ho;
		my @he;
		my @ho_ratio;
		my @he_ratio;
		
		my @models_ho;
		my @models_he;
		my $nbhe = 0;
		my $nbho = 0;
		my $max_dp = 0;
		my $max_ratio = 0;
		my $bit_models = 0 ;
		foreach my $pid (keys %{$v->{patients_calling}}){
			my @values;
			my $dp = $v->{patients_calling}->{$pid}->{dp};
			my $ratio = $v->{patients_calling}->{$pid}->{pc};
			my $alt = int($dp*($ratio/100));
			$max_dp = $dp if $dp > $max_dp;
			$max_ratio = int($ratio) if int($ratio) > $max_ratio;
			my $model= dejavu_duckdb::find_variant_model($vectors, int($i), $pid);
			$bit_models = $bit_models | $model;

			push(@values,$dp);
			push(@values,$alt);
			push(@values,$model);
			if ( $v->{patients_calling}->{$pid}->{gt} eq 'he'){
				push(@he,$pid);		
				push(@he_ratio,@values);
				$nbhe ++;
			}
			else {
				push(@ho,$pid);		
				push(@ho_ratio,@values);
				$nbho ++;
			}
		}
		my  $rocksdb_id = chunks::return_rocks_id_from_genbo_id($v->{id});
		if ($rocksdb_id =~ /0000000000/){
			warn "********* ".$rocksdb_id." ********* ".$v->{id};
			die();
		}

		my ($s,$all) = split("!",$rocksdb_id);
		$s = $s +0;
		my $vhh;
		$vhh->{chromosome} = $chr->ucsc_name;
		$vhh->{allele} = $all;
		$vhh->{start} = int($s);
		$vhh->{end} = int($s)+length($b);
		$vhh->{rocksid} = $rocksdb_id;
		$lift->add_region($vhh);
		$vhh->{max_ratio} = $max_ratio;
		$vhh->{max_dp} = $max_dp;
		$vhh->{transmissions} = $bit_models;
		$vhh->{patients} =  pack("w*",@he,@ho);;
		$vhh->{he} = scalar(@he);
		$vhh->{ho} = scalar(@ho);
		$vhh->{value} = pack("w*",@he_ratio,@ho_ratio);;
		
		$snps->{$rocksdb_id} = dclone($vhh);
		$cpt++;
		if ($cpt %100000 ==0){
			warn $chr->name."::"."::".$cpt. "time ".	abs(time - $time);
			$time = time;
		}
	}
	my $hash = $lift->liftOver_regions($project->name.".".$chr->name) ;
	my $tab;
	foreach my $a (values %$hash){
		foreach my $vl (@$a){
			my $rid = $vl->{rocksid};
			$snps->{$rid}->{LIFT} =  delete $vl->{LIFT};
		}
	}
	my $file = dejavu_duckdb::save_csv($chr,$snps,$dir_tmp);
	$final_polyviewer_all->close();
	return $file;
}
	


	
	




exit(0);








