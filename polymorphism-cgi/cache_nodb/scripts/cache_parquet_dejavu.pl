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
use lib "$RealBin/../../../GenBo/lib/obj-nodb/packages";
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
use lib "$RealBin/../../../polypipeline/scripts/scripts_pipeline/dejavu/utility/";
use lib "$RealBin/../../../polypipeline/scripts/scripts_pipeline/dejavu/hg38/variant/duckdb/";
use liftOverRegions;
use chunks;
use Sereal qw(sereal_encode_with_object sereal_decode_with_object write_sereal_file read_sereal_file);
#use GenBoNoSqlDuckDBDejaVu;
use Text::CSV;
use MIME::Base64;
use dejavu_duckdb;
use MCE::Loop;
use MCE::Flow;
use Storable qw(dclone);
my $fork = 1;
my $project_name;
my $force;
GetOptions(
	'fork=s' => \$fork,
	'project=s' => \$project_name,
	'force=s' => \$force,
);




my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => "$project_name");
$project->genome_version;
 
my $dir_tmp = $project->lmdb_pipeline_dir();
my $data_lift;	

my @chromosomes = shuffle (@{$project->getChromosomes});
#my $dir_parquet = "/data-beegfs/parquet_HG19_HG38/";

my $pm2 = new Parallel::ForkManager($fork);

my $dir_parquet = $buffer->dejavu_parquet_dir();
my $parquet_file = $dir_parquet."/".$project->name.".".$project->id.".parquet";

my $can_dejavu = 1;
foreach my $pname (@{$buffer->getQuery->listProjectsWithoutDejaVu()}) {
	next if $pname ne $project_name;
	$can_dejavu = undef;
	last;
}

if (not $can_dejavu) { $parquet_file .= '.no_dejavu'; }

exit(0) if -e $parquet_file and not $force;


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
   	my $file  = HG19_HG38($chr,$lift);
    push(@$files,$file);
   	 
   }
	 MCE->gather($files);
} @chromosomes;	


#my $final;
# $pm2->run_on_finish(
#	sub {
#		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;
#		
#		unless (defined($hRes) or $exit_code > 0) {
#			print qq|No message received from child process $exit_code $pid!\n|;
#			return;
#		}
#	
#	}
#);
#
#
#foreach my $chr (@chromosomes){	
#	my $pid = $pm2->start and next;
#	foreach my $patient (@{$project->getPatients()}) { $patient->setOrigin($chr); }
#	my $lift = liftOverRegions->new(project=>$project,version=>$project->lift_genome_version);
#	 my $array = HG19_HG38($chr,$lift);
#	 $pm2->finish( 0, {lift_variants=>$array});
#}
#$pm2->wait_all_children();	
#warn "end";


my @files;
foreach my $t (@results) {
	push(@files,@$t);
}
@files = sort {
    ($b =~ /_X\.csv$/) <=> ($a =~ /_X\.csv$/)  # "_X" en premier
    || $a cmp $b                                # Trie alphabétique normal après
} @files;
my $filename = join(",",@files);


my $dbh = DBI->connect("dbi:ODBC:Driver=DuckDB;Database=:memory:", "", "", { RaiseError => 1 , AutoCommit => 1});
my $query = "
	COPY (
        SELECT * from read_csv_auto([$filename]) order by chr38,pos38,chr19,pos19,allele
    )
    TO '$parquet_file' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE,ROW_GROUP_SIZE 1000000);";
$dbh->do($query);
$filename =~s/,/ /g;
system("rm $filename");


print "\nDONE!\n";
exit(0);

sub open_csv {
	my ($project,$chr) = @_;
	my $dir_tmp_cvs = "/tmp/pipeline/";
	my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
	my $filename = "$dir_tmp_cvs/".$project->name."_".$chr->name.".csv";
	warn $filename;
	my $fh;		
	open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
	$csv->print($fh, ["project","pos38","pos19","allele","max_ratio","max_dp","transmissions","he","ho","value"]); 
	return($filename,$fh,$csv);
}



sub HG19_HG38 {
	my($chr,$lift) =@_;
	warn $chr->name." --- ";
	
	my @lVarIds = @{$chr->getListVarIds($chr->getVariantsVector())};
	
	my $snps;
	my $lt = time;
	
	my $t =time;
	my $snps;
	my $variations;
	my $hregions;
	my $cpt =0;
	
	my $lmdb_dir = $project->lmdb_cache_variations_dir();
	my $lmdb = GenBoNoSqlLmdb->new(
		dir         => $lmdb_dir,
		mode        => 'r',
		is_index    => 1,
		name        => $chr->name,
		is_compress => 1,
		vmtouch     => 1,
	);
	
	my $time = time;
	my $vectors = dejavu_duckdb::get_hash_model_variant($chr);
	foreach my $vid (@lVarIds){
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
			
			
			$max_dp = $dp if $dp > $max_dp;
			$max_ratio = int($ratio) if int($ratio) > $max_ratio;
			my $model= dejavu_duckdb::find_variant_model($vectors, int($index_lmdb), $pid);
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
		if ($rocksdb_id =~ /0000000000/){
			warn "********* ".$rocksdb_id." ********* ".$vid;
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
	return $file;
}




exit(0);
