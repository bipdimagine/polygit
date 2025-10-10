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
use liftOver;
use Time::HiRes;
use Sereal qw(sereal_encode_with_object sereal_decode_with_object write_sereal_file read_sereal_file);
use Sereal::Decoder qw(decode_sereal sereal_decode_with_object scalar_looks_like_sereal);
use Archive::Tar;
use List::Util qw(
  shuffle
);
use lib "$RealBin/../../../utility";
use chunks;
use File::Temp qw(tempdir);
use Compress::Snappy;
use MCE::Loop;
use MCE::Flow;
use GenBoNoSqlDuckDBDejaVu;
use MIME::Base64;



my $fork = 10;
my $project_name;
my $chr;
my $version;
GetOptions(
	'chr=s'		   => \$chr,
	'fork=s'		   => \$fork,
	'version=s'   => \$version,
);
my $buffer = new GBuffer;
die ("version ?") unless $version;
<<<<<<< HEAD
 
my $dir_final = "/data-beegfs/tmp/rocks-".$version.'.pn/';
#my $dir_final = '/data-isilon/DejaVu/'.uc($version).'/variations/';
=======

my $dir_final = $buffer->config_path("root","dejavu").'/'.uc($version).'/variations/rocks/';
>>>>>>> branch 'duckdb' of https://github.com/bipdimagine/polygit.git

 $| = 1;
my $rg38 = GenBoNoSqlRocksGenome->new(chunk_size=>10_000_000,dir=>$dir_final,mode=>"c",index=>"genomic",chromosome=>$chr,genome=>$version,pack=>"",description=>[]);
my $rg38_1 = GenBoNoSqlRocksGenome->new(chunk_size=>10_000_000,dir=>$dir_final,mode=>"c",index=>"genomic",chromosome=>$chr,genome=>$version,pack=>"",description=>[]);
my $regionss = $rg38->regions();
#save_regions($regionss->[0]);
#die();
MCE::Loop->init(
   max_workers => $fork, chunk_size => 'auto'
);
if (scalar(@$regionss) ==1){
		save_regions($regionss->[0]);
		exit(0);
}
mce_loop {
   my ($mce, $chunk_ref, $chunk_id) = @_;
   if (ref($chunk_ref) ne "ARRAY") {
   	save_regions($chunk_ref);
   }
   else {
   foreach my $region (@$chunk_ref){
		save_regions($region)
	}
   }
   	
} sort{$a->{start} <=> $b->{start}} @{$regionss};


exit(0);



sub save_regions {
	my ($region) = @_;
	my $dir = $buffer->dejavu_parquet_dir();
#	my $dir = "/data-beegfs/projects.parquet/";#$buffer->deja_vu_project_dir($version);
 	my %hash;
 	my $start = $region->{start};
 	my $end = $region->{end};
 	my %pr;
 	my $nb ;
 	my $file1  = "/tmp/pipeline/".$chr.".".$region->{id}.".csv";
 	
 	my $ff= "$dir/NGS*.parquet";
 	
 	my $start_string = $start;
 	my $end_string =$end;
 	my $sql_file = "/tmp/".$region->{id}."sql";
 	my $string_chr = $chr."!";
	my $sql = qq {PRAGMA threads=4;COPY (
          SELECT pos38,allele, STRING_AGG(project || ';' || he || ';' || ho, ';') AS value
          FROM read_parquet('$ff') where chr38 = '$chr' and  pos38 between '$start_string' and '$end_string' 
          GROUP BY pos38,allele order by pos38 
      ) TO '$file1' (FORMAT 'csv', COMPRESSION 'ZSTD');
 	};

	if ($version eq "HG19"){
		$sql = qq {PRAGMA threads=8;COPY (
          SELECT pos19,allele, STRING_AGG(project || ';' || he || ';' || ho, ';') AS value
          FROM read_parquet('$ff') where chr19 = '$chr' and  pos19 between '$start_string' and '$end_string'
          GROUP BY pos19,allele  
      ) TO '$file1' (FORMAT 'csv', COMPRESSION 'ZSTD');
 	}
	}

 open(my $fh, '>', $sql_file) or die "Impossible d'ouvrir $sql_file: $!";
print $fh $sql;
close($fh);
# Exécuter DuckDB avec le fichier SQL
system("/software/bin/duckdb < $sql_file");
unlink $sql_file;

#unlink $file1;
#return;
my $rocks =  $rg38->nosql_rocks_tmp($region);
open(CSV ,"zstdcat $file1 | ")  or die "Impossible d'ouvrir $file1 : $!";;
my $xx;
$rocks->put_batch_raw("xx","zz");



my $sql = q{
    SELECT p2.phenotype_id, p2.project_id
    FROM PolyPhenotypeDB.phenotype AS p1
    JOIN PolyPhenotypeDB.phenotype_project AS p2
      ON p1.phenotype_id = p2.phenotype_id
};

# Préparation et exécution
my $dbh = $buffer->dbh;
my $sth = $dbh->prepare($sql);
$sth->execute();

# Création de la table de hash : project_id ? [ phenotype_id list ]
my $project_to_pheno;
my $pheno;
while (my $row = $sth->fetchrow_hashref) {
    my $project_id   = $row->{project_id};
    my $phenotype_id = $row->{phenotype_id};
    
   # $project_to_pheno->{$project_id}->{$phenotype_id} = 0;
    push @{ $project_to_pheno->{$project_id} }, $phenotype_id;
}

$sth->finish;
$dbh->disconnect;



while(my $line = <CSV>){
	$xx++;
	chomp($line);
	my($a,$c,$b) = split(",",$line);
	next if $c eq "allele";
	my @z = split(";",$b);
		my $pheno_final;
	for (my $i = 0; $i < @z; $i += 3) {
	 unless (exists $project_to_pheno->{$z[$i]}) {
	 		$pheno_final->{-1}->{project} ++;
			$pheno_final->{-1}->{he} += $z[$i+1];
			$pheno_final->{-1}->{h0} += $z[$i+2];
			next;
	 	}
	
		foreach my $pheno (@{$project_to_pheno->{$z[$i]}} ){
			$pheno_final->{$pheno}->{project} ++;
			$pheno_final->{$pheno}->{he} += $z[$i+1];
			$pheno_final->{$pheno}->{h0} += $z[$i+2];
		}
	}
	my @vs ; 
	foreach my $k (keys %$pheno_final){
		$pheno_final->{$k}->{ho} += 0;
		$k +=0;
		push(@vs,($k,$pheno_final->{$k}->{project},$pheno_final->{$k}->{he},$pheno_final->{$k}->{ho}));
	}
	#warn  $b unless @vs;
	my $v = pack("w*",@z);
	my $v2 = pack("w*",@vs);
	my ($chr,$pos) = split("!",$a);
	$rocks->put_batch_raw(sprintf("%010d", $a)."!".$c,$v);
}
close (CSV);
$rocks->close();
warn "end $xx ".$region->{id};
warn $file1;
unlink $file1;
return;

		
}

 