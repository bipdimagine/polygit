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
use File::Slurp qw(write_file);
use List::Util qw/shuffle/;
use lib "$RealBin/../../../utility";
use liftOverRegions;
use chunks;
use Sereal qw(sereal_encode_with_object sereal_decode_with_object write_sereal_file read_sereal_file);
use GenBoNoSqlDuckDBDejaVu;
use Text::CSV;
use MIME::Base64;
use GenBoNoSqlRocks;

my $fork = 1;


my $buffer1 = new GBuffer;
my $dbh = $buffer1->dbh();
my $sql = qq{
   SELECT pr.name FROM PolyprojectNGS.run_machine as rm ,PolyprojectNGS.patient as p,PolyprojectNGS.projects as pr where rm.machine_id=40 and p.run_id=rm.run_id and pr.project_id=p.project_id group by p.project_id;
};
my $sth = $dbh->prepare($sql);

$sth->execute();
my $hash;
while (my ($project_name) = $sth->fetchrow_array) {
	my $buffer = new GBuffer();
my $project = $buffer->newProjectCache( -name => "$project_name");
my $parquet = $project->parquet_cache_variants();
next unless -e $parquet;

foreach my $patient (@{$project->getPatients}){
	
	next unless $patient->getRun->isRevio();
	warn $patient->name;
	warn  $patient->getRun->machine();
	my $pcol = "patient_".$patient->id."_type";
	my $sql = qq{

    copy (

        select variant_rocksdb_id

        from '$parquet'

        where $pcol <> 0 and variant_gnomad_ac < 1000 and variant_other_patients < 50

    ) to stdout (format csv, header false)

};

open(my $fh, "-|", "duckdb", "-c", $sql)

    or die "duckdb failed";

while (my $line = <$fh>) {

    chomp $line;

    $hash->{$line}++;

}

close($fh);
}
}
warn "start write  !!! ";

my $dir = "/data-pure/public-data/dejavu/HG38/revio/";
system("mkdir $dir") unless -e $dir;
my $dir_tmp = "/data-beegfs/tmp/";
my $filename = $dir_tmp."/revio.tmp.csv";

	my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
	
	my $fh;		
	open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
	$csv->print($fh, ["rocksdb_id","nb",]);
	
my $rocks =  GenBoNoSqlRocks->new(dir=>$dir,mode=>"c",name=>"rocks.dejavu",pipeline=>1);
my $nb = scalar keys %$hash;
my $aa =0;
foreach my $k (keys %$hash){
	warn $aa."/".$nb if $aa % 1000000 == 0;
	$aa++;
	next if $hash->{$k} == 1;
	$rocks->put_raw($k,$hash->{$k});
	$csv->print($fh,[$k,$hash->{$k}]);
}
$rocks->write_batch();
$rocks->close;
close($fh);
warn "end rocks";
my $parquet_file = $dir."/revio.parquet";
my $query = "
	COPY (
        SELECT * from read_csv_auto(['$filename']) order by rocksdb_id
    )
    TO '$parquet_file' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE);";
    my $exit_code = system("duckdb", ":memory:", "-c", $query);
die("duckdb failed with exit code: $exit_code") if $exit_code != 0;
warn $parquet_file;
warn "end";
exit(0);


