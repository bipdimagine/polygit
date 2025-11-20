#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use strict;
use DBD::mysql;
use Config::Std;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/obj-nodb/packages/";
use rocksid;
use Text::CSV;
use File::Temp;
use Getopt::Long;

my @versions =("HG38","HG19");

foreach my $version (@versions){

# Exemple de connexion DuckDB (à adapter si tu utilises un fichier .duckdb)
 my $config_dir =  "../../GenBoConfig/genbo/";
my $filename2 =  $config_dir."genbo.cfg";
		my %config1;
		die($filename2) unless -e $filename2;
		read_config $filename2 =>  %config1;
		my $config = $config1{polyprod};
	my $ip = $config->{ip};
	my $db_user_name = $config->{user};
	my $db_password = $config->{pw}."";
	my $port = $config->{port};
	$port = 3306 unless $port;


#directiory 
		
# Exemple de connexion DuckDB (à adapter si tu utilises un fichier .duckdb)

 
 
 my $filename1 = $config_dir."paths.cfg";
confess($filename1) unless -e $filename1;
read_config $filename1 => my %config2;

my $d1 = $config2{root}->{public_data}."/repository/$version/local_validation/";
my $parquet_file_final = $d1."local.parquet";


my $fht = File::Temp->new(
    	DIR    => $config2{root}->{tmp},     # Optionnel : répertoire pour le fichier
    	SUFFIX => '.csv'      # Optionnel : suffixe du fichier
	);
	

	
	my $dsn = "DBI:mysql::$ip;port=$port\n";

my $validation =  {
	"Pathogenic" => 5,
	"Likely pathogenic" => 4,
	"Uncertain significance" => 3,
	"Likely benign" => 2,
	"Benign" => 1,
	"False Positive" => -1,
	"ToDo" =>- 3,
	"Unknown" =>- 0,
};
my %rev_validation = reverse %$validation;

	 	my $dbh = DBI->connect($dsn, $db_user_name, "$db_password")|| die "Database connection not made: $DBI::errstr";
my $sql = qq{
   SELECT v.polyid, va.validation,acmg.term as term
FROM validation_ACMG.variations AS v  JOIN validation_ACMG.variations AS v2 ON v.version=\'$version\' and v.uniq_id = v2.uniq_id 
join  validation_ACMG.validations AS va ON v2.variation_id = va.variation_id,
validation_ACMG.acmg as acmg
WHERE va.validation_id = (
    SELECT MAX(va.validation_id)
    FROM validation_ACMG.validations AS va2
    WHERE va2.variation_id = v2.variation_id
) and va.validation= acmg.idacmg;
};
warn $sql;
my $sth = $dbh->prepare($sql);
$sth->execute();

my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
	my $filename = $fht->filename;
#$dir_tmp_cvs."/".$chr->name.".".$project->name.".csv";#$dir_tmp_cvs/".$project->name."_".$chr->name.".csv";
	my $fh;		
	open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
# my $diro = $project->rocks_directory();
	
	
	$csv->print($fh, ["rocksdb_id","local_value","local_validation"]); 

while (my $row = $sth->fetchrow_hashref) {
	my $id =$row->{polyid};
	$id =~s/_/-/g; 
	my $rid = rocksid::return_genomic_rocks_id_from_gnomad_id($id);
	my $v = $row->{validation};
	#
	$csv->print($fh,[$rid,$v,$row->{term}]);
	my $type = $rev_validation{$v};
	#warn $rid." ".$v." ".$type;
	#$duckdb->do("INSERT INTO temp_data VALUES (\'$rid\', $v, \'$type\')");
}
close($fh);

$sth->finish;
$dbh->disconnect;

my $query = qq{
	COPY (
        SELECT * from read_csv_auto([\'$filename\']) order by rocksdb_id
    )
    TO '$parquet_file_final' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE);};
 my $cmd = "duckdb  -c \"$query\"";
 system($cmd);
 warn $cmd;
 }
 