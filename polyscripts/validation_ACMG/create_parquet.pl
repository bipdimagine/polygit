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

 my $config_dir =  "../../GenBoConfig/genbo/";

		
# Exemple de connexion DuckDB (à adapter si tu utilises un fichier .duckdb)

 
 
 my $filename = $config_dir."paths.cfg";
confess($filename) unless -e $filename;
read_config $filename => my %config1;

my $d1 = $config1{root}->{public_data}."/repository/HG38/local_validation/";
my $parquet_file_final = $d1."local.parquet";
#perl monscript.pl | duckdb -c "SELECT * FROM read_csv_auto('/dev/stdin');"
 my $cmd = qq{perl $Bin/transform_in_parquet.pl | duckdb -c \"COPY (SELECT * FROM read_csv_auto('/dev/stdin')) TO \'$parquet_file_final\' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE);\"};
system("$cmd");
warn $d1;