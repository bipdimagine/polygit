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
my $file;
my $dir1;
GetOptions(
	'chr=s'		   => \$chr,
	'fork=s'		   => \$fork,
	'version=s'   => \$version,
	'file=s' => \$file,
	'dir=s' => \$dir1,
);
my $buffer = new GBuffer;
die ("version ?") unless $version;

my $dir_final = $dir1."/$file/";


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
my $x;
while (my $row = $sth->fetchrow_hashref) {
    my $project_id   = $row->{project_id};
    my $phenotype_id = $row->{phenotype_id};
    $pheno->{$row->{phenotype_id}}++;
   # $project_to_pheno->{$project_id}->{$phenotype_id} = 0;
    push @{ $project_to_pheno->{$project_id} }, $phenotype_id;
    $x++ if $phenotype_id ==22;
}

#warn Dumper $pheno;
#warn scalar keys (%$pheno);
#die();

$sth->finish;
$dbh->disconnect;
 
#my $dir_final = '/data-isilon/DejaVu/'.uc($version).'/variations/';

system("mkdir -p $dir_final" ) unless -e $dir_final;
my $dir_final_pheno = $dir1."/$file.phenotype/";
system("mkdir -p $dir_final_pheno" ) unless -e $dir_final_pheno;
warn $dir_final_pheno;
 $| = 1;
my $rg38 = GenBoNoSqlRocksGenome->new(chunk_size=>10_000_000,dir=>$dir_final,mode=>"c",index=>"genomic",chromosome=>$chr,genome=>$version,pack=>"",description=>[]);
my  $r_pheno = GenBoNoSqlRocks->new(mode=>"c",pipeline=>1,dir=>$dir_final_pheno, name=>$chr);

my $regionss = $rg38->regions();

#copie la liste des projets dans le dejavu au temps T
my $dir = $buffer->dejavu_parquet_dir();

opendir(my $dh, $dir) or die "Impossible d'ouvrir $dir: $!";

while (my $file = readdir($dh)) {
    next unless $file =~ /\.([0-9]+)\.parquet/;  # cherche le motif
    $r_pheno->put_batch_raw($1,1);
}

closedir($dh);
#save_regions($regionss->[0]);
#die();
my $jobs = 0;
MCE::Loop->init(
   max_workers => $fork, chunk_size => 'auto',
    gather => sub {
        my ($mce, $data) = @_;
        foreach my $id (keys %$data){
        	$r_pheno->put_batch_raw($id,$data->{$id});
		}
    }
);
push (@{$regionss},  {none=>1});


mce_loop {
   my ($mce, $chunk_ref, $chunk_id) = @_;
   	warn "start ".$chunk_id;
   if (ref($chunk_ref) ne "ARRAY") {
   	confess();
   }
   else {
   foreach my $region (@$chunk_ref){
   		next if exists $region->{none};
		my $hash = save_regions($region);
		MCE->gather($chunk_id,$hash);
	}
   }
   	
} sort{$a->{start} <=> $b->{start}} @{$regionss};

MCE::Loop->finish;
$r_pheno->close();
confess($jobs." ".scalar(@{$regionss})) if $jobs != scalar(@{$regionss});

confess() if ($jobs < scalar(@{$regionss}));
warn $dir_final;
warn $dir_final_pheno;
exit(0);



sub save_regions {
	my ($region) = @_;
	my $dir = $buffer->dejavu_parquet_dir();
 	my %hash;
 	my $start = $region->{start};
 	my $end = $region->{end};
 	my %pr;
 	my $nb ;
 	
 	my $ff= "$dir/NGS*.parquet";
 	
 	my $start_string = $start;
 	my $end_string =$end;
 	my $sql_file = "/tmp/".$region->{id}."sql";
 	my $string_chr = $chr."!";
	my $sql = qq {PRAGMA threads=5;
          SELECT pos38,allele, STRING_AGG(project || ';' || he || ';' || ho, ';') AS value
          FROM read_parquet('$ff') where chr38 = '$chr' and  pos38 between '$start_string' and '$end_string' 
          GROUP BY pos38,allele order by pos38 ;
 	};

	if ($version eq "HG19"){
		$sql = qq {PRAGMA threads=5;
          SELECT pos19,allele, STRING_AGG(project || ';' || he || ';' || ho, ';') AS value
          FROM read_parquet('$ff') where chr19 = '$chr' and  pos19 between '$start_string' and '$end_string'
          GROUP BY pos19,allele; }; 
 	}
	
	
my $duckdb = $buffer->software("duckdb");


my $cmd = qq{$duckdb -csv -noheader -c "$sql"};
my $rocks =  $rg38->nosql_rocks_tmp($region);
open(CSV ,"-|", $cmd)  or die "Impossible d'ouvrir  : $!";;
my $xx;
$rocks->put_batch_raw("xx","zz");

my $hh;
my $nb_rocks=0;

while(my $line = <CSV>){
	$xx++;
	$nb_rocks ++;
	chomp($line);
	my($a,$c,$b) = split(",",$line);
	next if $c eq "allele";
	my @z = split(";",$b);
		my $pheno_final;
		my $ptotal;
	for (my $i = 0; $i < @z; $i += 3) {
		$ptotal->{project} ++;
		$ptotal->{he} += $z[$i+1];
		$ptotal->{ho} += $z[$i+2];
	 unless (exists $project_to_pheno->{$z[$i]}) {
	 		$pheno_final->{0}->{project} ++;
			$pheno_final->{0}->{he} += $z[$i+1];
			$pheno_final->{0}->{h0} += $z[$i+2];
			next;
	 	}
	else {
		my $pheno = join (";",@{$project_to_pheno->{$z[$i]}});;
		$hh->{$z[$i]} = 1; 
		$pheno_final->{$pheno}->{project} ++;
		$pheno_final->{$pheno}->{he} += $z[$i+1];
		$pheno_final->{$pheno}->{h0} += $z[$i+2];	
	
	}
	}
	my @vs ; 
	push(@vs,(99,$ptotal->{project},$ptotal->{he},$ptotal->{ho}));
	foreach my $k (keys %$pheno_final){
		$pheno_final->{$k}->{ho} += 0;
		$k +=0;
		push(@vs,($k,$pheno_final->{$k}->{project},$pheno_final->{$k}->{he},$pheno_final->{$k}->{ho}));
	}
	#warn  $b unless @vs;
	
	my $v = pack("w*",@z);
	my $pheno = shift (@vs);
	my $v2 = pack("w*",@vs);
	my ($chr,$pos) = split("!",$a);
	my $rid = sprintf("%010d", $a)."!".$c;
	#warn $rid;
	#19_2110764_G_A
	

	
	$hh->{sprintf("%010d", $a)."!".$c} = $v2;
	
	$rocks->put_batch_raw(sprintf("%010d", $a)."!".$c,$v);
	$rocks->write_batch() if $nb_rocks % 10000 == 0;
}
close (CSV);
$rocks->write_batch();
$rocks->close();
warn "end $xx ".$region->{id}." ".scalar keys %$hh;
return ($hh); 
}

 