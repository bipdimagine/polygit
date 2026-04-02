#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use strict; 
use GenBoNoSqlRocks;
use Data::Printer;
use Getopt::Long;
use Carp;
use Data::Dumper;
use Getopt::Long; 
use GBuffer;
use GenBoProject;
use Parallel::ForkManager;
use GenBoNoSql;
 use Time::ETA;
use List::MoreUtils qw(any uniq natatime);
use Bio::DB::HTS::Tabix;
use POSIX;
my $allsnps;
use Getopt::Long;
use Date::Tiny;
use Bio::DB::HTS::Faidx;
use Text::CSV;
use GenBoNoSqlRocksGenome;
use MCE::Loop;
use feature 'state';

my $fork = 3;
my $infile;
my $outfile;
my $tmp_dir;
GetOptions(
	'in=s' => \$infile,
	'out=s' => \$outfile,
	'tmp_dir=s' => \$tmp_dir,
	'fork=s' => \$fork,
);

die if not $infile;
die if not $outfile;
die if not $tmp_dir;

print "\n----- NCBOOST ----\n";
print "1- parse files\n";

my @lChr = (1..22, 'X', 'Y');
my @l_tabix;

my $buffer_init = new GBuffer;
my $project_name = $buffer_init->getRandomProjectName('HG38_DRAGEN');
my $project_init = $buffer_init->newProject( -name => $project_name );
foreach my $chr_id (@lChr) {
	my $chr = $project_init->getChromosome($chr_id);
	my $start = 1;
	my $end = 10000000;
	while ($end < $chr->length()) {
		my $region = $chr_id.'-'.$start.'-'.$end;
		push(@l_tabix, $region);
		$start += 10000000;
		$end += 10000000;
	}
}
$project_init = undef;
$buffer_init = undef;

MCE::Loop->init(
   max_workers => $fork,
   chunk_size => 'auto',
   gather => sub {
        my ($data) = @_;
        print $data->{ok} ." DONE\n";
    }
);
mce_loop {
	my ($mce, $chunk_ref, $index) = @_;
	my $v1 = Bio::DB::HTS::Tabix->new( filename => $infile );
	foreach my $chr_part (@$chunk_ref) {
		print '-> chr'.$chr_part."\n";
		my $filename = $tmp_dir.'/part.'.$chr_part.'.csv';
		if (not -e $filename) {
			my $buffer = new GBuffer;
			my $project_name = $buffer->getRandomProjectName('HG38_DRAGEN');
			my $project = $buffer->newProject( -name => $project_name );
			my $fh;	
			my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
			open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
			$csv->print($fh, ["chr","pos","rocksdb_id","score","gnomad_ac","gnomad_ho","cadd"]); 
			my @ltmp = split('-',$chr_part);
			my $region = $ltmp[0].':'.$ltmp[1].'-'.$ltmp[2];
			my $iter = $v1->query($region);
			while ( my $row = $iter->next ) {
				my @lCol = split(' ', $row);
				my $this_chr_id = $lCol[0];
				my $pos = $lCol[1];
				my $ref = $lCol[2];
				my $alt = $lCol[3];
				my $score = $lCol[-1];
				$score = sprintf("%.3f", $score);

				my ($genomic_rocksid, $rocksid);
				if (length($ref) != length($alt)) {
					$pos++;
					my $var = $project->_newVariant($this_chr_id.'_'.$pos.'_'.$ref.'_'.$alt);
					$genomic_rocksid = $var->genomic_rocksdb_id();
					$rocksid = $var->rocksdb_id();
					$var = undef;
				}
				else {
					$genomic_rocksid = $this_chr_id.'!'.sprintf("%010d", $pos).'!'.$alt;
					$rocksid = sprintf("%010d", $pos).'!'.$alt;
				}
				
				my ($gac, $gho);
				my $h_gad = $project->getChromosome($this_chr_id)->rocksdb('gnomad')->value($rocksid);
				if ($h_gad) {
					$gac = $h_gad->{ac};
					$gho = $h_gad->{ho};
				}
				else {
					$gac = 0;
					$gho = 0;
				}
				my $cadd = 0;
				my $h_cadd = $project->getChromosome($this_chr_id)->rocksdb("cadd")->value($rocksid);
				$cadd = $h_cadd->{cadd_score} if $h_cadd;
				
				$csv->print($fh,[$this_chr_id, $pos, $genomic_rocksid, $score, $gac, $gho, $cadd]);
			}
			close($fh);
			$project = undef;
			$buffer = undef;
		}
	}
	my $hres;
	$hres->{ok} = 'part ';
	MCE->gather($hres);
} @l_tabix;	

MCE::Loop->finish();

print "\n----- NCBOOST ----\n";
print "2- parquet file\n";

my @lfiles_all;
foreach my $chr_id (@lChr) {
	foreach my $region (@l_tabix) {
		my @ltmp = split('-', $region);
		next if $ltmp[0] ne $chr_id;
		my $filename = $tmp_dir.'/part.'.$region.'.csv';
		push(@lfiles_all, $filename);
	}
	
}


my $filenames = join("','", @lfiles_all);
my $parquet_file_chr = '/data-pure/public-data/repository/HG38/ncboost/20260301/parquet/global/';
if (not -e $parquet_file_chr) {
	my $query = qq{
		PRAGMA threads=20;
		PRAGMA memory_limit='125GB';
		COPY (
		    SELECT *
		    FROM read_csv(
		        ['$filenames'],
		        columns = {
		            'chr': 'VARCHAR',
		            'pos': 'INTEGER',
		            'rocksdb_id': 'VARCHAR',
		            'score': 'DOUBLE',
		            'gnomad_ac': 'INTEGER',
		            'gnomad_ho': 'INTEGER',
		            'cadd': 'INTEGER',
		        },
		        delim=',',
		        header=true
		    )
		)
		TO '$parquet_file_chr'
		(FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE, PARTITION_BY (chr));
	};
	my $cmd = "duckdb -c \"$query\"";
	system($cmd);
}


my $filename_out = '/data-pure/public-data/repository/HG38/ncboost/20260301/parquet/dejavu_filter/';
foreach my $chr_id (@lChr) {
	print '-> '.$chr_id."\n";
	my $dir = "/data-pure/public-data/repository/HG38/ncboost/20260301/parquet/dejavu_filter/chr\=$chr_id/";
	next if -e $dir.'/data_0.parquet';
	mkdir ($dir); 
	my $query_2 = qq{
		PRAGMA threads=20;
        PRAGMA memory_limit='125GB';
        COPY (
            WITH positions AS (
                SELECT chr38, pos38,
                       CAST(SUM(he) AS INTEGER) AS he_sum,
                       CAST(SUM(ho) AS INTEGER) AS ho_sum
                FROM read_parquet('/data-pure/public-data/dejavu/projects_parquet/NGS*.parquet')
                WHERE chr38 = '$chr_id'
                GROUP BY chr38, pos38
            )
            SELECT n.*, p.he_sum, p.ho_sum
            FROM read_parquet('$parquet_file_chr/chr=$chr_id/*.parquet') n
            INNER JOIN positions p
                ON n.chr = p.chr38 AND n.pos = p.pos38
        )
        TO '$dir/data_0.parquet'
        (FORMAT PARQUET, COMPRESSION ZSTD);
	};
	warn $query_2;
	my $cmd2 = "duckdb -c \"$query_2\"";
	system($cmd2);
}

