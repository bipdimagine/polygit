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

my $fork = 3;
my $infile_snp;
my $infile_indels;
my $outfile;
my $tmp_dir;
GetOptions(
	'in_snp=s' => \$infile_snp,
	'in_indels=s' => \$infile_indels,
	'out=s' => \$outfile,
	'tmp_dir=s' => \$tmp_dir,
	'fork=s' => \$fork,
);

die if not $infile_snp;
die if not $infile_indels;
die if not $outfile;
die if not $tmp_dir;

print "\n----- CADD ----\n";
print "1- parse files\n";

my @lFiles;
push(@lFiles, $infile_snp);
push(@lFiles, $infile_indels);

my @lChr = (1..22, 'X', 'Y', 'M');
my $h_tabix_chr;

my $buffer_init = new GBuffer;
my $project_name = $buffer_init->getRandomProjectName('HG38_DRAGEN');
my $project_init = $buffer_init->newProject( -name => $project_name );
foreach my $chr_id (@lChr) {
	my $chr = $project_init->getChromosome($chr_id);
	my $start = 1;
	my $end = 5000000;
	while ($end < $chr->length()) {
		my $region = $chr_id.'-'.$start.'-'.$end;
		push(@{$h_tabix_chr->{$chr_id}}, $region);
		$start += 5000000;
		$end += 5000000;
	}
}
$project_init = undef;
$buffer_init = undef;

foreach my $chr_id (@lChr) {
	next if not exists $h_tabix_chr->{$chr_id};
	my @l_tabix = @{$h_tabix_chr->{$chr_id}};
	
	MCE::Loop->init(
	   max_workers => $fork,
	   chunk_size => 1,
	   gather => sub {
	        my ($data) = @_;
	        print $data->{ok} ." DONE\n";
	    }
	);
	mce_loop {
		my ($mce, $chunk_ref, $index) = @_;
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
				$csv->print($fh, ["chr","pos","rocksdb_id","score"]); 
				foreach my $file (@lFiles) {
					my $v1 = Bio::DB::HTS::Tabix->new( filename => $file );
					my @ltmp = split('-',$chr_part);
					my $region = $ltmp[0].':'.$ltmp[1].'-'.$ltmp[2];
					my $iter = $v1->query($region);
					next unless $iter;
					while ( my $row = $iter->next ) {
						my ($this_chr_id, $pos, $ref, $alt, $raw, $score) = split(' ', $row);
						my $rocksid;
						if (length($ref) != length($alt)) {
							$pos++;
							my $var = $project->_newVariant($this_chr_id.'_'.$pos.'_'.$ref.'_'.$alt);
							$rocksid = $var->genomic_rocksdb_id();
							$var = undef;
						}
						else {
							$rocksid = $this_chr_id.'!'.sprintf("%010d", $pos).'!'.$alt;
						}
						$csv->print($fh,['chr'.$this_chr_id, $pos, $rocksid, int($score)]);
					}
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
}

print "\n----- CADD ----\n";
print "2- parquet file\n";

my @lfiles_all;
foreach my $chr_id (@lChr) {
	my @lfiles;
	my @l_tabix = @{$h_tabix_chr->{$chr_id}};
	foreach my $region (@l_tabix) {
		my @ltmp = split('-', $region);
		next if $ltmp[0] ne $chr_id;
		my $filename = $tmp_dir.'/part.'.$region.'.csv';
		push(@lfiles, $filename);
		push(@lfiles_all, $filename);
	}
	
	my $filenames = join("','", @lfiles);
	my $parquet_file_chr = $tmp_dir.'/../parquet.TEST/cadd.'.$chr_id.'.parquet';
	my $query = qq{
	COPY (
        SELECT * from read_csv_auto([\'$filenames\']) order by chr,pos
    )
    TO '$parquet_file_chr' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE);};
	my $cmd = "duckdb -c \"$query\"";
	system($cmd);
}
#my $filenames = join("','", @lfiles_all);
#my $parquet_file_chr = $tmp_dir.'/../parquet.TEST/cadd.parquet';
#my $query = qq{
#COPY (
#       SELECT * from read_csv_auto([\'$filenames\']) order by chr,pos
#   )
#TO '$parquet_file_chr' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE);};
#my $cmd = "duckdb -c \"$query\"";
#system($cmd);
