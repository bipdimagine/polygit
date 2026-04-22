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

#die if not $infile;
#die if not $outfile;
#die if not $tmp_dir;

print "\n----- NCBOOST ----\n";
print "1- parse files\n";

my @lChr = (1..22, 'X', 'Y');
my @l_all_tabix;
foreach my $chr_id (@lChr) {
	my @l_tabix;
	my $buffer_init = new GBuffer;
	my $project_name = $buffer_init->getRandomProjectName('HG38_DRAGEN');
	my $project_init = $buffer_init->newProject( -name => $project_name );
	my $chr = $project_init->getChromosome($chr_id);
	my $start = 1;
	my $end = 50000;
	while ($end < $chr->length()) {
		my $region = $chr_id.'-'.$start.'-'.$end;
		push(@l_tabix, $region);
		push(@l_all_tabix, $region);
#		warn $region;
		$start += 50000;
		$end += 50000;
	}
	
	next;
	my $pm = new Parallel::ForkManager($fork);
	my $chunk_size = $fork;
	while (my @chunk = splice(@l_tabix, 0, $chunk_size)) {
	    my $pid = $pm->start() and next;
	
	    my $buffer = new GBuffer;
		my $project_name = $buffer->getRandomProjectName('HG38_DRAGEN');
		my $project = $buffer->newProject( -name => $project_name );
	    my $v1 = Bio::DB::HTS::Tabix->new(filename => $infile);
	
		my %chr_cache;
	    foreach my $chr_part (@chunk) {
			my $filename = $tmp_dir.'/part.'.$chr_part.'.csv';
			next if -e $filename;
	        my ($chr_id, $start, $end) = split('-', $chr_part);
			my $fh;	
			my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
			open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
			$csv->print($fh, ["chr","pos","rocksdb_id","score","annotation","gnomad_ac","gnomad_ho","cadd","dejavu", "dejavu_ho"]); 
			my @ltmp = split('-',$chr_part);
	        
	        my $entry = $chr_cache{$chr_id} //= {
			    chr  => $project->getChromosome($chr_id),
			    gad  => $project->getChromosome($chr_id)->rocksdb('gnomad'),
			    cadd => $project->getChromosome($chr_id)->rocksdb('cadd'),
			    dejavu => $project->getChromosome($chr_id)->rocks_dejavu(),
			};
			
			my $chr_obj = $entry->{chr};
			my $no_gad  = $entry->{gad};
			my $no_cadd = $entry->{cadd};
			my $no_dv = $entry->{dejavu};
	
	        my $iter = $v1->query("$chr_id:$start-$end");
	        while ( my $row = $iter->next ) {
				my @lCol = split(' ', $row);
				my $this_chr_id = $lCol[0];
				my $pos = $lCol[1];
				my $ref = $lCol[2];
				my $alt = $lCol[3];
				my $score = $lCol[-1];
				$score = sprintf("%.3f", $score);
				my ($genomic_rocksid, $rocksid, @lCons);
				$pos++;
				$genomic_rocksid = $this_chr_id.'!'.sprintf("%010d", $pos).'!'.$alt;
				$rocksid = sprintf("%010d", $pos).'!'.$alt;
				my ($gac, $gho);
				my $h_gad = $no_gad->value($rocksid);
				if ($h_gad) {
					$gac = $h_gad->{ac};
					$gho = $h_gad->{ho};
				}
				else {
					$gac = 0;
					$gho = 0;
				}
				my $cadd = 0;
				my $h_cadd = $no_cadd->value($rocksid);
				$cadd = $h_cadd->{cadd_score} if $h_cadd;
#				if ($chr_obj->intergenic_intspan->contains($pos)) {
#				    $csv->print($fh, [$this_chr_id, $pos, $genomic_rocksid, $score, 'intergenic', $gac, $gho, $cadd]);
#				    next;
#				}
#				else {
					my $var = $project->_newVariant($this_chr_id.'_'.$pos.'_'.$ref.'_'.$alt);
					my $h = $no_dv->dejavu($var->rocksdb_id());
					my $nb_he = 0;
					my $nb_ho = 0;
					my @lProjects;
					foreach my $pid (keys %{$h}) {
						$nb_he += $h->{$pid}->{'he'};
						$nb_ho += $h->{$pid}->{'ho'};
					}
					my $nb_all = $nb_he + $nb_ho;
					my $var_annot = $var->variationTypeInterface();
					foreach my $this_annot (split(',', $var_annot)) {
						$this_annot =~ s/ /_/g;
						$csv->print($fh,[$this_chr_id, $pos, $genomic_rocksid, $score, lc($this_annot), $gac, $gho, $cadd, $nb_all, $nb_ho]);
					}
					$var = undef;
#				}
			}
			close($fh);
		}
		foreach my $e (values %chr_cache) {
		    $e->{gad}->close();
		    $e->{cadd}->close();
		    $e->{dejavu}->close();
		}
		$project->disconnect();
		$project = undef;
		$buffer = undef;
		$pm->finish();
	}
	$pm->wait_all_children;
	
	$project_init = undef;
	$buffer_init = undef;
}



#foreach my $chr_part (@l_tabix) {
#	my $pid = $pm->start() and next;
#	my $v1 = Bio::DB::HTS::Tabix->new( filename => $infile );
#	print '-> chr'.$chr_part."\n";
#	my $filename = $tmp_dir.'/part.'.$chr_part.'.csv';
#	if (not -e $filename) {
#		my $buffer = new GBuffer;
#		my $project_name = $buffer->getRandomProjectName('HG38_DRAGEN');
#		my $project = $buffer->newProject( -name => $project_name );
#		my $fh;	
#		my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
#		open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
#		$csv->print($fh, ["chr","pos","rocksdb_id","score","annotation","gnomad_ac","gnomad_ho","cadd"]); 
#		my @ltmp = split('-',$chr_part);
#		my $region = $ltmp[0].':'.$ltmp[1].'-'.$ltmp[2];
#		my $iter = $v1->query($region);
#		my $no_gad = $project->getChromosome($ltmp[0])->rocksdb('gnomad');
#		my $no_cadd = $project->getChromosome($ltmp[0])->rocksdb('cadd');
#		while ( my $row = $iter->next ) {
#			my @lCol = split(' ', $row);
#			my $this_chr_id = $lCol[0];
#			my $pos = $lCol[1];
#			my $ref = $lCol[2];
#			my $alt = $lCol[3];
#			my $score = $lCol[-1];
#			$score = sprintf("%.3f", $score);
#
#			my ($genomic_rocksid, $rocksid, @lCons);
#			$pos++;
#			$genomic_rocksid = $this_chr_id.'!'.sprintf("%010d", $pos).'!'.$alt;
#			$rocksid = sprintf("%010d", $pos).'!'.$alt;
#			if ($project->getChromosome($this_chr_id)->intergenic_intspan->contains($pos)) {
#				push(@lCons, 'intergenic');
#			}
#			else {
#				my $var = $project->_newVariant($this_chr_id.'_'.$pos.'_'.$ref.'_'.$alt);
#				my $var_annot = $var->variationTypeInterface();
#				foreach my $this_annot (split(',', $var_annot)) {
#					$this_annot =~ s/ /_/g;
#					push(@lCons, lc($this_annot));
#				}
#				$var = undef;
#			}
#			
#			my ($gac, $gho);
#			my $h_gad = $no_gad->value($rocksid);
#			if ($h_gad) {
#				$gac = $h_gad->{ac};
#				$gho = $h_gad->{ho};
#			}
#			else {
#				$gac = 0;
#				$gho = 0;
#			}
#			my $cadd = 0;
#			my $h_cadd = $no_cadd->value($rocksid);
#			$cadd = $h_cadd->{cadd_score} if $h_cadd;
#			foreach my $cons (@lCons) {
#				$csv->print($fh,[$this_chr_id, $pos, $genomic_rocksid, $score, $cons, $gac, $gho, $cadd]);
#			}
#		}
#		close($fh);
#		$no_gad->close();
#		$no_cadd->close();
#		$project->disconnect();
#		$project = undef;
#		$buffer = undef;
#	}
#	$pm->finish();
#}
#$pm->wait_all_children;


print "\n----- NCBOOST ----\n";
print "2- parquet file\n";


my $batch = RocksDB::WriteBatch->new();




my $parquet_file_chr = '/data-pure/public-data/repository/HG38/ncboost/20260301/parquet-TEST/global/';
#foreach my $chr_id (@lChr) {
#	my @lfiles_all;
#	foreach my $region (@l_all_tabix) {
#		my @ltmp = split('-', $region);
#		next if $ltmp[0] ne $chr_id;
#		my $filename = $tmp_dir.'/chr'.$chr_id.'/part.'.$region.'.csv';
#		next if not -e $filename;
#		push(@lfiles_all, $filename);
#		warn $filename;
#	}
#	my $filenames = $tmp_dir.'/chr'.$chr_id.'/part.*.csv';
#	#my $filenames = join("','", @lfiles_all);
#	my $dir_chr = "$parquet_file_chr/chr\=$chr_id/";
#	mkdir $dir_chr if not -d $dir_chr;
#	my $query = qq{
#		PRAGMA threads=20;
#		PRAGMA memory_limit='125GB';
#		COPY (
#		    SELECT 
#				*
#		    FROM read_csv(
#		        ['$filenames'],
#		        columns = {
#		            'chr': 'VARCHAR',
#		            'pos': 'INTEGER',
#		            'rocksdb_id': 'VARCHAR',
#		            'score': 'DOUBLE',
#		            'annotation': 'VARCHAR',
#		            'gnomad_ac': 'INTEGER',
#		            'gnomad_ho': 'INTEGER',
#		            'cadd': 'INTEGER',
#		            'dejavu': 'INTEGER',
#		            'dejavu_ho': 'INTEGER',
#		        },
#		        delim=',',
#		        header=true,
#		        ignore_errors=true
#		    )
#		)
#		TO '$dir_chr/data_0.parquet'
#		(FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE);
#	};
#	my $cmd = "duckdb -c \"$query\"";
#	warn $cmd;
#	system($cmd);
#}
#
#die;

print "\n----- NCBOOST ----\n";
print "3- parquet filter file\n";

my $filename_out = '/data-pure/public-data/repository/HG38/ncboost/20260301/parquet-TEST/dejavu_filter/';
foreach my $chr_id (@lChr) {
#	print '-> '.$chr_id."\n";
#	my $dir = "/data-pure/public-data/repository/HG38/ncboost/20260301/parquet-TEST/dejavu_filter/chr\=$chr_id/";
#	next if -e $dir.'/data_0.parquet';
#	mkdir ($dir); 
#	my $query_2 = qq{
#		PRAGMA threads=20;
#        PRAGMA memory_limit='125GB';
#        COPY (
#            WITH positions AS (
#                SELECT chr38, pos38
#                FROM read_parquet('/data-pure/public-data/dejavu/projects_parquet/NGS*.parquet')
#                WHERE chr38 = '$chr_id'
#                GROUP BY chr38, pos38
#            )
#            SELECT n.*
#            FROM read_parquet('$parquet_file_chr/chr=$chr_id/*.parquet') n
#            INNER JOIN positions p
#                ON n.chr = p.chr38 AND n.pos = p.pos38
#            WHERE score>=0.8 and dejavu < 100 and dejavu_ho < 100 and gnomad_ac < 100
#        )
#        TO '$dir/data_0.parquet'
#        (FORMAT PARQUET, COMPRESSION ZSTD);
#	};
#	warn $query_2;
#	my $cmd2 = "duckdb -c \"$query_2\"";
#	system($cmd2);
}

print "\n----- NCBOOST ----\n";
print "4- rocks\n";
my $pm = new Parallel::ForkManager(25);
foreach my $chr_id (@lChr) {
	my $pid = $pm->start and next;
	#next if $chr_id ne '21';
	warn 'create rocks for chr'.$chr_id;
	my $parquet = $parquet_file_chr.'/chr='.$chr_id.'/data_0.parquet';
	next if not -e $parquet;
	my $rocks_path = '/data-pure/public-data/repository/HG38/ncboost/20260301/rocksdb/';
	my $pack = "C";
	my $no = GenBoNoSqlRocksAnnotation->new(dir=>$rocks_path,mode=>"c",name=>"$chr_id",version=>'20260301',pipeline=>1);
	my $sql = "PRAGMA threads=1;";
	$sql .= "select rocksdb_id, score from '$parquet';";
	my $b = new GBuffer;
	my $duckdb = $b->software('duckdb');
	open(my $fh, "-|", "$duckdb -csv -c \"$sql\"") or die "duckdb failed";
	while (my $line = <$fh>) {
	    chomp $line;
	    my ($rocksid, $score) = split(',', $line);
	    next if $rocksid eq 'rocksdb_id';
	    my ($chr_t, $pos_t, $alt_t) = split('!', $rocksid);
		$no->put_batch_raw($pos_t.'!'.$alt_t, $score);
	}
	close($fh);
	sleep(5);
	$no->write_batch();
	$no->close();
	$pm->finish();
}
$pm->wait_all_children();

