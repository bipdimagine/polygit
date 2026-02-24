#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use strict; 
use Data::Dumper;
use Parallel::ForkManager;
use JSON::XS;
use GBuffer;
use Getopt::Long; 
use Carp;
use Text::CSV;
use JSON;
use GenBoNoSqlRocksGenome;
use Bio::DB::HTS::Tabix;
use List::MoreUtils qw(natatime);

my $only_chr;
my $fork = 1;
my $genome_version = "HG38";
my ($infile);
GetOptions(
	'file=s' => \$infile,
	'fork=s' => \$fork,
	'chr=s'	 => \$only_chr,
);

#/data-pure/public-data/repository/HG38/promoterAI/file/promoterAI_tss500.tsv.gz

confess("\n\nERROR: $infile not exists. Die\n\n") if not -e $infile;
#confess("\n\nERROR: chr option mandatory. Die\n\n") if not $only_chr;

my $buffer = new GBuffer;
my $project_name = $buffer->getRandomProjectName('HG38_DRAGEN');
my $project = $buffer->newProject( -name => $project_name );

my $global_dir = "/data-pure/public-data/repository/HG38/promoterAI/1.0.test/";
my $parquet_file_final = $global_dir."/parquet/promoterAI.parquet";
my $parquet_file_final_f = $global_dir."/parquet/promoterAI.only_promoter.parquet";
my $dir_out = $global_dir."/rocksdb/";

my $global_dir_19 = "/data-pure/public-data/repository/HG19/promoterAI/1.0.test/";
my $parquet_file_final_19 = $global_dir_19."/parquet/promoterAI.parquet";
my $parquet_file_final_f_19 = $global_dir_19."/parquet/promoterAI.only_promoter.parquet";
my $dir_out_19 = $global_dir_19."/rocksdb/";



my $pack = "C";

my $i = 0;
my $chr_old = '';


my $pm = new Parallel::ForkManager($fork);
my (@lCsv_hg38, @lCsv_hg38_filtre);
my @lHandler_hg38;

$pm->run_on_finish(
	sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data) = @_;
		push(@lCsv_hg38, $data->{hg38}->{file});
		push(@lCsv_hg38_filtre, $data->{hg38}->{file_f});
	}
);

foreach my $chrID (1..22, 'X', 'Y') {
	
	#next if $chrID ne '21';
	#next if $chrID ne $only_chr;
	
	my $pid = $pm->start and next;
	my $hres;
	
	my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
	my $filename = $global_dir.'/file/part.'.$chrID.'.csv';
	my $fh;	
	open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
	$csv->print($fh, ["chr","pos","ref","alt","rocksdb_id","gene","geneid","trid","strand","tss_pos","promoterAI"]); 
	
	my $csv_f = Text::CSV->new({ binary => 1, eol => "\n" });
	my $filename_f = $global_dir.'/file/part_filter.'.$chrID.'.csv';
	my $fh_f;	
	open( $fh_f, ">", $filename_f) or die "Impossible d'ouvrir $filename_f: $!";
	$csv_f->print($fh_f, ["rocksdb_id","gene","geneid","promoterAI"]); 
	
	my $no = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out,mode=>"c",name=>$chrID,pipeline=>1);
	my $no_lift = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out,mode=>"c",name=>$chrID.'.lift',pipeline=>1);
	
	my $h_liftover;
	my $tabix  = Bio::DB::HTS::Tabix->new( filename => $infile );
	my $res = $tabix->query( 'chr'.$chrID );
	my $i = 0;
	print '-> '.$chrID."\n";
	my $huids;
	my $last_uid;
	while ( my $line = $res->next ) {
		my ($chr_id, $pos, $ref, $alt, $gene, $gid, $trid_f, $strand, $tss_pos, $score) = split(' ', $line);
		next if $chr_id eq 'chrom';
		$chr_id =~ s/chr//;
		my ($trid, $o) = split('\.', $trid_f);
		my $rocks_id = $no->return_rocks_id($pos,$ref,$alt);
		my $uid = $chr_id."!".$rocks_id;
		push(@{$huids->{$uid}}, "$trid;$strand;$tss_pos;$score");
		
		my $chr_id_2 = $chr_id;
		$chr_id_2 =~ s/chr//;
		
		my $gene_infos = $gene.", ".$gid."_".$chr_id_2.", ".$trid.", ".$strand.", ".$tss_pos.", ".$score;
		push(@{$h_liftover->{$uid}->{gene_infos}}, $gene_infos);
		push(@{$h_liftover->{$uid}->{to_put}}, "$trid;$strand;$tss_pos;$score");
		$h_liftover->{$uid}->{score} = $score;
		$h_liftover->{$uid}->{alt} = $alt;
		$h_liftover->{$uid}->{ref} = $ref;
		$h_liftover->{$uid}->{gid} = $gid;
		$h_liftover->{$uid}->{gene} = $gene;
		$h_liftover->{$uid}->{uid} = $uid;
		
		$csv->print($fh,['chr'.$chr_id, $pos, $ref, $alt, $uid, $gene, $gid.'_'.$chr_id_2, $trid, $strand, $tss_pos, $score]);
		$csv_f->print($fh_f,[$uid, $gene, $gid.'_'.$chr_id_2, $score]) if (abs($score) >= 0.2);
		
		if ($last_uid and $last_uid ne $uid) {
			$no_lift->put_batch_raw($uid, encode_json $h_liftover->{$uid});
			my $values = join(',', @{$huids->{$uid}});
			$no->put_batch_raw($uid, $values);
			delete $h_liftover->{$uid};
			delete $huids->{$uid};
		}
		
		$last_uid = $uid;
		$i++;
	}
	
#	foreach my $uid (keys %$h_liftover) {
#		$no_lift->put_batch_raw($uid, $h_liftover->{$uid});
#	}
#	
#	foreach my $uid (keys %$huids) {
#		my $values = join(',', @{$huids->{$uid}});
#		$no->put_batch_raw($uid, $values);
#	}
	$no->write_batch;
	$no->close();
	$no_lift->write_batch;
	$no_lift->close();
	close($fh);
	close($fh_f);
	
	if ($i > 0) {
		$hres->{hg38}->{file} = $filename;
		$hres->{hg38}->{file_f} = $filename_f;
 		$pm->finish(0, $hres);
	}
	else {
		$pm->finish(0);
	}
}
sleep(3);
$pm->wait_all_children();

my (@lCsv_hg19, @lCsv_hg19_filtre);
my @lHandler_hg19;
my $hres_toput_19;


foreach my $chrID (1..22, 'X', 'Y') {
	
#	next if $chrID ne '21';
	print '-> '.$chrID;
	my (@lRocksIds, $h_liftover_global);
	my $no_lift = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out,mode=>"r",name=>$chrID.'.lift');
	$no_lift->start_iter($chrID."!");
	while (my $json = $no_lift->next($chrID."!")){
		my $h = decode_json($json);
		my $uid = $h->{uid};
		delete $h->{uid};
		push(@lRocksIds, $uid);
		$h_liftover_global->{$uid} = $h;
	}
	$no_lift->close();
	print ' -> found lift: '.scalar(scalar(@lRocksIds))."\n";

	$fork = 10 if $fork > 10;	
	my $pm2 = new Parallel::ForkManager($fork);
	$pm2->run_on_finish(
		sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data) = @_;
			push(@lCsv_hg19, $data->{hg19}->{file});
			push(@lCsv_hg19_filtre, $data->{hg19}->{file_f});
			foreach my $chr_id (keys %{$data->{to_put}}) {
				foreach my $uid_19 (keys %{$data->{to_put}->{$chr_id}}) {
					$hres_toput_19->{$chr_id}->{$uid_19} = $data->{to_put}->{$chr_id}->{$uid_19};
				}
			}
		}
	);
	
	$project->getChromosomes();
	my $nb = 50000;
	my $iter = natatime($nb, @lRocksIds);
	my $max = int(scalar(@lRocksIds) / $nb);
	$project->disconnect();
	my $nb_part = 0;
	
	while( my @tmp = $iter->() ){
		$nb_part++;
		my $pid = $pm2->start and next;
#		$project->disconnect();
		my $buffer_local = new GBuffer;
		my $project_local = $buffer_local->newProject( -name => $project_name );

		my $hres;
		
		my $csv_19 = Text::CSV->new({ binary => 1, eol => "\n" });
		my $filename_19 = $global_dir_19.'/file/part.chr'.$chrID.'-'.$nb_part.'.csv';
		my $fh_19;	
		open( $fh_19, ">", $filename_19) or die "Impossible d'ouvrir $filename_19: $!";
		$csv_19->print($fh_19, ["chr","pos","ref","alt","rocksdb_id","gene","geneid","trid","strand","tss_pos","promoterAI"]); 
		
		my $csv_f_19 = Text::CSV->new({ binary => 1, eol => "\n" });
		my $filename_f_19 = $global_dir_19.'/file/part_filter.chr'.$chrID.'-'.$nb_part.'.csv';
		my $fh_f_19;	
		open( $fh_f_19, ">", $filename_f_19) or die "Impossible d'ouvrir $filename_f_19: $!";
		$csv_f_19->print($fh_f_19, ["rocksdb_id","gene","geneid","promoterAI"]); 
		
		my $h_var_id_38;
		my @lVar;
		foreach my $global_rocks (@tmp) {
			my @ltmp = split('!', $global_rocks);
			my $chr = $project_local->getChromosome($ltmp[0]);
			my $rocks_id = $ltmp[1].'!'.$ltmp[2];
			my $var_id = $chr->transform_rocksid_to_varid($rocks_id);
			push(@lVar, $project_local->_newVariant($var_id));
			$h_var_id_38->{$var_id} = $global_rocks;
		}
		
		my $lift = liftOver->new(project=>$project_local, version=>$project_local->lift_genome_version, can_except_errors=>1);
		$lift->lift_over_variants(\@lVar);
		my $h_dv_var_ids_erros_lift = $lift->liftover_variants_errors();
		
		my $huids;
		foreach my $var (@lVar) {
			my $var_id = $var->id();
			next if (exists $h_dv_var_ids_erros_lift->{$var_id});
			my $var_id_19 = $var->{lift_over_HG19}->{id};
			my $uid_38 = $var->genomic_rocksdb_id();
			my @ltmp = split('!', $uid_38);
			$ltmp[1] = sprintf("%010d", $var->{lift_over_HG19}->{position});;
			my $uid_19 = join('!', @ltmp);
			
			foreach my $gene_infos (@{$h_liftover_global->{$uid_38}->{gene_infos}}) {
				my $chr_id = $var->getChromosome->id();
				my $pos = $var->{lift_over_HG19}->{position};
				my $ref = $h_liftover_global->{$uid_38}->{ref};
				my $alt = $h_liftover_global->{$uid_38}->{alt};
				$gene_infos =~ s/"//g;
				$csv_19->print($fh_19,['chr'.$chr_id, $pos, $ref, $alt, $uid_19, $gene_infos]);
			}
	
			my $score = $h_liftover_global->{$uid_38}->{score};
			my $gene = $h_liftover_global->{$uid_38}->{gene};
			my $gid = $h_liftover_global->{$uid_38}->{gid};
			
			if (abs($score) >= 0.2) {
				$csv_f_19->print($fh_f_19,[$uid_19, $gene, $gid.'_'.$var->getChromosome->id(), $score]);
			}
			
			$hres->{to_put}->{$var->getChromosome->id()}->{$uid_19} = $h_liftover_global->{$uid_38}->{to_put};
			
		}
		
#		$project->disconnect();
		close $fh_19;
		close $fh_f_19;
		$hres->{hg19}->{file} = $filename_19;
		$hres->{hg19}->{file_f} = $filename_f_19;
		print "$nb_part/$max \n";
	 	$pm2->finish(0, $hres);
	}
	sleep(3); 
	$pm2->wait_all_children();
 	print 'END '.$chrID;
}
	


foreach my $chr (@{$project->getChromosomes()}) {
#	next if $chr->id() ne $only_chr;
	my $no_19 = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out_19,mode=>"c",name=>$chr->id,pipeline=>1);
	if (exists $hres_toput_19->{$chr->id()}) {
		foreach my $uid_19 (keys %{$hres_toput_19->{$chr->id()}}) {
			my $values = join(',', @{$hres_toput_19->{$chr->id()}->{$uid_19}});
			$no_19->put_batch_raw($uid_19, $values);
		}
	}
	$no_19->close();
}

my $filenames = join("','", @lCsv_hg38);
my $query = qq{
	COPY (
        SELECT * from read_csv_auto([\'$filenames\']) order by chr,pos
    )
    TO '$parquet_file_final' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE);};
my $cmd = "duckdb -c \"$query\"";
system($cmd);
warn $cmd;
warn 'remove tmp files 1';
foreach my $file (@lCsv_hg38) { system("rm $file") if -e $file; }

my $filenames_f = join("','", @lCsv_hg38_filtre);
my $query_f = qq{
	COPY (
        SELECT * from read_csv_auto([\'$filenames_f\']) order by rocksdb_id
    )
    TO '$parquet_file_final_f' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE);};
my $cmd_f = "duckdb -c \"$query_f\"";
system($cmd_f);
warn $cmd_f;
warn 'remove tmp files 2';
foreach my $file (@lCsv_hg38_filtre) { system("rm $file") if -e $file; }

my $filenames_19 = $global_dir_19.'file/part.*';
my $query2 = qq{
	COPY (
        SELECT * from read_csv_auto(\'$filenames_19\') order by chr,pos
    )
    TO '$parquet_file_final_19' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE);};
my $cmd2 = "duckdb -c \"$query2\"";
system($cmd2);
warn $cmd2;
warn 'remove tmp files 3';
foreach my $file (@lCsv_hg19) { system("rm $file") if -e $file; }

my $filenames_f_19 = $global_dir_19.'file/part_filter.*';
my $query_f_2 = qq{
	COPY (
        SELECT * from read_csv_auto(\'$filenames_f_19\') order by rocksdb_id
    )
    TO '$parquet_file_final_f_19' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE);};
my $cmd_f_2 = "duckdb -c \"$query_f_2\"";
system($cmd_f_2);
warn $cmd_f_2;
warn 'remove tmp files 4';
foreach my $file (@lCsv_hg19_filtre) { system("rm $file") if -e $file; }


warn 'All Done. ok';
 
 