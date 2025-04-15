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
use File::Basename;
use File::Slurp qw(write_file);
use GenBoNoSqlDejaVuCNV;
use lib "$RealBin/../../../utility";
use liftOverRegions;
use Text::CSV;
use MIME::Base64;

my $buffer = GBuffer->new();

my $project_name;

my $dir_tmp = "/data-beegfs/tmp/";
GetOptions(
	'project=s' => \$project_name,
);

read_cnv($project_name);



exit(0);


sub depth {
    my ($h) = @_;
    my $d = 0;
    while () {
        return $d if ref($h) ne 'HASH';
        ($h) = values %$h;
        $d++;
    }
}


sub read_cnv{
	my ($project_name) = @_;

	my $cnv_callers = {
    "wisecondor"          => 1 << 0,  # 2^0 = 1
    "canvas"        => 1 << 1,  # 2^1 = 2
    "manta"        => 1 << 2,  # 2^2 = 4
};
	my $buffer = GBuffer->new();
	my $project = $buffer->newProject( -name => $project_name );
	
	my $variationsDir = $project->getCNVDir();
	my $CNVFile = $variationsDir.$project->name."_dejavu.allCNV";
	die() unless -e $CNVFile;
	my $lift = liftOverRegions->new(project=>$project,version=>$project->lift_genome_version);
	my $hp;
	foreach my $p (@{$project->getPatients}){
		$hp->{$p->name}  ++;
	}
	my $hdejavu = retrieve($CNVFile) or die "Can't retrieve datas from ".$CNVFile." !\n";
	my $hcnvs;
	my $depth = depth($hdejavu);
	my $all_hash;
	if ($depth == 5 ){
		foreach my $type (keys %{$hdejavu}) {
			foreach my $num (keys %{$hdejavu->{$type}}){
				foreach my $id (keys %{$hdejavu->{$type}->{$num}}){
				$all_hash->{$id} = $hdejavu->{$type}->{$num}->{$id};
				}
		}
		}
	}
	elsif ($depth == 3){
		$all_hash = $hdejavu;
	}
	
		foreach my $id (keys %{$all_hash})
				{	
					my ($type,$chr,$start,$end) = split("_",$id);	
					my $hcnv ;
					$hcnv->{chromosome} = $project->getChromosome($chr)->ucsc_name;
					$hcnv->{type} = $type;
					$hcnv->{id} = $id;
					my ($t,$c,$s,$end) = split("_",$id);
					$hcnv->{start} = $start;
					$hcnv->{end} = $end;
					my $callers = 0;
					foreach my $patient (keys %{$all_hash->{$id}}) {
						next unless exists $hp->{$patient};
						my $bcaller =0;
						foreach my $caller (keys %{$all_hash->{$id}->{$patient}}){
							die("-$caller") unless exists $cnv_callers->{$caller};
							$bcaller = $bcaller | $cnv_callers->{$caller};
							$callers = $callers | $cnv_callers->{$caller};
						}
						my $pid = $project->getPatient($patient)->id;
						push(@{$hcnv->{patients}},$pid);
						my $pid = $project->getPatient($patient)->id;
						push(@{$hcnv->{patients_callers}},$bcaller);
						
					}
					$hcnv->{callers} = $callers;
					$lift->add_region_id($hcnv);
				}

	my $lift =  $lift->liftOver_regions_cnv($project->name);
	my $filename = save_csv($project,$lift,$dir_tmp);
		my $dbh = DBI->connect("dbi:ODBC:Driver=DuckDB;Database=:memory:", "", "", { RaiseError => 1 , AutoCommit => 1});
	my $parquet_file = "/data-beegfs/projects.cnv.parquet/".$project->name.".parquet";
	my $query = "
	COPY (
        SELECT * from read_csv_auto([$filename]) order by type,chr38,end38,chr19,end19
    )
    TO '$parquet_file' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE,ROW_GROUP_SIZE 1000000);";
    warn $query;
	$dbh->do($query);
	$dbh->disconnect;
		
}

sub save_csv {
	my ($project,$snps,$dir_tmp) = @_;
	
	my $filename = "$dir_tmp/".$project->name.".csv";
	
	my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
	
	my $fh;		
	open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
	$csv->print($fh, ["project","chr38","start38","end38","chr19","start19","end19","type","callers","patients","patients_callers","nb_patient"]); 
	$csv->print($fh, [0,"Z",-1,-1,"Z",-1,-1,"W",0,"w","w",0]); 
	
	
	
	
	foreach my $vhh (values %$snps){
		my $chr = $project->getChromosome($vhh->{chromosome});
		my $chr0 = $chr->name;
		my $pos0 = $vhh->{start};
		my $startlift =-1;
		my $endlift =-1;
		my $chrlift ="Z";
		my $mt;
		if ($chr->name eq "MT"){
			$mt=1 if $project->getChromosome("MT")->length() == 16571 or $project->current_genome_version eq "HG38";
		}
		if (exists $vhh->{LIFT} && $vhh->{LIFT}->{MULTI}  == 1 ) {
			 if ($project->isChromosomeName($vhh->{LIFT}->{chromosome})){
			 	
				$startlift = $vhh->{LIFT}->{start};
				$endlift = $vhh->{LIFT}->{end};
				$chrlift = $project->getChromosome($vhh->{LIFT}->{chromosome})->name;
				
			 }
			
		}
		
		my $start38 = $vhh->{start};
		my $end38 = $vhh->{end};
		my $chr38 = $chr->name;
		my $start19 = $startlift ;
		my $end19 = $endlift;
		my $chr19 = $chrlift ;
		if ($project->current_genome_version eq "HG19"){
			$start19 = $vhh->{start};
			$end19 = $vhh->{end};
			$chr19 = $chr->name;
			$start38 = $startlift ;
			$end38 = $endlift ;
			$chr38 = $chrlift ;
		}
		
		if ($mt == 1 ) {
				$start19 = $start38;
				$end19 = $end38;
		}
		
		my $encoded_patients_calling = join(",",@{$vhh->{patients_callers}}); 
		my $nb_patient = scalar(@{$vhh->{patients}});
		my $encoded_data_patients = join(",",@{$vhh->{patients}}); 
		
		$csv->print($fh, [$project->id,$chr38,$start38,$end38,$chr19,$start19,$end19,$vhh->{type},$vhh->{callers},$encoded_data_patients,$encoded_patients_calling,$nb_patient]);
	}
	close($fh);
	return "\'".$filename."\'";
	
}


#sub save2 {
#	my ($dir,$ids) = @_;
#	my $nodejavu = GenBoNoSqlDejaVu->new( dir => $dir, mode => "c" );
#
#	 
#foreach my $chr (keys %{$ids})
#{
#	next if $chr =~ /KI/;
#	next if $chr =~ /GL/;
#	next if length($chr) > 4;
#	$nodejavu->create_table($chr);
#	my $sth = $nodejavu->dbh($chr)->prepare(
#		'insert into  __DATA__(_key,_value,start,end,variation_type,ho,projects)  values(?,?,?,?,?,?,?) ;') or die $DBI::errstr;
#		$sth->execute();
#				my $tree;
#				
#				foreach my $id  (keys %{$ids->{$chr}})
#				{
#								my ( $t, $c, $d, $f ) = split( /_/,$id);
#								my $text;
#								$sth->execute($id,$nodejavu->encode($ids->{$c}->{$id}->{data}),$d,$f,$t,0,0);
#								my $id2 = $nodejavu->dbh($chr)->sqlite_last_insert_rowid();
#								#warn Dumper $total->{$id};
#	  							$nodejavu->sth_insert_position_cached($chr)->execute($id2,$d,$f) ;
#				}
#$nodejavu->dbh($chr)->do(qq{CREATE UNIQUE INDEX if not exists _key_idx  on __DATA__ (_key);});
#$nodejavu->dbh($chr)->do(qq{CREATE  INDEX if not exists _start_idx  on __DATA__ (start);});
#$nodejavu->dbh($chr)->do(qq{CREATE  INDEX if not exists _end_idx  on __DATA__ (end);});
#$nodejavu->dbh($chr)->do(qq{CREATE  INDEX if not exists _type_idx  on __DATA__ (variation_type);});
#
#}
#
#$nodejavu->close();
#}
