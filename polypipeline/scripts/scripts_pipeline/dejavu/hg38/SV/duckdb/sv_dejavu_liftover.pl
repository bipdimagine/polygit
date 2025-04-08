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

sub return_chr_name {
	my ($n) = @_;
	return "X" if $n eq "23";
	return "Y" if $n eq "24";
	return "MT" if $n eq "25";
	return $n;
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
	
	my $variationsDir = $project->getSVeqDir();();
	my $global;
	foreach my $patient (@{$project->getPatients}){
		my $CNVFile = $variationsDir.$patient->name.".allBND.store";
		unless (-e $CNVFile){
			$CNVFile = $variationsDir.$patient->name.".allBND";
		}
		next unless (-e $CNVFile);
		my $hdejavu = retrieve($CNVFile) or die "Can't retrieve datas from ".$CNVFile." !\n";
		foreach my $id (keys %$hdejavu){
			my ($c1,$pos1,$c2,$pos2) = split("_",$id);
			unless (exists $global->{$id}){
				$c1 = return_chr_name($c1);
				next unless $project->isChromosomeName($c1);
				$c2 = return_chr_name($c2);
				next unless $project->isChromosomeName($c2);
			$global->{$id}->{chromosome1} =$project->getChromosome($c1)->ucsc_name;
			$global->{$id}->{chromosome2} = $project->getChromosome($c2)->ucsc_name;;
			$global->{$id}->{position1} = $pos1;
			$global->{$id}->{position2} = $pos2;
			$global->{$id}->{type} = $hdejavu->{$id}->{TYPE};
			
			}
			push(@{$global->{$id}->{patients}},$patient->id);
			if (exists $hdejavu->{$id}->{INFOS}){
							push(@{$global->{$id}->{sr}},$hdejavu->{$id}->{INFOS}->{SR}->[0].":".$hdejavu->{$id}->{INFOS}->{SR}->[1]);
			push(@{$global->{$id}->{pr}},$hdejavu->{$id}->{INFOS}->{PR}->[0].":".$hdejavu->{$id}->{INFOS}->{PR}->[1]);
			}
			else {
				push(@{$global->{$id}->{sr}},"-1:-1");
				push(@{$global->{$id}->{pr}},"-1:-1");
			}

			
		#	my $x  = pack("w*",$patient->id,$hdejavu->{INFOS}->{pr}->[0],$hdejavu->{INFOS}->{pr}->[1],$hdejavu->{INFOS}->{sr}->[0],$hdejavu->{INFOS}->{sr}->[1]);
			#my $x = 0;#pack("w*",$patient->id,$hdejavu->{INFOS}->{pr}->[0],$hdejavu->{INFOS}->{pr}->[1],$hdejavu->{INFOS}->{sr}->[0],$hdejavu->{INFOS}->{sr}->[1])
		#	push(@{$global->{$id}->{patients}},$x);
			$global->{$id}->{nb_patients} ++;
			
		}
	}
	my $lift = liftOverRegions->new(project=>$project,version=>$project->lift_genome_version);
	
	foreach my $id (keys %$global){
		my $h = $global->{$id};
		$h->{id} = $id;
		
		$lift->add_region_bnd($h);
		
	}
	
	my $lifto =  $lift->liftOver_bnd($project->name);

	my $filename = save_csv_translocation($project,$lifto,$dir_tmp);
		my $dbh = DBI->connect("dbi:ODBC:Driver=DuckDB;Database=:memory:", "", "", { RaiseError => 1 , AutoCommit => 1});
	my $parquet_file = "/data-beegfs/projects.sv.parquet/".$project->name.".parquet";
	my $query = "
	COPY (
        SELECT * from read_csv_auto([$filename]) order by type
    )
    TO '$parquet_file' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE,ROW_GROUP_SIZE 1000000);";
    warn $query;
	$dbh->do($query);
	$dbh->disconnect;
		
}

sub save_csv_translocation {
	my ($project,$snps,$dir_tmp) = @_;
	
		my $cnv_callers = {
    "wisecondor"          => 1 << 0,  # 2^0 = 1
    "canvas"        => 1 << 1,  # 2^1 = 2
    "manta"        => 1 << 2,  # 2^2 = 4
};
	
	
	my $filename = "$dir_tmp/".$project->name.".csv";
	
	my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
	
	my $fh;		
	open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
	$csv->print($fh, ["project","type","chr1_38","position1_38","chr2_38","position2_38","chr1_19","position1_19","chr2_19","position2_19","callers","nb_patient","patients","sr","pr"]); 
	$csv->print($fh, [0,"NONE","N",-1,"N",-1,"N",-1,"N","-1",0,0,"N","N","N"]); ; 
	
	my $mt ;
	$mt = 1 if $project->getChromosome("MT")->length() > 16571 ;
	#or $project->current_genome_version eq "HG38";	
	
	foreach my $vhh (values %$snps){
			my $chr1 = $project->getChromosome($vhh->{chromosome1})->name;
			my $chr2 = $project->getChromosome($vhh->{chromosome2})->name;
			my $position1 =$vhh->{position1};
			my $position2 =$vhh->{position2};
			my $chrlift1 = "Z";
			my $chrlift2 = "Z";
			my $position1lift = -1;
			my $position2lift = -1;
			if ($vhh->{NB_LIFT} == 2 ){
				my @alift = sort{ $project->getChromosome($vhh->{chromosome1})->karyotypeId <=>  $project->getChromosome($vhh->{chromosome1})->karyotypeId } @{$vhh->{LIFT}};
				if ( $project->isChromosomeName($alift[0]->{chromosome}) && $project->isChromosomeName($alift[1]->{chromosome})){
					$chrlift1 = $project->getChromosome($alift[0]->{chromosome})->name;;
					$chrlift2 = $project->getChromosome($alift[1]->{chromosome})->name;;
					$position1lift = $alift[0]->{position};
					$position2lift = $alift[1]->{position};
				}
		}
		
		
		if ($chr1 eq "MT")  {
			if ($mt == 1){
					$position1 = $position1lift if $position1lift > -1;
			}
			else {
				$position1lift = $position1;
			}
		}
		
		if ($chr2 eq "MT" ) {
			if ($mt == 1){
					$position2 = $position2lift if $position2lift > -1;
			}
			else {
				$position2lift = $position2;
			}
		}
		
		
		my $type =  $vhh->{type};
		my $nb_patient = scalar(@{$vhh->{patients}});
		my $spatient = join(",",@{$vhh->{patients}});
		my $sr = join(",",@{$vhh->{sr}});
		my $pr = join(",",@{$vhh->{pr}});	
		if ($project->current_genome_version eq "HG19") {
			$csv->print($fh, [$project->id,$type,$chrlift1,$position1lift,$chrlift2,$position2lift,$chr1,$position1,$chr2,$position2,$cnv_callers->{manta},$nb_patient,$spatient,$sr,$pr]); 
		}
		else {
			$csv->print($fh, [$project->id,$type,$chr1,$position1,$chr2,$position2,$chrlift1,$position1lift,$chrlift2,$position2lift,$cnv_callers->{manta},$nb_patient,$spatient,$sr,$pr]); 
		}
		
	}
	close($fh);
	return "\'".$filename."\'";
	
}

