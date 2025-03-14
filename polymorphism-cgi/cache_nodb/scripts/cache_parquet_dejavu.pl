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
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use lib "$RealBin/../../../GenBo/lib/obj-nodb/packages";
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
use List::Util qw/shuffle/;
use lib "$RealBin/../../../polypipeline/scripts/scripts_pipeline/dejavu/utility/";
use liftOverRegions;
use chunks;
use Sereal qw(sereal_encode_with_object sereal_decode_with_object write_sereal_file read_sereal_file);
#use GenBoNoSqlDuckDBDejaVu;
use Text::CSV;
use MIME::Base64;
use Storable qw(dclone);
my $fork = 1;
my $project_name;
GetOptions(
	'fork=s' => \$fork,
	'project=s' => \$project_name,
);



my $h_models_ids;
$h_models_ids->{solo} = 'so';
$h_models_ids->{father} = 'fa';
$h_models_ids->{mother} = 'mo';
$h_models_ids->{both} = 'bo';
$h_models_ids->{is_parent} = 'pa';
$h_models_ids->{recessif} = 're';
$h_models_ids->{dominant} = 'do';
$h_models_ids->{denovo} = 'de';
$h_models_ids->{strict_denovo} = 'sd';
$h_models_ids->{error} = 'pb';


my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => "$project_name");
$project->genome_version;
 
# my $dir1 = $buffer->deja_vu_public_dir("HG19","projects_lite");
# my @files = `ls $dir1/NGS*.lite`;
# chomp(@files);
# foreach my $f (@files){
# 	my @fs = split("/",$f);
# 	my $pn = $fs[-1];
# 	$pn =~ s/\.lite//;
# 	print $pn."\n";
# 	
# }
# die();
my $data_lift;	

my @chromosomes = shuffle (@{$project->getChromosomes});
#my $dir_parquet = "/data-beegfs/parquet_HG19_HG38/";
my $dir_parquet = $project->deja_vu_public_projects_parquet().'/';

my $pm2 = new Parallel::ForkManager($fork);
my $parquet_file = $dir_parquet."/".$project->name.".".$project->id.".parquet";
exit(0) if -e $parquet_file;

$project->getPatients;
$project->preload_patients();
my $final;
 $pm2->run_on_finish(
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;
		
		unless (defined($hRes) or $exit_code > 0) {
			print qq|No message received from child process $exit_code $pid!\n|;
			return;
		}
	
	}
);


foreach my $chr (@chromosomes){	
	my $pid = $pm2->start and next;
	foreach my $patient (@{$project->getPatients()}) { $patient->setOrigin($chr); }
	my $lift = liftOverRegions->new(project=>$project,version=>$project->lift_genome_version);
	 my $array = HG19_HG38($chr,$lift);
	 $pm2->finish( 0, {lift_variants=>$array});
}
$pm2->wait_all_children();	
warn "end";


my $dir_tmp_cvs = "/tmp/pipeline/";
my $filename = "$dir_tmp_cvs/*_".$project->name.".csv";
my $dbh = DBI->connect("dbi:ODBC:Driver=DuckDB;Database=:memory:", "", "", { RaiseError => 1 , AutoCommit => 1});
my $query = "
	COPY (
        SELECT * from '$filename' order by pos38,pos19,allele,max_ratio,max_dp,transmissions
    )
    TO '$parquet_file' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE,ROW_GROUP_SIZE 1000000);";
$dbh->do($query);
system("rm $filename");
warn $parquet_file;
exit(0);

sub open_csv {
	my ($project,$chr) = @_;
	my $dir_tmp_cvs = "/tmp/pipeline/";
	my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
	my $filename = "$dir_tmp_cvs/".$project->name."_".$chr->name.".csv";
	warn $filename;
	my $fh;		
	open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
	$csv->print($fh, ["project","pos38","pos19","allele","max_ratio","max_dp","transmissions","he","ho","value"]); 
	return($filename,$fh,$csv);
}



sub HG19_HG38 {
	my($chr,$lift) =@_;
	warn $chr->name." --- ";
	
	my @lVarIds = @{$chr->getListVarObjects($chr->getVariantsVector())};
	
	my $snps;
	my $lt = time;
	
	my $t =time;
	my $snps;
	my $variations;
	my $hregions;
	my $cpt =0;
	
	
	my $time = time;
	my $vectors = get_hash_model_variant($chr);
	foreach my $v (@lVarIds){
		my $vid = $v->id;
		$cpt ++ ;
		
		my $vh;
		my ($c,$p,$a,$b) = split("_",$vid);
		my $vhh;
		$vhh->{chromosome} = $chr->ucsc_name;
		$vhh->{start} = $p;
		$vhh->{end} = $p+length($b);
		
		my $index_lmdb = $v->vector_id();
		my @ho;
		my @he;
		my @models_ho;
		my @models_he;
		my $nbhe = 0;
		my $nbho = 0;
		my $max_dp = 0;
		my $max_ratio = 0;
		foreach my $patient (@{$v->getPatients}){
			my @values;
			push(@values,$patient->id);
			my $dp = $v->getDP($patient);
			my $alt = $v->getNbAlleleAlt($patient);
			$dp = 1 if not $dp;
			$max_dp = $dp if $dp > $max_dp;
			my $ratio = ($alt/$dp)*100;
			$max_ratio = int($ratio) if int($ratio) > $max_ratio;	
			push(@values,$dp);
			push(@values,$alt);
			my $model = find_variant_model($vectors, $index_lmdb, $patient->id);
			if ($v->isHeterozygote($patient) == 1){
				push(@he,@values);		
				$nbhe ++;
				push(@models_he, $model);
			}
			else {
				push(@ho,@values);		
				$nbho ++;
				push(@models_ho, $model);
			}
		}
		my  $rocksdb_id = chunks::return_rocks_id_from_genbo_id($vid);
		if ($rocksdb_id =~ /0000000000/){
			warn "********* ".$rocksdb_id." ********* ".$vid;
			die();
		}
		my ($s,$all) = split("!",$rocksdb_id);
		$s = $s +0;
		$vhh->{chromosome} = $chr->ucsc_name;
		$vhh->{allele} = $all;
		$vhh->{start} = int($s);
		$vhh->{end} = int($s)+length($b);
		$vhh->{rocksid} = $rocksdb_id;
		$lift->add_region($vhh);
		my $value = compress1($project->id,\@he,\@ho);
		$vhh->{max_ratio} = $max_ratio;
		$vhh->{max_dp} = $max_dp;
		$vhh->{transmissions} = join(',',@models_he);
		$vhh->{transmissions} .= ',' if $vhh->{transmissions};
		$vhh->{transmissions} .= join(',',@models_ho);
		$vhh->{value} = $value;
		$snps->{$rocksdb_id} = dclone($vhh);
		
		if ($cpt %100000 ==0){
			warn $chr->name."::".scalar(@lVarIds)."::".$cpt. "time ".	abs(time - $time);
			$time = time;
		}
	}
	my $hash = $lift->liftOver_regions($project->name.".".$chr->name) ;
	my $tab;
	foreach my $a (values %$hash){
		foreach my $vl (@$a){
			my $rid = $vl->{rocksid};
			$snps->{$rid}->{LIFT} =  delete $vl->{LIFT};
		}
	}
	
	my ($filename,$fh,$csv) = open_csv($chr,$project);
	foreach my $vhh (values %$snps){
		my @vvs = unpack("w*",$vhh->{value});
		my $pid = shift @vvs;
		my $he = shift @vvs;
		my $ho = shift @vvs;
		my $value = pack("w*",@vvs);
		my ($pos19, $pos38);
		if ($project->genome_version eq "HG19"){
			$pos19 = $chr->name."!".sprintf("%010d", $vhh->{start});
			$pos38 ="0!0";
			if (exists $vhh->{LIFT}){
				$pos38 =$project->getChromosome($vhh->{LIFT}->{chromosome})->name."!".sprintf("%010d", $vhh->{LIFT}->{start}) if $project->isChromosomeName($vhh->{LIFT}->{chromosome});
			}
		}
		else {
			$pos38 = $chr->name."!".sprintf("%010d", $vhh->{start});
			$pos19 ="0!0";
			if (exists $vhh->{LIFT}){
				$pos19 =$project->getChromosome($vhh->{LIFT}->{chromosome})->name."!".sprintf("%010d", $vhh->{LIFT}->{start}) if $project->isChromosomeName($vhh->{LIFT}->{chromosome});
			}
		}
		if ($chr->name eq 'MT' ){
			if ($project->genome_version eq "HG19"){
				$pos19 = $pos38;
			}
			else {
				$pos38 = $pos19;
			}
		}
		if ($pos19 eq '0') {
			warn Dumper $vhh;
		}
		my $encoded_data = encode_base64($value,""); 
		my $toto = decode_base64($encoded_data);
		die ("\n\nPROBLEM with encode\n\n") if $toto ne $value;
		$csv->print($fh, [$project->id,$pos38,$pos19,$vhh->{allele},$vhh->{max_ratio},$vhh->{max_dp},$vhh->{transmissions},$he,$ho,$encoded_data]);
	}
	close($fh);
}

sub find_variant_model {
	my ($h_models, $vid, $patient_id) = @_;
	return $h_models_ids->{solo} if exists $h_models->{$patient_id}->{solo};
	if (exists $h_models->{$patient_id}->{dominant}) {
		return $h_models_ids->{dominant} if $h_models->{$patient_id}->{dominant}->contains($vid);
	}
	foreach my $model_name ('strict_denovo', 'denovo', 'recessif', 'father', 'mother', 'both') {
		return $h_models_ids->{$model_name} if $h_models->{$patient_id}->{$model_name}->contains($vid);
	}
	return $h_models_ids->{error};
}
	
sub get_hash_model_variant {
	my ($chr, $vector_id) = @_;
	my $hvector;
	foreach my $patient (@{$project->getPatients()}) {
		my $fam = $patient->getFamily();
		my $patient_id = $patient->id();
		if ($fam->isTrio()) {
			if ($fam->isDominant()) {
				$hvector->{$patient_id}->{dominant} = $fam->getVector_individual_dominant($chr, $patient);
			}
			$hvector->{$patient_id}->{strict_denovo} = $fam->getVector_individual_strict_denovo($chr,$patient);
			$hvector->{$patient_id}->{denovo} = $fam->getVector_individual_denovo($chr,$patient);
			$hvector->{$patient_id}->{recessif} = $fam->getVector_individual_recessive($chr,$patient);
			$hvector->{$patient_id}->{father} = $fam->getFatherVector($chr);
			$hvector->{$patient_id}->{mother} = $fam->getMotherVector($chr);
			$hvector->{$patient_id}->{both}  = $hvector->{father}->{$patient_id} & $hvector->{mother}->{$patient_id};
			$hvector->{$patient_id}->{both} -= $hvector->{recessif}->{$patient_id};
		}
		else {
			$hvector->{$patient_id}->{solo} = 1;
		}
	}
	return $hvector;
}	
exit(0);


sub compress1 {
	my ($pid,$list1,$list2) =@_;
	$list1 = [] if $list1 ==undef;
	$list2 = [] if $list2 ==undef;
	my $compressed = pack("w*",$pid, scalar(@$list1)/3, scalar(@$list2)/3, @$list1, @$list2);
	return $compressed;
}

sub decompress1{
	my ($c1) =@_;
	my ($a,$b,@t) = unpack("w*",$c1);
	my @decompressed_list1 = @t[0..$a-1];
	my @decompressed_list2 = @t[$a..$a+$b-1];
	warn $a." ".$b." ".join(";",@decompressed_list1)." ++ ".join(";",@decompressed_list2);
}




