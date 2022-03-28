package Cache_nodb_bitvector;
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use Parallel::ForkManager;
use Bio::DB::Sam;
use JSON;
use CacheGenesData_bitvector;
use Storable qw(store retrieve freeze thaw);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
#use KyotoCabinet;
use Cwd 'abs_path';
use Digest::MD5 qw(md5_hex);
use lib "$RealBin/../../GenBo/lib/obj-nodb/packages/";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;


sub createGlobalInfosKct {
	my $project = shift;
	my $hashInfos = $project->infosProject();
	my ($hash, $hash_captures, $hash_align, $hash_calling);
	$hash->{global_infos}->{name} = $hashInfos->{name};
	$hash->{global_infos}->{id} = $hashInfos->{id};
	$hash->{global_infos}->{description} = $hashInfos->{description};
	$hash->{global_infos}->{creation_date} = $hashInfos->{creation_date};
	$hash->{global_infos}->{project_type} = $hashInfos->{projectType};
	$hash->{global_infos}->{project_type_id} = $hashInfos->{projectTypeId};
	$hash->{global_infos}->{dbname} = $hashInfos->{dbname};
	$hash->{analyse}->{build} = $project->version();
	$hash->{analyse}->{exome} = $project->isExome();
	$hash->{analyse}->{genome} = $project->isGenome();
	$hash->{analyse}->{diagnostic} = $project->isDiagnostic();
	foreach my $capture (@{$project->getCaptures()}) {
		$hash_captures->{$capture->name()} = undef;
	}
	my @lCapturesName = sort keys %$hash_captures;
	$hash->{analyse}->{capture} = \@lCapturesName;
	my @lTypes = ('evs', '1000genomes', 'dbsnp', 'prediction_matrix', 'cosmic', 'exac');
	foreach my $type (@lTypes) {
		my $root_dir = $project->buffer->config->{public_data}->{$project->version()};
		my $extend =  $project->buffer->config->{kyoto}->{$type};
		my $dir = $root_dir.$extend;
		my $dir_abs_path = abs_path($dir);
		my @lTmp = split('/', $dir_abs_path);
		$hash->{analyse}->{version}->{$type} = $lTmp[-1];
	}
	foreach my $patient (@{$project->getPatients()}) {
		foreach my $method (@{$patient->callingMethods()}) { $hash_calling->{$method} = undef; }
		foreach my $method (@{$patient->alignmentMethods()}) { $hash_align->{$method} = undef; }
		foreach my $file (@{$patient->getVariationsFiles()}) {
			$hash->{check}->{vcf}->{variations}->{$patient->name()}->{$file} = md5_hex($file);
		}
		foreach my $file (@{$patient->getIndelsFiles()}) {
			$hash->{check}->{vcf}->{indels}->{$patient->name()}->{$file} = md5_hex($file);
		}
		#$hash->{check}->{sex}->{$patient->name()}->{exp} = $patient->sex();
		#$hash->{check}->{sex}->{$patient->name()}->{obs} = $patient->compute_sex();
	}
	my @lMethodsAlign = sort keys(%$hash_align);
	my @lMethodsCalling = sort keys(%$hash_calling);
	$hash->{analyse}->{alignment} = \@lMethodsAlign;
	$hash->{analyse}->{calling} = \@lMethodsCalling;
	$hash->{analyse}->{cache}->{dejavu} = strftime '%Y-%m-%d', localtime;
	unless (-d $project->getCacheBitVectorDir()) {
		my $cmd1 = "mkdir ".$project->getCacheBitVectorDir();
		my $cmd2 = "chmod 777 ".$project->getCacheBitVectorDir();
		`$cmd1`;
		`$cmd2`;
	}
	my $freeze_infos = $project->getCacheBitVectorDir().'/global_infos.freeze';
	`rm $freeze_infos` if (-e $freeze_infos);
	store($hash, $freeze_infos);
	`chmod 777 $freeze_infos`;
	return;
}

sub cache {
	my ($project, $fork, $chr_name, $interval) = @_;
	my $project_name = $project->name();
	my $now = strftime "%H:%M", localtime if ($project->cache_verbose());
	if ($project->cache_verbose()) {
		if ($interval) { warn "\n######## [START] [Limit:$interval] [$now] CACHE CHR $chr_name ########\n"; }
		else { warn "\n######## [START] [$now] CACHE CHR $chr_name ########\n"; }
	}
	CacheGenesData_bitvector::create_cache_genes($project, $chr_name, $fork, $interval);
	$now = strftime "%H:%M", localtime;
	warn "\n\n######## [END] [$now] CACHE CHR $chr_name ########\n" if ($project->cache_verbose());
	return $chr_name.' done';
}

sub cache_dejavu {
	my ($project, $fork, $chr_name, $interval) = @_;
	my $project_name = $project->name();
	my $now = strftime "%H:%M", localtime if ($project->cache_verbose());
	warn "\n######## [START] [$now] DEJAVU CHR $chr_name ########\n" if ($project->cache_verbose());
	CacheGenesData_bitvector::updateDejaVu($project, $chr_name, $fork, $interval);
	$now = strftime "%H:%M", localtime;
	warn "\n\n######## [END] [$now] DEJAVU CHR $chr_name ########\n" if ($project->cache_verbose());
	return 1;
}

sub cache_strictdenovo {
	my ($project, $fork, $chr_name) = @_;
	my $project_name = $project->name();
	my $now = strftime "%H:%M", localtime;
	warn "\n\n######## [START] [$now] STRICT_DENOVO CHR $chr_name ########\n" if ($project->cache_verbose());
	if (-e $project->getPedigreeFile()) {
		my $cmd = "$RealBin/../../polymorphism-cgi/json_output_nodb/interface_json.pl project=$project_name filter_chromosome=$chr_name prepare_strict_denovo=1 fork=$fork";
		$cmd .= " no_verbose=1" unless ($project->cache_verbose());;
		warn "\n\n### STRICT DENOVO: prepare strict-denovo model\n-> $cmd\n" if ($project->cache_verbose());
		`$cmd`;
	}
	$now = strftime "%H:%M", localtime;
	warn "\n\n######## [END] [$now] STRICT_DENOVO CHR $chr_name ########\n\n" if ($project->cache_verbose());
	return $chr_name.' done';
}
 
1;