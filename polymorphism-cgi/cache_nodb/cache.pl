#!/usr/bin/perl
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($RealBin);
use lib "$RealBin/../GenBo";
use lib "$RealBin/../GenBo/lib/GenBoDB";
use lib "$RealBin/../GenBo/lib/obj-nodb";
use lib "$RealBin/../GenBo/lib/obj-nodb/packages";
use lib "$RealBin";
use lib "$RealBin/../packages/cache";
use lib "$RealBin/../packages/cache/vector";
use lib "$RealBin/../packages/validation_variation";
#use Try::Tiny;
#use connect; 

use GBuffer;
#use GenBoStorable;
use Data::Dumper;
use CacheGenesData_nodb; 
use CacheGenesData_bitvector; 
use Cache_nodb_bitvector;
#use util_file;
use Sys::Hostname;
use Cache_nodb;
use Storable qw/freeze thaw nfreeze nstore_fd nstore retrieve/;
use image_coverage;
use Getopt::Long;
use preload_coverage;
use PBS::Client;
use Parallel::ForkManager;
use String::ProgressBar;
use POSIX qw(strftime);
 use Time::Elapsed qw( elapsed );
 use Cache_chromosome_bitvector;
 
my $fork  = 1;
my $verbose;
my ($project_name, $checkOnly, $type, $chr_name, $launch_qsub, $force, $interval_cache, $bds, $test);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'checkOnly=s'  => \$checkOnly,
	'type=s'       => \$type,
	'chr=s'        => \$chr_name,
	'qsub'         => \$launch_qsub,
	'force=s'      => \$force,
	'verbose=s'	   => \$verbose,
	'interval=s'   => \$interval_cache,
	'bds'   => \$bds,
	'test'   => \$test,
);


my $buffer = new GBuffer;
my $host = hostname;

$fork = 1 unless $fork;
my $project = $buffer->newProject( -name => $project_name );
die("\n\nERROR: unknown project $project_name. Die...\n\n") unless $project;
my $diag = $project->isDiagnostic;
if ($verbose) { $project->cache_verbose(1); }
else { $project->cache_verbose(0); }


##### LAUNCH QSUB #####
if ($test){
	Cache_chromosome_bitvector::cache_cnv($project->name,$fork);
	die();
	if ($chr_name && $type eq "part1"){
	my $fileout = 	$project->getCacheBitVectorDir().'/log/cache_'.$chr_name."_P1.log";
	warn $fileout;
	unlink $fileout if -e $fileout;
	Cache_chromosome_bitvector::storeIds($project->getChromosome($chr_name),$fork);
	Cache_chromosome_bitvector::store_annotations($project->name,$chr_name,$fork);
	warn "end **************************";
	system("date > $fileout && chmod a+rw $fileout");
	exit(0);
	#Cache_chromosome_bitvector::construct_marc_hash($project_name,$chr_name,$fork);
	}

}

if ($bds){

if ($chr_name && $type eq "part1"){
		my $fileout = 	$project->getCacheBitVectorDir().'/log/cache_'.$chr_name."_P1.log";
		unlink $fileout if -e $fileout;
		launchCacheVectorByChromosome($project, $fork,$chr_name);
		launchDenovoByChromosome($project, $fork,$chr_name);
		system("date > $fileout && chmod a+rw $fileout");
		exit(0);
}
elsif ($type eq "global"){
	my $fileout = 	$project->getCacheBitVectorDir()."/log/$project_name.global.log";
	unlink $fileout if -e $fileout;
	system("date > $fileout && chmod a+rw $fileout");
	Cache_nodb_bitvector::createGlobalInfosKct($project);
	CacheGenesData_bitvector::cache_lite_for_dejavu($project_name,$fork);
	eval {
			Cache_nodb::cache_coverage($project_name,$fork);
	};
		if ($@){
			warn "\n\n\n-@--@--@--@--@--@--@--@--@--@--@--@--@-\n";	
			warn "__  CACHE::ERROR cache_cnv  $project_name $@";
			warn "-@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@\n"	;
		exit(1);
		};	
		
		
	if ($diag) {
			eval{
		Cache_nodb::cache_polydiag($project_name,$fork);
		};
		if ($@){
		warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
		warn "__  CACHE::ERROR cache_polydiag  $project_name $@";
		warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
		exit(1);
	};
	eval{
	Cache_nodb::cache_cnv($project_name,$fork);
		};
		if ($@){
		warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
		warn "__  CACHE::ERROR cache_polydiag  $project_name $@";
		warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
	exit(1);
	};
	die("problem ") unless $project->isCacheVectorDone();
	
}
system("date > $fileout");
system("date > $fileout && chmod a+rw $fileout");
	exit(0);
}
else{
	exit(1);
}

}
if ($launch_qsub) {
	my ($listCmds, $list_ofile, $list_efile) = get_list_cmds();
	my $i = 0;
	foreach my $cmd (@$listCmds) {
		if ($project->buffer->config->{cache_vector}->{qsub_api} eq 'PBS') {
			my $pbs = PBS::Client->new;
			my $job = PBS::Client::Job->new(
				name => 'cache_'.$project_name.'_'.($i+1),
				ofile => $list_ofile->[$i],
				efile => $list_efile->[$i],
				nodes => 1,
				ppn => $fork,
				cmd => $cmd,
			);
			warn 'Launch -> '.$cmd;
			$pbs->qsub($job);
			sleep(3);
		}
		$i++;
	}
	exit(0);
}
my $total = time;

##### LAUNCH QSUB #####


if ($type eq "vector") {
	launchCacheVector();
	exit(0);
}

if ($type eq "vector_strictdenovo") {
	launchCacheVector_strictdenovo();
	exit(0);
}

if ($type eq "vector_dejavu") {

	launchCacheVector_update_dejavu();
	exit(0);
}

if ($type eq "polydiag"){
	if ($diag){
		eval{
			warn "start";
			CacheGenesData_bitvector::cache_lite_for_dejavu($project_name,$fork);
			Cache_nodb::cache_polydiag($project_name,$fork);
			warn "end";
		};
		if ($@){
			warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
			warn "__  CACHE::ERROR cache_polydiag  $project_name $@";
			warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
			exit(0);
		};
	}
	exit(0);
}

if ($type eq "special"){
	eval{
		my $status = Cache_nodb::forkCache($project,$fork);
	};
	if ($@){
		warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
		warn "__  CACHE::ERROR cache_polydiag  $project_name $@";
		warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
		exit(0);
	};
	exit(0);
}

if ($type eq "polyweb"){
	eval{
		my $status = Cache_nodb::forkCache($project,$fork);
	};
	if ($@){
		warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
		warn "__  CACHE::ERROR cache_polyweb  $project_name $@";
		warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
		exit(0);
	};
	exit(0);
}

if ($type eq "image"){
	if ($diag){
		eval{
			Cache_nodb::cache_images($project_name,$fork);
		};
		if ($@){
			warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
			warn "__  CACHE::ERROR cache_polydiag  $project_name $@";
			warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
			exit(0);
		};
	}
	exit(0);
}

if ($type eq "coverage"){
		eval{
		Cache_nodb::cache_coverage($project_name,$fork);
	};
	if ($@){
	warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
	warn "__  CACHE::ERROR cache_polydiag  $project_name $@";
	warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
	exit(0);
};

	exit(0);
}

if ($type eq "dejavu"){
		eval{
		#warn "start dejavu";
		CacheGenesData_bitvector::cache_lite_for_dejavu($project_name,$fork);
		#Cache_nodb::cache_dejavu($project_name,$fork);
	};
	if ($@){
	warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
	warn "__  CACHE::ERROR cache_polydiag  $project_name $@";
	warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
	exit(0);
};

	exit(0);
}

if ($type eq "cnv"){
	if ($diag){
		eval{
		Cache_nodb::cache_cnv($project_name,$fork);
	};
	if ($@){
	warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
	warn "__  CACHE::ERROR cache_polydiag  $project_name $@";
	warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
	exit(1);
};

	}
	exit(0);
}

launchCacheVector();
#exit(0);
#launchCacheVector_strictdenovo();


#	warn "===============================================\n";
#	warn " Start cache with fork = ".$fork."\n";
#	warn "===============================================\n";
#eval{	
#	my $status = Cache_nodb::forkCache($project,$fork);
#};
#		if ($@){
#			warn "\n\n\n-:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:-\n";	
#			warn "__  CACHE::ERROR forkCache  $project_name $@";
#			warn "-:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:--:-\n"	;
#		exit(0);
#		};

	eval {
			Cache_nodb::cache_coverage($project_name,$fork);
			
	};
		if ($@){
			warn "\n\n\n-@--@--@--@--@--@--@--@--@--@--@--@--@-\n";	
			warn "__  CACHE::ERROR cache_cnv  $project_name $@";
			warn "-@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@\n"	;
		exit(1);
		};	
		
	eval{
		CacheGenesData_bitvector::cache_lite_for_dejavu($project_name,$fork);		
		#Cache_nodb::cache_dejavu($project_name,$fork);
	};
	if ($@){
	warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
	warn "__  CACHE::ERROR cache_polydiag  $project_name $@";
	warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
	exit(1);
};
		
	if ($diag){
			eval{
		Cache_nodb::cache_polydiag($project_name,$fork);
		};
		if ($@){
		warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
		warn "__  CACHE::ERROR cache_polydiag  $project_name $@";
		warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
		exit(1);
	};
	eval{
	Cache_nodb::cache_cnv($project_name,$fork);
		};
		if ($@){
		warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
		warn "__  CACHE::ERROR cache_polydiag  $project_name $@";
		warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
	exit(1);
	};
	}
	
	
	
 warn "\n\n\n\n\n\n-*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*-\n";
 warn "-*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*-\n";
	warn "CACHE::OK  $project_name => " .elapsed(abs(time-$total))."\n";
 warn "-*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*-\n";
  warn "-*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*-\n";
#

#warn "===============================================\n";
#warn " Check parsing VCF for each patient";
#warn "===============================================\n";
#Cache_nodb::checkParseVcfByPatient($project);

exit(0);


sub get_list_cmds {
	return get_list_cmds_vector() if ($type =~ /vector/);
	my (@lChr, @lCmds, @lO, @lE);
	my $cmd = "/usr/bin/perl $RealBin/cache.pl -fork=$fork -project=$project_name -type=$type";
	push(@lCmds, $cmd);
	mkdir $project->getCacheDir().'/log' unless (-d $project->getCacheDir().'/log');
	push(@lO, $project->getCacheDir().'/log/cache_'.$chr_name.'_OK.log');
	push(@lE, $project->getCacheDir().'/log/cache_'.$chr_name.'_ERROR.log');
	return (\@lCmds, \@lO, \@lE);
}

sub get_list_cmds_vector {
	mkdir $project->getCacheBitVectorDir().'/log' unless (-d $project->getCacheBitVectorDir().'/log');
	my (@lChr, @lCmds, @lO, @lE);
	unless ($chr_name) { @lChr = (1..22, 'X', 'Y', 'MT'); }
	else { push(@lChr, $chr_name); }
	foreach my $chr_name (@lChr) {
		my $cmd = "/usr/bin/perl $RealBin/cache.pl -fork=$fork -project=$project_name -chr=$chr_name -type=$type -force=1";
		push(@lCmds, $cmd);
		push(@lO, $project->getCacheBitVectorDir().'/log/cache_'.$chr_name.'_OK.log');
		push(@lE, $project->getCacheBitVectorDir().'/log/cache_'.$chr_name.'_ERROR.log');
	}
	return (\@lCmds, \@lO, \@lE);
}

sub launchCacheVector_update_dejavu {
	my @lChr;
	unless ($chr_name) { @lChr = (1..22, 'X', 'Y', 'MT'); }
	else { push(@lChr, $chr_name); }
	foreach my $chr_name (@lChr) {
		my $interval;
		$interval = $project->buffer->config->{cache_vector}->{interval_genome} if ($project->isGenome);
		eval{
			Cache_nodb_bitvector::cache_dejavu($project, $fork, $chr_name, $interval);
		};
		if ($@){
			warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
			warn "__  CACHE::ERROR cache_vector_dejavu  $project_name $chr_name $@";
			warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
			exit(0);
		};
	}
	return;
}

sub launchCacheVector_strictdenovo {
	my $now = strftime "%H:%M:%S", localtime;
	warn "\n\n#### [START] [Fork:$fork] [$now] STRICT DENOVO $project_name ####\n\n";
	my @lChr;
	unless ($chr_name) { @lChr = (1..22, 'X', 'Y', 'MT'); }
	else { push(@lChr, $chr_name); }
	$project->buffer->dbh_reconnect();
	foreach my $chr_name (@lChr) {
		eval{
			Cache_nodb_bitvector::cache_strictdenovo($project, $fork, $chr_name);
		};
		if ($@){
			warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
			warn "__  CACHE::ERROR cache_vector_strictdenovo  $project_name $chr_name $@";
			warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
			exit(0);
		};
		my $now = strftime "%H:%M:%S", localtime;
		print "[$now] $chr_name Done !\n";
	}
	my $now = strftime "%H:%M:%S", localtime;
	warn "\n#### [END] [FORK:$fork] [$now] STRICT DENOVO $project_name #####\n\n" unless ($project->cache_verbose());
	return;
}

sub launchCacheVector {
	if ($project->isCacheVectorDone()) {
		if ($force) {
			my $dir = $project->getCacheBitVectorDir();
			rmdir($dir) unless ($chr_name);
		}
		else {
			warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
			warn "__  CACHE::ALREADY_DONE cache_vector  $project_name";
			warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n\n\n"	;
			exit(0);
		}
	}
	foreach my $patient (@{$project->getPatients()}) {
		$patient->alignmentMethods();
		$patient->callingMethods();
	}
	$project->getChromosomes();
	
	$project->getCaptures();
	my $project_name = $project->name();
	my $now = strftime "%H:%M:%S", localtime;
	if ($interval_cache) {
		warn "\n\n######## [START] [Fork:$fork] [Limit:$interval_cache] [$now] CACHE $project_name ########\n\n";
	}
	else {
		warn "\n\n######## [START] [Fork:$fork] [$now] CACHE $project_name ########\n\n";
	}
	my $patients = $project->getPatients();
	my $nb = scalar(@$patients);
	my $type = "Exome";
	$type ="Somatic" if $project->isSomaticStudy() ;
	$type ="target"  if $project->isDiagnostic() ; 
	$type ="Ciliome"  if $project->isCiliome() ; 
	$type ="Genome" if $project->isGenome() ;
	
	warn '# Project : '.$project->name()."\n";
	warn '# Type : '.$type."\n";
	warn '# Samples: '.join(",",map{$_->name}@$patients)."\n";
	warn '# Captures : '.join(",",map{$_->name}@{$project->getCaptures})."\n\n";
	
	my $t = time;
	if ($project->isGenome) {
		my @lChr;
		if ($chr_name) { push(@lChr, $chr_name); }
		else { @lChr = (1..22, 'X', 'Y', 'MT'); }
		foreach my $this_chr_name (@lChr) {
			my $interval = $project->buffer->config->{cache_vector}->{interval_genome};
			launchCacheVectorGenome($project, $this_chr_name, $interval);
		}
	}
	else { launchCacheVectorOthers($project); }
	$project->buffer->dbh_reconnect();
	Cache_nodb_bitvector::createGlobalInfosKct($project);
	$project->buffer->dbh_reconnect();
	unless ($chr_name) { launchCacheVector_strictdenovo(); }
	my $now = strftime "%H:%M:%S", localtime;
	my $eta = elapsed(abs($t-time))."\n";		
	warn "\n######## [END] [FORK:$fork] [$now] [$eta] CACHE $project_name ########\n\n";
	unless ($project->isCacheVectorDone()) {
		warn "\n\n### ERROR -> CACHE VECTOR !!! DIE !!! ###\n\n";
		die;
	}
	return;
}

sub launchCacheVectorGenome {
	my ($project, $chr_name, $interval) = @_;
	foreach my $patient (@{$project->getPatients()}) {
		$patient->alignmentMethods();
		$patient->callingMethods();
	}
	$project->getCaptures();
	my $project_name = $project->name();
		
	eval{
		Cache_nodb_bitvector::cache($project, $fork, $chr_name, $interval);
	};
	if ($@){
		warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
		warn "__  CACHE::ERROR cache_vector  $project_name $chr_name $@";
		warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
		exit(0);
	};
	return;
}

sub launchCacheVectorByChromosome {
	my ($project, $fork,$chr_name) = @_;
	my $interval = $project->buffer->config->{cache_vector}->{interval_genome};
	my $fileout = 	$project->getCacheBitVectorDir().'/log/cache_'.$chr_name."_P1.log";
	unlink $fileout if -e $fileout;
	foreach my $patient (@{$project->getPatients()}) {
		$patient->alignmentMethods();
		$patient->callingMethods();
	}
	$project->getCaptures();
	my $project_name = $project->name();
		
	eval{
		$interval = undef;

		Cache_nodb_bitvector::cache($project, $fork, $chr_name, $interval);
	};
	if ($@){
		warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
		warn "__  CACHE::ERROR cache_vector  $project_name $chr_name $@";
		warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
		exit(1);
	};
	system ("date > $fileout");
	return;
}
sub launchDenovoByChromosome {
	my ($project, $fork,$chr_name) = @_;
	my $interval = $project->buffer->config->{cache_vector}->{interval_genome};
	my $fileout = 	$project->getCacheBitVectorDir().'/log/cache_'.$chr_name."_P2.log";
	unlink $fileout if -e $fileout;
	foreach my $patient (@{$project->getPatients()}) {
		$patient->alignmentMethods();
		$patient->callingMethods();
	}
	$project->getCaptures();
	my $project_name = $project->name();
		
	eval{
		Cache_nodb_bitvector::cache_strictdenovo($project, $fork, $chr_name);
	};
	if ($@){
		warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
		warn "__  CACHE::ERROR cache_vector  $project_name $chr_name $@";
		warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
		exit(1);
	};
	warn $fileout;
	system ("date > $fileout");
	return;
}
sub launchCacheVectorOthers {
	my ($project) = @_;

	#$project->cache_verbose(0);
	my $i_progress = 1;
	my @lChr;
	if ($chr_name) {
		warn "# Chr$chr_name in progress only !\n\n";
		push(@lChr, $chr_name);
	}
	else { @lChr = (1..22, 'X', 'Y', 'MT'); }
	my $chrs = $project->getChromosomes();
	my %hreal_chr;
	foreach my $c (@{$chrs}){
		$hreal_chr{$c->name} ++;
		#warn $c->name();
	}
	my %chrs;
	my $n =0;
	foreach my $k (@lChr){
		$chrs{$k} = $n;
		$n++;
	}
	my $pm = new Parallel::ForkManager($fork);
		my $running_chr;
	$pm->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    		my $now = strftime "%H:%M:%S", localtime;
    		my $eta = $data->{eta};
    		my $chr_name = $data->{chr};
    			delete $chrs{$chr_name};
    			delete $running_chr->{$chr_name};
    			print "[$now] $chr_name Done : $eta!  ".scalar( keys%chrs)."/ 24 \n";
	    #}
    }
  );
		

	foreach my $chr_name (@lChr) {
		next unless exists $hreal_chr{$chr_name};
		$running_chr->{$chr_name} ++;
		my $pid = $pm->start and next;
		my $t = time;
		
		$project->buffer->dbh_reconnect();
		if ($interval_cache) {
			Cache_nodb_bitvector::cache($project, 1, $chr_name, $interval_cache);
		}
		else {
			Cache_nodb_bitvector::cache($project, 1, $chr_name);
		}
		$project->buffer->dbh_deconnect();
		
		my $eta = elapsed(abs($t-time));
		my $res;
		$res->{eta} = $eta;
		$res->{chr} = $chr_name;
		
		$pm->finish(0,$res),;
	}
	$pm->wait_all_children();
	
	confess(Dumper $running_chr) if scalar(keys %$running_chr);
	
	return;
}
