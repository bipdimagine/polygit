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

use GBuffer;
use Cache_PolyQuery;
use Cache_PolyDiag;
use Data::Dumper;
use Sys::Hostname;
use Storable qw/freeze thaw nfreeze nstore_fd nstore retrieve/;
use Getopt::Long;
use PBS::Client;
use Parallel::ForkManager;
use String::ProgressBar;
use POSIX qw(strftime);
use Time::Elapsed qw( elapsed );
 
my $fork = 1;
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

unless ($project_name) { die("\n\nERROR: -project option missing... Die...\n\n"); }
unless ($chr_name) { die("\n\nERROR: -chr option missing... Die...\n\n"); }
unless ($type) { die("\n\nERROR: -type option missing... Die...\n\n"); }


my $nb_step = 1;
my $buffer = new GBuffer;
my $host = hostname;

$fork = 1 unless $fork;
my $project = $buffer->newProject( -name => $project_name );
die("\n\nERROR: unknown project $project_name. Die...\n\n") unless $project;
my $diag = $project->isDiagnostic;
if ($verbose) { $project->cache_verbose(1); }
else { $project->cache_verbose(0); }

if ($type eq 'polyquery' or $type eq 'exome' or $type eq 'genome') {
	Cache_PolyQuery::create_global_infos_freeze($project);
		launch_cache_polyquery($project, $chr_name, $fork);
	launch_cache_coverage($project, $chr_name, $fork);
	warn "\n### END CACHE ".$project->name()." ###\n";
}
elsif ($type eq 'cache_polyquery' ) {
		launch_cache_polyquery($project, $chr_name, $fork);
		warn "\n### END CACHE ".$project->name()." ###\n";
}
elsif ($type eq 'coverage' ) {
	Cache_PolyQuery::create_global_infos_freeze($project);
	warn "\n### END CACHE ".$project->name()." ###\n";
}
elsif ($type eq 'global_infos' ) {
	Cache_PolyQuery::create_global_infos_freeze($project);
	warn "\n### END CACHE ".$project->name()." ###\n";
}
elsif ($type eq 'coverage' ) {
	Cache_PolyQuery::create_global_infos_freeze($project);
	warn "\n### END CACHE ".$project->name()." ###\n";
}
elsif ($type eq 'polydiag' or $type eq 'diag') {
	launch_cache_coverage($project, $chr_name, $fork);
	launch_cache_polydiag($project, $chr_name, $fork);
}

elsif ($type eq 'strict_denovo' or $type eq 'strict-denovo' or $type eq 'strictdenovo') {
	launch_cache_strictdenovo($project, $chr_name);
}

elsif ($type eq 'somatic_loh' or $type eq 'somatic-loh' or $type eq 'somaticloh' or $type eq 'loh') {
	launch_cache_somatic_loh($project, $chr_name);
}

elsif ($type eq 'global_infos' or $type eq 'global-infos' or $type eq 'globalinfos') {
	warn "\n";
	warn "### STEP $nb_step: Cache PolyQuery -> global infos.\n";
	Cache_PolyQuery::create_global_infos_freeze($project);
	$nb_step++;
}

elsif ($type eq 'check_polyquery') {
	warn "\n\n";
	warn "###################################################\n";
	warn $project->name()." - Checking ALL steps cache for polyquery\n";
	warn "###################################################\n";
	launch_check_cache_done($project, 'polyquery', $chr_name, 1, $fork);
	launch_check_cache_done($project, 'global_infos', $chr_name, 1);
	launch_check_cache_done($project, 'strict-denovo', $chr_name, 1);
	launch_check_cache_done($project, 'somatic_loh', $chr_name, 1);
	launch_check_cache_done($project, 'coverage', $chr_name, 1);
	$nb_step++;
}

elsif ($type eq 'check') {
	warn "\n\n";
	warn "###################################################\n";
	warn $project->name()." - Checking ALL steps cache\n";
	warn "###################################################\n";
	launch_check_cache_done($project, 'polyquery', $chr_name, 1, $fork);
	launch_check_cache_done($project, 'global_infos', $chr_name, 1);
	launch_check_cache_done($project, 'strict-denovo', $chr_name, 1);
	launch_check_cache_done($project, 'somatic_loh', $chr_name, 1);
	launch_check_cache_done($project, 'coverage', $chr_name, 1);
	launch_check_cache_done($project, 'polydiag', $chr_name, 1);
	$nb_step++;
}

elsif ($type eq 'all') {
	Cache_PolyQuery::create_global_infos_freeze($project);
	launch_cache_polyquery($project, $chr_name, $fork);
	launch_cache_coverage($project, $chr_name, $fork);
	launch_cache_dejavu($project, $chr_name, $fork);
	launch_cache_polydiag($project, $chr_name, $fork) if ($project->isDiagnostic());
}
else {
	die("\n\n -type $type not recognized... Die...\n\n");
}

warn "\n";
exit(0);



############### METHODS ###############



sub launch_cache_dejavu {
	my ($project, $chr_name, $fork) = @_;
	Cache_PolyQuery::cache_lite_for_dejavu($project->name(), $chr_name, $fork);
}

sub launch_cache_coverage {
	my ($project, $chr_name, $fork) = @_;
	eval {
		Cache_PolyDiag::cache_coverage($project->name(), $fork);
	};
	if ($@){
		warn "\n\n\n-@--@--@--@--@--@--@--@--@--@--@--@--@-\n";	
		warn "__  CACHE::ERROR cache_cnv  $project_name $@";
		warn "-@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@\n"	;
		exit(1);
	};
	launch_check_cache_done($project, 'coverage');
}

sub launch_cache_polydiag {
	my ($project, $chr_name, $fork) = @_;
	eval{
		Cache_PolyDiag::cache_polydiag($project->name(),$fork);
	};
	if ($@){
		warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
		warn "__  CACHE::ERROR cache_polydiag  $project_name $@";
		warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
		exit(1);
	};
	eval{
		Cache_PolyDiag::cache_cnv($project->name(),$fork);
	};
	if ($@){
		warn "\n\n\n_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n";	
		warn "__  CACHE::ERROR cache_polydiag  $project_name $@";
		warn "_!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!__!_\n"	;
	exit(1);
	};
	launch_check_cache_done($project, 'polydiag');
}

sub launch_cache_polyquery {
	my ($project, $chr_name, $fork) = @_;
	my @lChr;
	if ($chr_name eq 'all') { @lChr = @{$project->getChromosomes()}; }
	else { push(@lChr, $project->getChromosome($chr_name)); }
	foreach my $chr (@lChr) {
		my $this_chr_name = $chr->id();
		$nb_step = 1;
		warn "\n### [CHR$this_chr_name]\n";
		warn "# STEP $nb_step: Cache PolyQuery -> store ids:\n";
		my $fileout = $project->getCacheBitVectorDir().'/log/cache_'.$this_chr_name."_P1.log";
		warn 'LOG: '.$fileout."\n";
		unlink $fileout if -e $fileout;
		my ($nb_var) = Cache_PolyQuery::store_ids($chr, $fork);
		warn "-> found $nb_var variants in chr$this_chr_name\n";
		$nb_step++;
		if ($nb_var > 0) {
			warn "# STEP $nb_step: Cache PolyQuery -> annotate.\n";
			$chr->{cache_hash_get_gene_ids} = Cache_PolyQuery::store_annotations($project->name(), $this_chr_name, $fork);
			$nb_step++;
			$nb_step++;
		}
		system("date > $fileout && chmod a+rw $fileout") if -e $fileout;
	}
	launch_check_cache_done($project, 'global_infos');
	launch_check_cache_done($project, 'polyquery', $chr_name);
	launch_cache_strictdenovo($project, $chr_name);
	launch_cache_somatic_loh($project, $chr_name);
}

sub launch_cache_strictdenovo {
	my ($project, $chr_name) = @_;
	my ($status) = Cache_PolyQuery::cache_strictdenovo($project, $fork, $chr_name);
	launch_check_cache_done($project, 'strict-denovo', $chr_name) if ($status ne 'not_necessary');
}

sub launch_cache_somatic_loh {
	my ($project, $chr_name) = @_;
	my ($status) = Cache_PolyQuery::cache_somatic_loh($project, $fork, $chr_name);
	launch_check_cache_done($project, 'somatic_loh', $chr_name) if ($status ne 'not_necessary');
}

sub launch_check_cache_done {
	my ($project, $type, $chr_name, $only_check, $fork) = @_;
	Cache_PolyQuery::check_cache_done($project, $type, $chr_name, $only_check, $fork);
}
