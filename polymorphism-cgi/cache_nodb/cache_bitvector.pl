#!/usr/bin/perl
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($RealBin);
use lib "$RealBin/../GenBo";
use lib "$RealBin/../GenBo/lib/GenBoDB";
use lib "$RealBin/../GenBo/lib/obj-nodb";
use lib "$RealBin";
use lib "$RealBin/../packages/cache";
use GBuffer;
use Data::Dumper;
use CacheGenesData_nodb_bitvector;
use Cache_nodb_bitvector;
use Sys::Hostname;
use Storable qw/freeze thaw nfreeze nstore_fd nstore retrieve/;
use Getopt::Long;
use PBS::Client;

my $fork = 1;
my $interval = 1500000;
my ($project_name, $checkOnly, $chr_name, $noDiag, $onlyDiag, $launch_qsub, $only_strict_denovo, $only_deja_vu);

GetOptions(
	'fork=s' => \$fork,
	'qsub' => \$launch_qsub,
	'project=s' => \$project_name,
	'checkOnly=s' => \$checkOnly,
	'chr=s' => \$chr_name,
	'onlyDiag=s' => \$onlyDiag,
	'interval=s' => \$interval,
	'noDiag' => \$noDiag,
	'only_strict_denovo' => \$only_strict_denovo,
	'only_deja_vu' => \$only_deja_vu,
);

my $buffer = new GBuffer;
$fork = 1 unless ($fork);
$fork = 24 if ($fork > 24);
my $project = $buffer->newProject( -name => $project_name );
die("\n\nERROR: unknown project $project_name. Die...\n\n") unless $project;

my $only = 'none';
$only = 'strict_denovo' if ($only_strict_denovo);
$only = 'deja_vu'  if ($only_deja_vu);

if ($launch_qsub) {
	unless (-d $project->getCacheBitVectorDir().'/log') {
		mkdir $project->getCacheBitVectorDir().'/log';
	}
	foreach my $chr_name (1..22, 'X', 'Y', 'MT') {
		my $cmd = "/usr/bin/perl /data-xfs/dev/mbras/Workspace/polymorphism-cgi/cache_nodb/cache_bitvector.pl -fork=$fork -project=$project_name  -chr=$chr_name";
		if ($only eq 'strict_denovo')    { $cmd .= " -only_strict_denovo"; }
		elsif ($only eq 'deja_vu') { $cmd .= " -only_deja_vu"; }
		else {
			$cmd .= " -interval=$interval";
			$cmd .= " -noDiag" if ($noDiag);
		}
		warn 'LAUNCH: '.$cmd."\n";
		my $pbs = PBS::Client->new;
		my $job = PBS::Client::Job->new(
			name => 'cache_'.$project_name.'_'.$chr_name,
			ofile => $project->getCacheBitVectorDir().'/log/cache_'.$chr_name.'_OK.log',
			efile => $project->getCacheBitVectorDir().'/log/cache_'.$chr_name.'_ERROR.log',
			nodes => 1,
			ppn => $fork,
			cmd => $cmd,
		);
		$pbs->qsub($job);
		sleep(3);
	}
	exit(0);
}


unless ($checkOnly) {
	warn "\n\n";
	warn "===============================================\n";
	warn " START CACHE with FORK = $fork\n";
	warn "===============================================\n";
	Cache_nodb_bitvector::forkCache($project, $fork, $chr_name, $interval, $only);
#	if ($noDiag) {
#		warn "\n\n";
#		warn "===============================================\n";
#		warn " END CACHE with FORK".$fork."\n";
#		warn "===============================================\n";
#	}
#	elsif ($project->isDiagnostic) {
#		Cache_nodb::cache_cnv($project_name, $fork);
#		warn "\n\n";
#		warn "===============================================\n";
#		warn " COMPUTE IMAGES with FORK".$fork."\n";
#		warn "===============================================\n";
#		Cache_nodb::cache_images($project_name, $fork);
#		warn "\n\n";
#		warn "===============================================\n";
#		warn " END CACHE with FORK".$fork."\n";
#		warn "===============================================\n";
#	}
#	else {
		warn "\n\n";
		warn "===============================================\n";
		warn " END CACHE with FORK".$fork."\n";
		warn "===============================================\n";
#	}
#	exit(0);
#	if ($project->getVersion() eq 'HG19') {
#		launch_liftover_hg19_to_hg38($project_name);
#		Cache_nodb::forkCache($project, $fork, 'HG38');
#		Cache_nodb::createKctIds_oldbuild_to_newbuild($project,$fork, 'HG19', 'HG38');
#	}
#	elsif ($project->getVersion() eq 'HG38') {
#		launch_liftover_hg38_to_hg19($project_name);
#		Cache_nodb::forkCache($project,$fork,'HG19');
#		Cache_nodb::createKctIds_oldbuild_to_newbuild($project,$fork, 'HG38', 'HG19');
#	}
}
#
#warn "===============================================\n";
#warn " Check parsing VCF for each patient";
#warn "===============================================\n";
#Cache_nodb::checkParseVcfByPatient($project);

exit(0);


sub launch_liftover_hg19_to_hg38 {
	my $projectName = shift;
	my $cmd = $RealBin.'../GenBo/script/ngs_exome/liftover/liftOver_hg19_to_hg38.pl '.$projectName;
	`$cmd`;
}

sub launch_liftover_hg38_to_hg19 {
	my $projectName = shift;
	my $cmd = $RealBin.'../GenBo/script/ngs_exome/liftover/liftOver_hg38_to_hg19.pl '.$projectName;
	`$cmd`;
}