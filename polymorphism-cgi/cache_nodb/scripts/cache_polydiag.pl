#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../GenBo";
use lib "$RealBin/../../GenBo/lib/obj-nodb";
use lib "$RealBin/../../GenBo/lib/obj-nodb/packages/cache/polydiag";
use Data::Dumper;
use Parallel::ForkManager;
use Cwd 'abs_path';
use String::ProgressBar;
use Getopt::Long;
use Carp;
use GBuffer;
require  "$RealBin/../../GenBo/lib/obj-nodb/packages/cache/polydiag/polydiag.pm";
require  "$RealBin/../../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";


my $fork = 1;
my ($project_name, $chr_name,$patient_name);
my $version;
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'patient=s'    => \$patient_name,
	'version=s'    => \$version,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }

$| = 1;



	my $buffer = new GBuffer();
	$buffer->vmtouch(1);
	my $project = $buffer->newProjectCache( -name => $project_name, -verbose => 1 , -version=>$version);
	$project->preload_patients();
	my $patients = $project->getPatients();
	if ($patient_name) {
		
		$patients = [$project->getPatient($patient_name)];
		#$project->get_only_list_patients($patient_name);
	}
	
	die() unless $project->getPatients();
	my $tbundle;

	foreach my $t ( @{ $project->bundle_transcripts() } ) {
		$tbundle->{$t}++;
	}
	my $nbt = scalar( @{ $patients } );
	my $pr  = String::ProgressBar->new( max => $nbt );
	my $nbb = 1;
	warn "###################################################\n";
	warn "Prepare Cache for  Polydiag !!!\n";
	warn "###################################################\n";
	$| = 1;
	my $error;
	my $pm = new Parallel::ForkManager($fork);
	my $error;
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference ) = @_;
			  if($exit_code ne 0){
			  	$error =1;
			  }
			$pr->update( $nbb++ );
			$pr->write();
			$error = 1 unless defined $data_structure_reference;
			die("\npb cache_polydiag") unless defined $data_structure_reference;
		}
	);

	my $project_name = $project->name();

	$pr->write();
	foreach my $p ( @{ $patients } ) {
		my $pname = $p->name();

		#next if $pname ne "AS1502552";
		my $pid = $pm->start and next;
		#my $resp= {};

		my $resp = polydiag::run_cache_polydiag_vector( $project, $p, $tbundle,$version );

		run_cache_web_polydiag($project,$p);
		$resp = $resp + 0;
		$pm->finish( 0, \$resp );
	}

	$pm->wait_all_children;
	$project->buffer->dbh->disconnect();
	
	#warn Dumper $project
	die() if $error;

sub run_cache_web_polydiag {
	my ($project,$patient) = @_;
	my $arg1 = "project=".$project->name." patients=".$patient->name;
	my $args = " panel= edit_mode=1 never=1 this=6  allele_quality=5 report_mode=1 transcripts=all";
	my $no_cache = $patient->get_lmdb_cache_polydiag("c");
	$no_cache->put("date",time);
	$no_cache->close();
	$no_cache = undef;
	my $polydiag = $RealBin."/../../validation_variation/patient_report.pl";
	for (my $f=2;$f<5;$f++){
		for (my $imp=2;$imp<4;$imp++){
			warn "$polydiag ".$arg1." $args impact=".$imp." frequence=".$f." pipeline=1 fork=1";
			system("$polydiag ".$arg1." $args impact=".$imp." frequence=".$f." pipeline=1 fork=1");
		}
	}
	
}

 exit(0);