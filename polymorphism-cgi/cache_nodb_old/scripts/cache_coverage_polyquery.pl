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
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use Cache_Commons;
use Cache_PolyDiag;



my $fork = 1;
my ($project_name, $chr_name);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }



warn "\n\n";
warn "###################################################\n";
warn "Prepare Cache for  COVERAGE !!!\n";
warn "###################################################\n";
Cache_PolyDiag::cache_delete_coverage_file( $project_name, $fork );
Cache_PolyDiag::cache_coverage_primers( $project_name, $fork );
Cache_PolyDiag::cache_coverage_list_primers( $project_name, $fork );
my $buffer1 = new GBuffer;
my $projectP =
  $buffer1->newProject( -name => $project_name, -verbose => 1 );
my $pm2           = new Parallel::ForkManager($fork);
my @patients_name = map { $_->name() } @{ $projectP->getPatients };
my $nbt           = scalar(@patients_name);
exit() if $nbt eq 0;
my $pr            = String::ProgressBar->new( "max" => $nbt );
my $nbb           = 1;
$pm2->run_on_finish(
	sub {
		my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) =
		  @_;
		$pr->update( $nbb++ );
		$pr->write();
	}
);
my $diag;
$diag = 1 if $projectP->isDiagnostic();
unless ($diag) {
	foreach my $patient_name (@patients_name) {
		my $pid = $pm2->start() and next;
		Cache_PolyDiag::compute_coverage_exome( $project_name, $patient_name );
		$pm2->finish(0);
	}
	$pm2->wait_all_children;
	exit();
}
my $nb = 0;
warn "coverage by patient  *******************\n";
$pr = String::ProgressBar->new( max => scalar(@patients_name) );
$nbb = 0;
$pr->update(0);
$pr->write();
foreach my $patient_name (@patients_name) {
	my $pid = $pm2->start() and next;
	Cache_PolyDiag::compute_coverage_diagnostic1( $project_name, $patient_name );
	$pm2->finish(0);
}
$pm2->wait_all_children;
warn "now Primers *******************";
$pr = String::ProgressBar->new( max => $nbt );
$nbb = 0;
$pr->update(0);
$pr->write();
foreach my $patient_name (@patients_name) {
	my $pid = $pm2->start() and next;
	Cache_PolyDiag::compute_coverage_diagnostic2( $project_name, $patient_name );
	$pm2->finish(0);
}
$pm2->wait_all_children;
warn "and finally ..  exons *****************************";
$pr = String::ProgressBar->new( max => $nbt );
$nbb = 0;
$pr->update(0);
$pr->write();
my @utrs = (1);

foreach my $patient_name (@patients_name) {
	$nb++;
	my $pid = $pm2->start() and next;
	Cache_PolyDiag::compute_coverage_diagnostic3( $project_name, $patient_name, 0 );
	$pm2->finish(0);
}
$pm2->wait_all_children;
print "\n";
