#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../GenBo";
use lib "$RealBin/../../GenBo/lib/obj-nodb";
use lib "$RealBin/../../GenBo/lib/obj-nodb/packages/cache/polydiag";
use Data::Dumper;
use Parallel::ForkManager;
use Storable qw(store retrieve freeze thaw);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Cwd 'abs_path';
use Digest::MD5 qw(md5_hex);
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
require "$RealBin/../../GenBo/lib/obj-nodb/packages/cache/polydiag/polydiag.pm";
require "$RealBin/../../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";
use Text::Table;


my $fork = 1;
my ($project_name, $chr_name);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }

$| = 1;


warn "\n### CACHE: polyquery coverage step\n";
warn "\n\n";
warn "###################################################\n";
warn "Prepare Cache for  COVERAGE !!!\n";
warn "###################################################\n";
print "\n------------------------------\n";
my @text_table;
my $tb = Text::Table->new(
        "Steps", "Status"
    );

my @text_table = (["delete","running"],["prepare primers","waiting"],["Coverage Patients","waiting"],["Coverage Primers","waiting"],["Coverage Exons","waiting"],["CNV","waiting"],["transcripts low","waiting"]);

#push(@text_table, ["delete","running"]);


$tb->load(@text_table);
print $tb;

my $steps =0;

polydiag::cache_delete_coverage_file( $project_name, $fork );
change_table_status($steps++);

#polydiag::cache_coverage_primers( $project_name, $fork );

#polydiag::cache_coverage_list_primers( $project_name, $fork );
#die();
change_table_status($steps++);
my $buffer1 = new GBuffer;

$buffer1->vmtouch(1);
my $projectP =
  $buffer1->newProject( -name => $project_name, -verbose => 1 );
my $diag;
$diag = 1 if $projectP->isDiagnostic();

my @patients_name = map { $_->name() } @{ $projectP->getPatients };
my $nbt           = scalar(@patients_name);
return if $nbt eq 0;
my $pr            = String::ProgressBar->new( "max" => $nbt );
my $nbb           = 1;
my $error;
my $pm2           = new Parallel::ForkManager($fork);
$pm2->run_on_finish(
	sub {
		my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) = @_;
		if ($exit_code ne 0){
		 $error ++;	
		}
		$pr->update( $nbb++ );
		$pr->write();
	}
);


unless ($diag) {
	
	foreach my $patient_name (@patients_name) {
		my $pid = $pm2->start() and next;
		polydiag::compute_coverage_exome( $project_name, $patient_name );
	    
		$pm2->finish(0);
	}
	$pm2->wait_all_children;
	die() if $error > 0;
	
	exit(0);
}
my $nb = 0;
$pr = String::ProgressBar->new( max => scalar(@patients_name) );
$nbb = 0;
$pr->update(0);
$pr->write();
foreach my $patient_name (@patients_name) {
	my $pid = $pm2->start() and next;
	polydiag::compute_coverage_diagnostic1( $project_name, $patient_name );

	$pm2->finish(0);
}
$pm2->wait_all_children;
die() if $error > 0;
$pr->update( $nbb++ );
	$pr->write();
change_table_status($steps++);
warn "end 2";
$pr = String::ProgressBar->new( max => $nbt );
$nbb = 0;
$pr->update(0);
$pr->write();
foreach my $patient_name (@patients_name) {
	my $pid = $pm2->start() and next;
	polydiag::compute_coverage_diagnostic2( $project_name, $patient_name );
	$pm2->finish(0);
}
$pm2->wait_all_children;

warn "end ---";
die() if $error > 0;
$pr->update( $nbb++ );
$pr->write();
change_table_status($steps++);

$pr = String::ProgressBar->new( max => $nbt );
$nbb = 0;
$pr->update(0);
$pr->write();
my @utrs = (1);

foreach my $patient_name (@patients_name) {
	$nb++;
	my $pid = $pm2->start() and next;
	polydiag::compute_coverage_diagnostic3( $project_name, $patient_name, 0 );
	$pm2->finish(0);
}
$pm2->wait_all_children;
$pr->update( $nbb++ );
$pr->write();
print "\n";
die() if $error > 0;
change_table_status($steps++);




polydiag::cache_cnv( $project_name, $fork );
change_table_status($steps++);
polydiag::compute_coverage_diagnostic4( $project_name, $fork );
change_table_status($steps++);
my $no_cache = $projectP->get_lmdb_cache_cnv("c");
$no_cache->put("date",time);
$no_cache->close();
my $polydiag = $RealBin."/../../validation_variation/table_images_cnv.pl project=".$project_name." panel= span=20 limit=30 patients=all transcripts=all cnv=1 order=name pipeline=1";
system($polydiag);
exit(0);


sub change_table_status {
	my ($id) = @_;
	$text_table[$id][1] = "OK";
	if (defined $text_table[$id+1]){
			$text_table[$id+1][1] = "running";
	}
	
 print "\n------------------------------\n";
my $tb = Text::Table->new(
        "Steps", "Status"
    );

$tb->load(@text_table);
print $tb;
 print "\n------------------------------\n";
}


