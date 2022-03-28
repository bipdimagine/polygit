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



my $fork = 1;
my ($project_name, $chr_name);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }


warn "\n### Cache For Deja Vu\n";
my $buffer1 = new GBuffer;
my $project = $buffer1->newProject( -name => $project_name, -verbose => 1 );
my $root_dir = $project->deja_vu_lite_dir() . "/projects/";
mkdir $root_dir unless -e $root_dir;
#exit(0) if -e $root_dir . "/" . $project_name . ".lite";
unlink $root_dir . "/" . $project_name . ".lite" if -e $root_dir . "/" . $project_name . ".lite";
my @chr_names = map { $_->name } @{ $project->getChromosomes };
my @patient_names = sort { $a cmp $b } map { $_->name } @{ $project->getPatients };
my $dir_out = $project->getCacheBitVectorDir() . "/lmdb_cache";
my $hpatients;
for ( my $i = 0 ; $i < @patient_names ; $i++ ) {
	$hpatients->{ $patient_names[$i] } = $i;
}
my $no = GenBoNoSql->new( dir => $root_dir, mode => "c" );
$no->put( $project_name, "patients", $hpatients );
my $atleast;
foreach my $chr ( @{ $project->getChromosomes } ) {
	warn ' -> dejavu chr'.$chr->id();
	my $fileout = $dir_out . "/" . $chr->name . ".dv.freeze";
	warn "miss $fileout " unless -e $fileout;
	next unless -e $fileout;
	$atleast++;
	my $h = retrieve $fileout;
	$no->put( $project_name, $chr->name, $h );
}
confess() unless $atleast;

