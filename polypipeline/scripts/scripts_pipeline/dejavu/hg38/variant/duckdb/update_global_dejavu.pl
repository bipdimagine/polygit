#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use lib "$RealBin/../../../../../../../GenBo/lib/obj-nodb/";
use Carp;
use Getopt::Long;
use GBuffer;


my $fork = 4;
my $version = 'HG38';
GetOptions(
	'fork=s' => \$fork,
	'version=s' => \$version,
);


print "#### UPDATE DEJAVU ROCKS (fork $fork) ####\n\n";

print "# Step 1: check parquets project files\n";
my $buffer = new GBuffer;
my $dir = $buffer->dejavu_parquet_dir();
print "-> dir: $dir\n";

my $cmd1 = $RealBin.'/change_statut_project_parquet.pl';
print "-> cmd: $cmd1\n";
system($cmd1);


print "\n\n# Step 2: update dejavu rocks by chromosome\n";
my $pm = new Parallel::ForkManager($fork);
foreach my $chr_id (1..22, 'X', 'Y', 'MT') {
	my $pid = $pm->start and next;
	#launch chr with fork 10
	my $cmd2 = $RealBin."/dejavu_rocks.pl -version=$version -fork=10 -chr=$chr_id";
	print "-> launch: $cmd2\n";
	system($cmd2);
	print "-> chr$chr_id done!\n";
	sleep(3);
}
$pm->wait_all_children();
print "> all done\n";


