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
use Bit::Vector;
use Bit::Vector::Overload;
use Carp;
use GBuffer;
use Cache_Commons;



my $fork = 1;
my ($project_name, $chr_name, $step_name);
GetOptions(
	'fork=s'     => \$fork,
	'project=s'  => \$project_name,
	'chr=s'      => \$chr_name,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }

my @lSteps = ('store_ids', 'store_annotations', 'strict_denovo', 'loh', 'global_infos');
foreach my $step_name (@lSteps) {
	my $cmd = "$RealBin/cache_check_step.pl -project=$project_name -chr=$chr_name -step=$step_name";
	`$cmd`;
}
