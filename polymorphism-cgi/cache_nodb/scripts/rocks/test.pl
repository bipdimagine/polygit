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
use lib "$RealBin/../";
use lib "$RealBin/../../../GenBo/lib/obj-nodb//packages";
use lib "$RealBin/../../../GenBo/lib/obj-nodb//polyviewer";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime); 
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp; 
use GBuffer;
use Cache_Commons;
use Sys::Hostname;
use PolyviewerVariant;
 my $host = hostname();


warn "*_*_*_*_*_".$host."*_*_*_*_*_";


my $fork = 1;
my ($project_name, $chr_name, $annot_version);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'annot_version=s'    => \$annot_version,
	'chr=s'        => \$chr_name,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }



my $buffer  = new GBuffer;
my $project = $buffer->newProjectCache( -name => $project_name );



my $chr     = $project->getChromosome($chr_name);
my $no      = $chr->rocks_polyviewer_variants("r");
#Y_22298473_A_dup-5095
my $o = $no->testPolyviewerVariant("Y_22298473_A_dup-5095");
warn Dumper $o;
warn $o;
warn $o->chromosome;
warn "is ? ".$o->isCnv;
