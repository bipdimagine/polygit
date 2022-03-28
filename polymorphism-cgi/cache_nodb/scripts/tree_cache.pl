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
use Set::IntSpan::Fast::XS;
use GenBoNoSqlLmdbIntervalTree;

#use Cache_Commons;

my $fork = 1;
my ($project_name, $chr_name);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'chr=s'    => \$chr_name,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
my $t = `hostname`;

my $nbErrors = 0;
my $buffer = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'familial' );
warn $project->dir_lite_cache_tree_objects();
 my $no = GenBoNoSqlIntervalTree->new(dir=>$project->dir_lite_cache_tree_objects(),mode=>"c");
my $no_lmdb = GenBoNoSqlLmdbIntervalTree->new(dir=>$project->dir_lite_cache_tree_objects(),mode=>"c",name=>$chr_name,is_compress=>1);
foreach my $chr (@{$project->getChromosomes}) {
	my $tree;
	my $no3 = $chr->lmdb_score_impact("r");
	my @keys = grep {$_ =~/_all/} @{$no3->get_keys()};
	foreach my $k (@keys){
		next unless $k =~/ENST/;
		my $v = $no3->get($k);
		my @pos =  $v->Index_List_Read();
		$k =~s/_all//;
		push(@$tree,[$k,$pos[0],$pos[-1]+1]);
	}
		$no_lmdb->put("transcript_vector",$tree);
		$no->put("transcripts_vector",$chr->name,$tree);
}

warn $project->dir_lite_cache_tree_objects();
$no->close();
$no_lmdb->close();


