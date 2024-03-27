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
my $project = $buffer->newProject( -name => $project_name );


my $chr     = $project->getChromosome($chr_name);

my $vs = $chr->getDeletions();
my $n = 0;

foreach my $v (sort {$a->start <=> $b->start}@$vs){
	$n++;
	warn $v->start." ".$v->name." ".$v->id." ".$v->rocksdb_id." ".$v->vcf_id;
	next if $v->id ne "2_128176287_CCCAGGGC_C";	
	warn $v->rocksdb_id;
	warn $v->getChromosome->sequence($v->start,$v->start+6);
	warn "+++".$v->name." ".$v->id ;
	warn Dumper $v->sequencing_details();
	
	die();
	
}

warn $n;