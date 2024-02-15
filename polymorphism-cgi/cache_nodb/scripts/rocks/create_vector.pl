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




warn "\n### CACHE: store annotations step\n";
my $nbErrors = 0;
my $buffer  = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProject( -name => $project_name );
if ($annot_version) {
	$project->changeAnnotationVersion($annot_version);
}
my $chr     = $project->getChromosome($chr_name);
my $no      = $chr->rocks_polyviewer_variants("r");

my $intspan_genes_categories  = {};
my $intspan_global_categories = {};
my $categories = Cache_Commons::categories();

foreach my $c ( keys %{ $categories->{global}->{frequency} } ) {
	$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
}

foreach my $c ( keys %{ $categories->{global}->{variation_type} } ) {
	$intspan_global_categories->{$c} = Set::IntSpan::Fast::XS->new();
}
my $max;
die();
#for (my $i =0 ;$i<=$max;$i++){
#
#
#foreach my $p (@{$variation->getPatients}){
#			my $type = $variation->getSequencingGenotype($p);
#			$hpatients->{$p->name}->{all} = Set::IntSpan::Fast::XS->new() unless exists $hpatients->{$p->name}->{all};
#			$hpatients->{$p->name}->{$type} =Set::IntSpan::Fast::XS->new() unless $hpatients->{$p->name}->{$type};
#			$hpatients->{$p->name}->{all}->add($lmdb_index);
#			$hpatients->{$p->name}->{$type}->add($lmdb_index);
#		}
#}
#
#}