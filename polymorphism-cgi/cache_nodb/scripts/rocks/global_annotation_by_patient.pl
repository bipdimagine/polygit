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
#use lib "/software/polyweb/poly-disk/poly-src/GenBo/lib/obj-nodb/";
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use GenBoNoSqlRocksGenome;
use GenBoNoSqlRocksVariation;
use List::Util qw(shuffle);
use Sys::Hostname;
use GenBoNoSqlRocksPolyviewerVariant;
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/packages/";
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/polyviewer/";
use PolyviewerVariant;
use Deep::Hash::Utils qw(reach slurp nest deepvalue);
use Carp;

my ($project_name, $chr_name, $no_verbose, $skip_pseudo_autosomal,$version,$annot_version);

GetOptions(
	'project=s'    => \$project_name,
	'annot_version=s'    => \$annot_version,
	'chr=s'        => \$chr_name,
	'no_verbose=s' => \$no_verbose,
	'skip_pseudo_autosomal=s' => \$skip_pseudo_autosomal,
	'version=s' => \$version,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }

my $buffer = new GBuffer;
$buffer->vmtouch(1);
#my $color = $colors[ rand @colors ];
my $project = $buffer->newProjectCache( -name => $project_name );




my $dir_pipeline = $project->rocks_pipeline_directory("patients");

my $no_p = {};


foreach my $patient (@{$project->getPatients}){
my $final_polyviewer_all = GenBoNoSqlRocks->new(dir=>$project->rocks_directory."/patients/",mode=>"c",name=>$patient->name);

foreach my $chr (@{$project->getChromosomes} ){
		warn $chr->name." ".$patient->name;
		my $no =  GenBoNoSqlRocks->new(dir=>$dir_pipeline."/".$patient->name,mode=>"r",name=>$chr->name);
		my $iter = $no->rocks->new_iterator->seek_to_first;
		while (my ($key, $value) = $iter->each) {
    		$final_polyviewer_all->put_batch_raw($key,$value);
		}
		$final_polyviewer_all->write_batch();
	}
	$final_polyviewer_all->close();
}


