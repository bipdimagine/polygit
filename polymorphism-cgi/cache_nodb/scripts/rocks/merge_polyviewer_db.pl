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
use Deep::Hash::Utils qw(reach slurp nest deepvalue);
use Carp;
use JSON::XS;

my ($project_name, $chr_name, $no_verbose, $skip_pseudo_autosomal,$version,$annot_version);

GetOptions(
	'project=s'    => \$project_name,
	'annot_version=s'    => \$annot_version,
	'chr=s'        => \$chr_name,
	'no_verbose=s' => \$no_verbose,
	'skip_pseudo_autosomal=s' => \$skip_pseudo_autosomal,
	'version=s' => \$version,
);
my $coder = JSON::XS->new;
my $enabled = $coder->get_convert_blessed;
warn $enabled;
unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }

my $buffer = new GBuffer;

$buffer->vmtouch(1);
#my $color = $colors[ rand @colors ];
my $project = $buffer->newProject( -name => $project_name );
my $final_polyviewer_json = GenBoNoSqlRocks->new(dir=>$project->rocks_directory("polyviewer_global_json"),mode=>"c",name=>$project->name);
my $final_polyviewer_all = GenBoNoSqlRocks->new(dir=>$project->rocks_directory("polyviewer_global"),mode=>"c",name=>$project->name);
my $nb =0;
foreach my $chr (@{$project->getChromosomes} ){
	
	warn "!! start ".$chr->name;
	my $final_polyviewer = GenBoNoSqlRocks->new(dir=>$project->rocks_directory("polyviewer_raw"),mode=>"r",name=>$chr->name.".polyviewer_variant");

		my $iter = $final_polyviewer->rocks->new_iterator->seek_to_first;
		my $nb = 0;
		warn "start chromosome ".$chr->name;
		while (my ($var_id, $value) = $iter->each) {
			
			$nb ++;
			my $c = $chr->name."!";
			next unless $var_id =~/^$c/;
			if  (ref ($var_id) =~ /HASH/) {
				warn Dumper $var_id;
				
			}
			warn $var_id." ".$nb if $nb%1000 == 0;
			my $v = $final_polyviewer_all->decode($value);
			my $test;
			foreach my $k (keys %$v){
				
				$test->{$k} = delete $v->{$k};
			}
		
			my $json  = $coder->encode($test); 
			$final_polyviewer_json->put_batch_raw($var_id,$json);
			$final_polyviewer_all->put_batch_raw($var_id,$value);
			
			
		#$finalrg->write_batch();
	}
	$final_polyviewer->close();
	$final_polyviewer_all->write_batch();
	$final_polyviewer_json->write_batch();
}
	$final_polyviewer_all->write_batch();
	$final_polyviewer_all->close();
	
	$final_polyviewer_json->write_batch();
	$final_polyviewer_json->close();
	