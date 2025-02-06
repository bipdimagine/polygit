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
use lib "$RealBin/../../../../../../../GenBo/lib/obj-nodb/";
use lib "$RealBin/../../../../../../../GenBo/lib/obj-nodb/packages";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use DBI;
use Compress::Snappy;
use Storable qw/thaw freeze/;
use GenBoNoSqlRocksGenome;
use GenBoNoSqlRocksAnnotation;
use File::Slurp qw(write_file);
use liftOver;
use Time::HiRes;
use Sereal qw(sereal_encode_with_object sereal_decode_with_object write_sereal_file read_sereal_file);
use Sereal::Decoder qw(decode_sereal sereal_decode_with_object scalar_looks_like_sereal);
use Archive::Tar;
use List::Util qw(
  shuffle
);
use lib "$RealBin/../../../utility";
use chunks;
use File::Temp qw(tempdir);
use Compress::Snappy;
use MCE::Loop;
use MCE::Shared;
use Term::StatusBar ;


# Créer un dataset avec des données


my $fork = 10;
my $project_name;
my $chr;
my $version;
my $force;
GetOptions(
	'chr=s'		   => \$chr,
	'fork=s'		   => \$fork,
	'version=s'   => \$version,
	'project=s'   => \$project_name,
	'force=s'   => \$force,
);

my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => "$project_name");
my $files =  chunks::get_index($project,$version);
my $status_bar = new Term::StatusBar (
                    label => 'jobs Done : ',
                   showTime=>1,
                   subTextAlign =>"center",
                    totalItems => scalar(@$files),  ## Equiv to $status->setItems(10)
                   
 );
system("clear");
$status_bar->start();
my $done=0;
sub make_progress {
   my ($total_size) = @_;
   return sub {
      my ($completed_size,@z) = @_;
      for (my $i =0;$i<abs($completed_size-$done);$i++){
      	$status_bar->update();
      }
      $done= $completed_size;
      
   };
}
my $shared_hash =  MCE::Shared->hash( );;
MCE::Loop->init(
   max_workers => 'auto', chunk_size => 'auto',progress => make_progress( scalar @$files )
);


my $decoder = Sereal::Decoder->new();
mce_loop {
   my ($mce, $files, $chunk_id) = @_;
   $shared_hash->set($chunk_id,1);
   foreach my $name (@$files) {
	my $ok = chunks::delete_project($project,"1.22500001.25000001",$version);
	die() unless $force; 
	
#	next unless @f;
	$shared_hash->del($chunk_id);
	 }
} @$files;

die() if $shared_hash->keys();
chunks::del_index($project,$version) if @$files;
exit(0);
	

	
	 
	 

