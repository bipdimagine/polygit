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
use Term::StatusBar;








my $fork = 10;
my $project_name;
my $chr;
my $version;
GetOptions(
	'chr=s'		   => \$chr,
	'fork=s'		   => \$fork,
	'version=s'   => \$version,
);
my $buffer = new GBuffer;
my $decoder = Sereal::Decoder->new();
my $dir38_sereal = $buffer->deja_vu_public_dir($version,"projects") ;
$dir38_sereal = "/data-isilon/DejaVu/$version/variations/projects.lmdb/";
my ($chunks_lift,$tree_lift) = chunks::chunks_and_tree($buffer,$version);

	
my @files;
		foreach my $chr ( keys %$chunks_lift){
			foreach my $r (@{$chunks_lift->{$chr}}){
				my $tar_file =  $dir38_sereal.$r->{id};
				next unless -e $tar_file;
				push(@files,$r->{id});
				}
			
		}
	
my $status_bar = new Term::StatusBar (
                    label => 'jobs Done : ',
                   showTime=>1,
                   subTextAlign =>"center",
                    totalItems => scalar(@files),  ## Equiv to $status->setItems(10)
                   
 );
system("clear");
$status_bar->start();
	my $done =0;
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
MCE::Loop->init(
   max_workers => 'auto', chunk_size => 'auto',progress => make_progress( scalar @files )
);

my @results = mce_loop {
   my ($mce, $files, $chunk_id) = @_;
   my $local_project;
   foreach my $tar_file (@$files){
   	 my $lmdb = GenBoNoSqlLmdb->new(dir=>$dir38_sereal,mode=>"r",name=>$tar_file);
		foreach my $p (@{$lmdb->get_keys}){
			push(@{$local_project->{$p}},$tar_file);
		}
	
	 }
	 MCE->gather($local_project);
} @files;	
	
my %projects;
 my $rocks = GenBoNoSqlLmdb->new(dir=>"$dir38_sereal",mode=>"c",name=>"projects.idx",sereal=>1);
foreach my $worker_result (@results) {
    next unless $worker_result;  # Ignorer si le worker n'a rien produit
    foreach my $p (keys %$worker_result) {
    	 $rocks->put($p,$worker_result->{$p});
    }
}
$rocks->close();

exit(0);

