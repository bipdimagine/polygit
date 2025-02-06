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
use lib "$RealBin/../../../../../../GenBo/lib/obj-nodb/";
use lib "$RealBin/../../../../../../GenBo/lib/obj-nodb/packages";
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
use lib "$RealBin/../../utility";
use chunks;
use File::Temp qw(tempdir);
use Compress::Snappy;








my $fork = 10;
my $project_name;
my $chr;
GetOptions(
	'chr=s'		   => \$chr,
	'fork=s'		   => \$fork,
);

my $dir38 = "/data-isilon/DejaVu/HG38/variations/projects/rocks.v2/";
#my $dir38_sereal = "/data-isilon/DejaVu/HG38/variations/projects/";
my $dir38_rocks = chunks::hg38_rocks_projects; 

#$dir38 = "/data-beegfs/tmp/projects/";
my $dir_final="/data-beegfs/tmp/rocks";





my $pm = new Parallel::ForkManager($fork);
 
my %jobs;
my $jobid = time;
my $chromosomes = chunks::chromosomes;
foreach my $chr (@{$chromosomes}){
	my $rg38 = GenBoNoSqlRocksGenome->new(chunk_size=>2_500_000,dir=>$dir38_rocks,mode=>"c",index=>"genomic",chromosome=>$chr,genome=>"HG38",pack=>"",description=>[]);
	$rg38->regions();

 foreach my $region (@{$rg38->regions}) {
	
	 $jobid ++;
   	$jobs{$jobid} ++;
   	my $pid = $pm->start and next;
	my $no =  $rg38->nosql_rocks($region);
	my $decoder = Sereal::Decoder->new();
	my $nb;
	my $nbv = 0;
	my $uid = 1;
	my $t = {};
	my $ti = time;
	my $at;
	my $at2;
	my $tar_file =  $dir38_sereal."/$chr/".$region->{id}.".tar";
	my $tar = Archive::Tar->new;
	$pm->finish( 0, {data=>[],job=>$jobid} ) unless -e $tar_file;
	 confess($tar_file)  unless -e $tar_file;
    $tar->read($tar_file);
	foreach my $project_name ($tar->list_files) {
		my $t3 =time;
		my $data = $tar->get_content( $project_name );
		my $array = $decoder->decode($data);
		
		$no->put($project_name,$array);
		}
		warn "end $tar_file";
		$no->write_batch;
		$no->rocks->compact_range();
		$pm->finish( 0, {data=>[],job=>$jobid} );
	}
	

 }
  $pm->wait_all_children();
