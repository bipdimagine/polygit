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
use Getopt::Long;
use GBuffer;
use Storable qw/thaw freeze/;
use File::Slurp qw(write_file);
use Time::HiRes;
use List::Util qw(
  shuffle
);
use lib "$RealBin/../../../utility";
use File::Temp qw(tempdir);
use Compress::Snappy;



my $file;
my $dir1;
GetOptions(
	'dir=s' => \$dir1,
);
my $buffer = new GBuffer;
my $dir = $buffer->dejavu_parquet_dir();

opendir(my $dh, $dir) or die "Impossible d'ouvrir '$dir': $!";
my @parquets = grep { -f "$dir/$_" && /\.parquet$/ } readdir($dh);
my  $no = GenBoNoSqlRocks->new(mode=>"c",dir=>$dir1, name=>"projects");
foreach my $l (@parquets){
	my @t = split(/\./,$l);
	warn $t[0];
	$no->put_raw($t[0],1);
	$no->put_raw($t[1],1);
	
}
$no->close();



