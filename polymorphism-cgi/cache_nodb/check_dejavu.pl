#!/usr/bin/perl

use FindBin qw($RealBin);
use FindBin qw($Bin);
use FindBin qw($RealBin);
use lib "$RealBin/../GenBo";
use lib "$RealBin/../GenBo/lib/GenBoDB";
use lib "$RealBin/../GenBo/lib/obj-nodb";
use lib "$RealBin/../GenBo/lib/obj-nodb/packages";
use GenBoNoSql;
use lib "$RealBin/packages";
use Carp;
use strict;
use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Carp;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use Vcf;
use Storable qw/thaw freeze/;
use JSON::XS;
use Parallel::ForkManager;
use Logfile::Rotate;
use String::ProgressBar;
use File::Temp;
use Digest::MD5::File qw( file_md5_hex );
use File::stat;
use Devel::Size qw(size total_size);
#use NoSQL;
use Tie::LevelDB;
use Compress::Snappy;
use DBI;

#use  List::Compare;
use GenBoNoSqlText;
use GenBoNoSqlDejaVu;
use Storable qw/thaw freeze/;
use List::MoreUtils qw(natatime);
use Logfile::Rotate;
use File::Util;
use File::Temp qw/ tempfile tempdir /;
use Digest::MD5 qw(md5 md5_hex md5_base64);

#use RocksDB;
use GenBoNoSqlLmdb;
use LMDB_File qw(:flags :cursor_op :error);
use List::MoreUtils qw{ natatime uniq};
use List::Util qw(max);

my $fork;
my $force;
my $verbose = 1;
my $chr;
GetOptions(
	'fork=i'    => \$fork,
	'verbose=i' => \$verbose,
	'force=i'   => \$force,
	'chr=s'     => \$chr,
);
$| = 1;

 my @chromosomes = (1..22,'X','Y','MT');
 my $dir_tmp = "/data-isilon/tmp/dejavu/";
 my $buffer    = GBuffer->new();
 my $dir_prod = $buffer->config_path("deja_vu")};
 my $dir_backup =  $dir_prod."/backup";
 

 system("mkdir $dir_backup && chmod a+rwx $dir_backup ") unless -e $dir_backup;
 my $error;
 foreach my $chr (@chromosomes){
 	
 	my $sqlite_tmp = $dir_tmp."/$chr.dejavu.lite";
 	my $sqlite_prod = $dir_prod."$chr.dejavu.lite";
 	warn "test : ".$sqlite_tmp;
 	 unless (-e $sqlite_tmp) {
 	 	$error .= "$sqlite_tmp";
 	 }
 	my $size1 = -s $sqlite_tmp;
 	my $size2 = -s $sqlite_prod;
 	warn $size1." ".$size2;
 	if ($size1 < (0.8  * $size2) && $size2 > 10000){
 		$error .= "size $sqlite_tmp ";
 	}
 }
 if ($error){
 	die();
 } 
 warn "tar";
   my $tar = $dir_prod."/dejavu.tar";
  system("cd $dir_prod;tar -cvf $tar *.*lite 2>/dev/null");
  #warn "logrotate ";
 logrotate_File($tar,$dir_backup) if -e $tar; 
 unlink $tar;
 foreach my $chr (@chromosomes){
 	my $sqlite_tmp = $dir_tmp."/$chr.dejavu.lite";
 	my $sqlite_prod = $dir_prod."/$chr.dejavu.lite";
 	system("mv $sqlite_tmp $dir_prod/");
 }

 system("mv $dir_tmp/projects.dejavu.lite $dir_prod/");
 
 sub logrotate_File {
	my ($file,$dir_backup) = @_;
	
	if (-e $file) {
		my $merge_rotate = new Logfile::Rotate (
				File  => $file,
				Count => 3,
				Gzip  => '/software/bin/bgzip.fork.sh',
				Dir	  => $dir_backup,
				Flock => 'no' );		
		$merge_rotate->rotate();
	}
	
}
