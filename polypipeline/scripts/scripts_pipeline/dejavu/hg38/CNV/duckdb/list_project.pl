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
use File::Basename;
use File::Slurp qw(write_file);
use GenBoNoSqlDejaVuCNV;
use lib "$RealBin/../../../utility";
use liftOverRegions;
use Text::CSV;
use MIME::Base64;

my $buffer = GBuffer->new();
my @releases = ("HG38");



my $dir_HG19 = $buffer->deja_vu_public_dir("HG19","CNV");

	my $dir_tmp = "/data-beegfs/tmp/";
my $dir_HG38 = "/data-beegfs/dejavu_cnv/";

my @HG19_files = `ls $dir_HG19/projects/*.dejavu`;

chomp(@HG19_files);
my %HG19_hash_files = @HG19_files;

my $id_HG19;
my $n;

my $total;
foreach my $CNVFile (keys %HG19_hash_files){
	
	my $filename = basename($CNVFile);
	my ($projectname,$rien) = split(/\./,$filename);
	print $projectname."\n";
	
}

my $dir_HG38 = $buffer->deja_vu_public_dir("HG38","CNV");
my @HG38_files = `ls $dir_HG38/projects/*.dejavu`;

chomp(@HG19_files);
my %HG38_hash_files = @HG38_files;

my $id_HG19;
my $n;

my $total;
foreach my $CNVFile (keys %HG38_hash_files){
	
	my $filename = basename($CNVFile);
	my ($projectname,$rien) = split(/\./,$filename);
	print $projectname."\n";
	
}
exit(0);