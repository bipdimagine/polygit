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

my $fork = 1;
my ($project_name, $chr_name);
my $version;
GetOptions(
	'project=s'    => \$project_name,
	'version=s'    => \$version,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }

my $chr = "Y";
warn "\n### Cache For Deja Vu\n";
my $buffer1 = new GBuffer;
$buffer1->vmtouch(1);
my $project = $buffer1->newProject( -name => $project_name, -verbose => 1 ,-version=>$version);


my $dir  = $project->deja_vu_lite_dir;
my $dbh = DBI->connect( "dbi:SQLite:dbname=$dir"."/$chr.dejavu.lite", "", "",{ sqlite_use_immediate_transaction => 0, } );
warn "$dir"."/$chr.dejavu.lite";
warn $dbh;

 my $sth = $dbh->prepare( 
      'SELECT _key,_value from __DATA__ '
      );			
 warn $sth;     
 my $rc = $sth->execute();
 my $fp = GenBoNoSqlRocks->new(dir=>$project->rocks_pipeline_directory("dejavu"),mode=>"c",name=>$chr);
 
 warn $project->rocks_pipeline_directory("dejavu");
 my $nb =0;
  while (my @s = $sth->fetchrow()) {
  	warn $nb if $nb%1000000 ==0;
  	 $fp->write_batch() if $nb%1000000 ==0; 
  	$nb ++;
  	  my $obj = thaw (decompress($s[1]));
  	  my $rid = $fp->return_rocks_id_from_genbo_id($s[0]);
  	  $fp->put_batch_raw($rid,$obj->{data});
 
   }
   warn "end";
   $fp->write_batch();
   $fp->close();
   