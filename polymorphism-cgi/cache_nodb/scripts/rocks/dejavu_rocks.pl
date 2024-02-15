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
use GenBoNoSqlRocksGenome;

my $fork = 1;
my ($project_name, $chr_name);
my $version;
my $chr;
GetOptions(
	'project=s'    => \$project_name,
	'version=s'    => \$version,
	'chr=s'		   => \$chr,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }

warn "\n### Cache For Deja Vu\n";
my $buffer1 = new GBuffer;
$buffer1->vmtouch(1);
my $project = $buffer1->newProject( -name => $project_name, -verbose => 1 ,-version=>$version);


my $dir  = $project->deja_vu_lite_dir;
my $dbh = DBI->connect( "dbi:SQLite:dbname=$dir"."/$chr.dejavu.lite", "", "",{ sqlite_use_immediate_transaction => 0, } );
warn "$dir"."/1.dejavu.lite";
warn $dbh;

 my $sth = $dbh->prepare( 
      'SELECT _key,_value from __DATA__ '
      );			
 warn $sth;     
 my $rc = $sth->execute();
 my $rg = GenBoNoSqlRocksGenome->new(dir=>$project->rocks_pipeline_directory("dejavu_genomic"),mode=>"c",index=>"genomic",chromosome=>$chr,genome=>"HG19",pack=>"",description=>[]);
 warn $project->rocks_pipeline_directory("dejavu");
 my $nb =0;
 
   warn "- end -";
 #  $fp->write_batch();
   $rg->close();
   my $fork =10;
   my $pm = new Parallel::ForkManager($fork);
   
   foreach my $r (@{$rg->regions}) {
   	my $pid = $pm->start and next;
   	warn "start";
   	warn Dumper $r->{start};
   	by_region($rg,$r);
   	warn "end";
   	$pm->finish( 0, {} );
}
$pm->wait_all_children();
   exit(0);
   
   sub by_region {
   	my ($rg,$region) = @_;
   	my $no = $rg->nosql_rocks($region);
   	my $dbh = DBI->connect( "dbi:SQLite:dbname=$dir"."/$chr.dejavu.lite", "", "",{ sqlite_use_immediate_transaction => 0, } );
	my $start = $region->{start};
	my $end = $region->{end};
	warn $start." ".$end;
	warn  'SELECT _key,_value from __DATA__ where start > $start and end <= $end';
 	my $sth = $dbh->prepare( 
      "SELECT _key,_value from __DATA__ where start > $start and end <= $end"
      );			
 warn $sth;     
 my $rc = $sth->execute();
  while (my @s = $sth->fetchrow()) {
  	warn $nb if $nb%1000000 ==0;
  	# $no->write_batch() if $nb%2000000 ==0; 
  
  	$nb ++;
  	  my $obj = thaw (decompress($s[1]));
  	  my $rid = $no->return_rocks_id_from_genbo_id($s[0]);
  	  next unless $rid;
  	 $no->put_raw($rid,$obj->{data});
  		#$fp->put_raw($rid,$obj->{data});
 
   }
   	warn "END ".$start." ".$end;
 # $no->write_batch();
  $no->close();
 
   }