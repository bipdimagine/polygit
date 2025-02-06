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
use MCE::Flow;





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
die ("version ?") unless $version;
my $dir38_sereal = $buffer->deja_vu_public_dir($version,"projects") ;

#$buffer->config_path("root","dejavu_projects")."/".$version. "/"."/projects.tar/";
#$dir38 = "/data-beegfs/tmp/projects/";

my $dir_final= $buffer->deja_vu_public_dir($version,"variations");
 $dir_final ="/data-beegfs/tmp/$version";
my $rg38 = GenBoNoSqlRocksGenome->new(chunk_size=>2_500_000,dir=>$dir_final,mode=>"c",index=>"genomic",chromosome=>$chr,genome=>$version,pack=>"",description=>[]);

$rg38->regions();

my $regionss = $rg38->regions();
MCE::Loop->init(
   max_workers => 20, chunk_size => 'auto'
);
if (scalar(@$regionss) ==1){
		save_regions($regionss->[0]);
		exit(0);
}
mce_loop {
   my ($mce, $chunk_ref, $chunk_id) = @_;
   if (ref($chunk_ref) ne "ARRAY") {
   	save_regions($chunk_ref);
   }
   else {
   foreach my $region (@$chunk_ref){
		save_regions($region)
	}
   }
   	
} sort{$a->{start} <=> $b->{start}} @{$regionss};






sub save_regions {
	my ($region) = @_;
	warn Dumper $region;
	my $dir = $buffer->deja_vu_project_dir($version);
	my $t;
	return unless -e $dir.$region->{id};
	
	system("vmtouch -t ".$dir.$region->{id}." -q");
 	my $lmdb = GenBoNoSqlLmdb->new(dir=>"$dir",mode=>"r",name=>$region->{id});
 	my %hash;
 
 	my %pr;
 	my $nb ;
	while (my($project_name,$array) = $lmdb->next_key_value){
		last unless $project_name;
		$nb ++;


#		warn scalar(@$array);
		foreach my $a (@$array) {
			my ($rid) = $a->{id};
			my ($value) = $a->{value};
			$pr{$project_name} ++ if $value == 21;
			#warn $rid." ".$value.' '.$project_name ;
			confess() if $value =~ /!-x-!/;
			push(@{$t->{$rid}},$value);
			}
		}
		my $no =  $rg38->nosql_rocks_tmp($region);
		$no->put_batch_raw("time",time);
		foreach my $key (keys %{$t}){
			$no->put_batch_raw($key,join("!-x-!",@{$t->{$key}}) );
		}
		#warn scalar(keys %{$t});
		$no->write_batch();
		$no->close();
}

 