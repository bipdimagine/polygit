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
use File::Slurp qw(write_file);

my $fork = 1;
my ($project_name, $chr_name);
my $version;
my $chr;
my $file;
GetOptions(
	'chr=s'		   => \$chr,
	'fork=s'		   => \$fork,
	'file=s' => \$file,
);

warn "\n### Cache For Deja Vu\n";
my $buffer = new GBuffer;
my $ucsc_chr = "chr".$chr;
$ucsc_chr = "chrM" if $ucsc_chr eq "chrMT";

my $dir38 = my $dir38 = $buffer->config->{deja_vu}->{path_rocks}."/HG38/".$buffer->config->{deja_vu}->{variations}."/rocks/";
$dir38 =~ s/HG19/HG38/;
my $rg38 = GenBoNoSqlRocksGenome->new(dir=>$dir38,mode=>"w",index=>"genomic",chromosome=>$chr,genome=>"HG38",pack=>"",description=>[]);


  my $v1 = Bio::DB::HTS::Tabix->new( filename => "$file");
   foreach my $region (@{$rg38->regions}) {
	my $no38 =  $rg38->nosql_rocks($region);
	$no38->put_raw("date",time);
	my $start = $region->{start};
	my $end = $region->{end};
	my $iter = $v1->query($ucsc_chr.":".$start."-".$end);
	
	if ($iter){
	my $xs;
	while (my $line = $iter->next) {
		chomp($line);
   		my ($chr,$pos38,$end,$id,$data) = split(" ",$line);
 		my ($a,$b,$c,$d) = split("_",$id);
   		my $id =  join("_",$a,$pos38,$c,$d);
   		my $rockid = $no38->return_rocks_id_from_genbo_id($id);
   		$no38->put_raw($rockid,$data);
	}
	
	warn "close ==> ".$ucsc_chr.":".$start."-".$end;
   }
   }
   $rg38->close();
   unlink $file;
   exit(0);



