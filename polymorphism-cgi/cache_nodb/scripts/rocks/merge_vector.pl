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
use lib "$RealBin/../../../GenBo/lib/obj-nodb/packages";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use Set::IntSpan::Fast::XS;
use GBuffer;
use Bit::Vector;
#use Cache_Commons;

my $fork = 1;
my $force;
my ($project_name, $chr_name,$annot_version);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'chr=s'        => \$chr_name,
	'annot_version=s'    => \$annot_version,
	'force=s'  => \$force,
);
my $hrun;
unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
#unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }
my $t = `hostname`;
my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'familial' );
my $nbErrors = 0;
my $values;
my $size = 0 ;
my @order;
foreach my $chr (@{$project->getChromosomes}){
	my $rocks1 = $chr->rocks_vector("r");
	my $iter = $rocks1->rocks->new_iterator->seek_to_first;
warn $chr->name;


while (my ($key, $value) = $iter->each) {
			next if $key =~ /^ENSG/;
			next if $key =~ /^intergenic/;
    		#$values->{vector}->{$key}->{$chr->name} = $rocks1->decode($value);
    		push(@{$values->{vector}->{$key}},$rocks1->decode($value));
    		$values->{size}->{$chr->name} = $rocks1->size;
    		my $start = $size;
    		$size+=$rocks1->size;
    		push(@order,{size=>$rocks1->size, end=>$size-1,start=>$start,name=>$chr->name});
 }
}

my $genomic_vector ;
my $rocks_genomic = GenBoNoSqlRocksVector->new(dir=>"/data-beegfs/test-cache//vector/",mode=>"c",name=>$project->name.".genomic.vector");
$rocks_genomic->size($size);
$rocks_genomic->put_batch("chromosome",\@order);
foreach my $key (keys %{$values->{vector}}){
	my @array =();
	warn $key;
	#foreach my $a (@order){
	#	push(@array,$values->{vector}->{$key}->{$a->{name}});
	#}
	my $v = Bit::Vector->Concat_List(@{$values->{vector}->{$key}});
	warn $v->Size;
	 $rocks_genomic->put_batch($key, $v);
 	# $rocks_genomic->put_batch($key, Bit::Vector->Concat_List(@array));
 	 $rocks_genomic->write_batch;
 	#warn $vector->Size;
}
$rocks_genomic->write_batch;
$rocks_genomic->close;
warn Dumper keys %$genomic_vector;

my $t = time;
 $genomic_vector->{"9802+all"}&=$genomic_vector->{"gnomad_ho_ac_100"};
 warn time -$t;
sleep(1);
warn $genomic_vector->{"9802+all"};
die();