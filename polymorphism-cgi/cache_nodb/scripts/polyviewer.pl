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
use Set::IntSpan::Fast::XS;
use Sys::Hostname;

 my $host = hostname();



warn "*_*_*_*_*_".$host."*_*_*_*_*_";
#use Cache_Commons;

my $fork = 1;
my $force;
my ($project_name, $chr_name, $annot_version);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'chr=s'        => \$chr_name,
);
my $buffer = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProjectCache( -name => $project_name );
my $chromosomes = $project->get_list_chromosomes($chr_name);

$project->preload_patients();
$project->buffer->dbh_deconnect();

#my $chr = $project->get_list_chromosomes($chr_name);


foreach my $chr (@$chromosomes){
	my $fileout = $chr->lmdb_cache_dir()."/lmdb.ok";
	my $toto = "/tmp/toto.txt";
		
	my $cmd1 = qq{$RealBin/update_vector_score.pl -project=$project_name -fork=$fork -chr=}.$chr->name;
	system("$cmd1 && touch $toto");
	die("problem update_vector_score") unless -e $toto;
	unlink $toto;
	my $cmd2 = qq{$RealBin/update_hash_variant_chromosome.pl -project=$project_name -fork=$fork -chr=}.$chr->name;
	
	system("$cmd2  && touch $toto");
	die("problem update_vector ".$chr->name ) unless -e "$toto";
	unlink $toto;
	die("problem update_vector 2 ".$chr->name ) unless -e $fileout;

}

exit(0);