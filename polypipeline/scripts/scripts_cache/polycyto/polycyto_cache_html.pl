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
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
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
use File::Temp qw/ tempfile tempdir /;
#use Mojo::DOM58; 
 
 my $host = hostname();


warn "*_*_*_*_*_".$host."*_*_*_*_*_";
#use Cache_Commons;

my $project_name;
my $patient_name;
my $set;
my $name;
my $run;
GetOptions(
	'project=s'    => \$project_name,
	'patient_name=s'    => \$patient_name,
	#'fork=s' => \$ppn,
); 

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name 			=> $project_name);
my $patient;
$patient = $project->getPatient($patient_name);
my $bin_dir = $RealBin."/../../../../polymorphism-cgi/validation_variation/";

	print $patient->name." ".$patient->vstatus."\n";

	if($project->isDefidiag && $patient->isParent() && $patient->status ne 1) {
		die("probleme status :".$patient->status." on parent ".$patient->name." ".$project->name);
	}
		die($patient->name) if $patient->vstatus == 0;

$patient_name = $patient->name();

my $cmd2 = qq{$bin_dir/../manta/PolyCytoRetrieveCNV.pl project=$project_name patient=$patient_name minlength=10000 maxlength=nomax transmission=all maxfreq=nomax dejavu=10 genotype=both select_best=1 chrom=all cytoband=all genes=all print=0 omim=0};
warn $cmd2;
system($cmd2);
#my $cmd2 = qq{$bin_dir/../manta/PolyCytoRetrieveCNV.pl projectname=$project_name filename=$patient_name minlength=10000 maxlength=nomax transmission=all maxfreq=nomax dejavu=10 genotype=both select_best=1 chrom=all cytoband=all genes=all print=0 omim=0 force=1};
#warn $cmd2;
#system($cmd2.">/dev/null");
exit(0);

