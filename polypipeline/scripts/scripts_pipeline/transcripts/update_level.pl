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
use Statistics::Descriptive::Smoother;
use File::Temp;

use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
#use lib "/software/polyweb/poly-disk/poly-src/GenBo/lib/obj-nodb/packages";
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/packages";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use Set::IntSpan::Fast::XS;
use image_coverage;
 use List::Util qw( min sum max);
use polyweb_dude;
use table_dude;
 
#use Cache_Commons;

my $fork = 1;
my $force;
my ($project_name, $patient_name);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'patient=s'    => \$patient_name,
	'force=s'  => \$force,
	
);


unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
#unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }
my $t = `hostname`;
my $nbErrors = 0;
my $buffer = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProjectCache( -name => $project_name);
 warn $project->name;
my $level_value = {"high" => 3, "medium" =>2 ,"low"=>1 };

my $lists = $project-> getListTranscripts();
$patient_name = "all" unless $patient_name;
#$project->get_only_list_patients($patient_name);
$project->getCaptures();
$project->getRuns();
#if ($patient_name eq "all" ){
	$project->preload_patients();
foreach my $patient (@{$project->getPatients}){
	if ( $patient_name ne "all"){
		next if $patient->name ne "$patient_name" ;
	}
	 if ($patient->isChild() or $patient->isIll){
		#my $fork =20;
 		my $trs = table_dude::get_transcripts($patient,$fork);
 	
 		my $hRes = table_dude::getListGenes($patient,['high','medium','low'],$fork);
 		 warn "save images ==>".scalar(@$hRes);
 		my $h = table_dude::get_images($patient,$hRes,$fork);
 		warn "end !!!!";
	 }
}	

exit(0);
 
