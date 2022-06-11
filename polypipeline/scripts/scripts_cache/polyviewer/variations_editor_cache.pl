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

my ($project_name);
my $patient_name;
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

	if($project->isDefidiag && $patient->isParent() && $patient->status ne 1) {
		die("probleme status :".$patient->status." on parent ".$patient->name." ".$project->name);
	}
		die($patient->name." status == 0 ") if $patient->status == 0;

$patient_name = $patient->name();
my $phenotypes = $project->getPhenotypes();
	my $ph ="";
	if (@$phenotypes){
		$ph = "phenotype='".$phenotypes->[0]->name."'";
	}
	
my $no_cache = $patient->get_lmdb_cache("w");
$no_cache->put("start_polyviewer_".time,time);
$no_cache->close();
my $cmd = qq{$bin_dir/variations_editor.pl  ac=5 ach=2 allele_quality=- annot="splicing+essential_splicing+nonsynonymous+stop+phase+maturemirna+frameshift+non-frameshift+predicted_splice_site" denovo=1 dv=3 dv_ho=2 edit_mode=1 in_this_run=6 keep_pathogenic=1 never=1 patients=$patient_name $ph project=$project_name recessive=1 report_mode=1 strict_denovo=1 user_name= xor=1};
warn $cmd;
system($cmd." >/dev/null");
my $cmd2 = qq{$bin_dir/variations_editor.pl  ac=5 ach=2 allele_quality=- annot="splicing+essential_splicing+nonsynonymous+stop+phase+maturemirna+frameshift+non-frameshift+predicted_splice_site" denovo=1 dv=4 dv_ho=2 edit_mode=1 in_this_run=6 keep_pathogenic=1 never=1 patients=$patient_name $ph project=$project_name recessive=1 report_mode=1 strict_denovo=1 user_name= xor=1};
system($cmd2." >/dev/null");
exit(0);

