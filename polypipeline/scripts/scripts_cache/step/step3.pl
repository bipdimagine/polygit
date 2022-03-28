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
my $set;
my $name;
my $run;
GetOptions(
	'project=s'    => \$project_name,
	'patient_name=s'    => \$patient_name,
	'set=s' => \$set,
	'name=s' => \$name,
	#'fork=s' => \$ppn,
); 
$set = "set".$set unless $set=~/set/;
my @list;
if($name){
 @list = `cat ../../../../defidiag/project/$name/$set.txt`;
 chomp(@list);
}
elsif($project_name){
	push(@list,split(",",$project_name));
}
die() unless @list;
foreach my $project_name (@list){
my $buffer = GBuffer->new();
die() unless $project_name;
my $project = $buffer->newProject( -name 			=> $project_name);
my $patients;
if ($patient_name && !$set){
	push(@$patients,$project->getPatient($patient_name));
}
else {
	$patients = $project->getPatients();
}
my $bin_dir = $RealBin."/../../../../polymorphism-cgi/validation_variation/";
foreach my $patient (@{$patients}){
	next if $patient->status ne 2;
	#next if $patient->sex ne 1;
	#next unless $patient->getFamily()->isTrio;

	#print $project_name."\n";
	#next;
	$patient_name = $patient->name;
	
my $cmd  = qq{$RealBin/../polyviewer/variations_editor_cache.pl -project=$project_name -patient=$patient_name :::20};
my $cmd2 =  qq{$RealBin/../polycyto/polycyto_cache_html.pl -project=$project_name -patient=$patient_name :::2};
print $cmd."\n";
print $cmd2."\n";

}
}
exit(0);

