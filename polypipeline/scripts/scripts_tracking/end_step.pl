#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";

use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use File::Basename;
use  File::Temp;
use Sys::Hostname;
use JSON::XS;
use Fcntl ':flock';
use Try::Tiny;

$| =1;
my $analysis_id =0;
my $step_id = "toto" ;
my $status = "titi";
my $cmd ;
my $project_name;
my $patient_name;
my $step_name;
my $prog_name ="";
my $cmd;

GetOptions(
	'project=s' => \$project_name,
	'patient=s' => \$patient_name,
	'step=s'  => \$step_name,
	'prog=s' => \$prog_name,
	'status=s' => \$status,
	'cmd=s' => \$cmd,
	);
	
my $buffer = GBuffer->new();	
my $project = $buffer->newProject( -name => $project_name );

foreach my $pname (split(",",$patient_name)){ 

my $patient = $project->getPatientOrControl($pname);

my $file = $patient->trackingFile();

my $repeat =0;
while (-e "$file.lock"){
	my $time = int(rand(10));
	sleep($time);
	$repeat ++;
	last  if $repeat == 100;
}

warn ("problem lock $file.lock") if (-e "$file.lock");
exit(0) if (-e "$file.lock");
system("touch $file.lock") ;


warn $file;
if (open my $fh, '+<', $file) {
	flock($fh, LOCK_EX) or die "Cannot lock $file - $!\n";
 	my $tracking= {};
 	try {
 	 	$tracking = decode_json <$fh>;
	}
	catch {
		close $fh;
		my $dir = $project->getPipelineTrackingDir()."/backup";

		system ("mv $file $dir/".$patient->name.".bad.".time);

		system ("mv $file $dir/$file.bad.".time);

		open (my $fh, '<', $file);
		flock($fh, LOCK_EX) or die "Cannot lock $file - $!\n";
	};
 	unless (exists  $tracking->{$step_name}){
 		my $s;
 		my $version = "";
 		$version  = $buffer->software_version($prog_name) if $prog_name;
		$s->{software} = $prog_name;
		$s->{tracking}->{start} = time;
		$s->{tracking}->{machine} = hostname;
		$s->{tracking}->{status} = "start";
		$s->{tracking}->{cmd} = $cmd;
		
		
		push(@{$tracking->{$step_name}},$s);
 	}
	die($step_name." not found ****") if exists $tracking->{$step_name}->[-1]->{tracking}->{end_time};
	$tracking->{$step_name}->[-1]->{tracking}->{end_time} = time;
	$tracking->{$step_name}->[-1]->{tracking}->{status} = $status;
 	seek $fh, 0, 0;

	print $fh encode_json $tracking;
	close $fh;
	system("chmod a+rw $file");
	warn "ok";
	unlink "$file.lock";
	next;
	}

	unlink "$file.lock";
	die();
}

exit(0);

	
	