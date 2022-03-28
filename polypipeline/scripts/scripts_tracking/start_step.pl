#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";

use Getopt::Long;
use Carp;
use GBuffer;
use Term::ANSIColor;
use Sys::Hostname;
use JSON::XS;
use Fcntl ':flock';
use Try::Tiny;

my $analysis_id =0;
my $step_id = "toto" ;
my $status = "titi";
my $cmd ;
my $project_name;
my $patient_name;
my $step_name;
my $prog_name;

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
try{


confess("\n\nERROR: $patient_name not found... Die...\n\n") unless ($patient);

$| =1;
my $version  = 	 $buffer->software_version($prog_name);

my $new ;
unless( -e $file) {
system("touch $file;echo $step_name>$file");
$new =1;
}
$new = 1 if -z $file;

my $repeat =0;
my $limit = 10 + int(rand(20));
while (-e "$file.lock"){
	my $time = int(rand(10));
	sleep($time);
	$repeat ++;
	warn $repeat."-->".$limit;
	last  if $repeat == $limit;
}

warn ("problem lock ==> $file.lock $repeat ") if (-e "$file.lock");
unlink if (-e "$file.lock");
system("touch $file.lock") ;
if (open (my $fh, '+<', $file)) {
	flock($fh, LOCK_EX) or die "Cannot lock $file - $!\n";
	my $tracking = {};
	try {
 	 	$tracking = decode_json <$fh> unless $new;
	}
	catch {
		close $fh;
		my $dir = $project->getPipelineTrackingDir()."/backup";
		system ("mv $file $dir/".$patient->name.".bad.".time);
		open (my $fh, '<', $file);
		flock($fh, LOCK_EX) or die "Cannot lock $file - $!\n";
	};
	my $s;
	$s->{software} = $version;
	$s->{tracking}->{start} = time;
	$s->{tracking}->{machine} = hostname;
	$s->{tracking}->{status} = "running";
	$s->{tracking}->{cmd} = $cmd;

	push(@{$tracking->{$step_name}},$s);
	seek $fh, 0, 0 unless $new;

	print $fh encode_json $tracking;
	close $fh;
	system("chmod a+rw $file");
	unlink "$file.lock";
	next;
}
unlink "$file.lock";
die();
}
catch {
	unlink "$file.lock" if (-e "$file.lock");
}
}
exit(0);
