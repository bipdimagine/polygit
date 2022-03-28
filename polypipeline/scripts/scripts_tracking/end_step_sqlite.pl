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
use DBI;

$| =1;
my $analysis_id =0;
my $step_id = "toto" ;
my $status = "titi";
my $cmd ;
my $project_name;
my $patient_name;
my $step_name;
my $prog_name ="";
my $run_id =1;

GetOptions(
	'project=s' => \$project_name,
	'patient=s' => \$patient_name,
	'step=s'  => \$step_name,
	'prog=s' => \$prog_name,
	'status=s' => \$status,
	'cmd=s' => \$cmd,
	'run_id=s' => \$run_id,
	);
	
my $buffer = GBuffer->new();	
my $project = $buffer->newProject( -name => $project_name );
my $version  = 	 $buffer->software_version($prog_name);

foreach my $patient (@{$project->get_list__controls_patients($patient_name)}){ 

my $file = $patient->trackingFile();
my $pid = $patient->id;
my $dbh = DBI->connect("dbi:SQLite:dbname=$file.lite","","");
$dbh->{AutoCommit} = 0;
my $t = time;
my $sql = qq{
	UPDATE steps SET status="$status",end=$t WHERE run_id="$run_id"  and step_name="$step_name"
};
$dbh->do($sql) or die();
$dbh->commit();
}

exit(0);