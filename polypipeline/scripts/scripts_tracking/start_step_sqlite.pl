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
use SQL::Abstract;
use UUID 'uuid';


$| =1;
my $analysis_id =0;
my $step_id = "toto" ;
my $status;
my $cmd ;
my $project_name;
my $patient_name;
my $step_name;
my $prog_name ="";
my $run_id;

GetOptions(
	'project=s' => \$project_name,
	'patient=s' => \$patient_name,
	'step=s'  => \$step_name,
	'prog=s' => \$prog_name,
	'status=s' => \$status,
	'cmd=s' => \$cmd,
	'run_id=s' => \$run_id,
	'status=s' => \$status,
	);
die() unless $run_id;	
my $buffer = GBuffer->new();	
my $project = $buffer->newProject( -name => $project_name );
my $version ={};
if ($prog_name){
try {	
 $version  = 	 $buffer->software_version($prog_name);
}
}
 my $version_json = encode_json $version;
 
foreach my $patient (@{$project->get_list__controls_patients($patient_name)}){ 

#my $patient = $project->getPatientOrControl($pname);
my $file = $patient->trackingFile();
my $dbh = DBI->connect("dbi:SQLite:dbname=$file.lite","","");
$dbh->{AutoCommit} = 0;
#my $sql_create = qq{CREATE TABLE IF NOT EXISTS `steps` (
#  `step_id` INTEGER PRIMARY KEY,
#  `run_id` TEXT ,
#  `project_id` TEXT ,
#  `software` TEXT ,
#  `software_version` TEXT ,
#  `machine` TEXT ,
#  `start` INTEGER ,
#  `status` TEXT ,
#  `cmd` TEXT ,
#  `step_name` TEXT ,
#  `end` INTEGER ,
#  UNIQUE(run_id,step_name)
#);
#};
#$dbh->do($sql_create);

my $sql = qq{INSERT INTO steps
                (run_id,step_name, software,software_version, machine, status, cmd,start)
                VALUES (?, ?, ?,?,?,?,?,?)};
 my $sth = $dbh->prepare($sql);               

my @values =($run_id,$step_name,$prog_name,$version_json,hostname,$status,$cmd,time);

$sth->execute(@values);
$dbh->commit();
}

exit(0);