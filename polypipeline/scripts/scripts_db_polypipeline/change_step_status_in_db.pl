#!/usr/bin/perl

use strict;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../../../../../GenBo/lib/";
use lib "$Bin/../../../../../GenBo";
use lib "$Bin/../../../../../GenBo/lib/GenBoDB";
use lib "$Bin/../../../../../GenBo/lib/obj-nodb";
use lib "$Bin/../../../../../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../../../../../GenBo/lib/kyoto";

use lib "$Bin/../../../GenBo/lib/obj-nodb/";

use Getopt::Long;
use GBuffer;


my $analysis_id =0;
my $step_name = "toto" ;
my $status = "titi";
my $patient ;
my $startdate  ;
my $enddate  ;

GetOptions(
	'analysis_id=s' => \$analysis_id,
	'step_name=s' => \$step_name,
	'status=s' => \$status,
	'patient=s' => \$patient,
	'start=s' => \$startdate,
	'end=s' => \$enddate
	);
	
change_step_status($analysis_id, $step_name, $status, $patient, $startdate, $enddate);
	
sub change_step_status {
	my ($analysis_id, $step_name, $new_status, $patient, $startdate, $enddate) = @_ ;
	my $buffer = GBuffer -> new();
	my $dbh = $buffer->dbh();
	$dbh->{AutoCommit} = 0;
	$dbh->do("use Polypipeline;");
	
	$ENV{'DATABASE'} = "";
	my $f = $patient."_";
	$step_name =~ s/$f//g ;
#	$step_name =~ s/^_//g ;
	my $sql= qq{UPDATE Polypipeline.`Analysis_Steps` SET status ="$new_status", start="$startdate", end="$enddate" WHERE analysis_id="$analysis_id" AND  step_name ="$step_name" AND patient = "$patient";};
	warn $sql ;
	$dbh->do($sql);
	$dbh->commit;
	$dbh->disconnect;
}


exit(0);






	

	