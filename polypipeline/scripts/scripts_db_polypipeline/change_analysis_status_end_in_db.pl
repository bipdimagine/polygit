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
my $status = "complete";



GetOptions(
	'analysis_id=s' => \$analysis_id,
	'status=s' => \$status
	);
	
change_step_status($analysis_id);
	
sub change_step_status {
	my ($analysis_id) = @_ ;
	my $buffer = GBuffer -> new();
	my $dbh = $buffer->dbh();
	$dbh->{AutoCommit} = 0;
	$dbh->do("use Polypipeline;");
	
	$ENV{'DATABASE'} = "";
	my $sql= qq{UPDATE Polypipeline.`Analysis` SET status ="$status", date=NOW() WHERE analysis_id="$analysis_id";};

	$dbh->do($sql);
	$dbh->commit;
	$dbh->disconnect;
}


exit(0);






	

	