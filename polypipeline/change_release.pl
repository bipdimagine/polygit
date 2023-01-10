#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo/lib/obj-nodb/";
use Getopt::Long;
use GBuffer;
use GenBoProject;
my $nobackup;
my $release ;
my $project_name;
GetOptions(
	'project=s' => \$project_name,
	'release=s' => \$release,
);

#$steps_name = "all" unless $steps_name;

my $report;
my $buffer = GBuffer->new();

my $dbh = $buffer->dbh;
$dbh->{AutoCommit} = 0;
my $project = $buffer->newProject( -name => $project_name );
update_release($dbh,$project->id);
warn "change project $project_name => $release \n" ;
$dbh->commit();
sub update_release{
	my ($dbh,$project_id) = @_;
	my $query = qq{
		UPDATE PolyprojectNGS.project_release AS t1 SET t1.release_id = (SELECT t2.release_id FROM PolyprojectNGS.releases AS t2 WHERE t2.name = "$release") where project_id = $project_id;
	};
	$dbh->do($query) ;
}

