#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use GBuffer;
use GenBoProject;
use Getopt::Long;
use Data::Dumper;
#use illumina_util;
use IO::Prompt;
use Sys::Hostname;
 use File::Find::Rule ;
use Text::Table;
use Term::Twiddle;
my $project_name;
my $project_name_origin;
my $filename;

my $name;
my $patients_name;
my $user_file;
my $steps_name;
my $force;
my $type;
my $fastq_ext;
GetOptions(
	'project=s' => \$project_name,
	'file=s' => \$user_file,
);

my @users = `cat $user_file`;
chomp (@users);
die() unless @users;

my $buffer = GBuffer->new();
my $dbh = $buffer->dbh;
$dbh->{AutoCommit} = 0;
my $project = $buffer->newProject( -name => $project_name );
my $query = $buffer->getQuery();

warn Dumper @users;
#test  

foreach my $u (@users) { 
	my $uid = $query->getUserId($u);
	warn $u." ==> ".$uid;
	die("not found this user => $u") unless $uid;
	#warn $uid;
	#next;
	addUser($buffer->dbh,$project->id,$uid);
}

$dbh->commit();
sub addUser{
	my ($dbh,$project_id,$user_id) = @_;
	my $query = qq{
		insert into PolyprojectNGS.user_projects  (PROJECT_ID,USER_ID) values ($project_id,$user_id);
	};
	warn $query;
	$dbh->do($query) ;
}

sub update_defidiag{
	my ($dbh,$project_id,$user_id) = @_;
	my $query = qq{
		update into PolyprojectNGS.projects  (validation_db) values ('defidiag') where project_id=$project_id;
	};
	$dbh->do($query) ;
}

sub update_phenotype {
	my ($dbh,$project_id,$user_id) = @_;
	my $query = qq{
			insert into PolyPhenotypeDB.phenotype_project  (PROJECT_ID,phenotype_id) values ($project_id,1);
	};
	$dbh->do($query) ;
}
