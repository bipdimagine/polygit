#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
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
	'user=s' => \$user_file,
);

my @users = `cat ./users/$user_file.txt`;
chomp (@users);
die() unless @users;

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $query = $buffer->getQuery();

#test 

foreach my $u (@users) { 
	my $uid = $query->getUserId($u);
	warn $uid;
	#next;
	addUser($buffer->dbh,$project->id,$uid);
}

sub addUser{
	my ($dbh,$project_id,$user_id) = @_;
	my $query = qq{
		insert into PolyprojectNGS.user_projects  (PROJECT_ID,USER_ID) values ($project_id,$user_id);
	};
	$dbh->do($query) ;
}

