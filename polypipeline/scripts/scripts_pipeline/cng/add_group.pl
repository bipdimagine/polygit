#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use FindBin qw($RealBin);
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
my $group;
my $set;
my $name;
my $solo;
GetOptions(
	'project=s' => \$project_name,
	'set=s' => \$set,
	'name=s' => \$name,
	'group=s' => \$group,
	'solo=s' => \$solo,
);
my @list;
my $dir2 = "";
$set = "set".$set unless $set=~/set/;
if($solo){
	$dir2 ="solo";
}
if($project_name) {
	push(@list,$project_name);
}
elsif ($set && $name)  {
	$set = "set".$set unless($set=~/set/);
 @list = `cat $RealBin/../../../../defidiag/project/$name/$set.txt`;
}
else {
	die("cmd -name= -set  or -project    (-group if not the same as -name)");
}
chomp(@list);
$group = $name unless $group;
die() unless(@list);

foreach my $project_name (@list){
my $buffer = GBuffer->new();
my $dbh = $buffer->dbh;
$dbh->{AutoCommit} = 0;
my $project = $buffer->newProject( -name => $project_name );
my $query = $buffer->getQuery();
warn "add $project_name to grpoup $group"; 
#test 
my $gid = getGroupId($buffer->dbh,$group);
die() unless $gid;
addGroup($buffer->dbh,$project->id,$gid);
$dbh->commit();
}
#$dbh->commit();
exit(0);
sub addUser{
	my ($dbh,$project_id,$user_id) = @_;
	my $query = qq{
		insert into PolyprojectNGS.user_projects  (PROJECT_ID,USER_ID) values ($project_id,$user_id);
	};
	$dbh->do($query) ;
}

sub addGroup {
        my ($dbh,$projectid,$groupid) = @_;
        my $sql = qq{    
                insert into PolyprojectNGS.ugroup_projects (ugroup_id,project_id) values ($groupid,$projectid);
        };
        return ($dbh->do($sql));
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
sub getGroupId {
my ($dbh,$name) = @_;
my $query = qq{
SELECT
*
FROM bipd_users.UGROUP
where NAME='$name';
};
my $sth = $dbh->prepare($query);
$sth->execute();
my $s = $sth->fetchrow_hashref();
return unless $s;
return $s->{UGROUP_ID};
}

