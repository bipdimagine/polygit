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
);


my $buffer = GBuffer->new();
my $dbh = $buffer->dbh;
$dbh->{AutoCommit} = 0;
my $project = $buffer->newProject( -name => $project_name );
system("./add_methods.pl -project=$project_name");
my $fams = $project->getFamilies();
shift(@$fams);
update_defidiag($dbh,$project->id);
update_phenotype($dbh,$project->id);
foreach my $fam (@$fams){
	
	my $pid = create_project($dbh,$fam->name);
	update_defidiag($dbh,$pid);
	update_phenotype($dbh,$pid);
	update_release($dbh,$pid);
	update_dbtable($dbh,$pid);
	my $patients = $fam->getMembers();
	foreach my $p (@$patients){
		warn $pid." ".$p->name();
		update_target_project($dbh,$pid,$p->id);
	}
	
}

$dbh->commit();

#test 

#my $res = queryPolyproject::getMethodFromName($buffer->dbh,@fieldM[$i]);

sub create_project {
my ($dbh,$description) = @_;
my $query = qq{CALL PolyprojectNGS.new_project("NGS",3,"$description");};
my $sth = $buffer->dbh()->prepare( $query );
$sth->execute();
my $res = $sth->fetchall_arrayref({});
sendError("No project created") if ( scalar(@$res) == 0 );
my $pid = $res->[0]->{project_id};
return $pid;
} 



sub update_defidiag{
	my ($dbh,$project_id) = @_;
	my $query = qq{
		update  PolyprojectNGS.projects  set validation_db = 'defidiag' where project_id=$project_id;
	};
	$dbh->do($query) ;
}

sub update_phenotype {
	my ($dbh,$project_id) = @_;
	my $query = qq{
			insert into PolyPhenotypeDB.phenotype_project  (PROJECT_ID,phenotype_id) values ($project_id,1);
	};
	$dbh->do($query) ;
}

sub update_release{
	my ($dbh,$project_id,$release) = @_;
	my $release_id = 929;
	
	my $query = qq{
		insert into PolyprojectNGS.project_release  (project_id,release_id,`default`) values ($project_id,$release_id,1);
	};

	$dbh->do($query) ;
}
sub update_dbtable{
	my ($dbh,$project_id,$release) = @_;
	my $release_id = 929;
	
	my $query = qq{
		insert into PolyprojectNGS.databases_projects  (project_id,db_id) values ($project_id,4);
	};

	$dbh->do($query) ;
}
sub update_target_project{
	my ($dbh,$project_id,$patient_id) = @_;
	
	#my $release_id = 929;
		my $query = qq{
		update  PolyprojectNGS.patient set  project_id_dest= $project_id where patient_id=$patient_id;
	};

	$dbh->do($query) ;
}

