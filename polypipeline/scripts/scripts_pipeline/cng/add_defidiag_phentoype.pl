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

update_defidiag($dbh,$project->id);
update_phenotype($dbh,$project->id);
warn "WARNING UPDATE defidiag and phenotype";
$dbh->commit();
	system("add_calling_method.sh -project=$project_name -method=canvas,manta,wisecondor");

sub update_defidiag{
	my ($dbh,$project_id) = @_;
	my $query = qq{
		update  PolyprojectNGS.projects  set validation_db ='defidiag'  where project_id=$project_id;
	};
	warn $query;
	$dbh->do($query) ;
}

sub update_phenotype {
	my ($dbh,$project_id) = @_;
	my $query = qq{
			insert into PolyPhenotypeDB.phenotype_project  (PROJECT_ID,phenotype_id) values ($project_id,1);
	};
	$dbh->do($query) ;
}
