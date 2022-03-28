#!/usr/bin/perl
use CGI qw/:standard :html3/;
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use colored;
use lib $Bin;
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use List::MoreUtils qw(part);
use File::Temp;
use JSON::XS;
use check_utils; 
use Statistics::Zscore; 
use IO::Prompt;
my $runs_id;
my $end_ext = "uni";
my $details;
my $fork = 1;
my $vcf_file;
my $log_file;
my $fast;
my $compute;
my $projectsName;
GetOptions(
	'runs=s'		=> \$runs_id,
	'run=s'		=> \$runs_id,
	'compute=s'		=> \$compute,
	'project=s' =>\$projectsName,
	'projects=s' =>\$projectsName,
);

my $date = `date`;

chomp($date);
if ($log_file){
	open (STDOUT,">>".$log_file);
}

colored::stabilo('blue',"START QUALITY CHECK FOR RUN $runs_id");
die("\n\nERROR: -project option mandatory. Exit...\n\n") unless ($runs_id);

my $buffer = GBuffer->new();
my $dbh = $buffer->dbh();
my $hprojects;
if ($projectsName){
	map{$hprojects->{$_} ++} split(",",$projectsName);
}
else {
foreach my $run_id (split(",",$runs_id)){
my $sql = qq{SELECT pj.name as name  FROM PolyprojectNGS.patient p, PolyprojectNGS.projects pj  where run_id =$run_id and pj.project_id = p.project_id group by pj.name};
my $sth = $dbh->prepare($sql);
$sth->execute();
	#my $s = $sth->fetchall_arrayref();	
	while ( my $v  = $sth->fetchrow_hashref() ){
		$hprojects->{$v->{name}} ++;
		
	}
}
}

if ($compute){
foreach my $name (keys %$hprojects){
		colored::stabilo('cyan',"START QUALITY CHECK FOR project $name");
		system ("$Bin/quality_check_project.pl -project=$name")
}
}
foreach my $name (keys %$hprojects){
		colored::stabilo('cyan',"START QUALITY CHECK FOR project $name");
		my $buffer1 = GBuffer->new();
		 my $project = $buffer1->newProject( -name => $name );
		 my $restart;
		 if (-e $project->getPedigreeFile()){
		 	die() unless check_utils::printErrorPedigree($project);
		 	my $mtime_ped = (stat $project->getPedigreeFile())[9];
		 	my $mtime_stat = (stat $project->getStatsFile())[9];
		 	$restart =1 if  ($mtime_ped > $mtime_stat);
		 }
		 $restart=1 unless -e $project->getStatsFile();
		 if ($restart)
		{
		 	colored::stabilo('cyan',"Change in pedigree detected restart quality for this project");
		 	system ("$Bin/quality_check_project.pl -project=$name");
		 }
		 	 if (-e $project->getPedigreeFile()){
		 		die() unless check_utils::printErrorPedigree($project);
		 	 }
		 	 
		 my $statistics = $project->stats();
		 check_utils::printTable($statistics,$project->name());
		 prompt("next project <return>",-d);
}
exit(0);