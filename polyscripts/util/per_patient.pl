#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
use GBuffer;
use Set::IntSpan::Fast::XS;
use Getopt::Long;
my $fork = 1;
my $cmd;
my ($project_name, $patient_name);
my $file;
my $run = 0;
my $url;
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'cmd=s'  => \$cmd,
	'file=s'  => \$file,
	'run=s'  => \$run,
	'url=s'  => \$url,
);


#unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }
my $t = `hostname`;
my $nbErrors = 0;
my @pp;
if($file){
 @pp = `cat $file`;
chomp(@pp);
}
else {
	push(@pp,split(",",$project_name));
}
foreach my $project_name (@pp){
my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name);
if ($run == 1){
	
foreach my $patient (@{$project->getPatients}){
	system($cmd." -patient=".$patient->name()." -project=$project_name -fork=$fork") ;
	
}
}
else{
foreach my $patient (@{$project->getPatients}){
	if ($url){
		print $cmd." patients=".$patient->name()." project=$project_name >/dev/null\n";
	}
	else {print $cmd." -patient=".$patient->name()." -project=$project_name -fork=$fork\n";
	}
	
}
}
}