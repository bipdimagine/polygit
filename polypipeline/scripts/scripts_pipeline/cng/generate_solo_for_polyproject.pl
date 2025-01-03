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
my $list;
my $group;
my $set;
my $name;

GetOptions(
	'set=s' => \$set,
	'name=s' => \$name,
	'project=s' =>\$project_name
);

my @projects;
if($set){
	 @projects = `cat /software/polyweb/poly-disk/poly-src/defidiag/project/$name/set$set.txt`;
	chomp(@projects);
}
elsif ($project_name){
	push(@projects,$project_name);
}
die() unless @projects;
my $nb =0;
print "Patient\tFamily\tBC\tSex\tStatus\n";
foreach my $project_name (@projects){
	next unless $project_name;
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );

foreach my $patient (@{$project->getPatients}){
	next unless $patient->isChild();
	print $patient->name()."\t".$patient->getFamily->name."\t".$patient->barcode()."\t".$patient->sex."\t".$patient->status."\n";
	$nb ++;
}
}
warn $nb;
