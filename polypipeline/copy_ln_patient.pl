#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo/lib/";
use lib "$Bin/..//GenBo/lib/obj-nodb/";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use lib "$Bin/..//GenBo/lib/GenBoDB/writeDB";
use lib "$Bin/packages";
use GBuffer;
use GenBoProject;
use Getopt::Long;
use Data::Dumper;
#use illumina_util;
use IO::Prompt;
use Sys::Hostname;
use File::Find::Rule;
 use File::Basename;

use colored; 
use file_util;
use check_utils;
use Text::Table;
use Term::Twiddle;
my $project_name;
my $project_name2;

my $filename;

my $name;
my $patient_name;
my $patient_name2;
my $steps_name;
my $force;
my $type;
my $fastq_ext;
GetOptions(
	'project=s' => \$project_name,
	'patient=s' => \$patient_name,
	'project2=s' => \$project_name2,
	'patient2=s' => \$patient_name2,
);
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $p = $project->getPatient($patient_name);
my $project2 = $buffer->newProject( -name => $project_name2 );
my $p2 = $project2->getPatient($patient_name2);
#$project2->makePath();
#$project2->getCoverageDir();
	my $dir1 = $project->getProjectPath();
  my @files = File::Find::Rule->file()
                                  ->name( "$patient_name.*" )
                                  ->in( $dir1);
  push(@files, File::Find::Rule->file()
                                  ->name( "$patient_name"."_aberrations.bed.gz" )
                                  ->in( $dir1));        
   push(@files, File::Find::Rule->file()
                                  ->name( "$patient_name"."_bins.bed.gz" )
                                  ->in( $dir1));       
foreach my $f (@files){
	my $target_file = $f;
	
	$target_file =~ s/$project_name/$project_name2/;
	$target_file =~ s/$patient_name/$patient_name2/;
	my $dir = dirname($target_file);
	system("mkdir -p $dir && chmod g+rwx $dir ") unless -e $dir;
	my $cmd = "ln -f -s $f $target_file";
	system($cmd);
}
                                