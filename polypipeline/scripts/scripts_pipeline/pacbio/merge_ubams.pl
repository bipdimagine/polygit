#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
#use Set::IntSpan;
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use Carp;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use String::ProgressBar;
use List::Util qw(sum);



 my $project_name;
 my $fork;
 my $callable_intspan_file;
 my $patient_name;
 #my $low_calling;
 my $method;
 
 
my $fork = 5;
GetOptions(
	'project=s'   => \$project_name,
	"patient=s" => \$patient_name,
);


my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatient($patient_name);

warn Dumper $patient->uBams_revio();

my $ubam = $patient->uBam_filename();
warn $ubam;

my $samtools = $buffer->software("samtools");
my $find = "SM:".$patient->name;
foreach my $bam (@{$patient->uBams_revio()}){
	
	open(my $fh, "samtools view -H $bam |") or die $!;

while (<$fh>) {

    next unless /^\@RG/;

    if (/\tSM:([^\t]+)/) {
    	my $sm = $1;
    	my $pname = $patient->name;
		if ($sm ne $patient->name ){
		print colored("$bam - $sm - $pname", 'red'), "\n";
		die();
		}

    }

}

close($fh);
	
}
warn $ubam;
my $cmd  = qq{$samtools merge $ubam -@ 10 }.join(" ",@{$patient->uBams_revio()});
system($cmd) unless -e $ubam;
  #system("$cmd");

 




