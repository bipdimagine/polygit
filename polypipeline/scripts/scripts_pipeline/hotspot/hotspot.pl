#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages/";
use  File::Temp;
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;

use Storable qw(store retrieve freeze);
 use List::MoreUtils qw(natatime);
 
 my $project_name;
my $final_vcf;
my $log_file;
my $patient_name;
my $fork;

GetOptions(
	'project=s'   => \$project_name,
	"patient=s"=>\$patient_name,
);
my $date = `date`;
chomp($date);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $rocks = $project->hotspot_rocks("c");
foreach my $patient (@{$project->getPatients()}){
my $hotspot = $patient->hotspot();
$rocks->put($patient->id,$hotspot);
}
$rocks->close();
exit(0);