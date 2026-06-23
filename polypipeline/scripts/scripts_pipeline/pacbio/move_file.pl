#! /usr/bin/env perl
use FindBin qw($Bin);
use strict;
use FindBin qw($Bin);
#use lib "/software/polyweb/poly-disk/poly-src/polygit/GenBo/lib/obj-nodb/";
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages/";
use Getopt::Long;
use GBuffer;
use autodie qw(system);
 
 my $project_name;
my $final_vcf;
my $log_file;
my $patient_name;
my $fork;

GetOptions(
	'project=s'   => \$project_name,
	"patient=s"=>\$patient_name,
	"fork=s" =>\$fork,
);
my $date = `date`;
chomp($date);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $bcftools = $buffer->software("bcftools");
my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");
my $samtools = $buffer->software("tabix");

#
#
# move_bam 
#
#

my $dir_out= $project->getCallingPipelineDir($patient->name);
my $bam_out = $dir_out."/".$patient->name.".bam";
die() unless -e $bam_out;
$fork = 64 if $fork > 64;
my $cram_out = $patient->getCramFileName("pbmm2");
my $cmd = "samtools -C -T $ref -o $cram_out $bam_out -@ $fork && samtools index  $cram_out -@ $fork";

#
# move_bam 
#


