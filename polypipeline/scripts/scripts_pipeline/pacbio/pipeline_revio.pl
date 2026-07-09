#! /usr/bin/env perl
use FindBin qw($Bin);
use strict;
use FindBin qw($Bin);
#use lib "/software/polyweb/poly-disk/poly-src/polygit/GenBo/lib/obj-nodb/";
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages/";
use Data::Dumper;
use Getopt::Long;
use GBuffer;
use autodie qw(system);
 
 my $project_name;
my $final_vcf;
my $log_file;
my $patient_name;
my $fork;
my $force;
GetOptions(
	'project=s'   => \$project_name,
	"log=s" =>\$log_file,
	"vcf=s" => \$final_vcf,
	"patient=s"=>\$patient_name,
	"fork=s" =>\$fork,
	"force=s"=>\$force,
);
my $date = `date`;
chomp($date);
$fork =256 unless $fork;
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatient($patient_name);
my $stforce ="";
$stforce = "force=1" if $force;
my $cmd_align = qq{$Bin/pbmm2.pl -project=$project_name -patient=$patient_name -fork=$fork $stforce};
warn $cmd_align;
my $dir_out= $project->getAlignmentPipelineDir($patient->name);
my $bam_out = $dir_out."/".$patient->name.".bam";
system($cmd_align);
warn "align ok ";
my $cmd_cram = qq{$Bin/bam2cram.pl -project=$project_name -patient=$patient_name -fork=$fork $stforce};
system($cmd_cram);
warn "align ok ";

my $cmd_deepvariant = qq{$Bin/deepvariant.pl -project=$project_name -patient=$patient_name -fork=$fork $stforce};
system($cmd_deepvariant);
warn "deepvariant ok ";
my $cmd_sawfish = qq{$Bin/sawfish.pl -project=$project_name -patient=$patient_name -fork=$fork $stforce};
system($cmd_sawfish);
warn "sawfish ok";

my $cmd_coverage = qq{perl $Bin/../coverage_genome.pl -project=$project_name -patient=$patient_name -fork=$fork -$stforce};
system($cmd_coverage);
my $cmd_wsiecondor = qq{$Bin/wisecondor.pl -project=$project_name -patient=$patient_name -fork=$fork $stforce};
system($cmd_wsiecondor);
my $cmd_wsiecondor2 = qq{$Bin/calling_wisecondor.pl -project=$project_name -patient=$patient_name -fork=$fork $stforce};
system($cmd_wsiecondor2);

