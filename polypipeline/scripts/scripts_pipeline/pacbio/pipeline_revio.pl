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

GetOptions(
	'project=s'   => \$project_name,
	"log=s" =>\$log_file,
	"vcf=s" => \$final_vcf,
	"patient=s"=>\$patient_name,
	"fork=s" =>\$fork,
);
my $date = `date`;
chomp($date);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatient($patient_name);

my $cmd_align = qq{$Bin/pbmm2.pl -project=$project_name -patient=$patient_name -fork=$fork};
system($cmd_align);
warn "align ok ";
my $cmd_deepvariant = qq{$Bin/deepvariant.pl -project=$project_name -patient=$patient_name -fork=$fork};
system($cmd_deepvariant);
warn "deepvariant ok ";
my $cmd_sawfish = qq{$Bin/sawfish.pl -project=$project_name -patient=$patient_name -fork=$fork};
system($cmd_sawfish);
warn "sawfish ok";
