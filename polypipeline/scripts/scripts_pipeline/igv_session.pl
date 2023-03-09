#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/";
use lib "$Bin/../../../GenBo/lib/GenBoDB";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../../../GenBo/lib/kyoto/";
use lib "$Bin/../../../GenBo/lib/GenBoDB/writeDB";
use lib "$Bin/../../packages";
use PBS::Client;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
use MooseX::Method::Signatures;
#use bds_steps;   
use GBuffer;
use GenBoProject;


  
my $bin_cecile=qq{$Bin/scripts/scripts_db_polypipeline};
my $bin_script_pipeline = qq{$Bin/scripts/scripts_pipeline};


my $projectName;
my $patients_name;
my $step;
my $bcl_dir;
my $run;
my $feature_ref;
my $no_exec;
my $aggr_name;
my $lane;


my $limit;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
);




my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
#$patients_name = "all" unless $patients_name;
$patients_name = $project->get_only_list_patients($patients_name);
my $out = $project->getProjectRootPath();
my $file = $out."/".$project->name()."_igv_session.xml";

open (IGV, ">$file");
print IGV '<Session genome="hg38" hasGeneTrack="true" hasSequenceTrack="true" locus="All" version="8">\n';
print IGV "<Resources>\n";

foreach my $patient (@{$patients_name}) {
	my $name=$patient->name();
	print IGV '<Resource path="http://www.polyweb.fr/NGS/'.$project->name."/HG38/align/hisat2/".$patient->name.".bam\"/>\n";
}
print IGV "</Resources>\n";
print IGV "</Session>\n";
