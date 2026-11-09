#!/usr/bin/env perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/";
use lib "$Bin/../../../GenBo/lib/GenBoDB";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../../../GenBo/lib/kyoto/";
use lib "$Bin/../../../GenBo/lib/GenBoDB/writeDB";
use lib "$Bin/../../packages";
#use PBS::Client;
use Getopt::Long;
use Data::Dumper;
#use IO::Prompt;
#use Sys::Hostname;
#use Parallel::ForkManager;
#use Term::ANSIColor;
#use Moose;
#use MooseX::Method::Signatures;
#use bds_steps;   
use GBuffer;
use Carp;


my $projectName;
my $patients_name;
my $version;

GetOptions(
	'project=s'		=> \$projectName,
	'patients=s'	=> \$patients_name,
	'version=s'		=> \$version,
);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName , -version => $version);
$patients_name = $project->get_only_list_patients($patients_name);
my $out = $project->getProjectRootPath();
my $file = $out."/".$project->name()."_igv_session.xml";

my $release = $project->annotation_genome_version;
my $genome;
$genome = 'hg38_1kg' if ($release =~ /^HG38/);
$genome = 'hg19' if ($release =~ /^HG19/);
$genome = 'mm39' if ($release =~ /^MM39/);
$genome = 'mm10' if ($release =~ /^MM10/);
confess ("No session genome ") unless ($genome);

open (IGV, ">$file");
print IGV qq{<Session genome="$genome" hasGeneTrack="true" locus="All" version="8">}."\n";
print IGV "<Resources>\n";

foreach my $patient (@{$patients_name}) {
	my $name=$patient->name();
	my $bam =$patient->getBamFile();
	my $bw = $bam =~ s/.(cr|b)am$/\.bw/r;
	$bw = $patient->getBigWigFile unless (-e $bw);
	$bam =~ s/(\/data-bipd)?\/data-(isilon|pure)\/sequencing\/ngs/https:\/\/www.polyweb.fr\/NGS\//;
	print IGV '<Resource path="'.$bam.'"'."/>\n";
	
	if (-e $bw) {
		$bw =~ s/\/data-(isilon|pure)\/sequencing\/ngs/https:\/\/www.polyweb.fr\/NGS\//;
		print IGV '<Resource path="'.$bw.'"'."/>\n";
	}
}
print IGV "</Resources>\n";
print IGV "</Session>\n";

print "IGV session: ".$file."\n";
