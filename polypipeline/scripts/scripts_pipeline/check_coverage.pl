#!/usr/bin/perl
use CGI qw/:standard :html3/;
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";

use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use JSON;
use check_utils;

my $projectName;
GetOptions(
	'project=s'   => \$projectName,
);
die("\n\nERROR: -project option mandatory. Exit...\n\n") unless ($projectName);

my @lHash;
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
foreach my $pat (@{$project->getPatients}){
	my $patName = $pat->name();
	my $hashCov = $pat->coverage();
	my $hash;
	$hash->{'patient'} = $patName;
	$hash->{'cov_1x'} = $hashCov->{'1x'};
	$hash->{'cov_5x'} = $hashCov->{'5x'};
	$hash->{'cov_15x'} = $hashCov->{'15x'};
	$hash->{'cov_30x'} = $hashCov->{'30x'};
	$hash->{'cov_mean'} = $hashCov->{'mean'};
	push(@lHash, $hash);
}
my $jsonFile = $project->getProjectPath().'../'.$projectName.'.resume.json';
my $jsonRes = check_utils::createJson(\@lHash, $jsonFile);
open(JSON_FILE, ">$jsonFile");
print JSON_FILE $jsonRes;
close(JSON_FILE);
exit(0);