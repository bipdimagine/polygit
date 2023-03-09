#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/obj-nodb";
use GBuffer;
use JSON;

my $cgi = new CGI();
my $type_project = $cgi->param('type_project');

my $bufferTMP = GBuffer->new();
my @lProj_names = @{$bufferTMP->listProjectsExomeForDejaVu()};

my $url_bam;
my $is_ok;
while (!$is_ok) {
	my $number_rand = rand(scalar(@lProj_names) - 1);
	my $project_name = $lProj_names[$number_rand];
	my $buffer = GBuffer->new();
	my $project = $buffer->newProject( -name => $project_name );
	
	next if $type_project eq 'exome' and not $project->isExome();
	next if $type_project eq 'genome' and not $project->isGenome();
	
	foreach my $patient (@{$project->getPatients()}) {
		next if $patient->isIll();
		next unless -e $patient->getBamFile();
		$url_bam = $patient->bamUrl();
		$is_ok = 1;
		last;
	}
	$project = undef;
	$buffer = undef;
}

my $hash;
$hash->{'url_bam'} = $url_bam;
print $cgi->header('text/json-comment-filtered');
print encode_json $hash;
exit(0);
