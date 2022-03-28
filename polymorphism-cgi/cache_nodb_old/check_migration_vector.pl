#!/usr/bin/perl
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($RealBin);
use lib "$RealBin/../GenBo";
use lib "$RealBin/../GenBo/lib/obj-nodb";
use GBuffer;
use Data::Dumper;
use Getopt::Long;

my ($fork, $year);
GetOptions(
	'fork=s' => \$fork,
	'year=s' => \$year,
);

unless ($fork) {
	warn "\n\nERROR: -fork option mandatory... Die...\n\n";
	die;
}
unless ($year) {
	warn "\n\nERROR: -year option mandatory... Die...\n\n";
	die;
}

my @lProjName;
my $buffer = GBuffer->new();
my @projectAll = @{$buffer->listProjects()};

my $path_vector = '/data-beegfs/polycache/';
my @lRef = ('EBV', 'EBV95-8', 'HG18', 'HG19', 'HG19b', 'HG38', 'MM38', 'MT', 'RN6');

my $to_do_all = 0;
my $to_do_error = 0;
my $to_do_not_run = 0;
my $to_do_ok = 0;

foreach my $project_name (sort @projectAll) {
	next unless ($project_name =~ /NGS$year/);
	$to_do_all++;
	my $exist_dir = 0;
	my $in_error = 0;
	my $ref_found;
	foreach my $ref (@lRef) {
		if (-d $path_vector.'/'.$ref.'/'.$project_name) {
			$exist_dir = 1;
			$ref_found = $ref;
		}
		if (-d $path_vector.'/'.$ref.'/'.$project_name.'_ERROR_1') {
			$to_do_error++;
			$in_error = 1;
		}
	}
	if ($in_error == 1) { next; }
	if ($exist_dir == 0) {
		$to_do_not_run++;
		next;
	}
	my $cmd = "/software/polyweb/poly-src/polymorphism-cgi/cache_nodb/cache_v2.pl -fork=$fork -chr=all -type=check_polyquery -project=$project_name";
	`$cmd`;
	if (-d $path_vector.'/'.$ref_found.'/'.$project_name) {
		$to_do_ok++;
	}
	elsif (-d $path_vector.'/'.$ref_found.'/'.$project_name.'_ERROR_1') {
		$to_do_error++;
	}
}


warn "\n\n";
warn 'Nb ALL:     '.$to_do_all."\n";
warn 'Nb ERROR:   '.$to_do_error."\n";
warn 'Nb NOT RUN: '.$to_do_not_run."\n";
warn 'Nb OK:      '.$to_do_ok."\n";
warn "\n";


