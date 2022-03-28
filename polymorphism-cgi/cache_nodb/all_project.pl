#!/usr/bin/perl

use FindBin qw($Bin);
use FindBin qw($RealBin);
use lib "$RealBin/../GenBo";
use lib "$RealBin/../GenBo/lib/GenBoDB";
use lib "$RealBin/../GenBo/lib/obj-nodb";
use lib "$RealBin";
use lib "$RealBin/../packages/cache";
use Carp;
use strict;
use Data::Dumper;
use GBuffer; 
use Getopt::Long;
use Carp;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use Vcf;
use KyotoCabinet;
use Storable qw/thaw freeze/;
use threads;
use Thread::Queue;
use JSON::XS;
use Parallel::ForkManager;
use Try::Tiny;

my ($fork, $projectName, $force,$year);
my $cmd;
GetOptions(
	'fork=i'  	=> \$fork,
	'project=s' => \$projectName,
	'force=s' 	=> \$force,
	'year=s' => \$year,
);
unless ($projectName) { die("\n\nERROR !! No project name specified (you can use 'all') !! Exit...\n\n"); }
unless ($fork) { $fork = 1; }


my @projectToDo;
my @projectAll;
my $buffer1 = GBuffer->new();	
unless ($projectName eq 'all') { @projectAll = split(',', $projectName); }
else { @projectAll = @{$buffer1->listProjects()}; }
open (ERROR,">$year.error.chr.txt");
open (ERROR_project,">$year.project.error.chr.txt");
open (OK,">$year.ok.chr.txt");

#my @chromosomes = (1..22,'X','Y','MT');

my $nb_node = 10;
$nb_node = 11 if $year eq "2011";
$nb_node = 12 if $year eq "2012";
$nb_node = 13 if $year eq "2013";
$nb_node = 14 if $year eq "2014";
$nb_node = 15 if $year eq "2015";
my $pm = new Parallel::ForkManager($fork);
my $nbp =0;
my @yproject = grep{$_ =~ /$year/}@projectAll ;
my $nb_project = scalar(@yproject);
foreach my $project_name (@projectAll){

	next if $project_name !~/NGS/;
	#next if $project_name eq "NGS2015_0775" or $project_name eq "NGS2015_0629" or $project_name eq "NGS2015_0608";
	if ($year) {
		next unless $project_name =~ /$year/;
		
	}
		
	my $buffer = GBuffer->new();
	my $project = $buffer->newProject(-name=>$project_name);
	next if $project->version() eq "HG18";
	next unless $project->isDiagnostic;
#	next unless $project->name =~ /NGS2015_09/;
	$nbp++;
	print "running $project_name : $nbp/".$nb_project."\n";
	my $cmd = qq{ run "perl $Bin/cache.pl -fork=$fork  -project=}.$project->name()." -type=polydiag\"";
	#my $cmd2 = qq{ $Bin/cache.pl -fork=$fork  -project=}.$project->name();
#	warn $cmd;
	system ("$cmd");
	#my $pname = $project->name();
#	try {
#		
#		my @toto = `perl $Bin/../json_output_nodb/./genes_json.pl project=$pname stat=all filter_type_variation=dbsnp+dbsnp_1p+dbsnp_none+1000genomes+1000genomes_1p+evs+evs_1p+silent+intergenic+intronic+pseudogene+filter_composite+filter_recessif+filter_family+filter_denovo+filter_strict_denovo+loh+dbl_evt+filter_dominant+polyphen_sift`;
#	}
#	catch {
#		warn "ERROR $pname";
#		print ERROR_project "$pname \n";
#		return;
#	}

		warn "end $project_name \n";

}

warn $nbp;
exit(0);