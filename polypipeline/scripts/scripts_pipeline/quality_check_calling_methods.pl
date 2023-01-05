#!/usr/bin/perl
use CGI qw/:standard :html3/;
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use colored;
use lib $Bin;
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use List::MoreUtils qw(part);
use File::Temp;
use JSON::XS;
use check_utils; 
use Statistics::Zscore; 
require "$Bin/../../../GenBo/lib/obj-nodb//packages/cache/polydiag/utility.pm";
require "$Bin/../../../polypipeline/scripts/scripts_pipeline/lib/quality_check.pm";

my $projectName;
my $end_ext = "uni";
my $details;
my $fork = 1;
my $vcf_file;
my $log_file;
my $noclean;
my $chr_names;
my $family_test;
GetOptions(
	'project=s'		=> \$projectName,
);

print "Quality Check - Calling Methods found in Cache\n\n";

my $buffer1 = new GBuffer;
my $project1 = $buffer1->newProjectCache( -name => $projectName );
quality_check::check_calling_methods_in_cache($project1);