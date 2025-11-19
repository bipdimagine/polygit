#!/usr/bin/perl

use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
#use lib "/data-isilon/bipd-src/pnitschk/git-master/polygit/GenBo/lib/obj-nodb";

use GBuffer;
use Set::IntSpan::Fast::XS;
use Getopt::Long;
use Digest::MD5 qw(md5 md5_hex md5_base64);
my $fork = 1;
Getopt::Long::Configure("pass_through");
my $type;
my ($project_name, $patient_name,$vid);
my $file;

GetOptions(
	'project=s'    => \$project_name,
	'type=s'    => \$type,
);
my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => $project_name);
my @extra_args = @ARGV;
warn Dumper @ARGV;

my $path;
 if( $project->gencode_version > 42 ){
 $path = "/data-pure/software/poly-src/duckdb/polygit/";
warn "coucou";

}

else {
	$path = "/software/polyweb/poly-disk/poly-src/polygit/";
}

if ($type eq "pipeline"){
	my $cmd = "$path//polypipeline/bds_pipeline.pl -project=$project_name ".join(" ",@ARGV);
	system($cmd);
}
else ($type eq "cache"){
	my $cmd = "$path//polypipeline/bds_cache.pl -project=$project_name ".join(" ",@ARGV);
	system($cmd);
	}

