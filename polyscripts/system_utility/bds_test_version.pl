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

my @extra_args = @ARGV;


my $path;
 $path = "/data-pure/software/poly-src/duckdb/polygit/";


if ($type eq "pipeline" or $type eq "cache"){
	my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => $project_name);
 if( $project->gencode_version < 43 ){
	$path = "/software/polyweb/poly-disk/poly-src/polygit/";

}


}
warn $type;

if ($type eq "pipeline"){
	my $cmd = "$path//polypipeline/bds_pipeline.pl -project=$project_name ".join(" ",@ARGV);
	system($cmd);
}
elsif ($type eq "cache"){
	my $cmd = "$path//polypipeline/bds_cache.pl -project=$project_name ".join(" ",@ARGV);
	system($cmd);
	}
elsif ($type eq "dragen_pipeline") {
	
	my $cmd = "$path/polypipeline/dragen/dragen_pipeline.pl -project=$project_name ".join(" ",@ARGV);
	system($cmd);
	 
}	
elsif ($type eq "dragen_demultiplex") {
	
	my $cmd = "$path/polypipeline/dragen/scripts/dragen_demultiplex.pl -project=$project_name ".join(" ",@ARGV);
	system($cmd);
	 
}	

elsif ($type eq "run_cluster") {
	my $cmd = "$path/polypipeline/polyscripts/system_utility/run_cluster.pl -project=$project_name ".join(" ",@ARGV);
	system($cmd);
	 
}	
elsif ($type eq "launch_dejavu") {
	my $cmd = "$path/polypipeline/scripts/scripts_pipeline/dejavu/hg38/variant/duckdb/launch_dejavu_rocks_duck.pl -project=$project_name ".join(" ",@ARGV);
	system($cmd);
	 
}
else {
	die("unkown option $type");
}