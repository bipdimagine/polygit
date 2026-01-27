#!/usr/bin/perl

use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
#use lib "/data-isilon/bipd-src/pnitschk/git-master/polygit/GenBo/lib/obj-nodb";

use GBuffer;
use Set::IntSpan::Fast::XS;
use Getopt::Long;
my $fork = 1;
my $cmd;
my $file;


my $buffer = new GBuffer;
my $genome = "HG38";
$buffer->annotation_genome_version($genome);
my $version;# = 21;
unless ($version){
	
	my @versions = sort {$a <=> $b} keys %{$buffer->public_data};
	$version = $versions[-1];
	warn "backup latest version $version";
}

my @dir ;

foreach my $database( keys %{$buffer->public_data()->{$version}}){
	my $database_config_name = $database;
	my $root =  "/$genome/$database/".$buffer->public_data->{$version}->{$database_config_name}->{config}->{version}."/".$buffer->public_data->{$version}->{$database_config_name}->{config}->{dir};
	if (exists $buffer->public_data->{$version}->{$database_config_name}->{config}->{semantic}){
			warn  "/semantic/$database/".$buffer->public_data->{$version}->{$database_config_name}->{config}->{version}."/".$buffer->public_data->{$version}->{$database_config_name}->{config}->{dir};
			next;
	}
	next unless -e $buffer->public_data_root.$root;
	next if $database =~/hgmd/; 
	push (@dir,".".$root);
}

my $cmd = "cd ".$buffer->public_data_root." && tar -cvf /data-pure/public-data/$genome.$version.tar ".join(" ",@dir);
warn $cmd;
die();
system($cmd);