#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use Parallel::ForkManager;
use Storable qw(store retrieve freeze thaw);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Cwd 'abs_path';
use Digest::MD5 qw(md5_hex);
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use Cache_Commons;



my $fork = 1;
my ($project_name, $chr_name);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }

warn "\n### CACHE: create local config\n";
my $buffer = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProject( -name => $project_name );
my $hash;
$hash->{ensembl_versions} = $project->buffer->config->{ensembl_versions};
$hash->{dbnsfp} = $project->buffer->config->{dbnsfp};
$hash->{ensembl} = $project->buffer->config->{ensembl};

my $file = $project->local_config_file();
open (FILE, '>'.$file);
foreach my $cat (sort keys %{$hash}) {
	print FILE "[$cat]\n";
	foreach my $key (sort keys %{$hash->{$cat}}) {
		my $value = $hash->{$cat}->{$key};
		print FILE $key.':'.$value."\n";
	}
	print FILE "\n";
}
close (FILE);
my $cmd = 'chmod 777 '.$file;
`$cmd`;
