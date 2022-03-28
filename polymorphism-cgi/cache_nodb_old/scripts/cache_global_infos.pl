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



warn "\n### CACHE: global infos step\n";
my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name );
my $hashInfos = $project->infosProject();
my ($hash, $hash_captures, $hash_align, $hash_calling);
$hash->{global_infos}->{name} = $hashInfos->{name};
$hash->{global_infos}->{id} = $hashInfos->{id};
$hash->{global_infos}->{description} = $hashInfos->{description};
$hash->{global_infos}->{creation_date} = $hashInfos->{creation_date};
$hash->{global_infos}->{project_type} = $hashInfos->{projectType};
$hash->{global_infos}->{project_type_id} = $hashInfos->{projectTypeId};
$hash->{global_infos}->{dbname} = $hashInfos->{dbname};
$hash->{analyse}->{build} = $project->version();
$hash->{analyse}->{exome} = $project->isExome();
$hash->{analyse}->{genome} = $project->isGenome();
$hash->{analyse}->{diagnostic} = $project->isDiagnostic();
foreach my $capture (@{$project->getCaptures()}) {
	$hash_captures->{$capture->name()} = undef;
}
my @lCapturesName = sort keys %$hash_captures;
$hash->{analyse}->{capture} = \@lCapturesName;
my @lTypes = ('evs', '1000genomes', 'dbsnp', 'prediction_matrix', 'cosmic', 'exac');
foreach my $type (@lTypes) {
	my $root_dir = $project->buffer->config->{public_data}->{$project->version()};
	my $extend =  $project->buffer->config->{kyoto}->{$type};
	my $dir = $root_dir.$extend;
	my $dir_abs_path = abs_path($dir);
	my @lTmp = split('/', $dir_abs_path);
	$hash->{analyse}->{version}->{$type} = $lTmp[-1];
}
foreach my $patient (@{$project->getPatients()}) {
	foreach my $method (@{$patient->callingMethods()}) { $hash_calling->{$method} = undef; }
	foreach my $method (@{$patient->alignmentMethods()}) { $hash_align->{$method} = undef; }
	foreach my $file (@{$patient->getVariationsFiles()}) {
		$hash->{check}->{vcf}->{variations}->{$patient->name()}->{$file} = md5_hex($file);
	}
	foreach my $file (@{$patient->getIndelsFiles()}) {
		$hash->{check}->{vcf}->{indels}->{$patient->name()}->{$file} = md5_hex($file);
	}
}
my @lMethodsAlign = sort keys(%$hash_align);
my @lMethodsCalling = sort keys(%$hash_calling);
$hash->{analyse}->{alignment} = \@lMethodsAlign;
$hash->{analyse}->{calling} = \@lMethodsCalling;
$hash->{analyse}->{cache}->{dejavu} = strftime '%Y-%m-%d', localtime;
unless (-d $project->getCacheBitVectorDir()) {
	my $cmd1 = "mkdir ".$project->getCacheBitVectorDir();
	my $cmd2 = "chmod 777 ".$project->getCacheBitVectorDir();
	`$cmd1`;
	`$cmd2`;
}
my $freeze_infos = $project->getCacheBitVectorDir().'/global_infos.freeze';
`rm $freeze_infos` if (-e $freeze_infos);
store($hash, $freeze_infos);
`chmod 777 $freeze_infos`;
