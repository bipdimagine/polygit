#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/";
use lib "$Bin/../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../GenBo/lib/obj-nodb/packages/";
use lib "$Bin/../../GenBo/lib/GenBoDB/writeDB/";
use lib "$Bin/packages";use Data::Dumper;
use Getopt::Long;
use Carp;
use IO::Prompt;
use colored;
use Cwd 'abs_path';
use Net::SFTP;
use Archive::Tar;

use GBuffer;
my $buffer = new GBuffer;


my $project_names;
my $patient_names;
my $file_types;
my $archive_dir = ".";
my $no_die;
my $archive_name;

GetOptions(
	'project=s'		=> \$project_names,		# Projet(s)
	'patients=s'	=> \$patient_names,		# Patient(s)
	'files=s'		=> \$file_types,		# fastq, bam, vcf or all
	'archive_dir=s'	=> \$archive_dir,
	'archive_name=s'=> \$archive_name,
	'nodie!'		=> \$no_die,			
);








