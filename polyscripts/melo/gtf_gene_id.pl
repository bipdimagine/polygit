#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/";
use lib "$Bin/../../GenBo/lib/GenBoDB";
use lib "$Bin/../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../polypipeline/packages";
use Logfile::Rotate;
use Cwd;
use PBS::Client;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
use MooseX::Method::Signatures;
#use bds_steps;   
use file_util;
use Class::Inspector;
use Digest::MD5::File ;
use GBuffer;
use GenBoProject;
use colored; 
use Config::Std;
use Text::Table;
use Time::Local 'timelocal';
use Cwd 'abs_path';
use File::Temp qw/ tempfile tempdir /;
use Term::Menus;
use Proc::Simple;
use Storable;
use JSON::XS;
use XML::Simple qw(:strict);

# Correct the gtf obtained from the convertion from the gff
# by adding the transcript_id field (eq to the gene_id value) when missing

my $gtf;
my $output;

GetOptions(
	'gtf=s'		=> \$gtf,
	'out=s'		=> \$output,
);

die ("'$gtf' does not exists") unless (-e $gtf);
$output = $gtf =~ s/\.gtf/_with_gene_id\.gtf/r unless ($output);
warn abs_path($output);

open(GTF, '<', $gtf) or die ("Can not open '$gtf': $!");
open(OUT, '>', $output) or die ("Can not open '$output': $!");
while(my $line = <GTF>) {
	chomp $line;
	print OUT $line;
	if ($line =~ /^#/) {
		print $line."\n";
		next;
	};
	my @fields = split("\t", $line);
#	warn Dumper \@fields;
	my @attributes = split(/; ?/, $fields[-1]);
#	warn Dumper \@attributes;
	my @transcript_id = grep(/^transcript_id/, @attributes);
#	warn @transcript_id;
	my $gene_id = grep(/^gene_id/,@attributes);
	unless ($gene_id) {
		print OUT '; ' unless ($line =~ /; ?$/);
		print OUT ' ' unless ($line =~ / $/);
		print OUT $transcript_id[0] =~ s/^transcript_id/gene_id/r;
	}
	print OUT "\n";
}
close(GTF);
close(OUT);
