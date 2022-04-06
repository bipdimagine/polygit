#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
#use Set::IntSpan;
use Carp;
use strict;
use Set::IntSpan::Fast::XS ;
use Data::Dumper;
use GBuffer; 
use Getopt::Long;
use Carp;
 use JSON::XS;
 use List::MoreUtils qw(natatime uniq);
 use Tabix;
 
 my $buffer = new GBuffer;
my $project_name= "NGS2017_1534";
my $fork;
my $dirin;
my $patient_name;
my $bamin;
my $bamout;
my $hg38;
GetOptions(
	'project=s' => \$project_name,
	'patient=s' => \$patient_name,
	'bamin=s'=> \$bamin,
	'bamout=s' => \$bamout,
	'fork=s' => \$fork,
);
my $size = -s $bamin;
my $tmp = "/tmp/";
$tmp = "/data-beegfs/tmp/" if ($size/1024**2) > 80;



my $bamsormadup = $buffer->software("bamsormadup");
my $bamcat = $buffer->software("bamcat");
my $samtools = $buffer->software("samtools");
my $cmd = qq{$bamcat $bamin | $bamsormadup threads=$fork  inputformat=bam  > $bamout && $samtools index $bamout };
system ($cmd);