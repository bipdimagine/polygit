#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
#use Set::IntSpan;
use GenBoNoSqlLmdb;
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
GetOptions(
	'project=s' => \$project_name,
	'patient=s' => \$patient_name,
	'bamin=s'=> \$bamin,
	'bamout' => \$bamout,
	'fork=s' => \$fork,
);
my $size = -s $bamin;
my $tmp = "/tmp/";
$tmp = "/data-beegfs/tmp/" if ($size/1024**2) > 80;
my $gatk = $buffer->software("gatk4");
warn $tmp;
my $cmd = qq{$gatk  --conf 'spark.local.dir=$tmp' --spark-runner LOCAL SortSamSpark -I $bamin -O $bamout};
system ($cmd);