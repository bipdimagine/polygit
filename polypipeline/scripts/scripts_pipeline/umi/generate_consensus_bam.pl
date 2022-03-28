#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages/";

use GBuffer;
use GenBoProject;
use Getopt::Long;
use Carp;
use strict;
use Set::IntSpan::Fast::XS;
use Data::Dumper;
use Data::Printer;
use Storable qw(store retrieve freeze);
use file_util;
use Term::Menus;
use IO::Prompt;
use GenBoNoSqlLmdb;
use Bio::Seq::Quality;
use Bio::SeqIO::fastq;
use GenOO::Data::File::FASTQ;
use Parallel::ForkManager;
use JSON::XS;

my $buffer = GBuffer->new();
my $patient_name;
my $project_name;
my $fork;
my $bamin;
my $bamout;

GetOptions(
	'project=s' => \$project_name,
	'patient=s' => \$patient_name,
	'fork=s' => \$fork,
	'bamin=s' => \$bamin,
	'bamout=s' => \$bamout,
);



	
#my $patient_name = "2006N080737_LEBMI";
my $project = $buffer->newProject(-name=>"$project_name");
my $patient = $project->getPatient($patient_name);


my $method = $patient->alignmentMethod();
my $dir = $project->getAlignmentPipelineDir($method."/".$patient->name);

my $fgbio =  $buffer->software('java')." -Djava.io.tmpdir=/data-beegfs/tmp -jar ".$buffer->software('fgbio');
die() unless -e $bamin;
my $ref =  $project->genomeFasta();
my $cmd = qq{$fgbio GroupReadsByUmi -s adjacency -o /dev/stdout -i $bamin | $fgbio CallMolecularConsensusReads  -M 1 --read-group-id=$patient_name  -o /dev/stdout -i /dev/stdin |  $fgbio FilterConsensusReads  -M 1  -N 30 -r $ref  -o $bamout -i /dev/stdin };
warn $cmd;
system($cmd."&& touch $dir/consensus.ok");

die() unless -e "$dir/consensus.ok";

unlink "$dir/consensus.ok";

exit(0);
