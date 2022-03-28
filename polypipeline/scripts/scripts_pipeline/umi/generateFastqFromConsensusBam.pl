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
my $out;

GetOptions(
	'project=s' => \$project_name,
	'patient=s' => \$patient_name,
	'fork=s' => \$fork,
	'bamin=s' => \$bamin,
	'bamout=s' => \$bamout,
	'out=s' => \$out,
);



	
#my $patient_name = "2006N080737_LEBMI";
my $project = $buffer->newProject(-name=>"$project_name");
my $patient = $project->getPatient($patient_name);

my $method = $patient->alignmentMethod();
my $dir = $project->getAlignmentPipelineDir($method."/".$patient->name);
my $file = "$dir/list_combined.txt";
if($out){
	print $file;
	exit(0);
}



my $pm = new Parallel::ForkManager(int($fork/4));

my $fgbio =  $buffer->software('java')." -Djava.io.tmpdir=/tmp -jar ".$buffer->software('fgbio');
my $picard =  $buffer->software('java')." -Djava.io.tmpdir=/tmp -jar ".$buffer->software('picard');
my $pigz =  $buffer->software('pigz');
my $bwa = $buffer->software("bwa");
my $sambamba = $buffer->software("sambamba");
my $samtools = $buffer->software("samtools");
my $index = $project->getGenomeIndex("bwa");#
$index .= "/all.fa";
die() unless -e $bamin;
my $ref =  $project->genomeFasta();
my $f1 = $dir."/".$patient->name."_R1.fastq";
my $f2 = $dir."/".$patient->name."_R2.fastq";

my @chrs = `samtools idxstats $bamin | cut -f1`;
chomp(@chrs);


foreach my $chr (@chrs){
	#next if $chr ne "chr21";
	my $pid = $pm->start and next;
	my $rname = $chr;
	$chr = "star" if $chr eq "*";
	my $bam3 = $dir."/$patient_name.$chr.align.bam";
	warn $bam3;
	$pm->finish() if -e $bam3;
	
	my $bam1 = $dir."/$patient_name.$chr.bam";
	my $f1 = $dir."/".$patient->name."._R1.$chr.fastq";
	my $f2 = $dir."/".$patient->name."_R2.$chr.fastq";
	my $cmd = qq{sambamba slice $bamin "$rname" > $bam1};
	warn $cmd;
	system($cmd) ;#unless -e $bam1;
	
	my $bam2 = $dir."/$patient_name.$chr.consensus.bam";
	my $cmd2 = qq{$fgbio GroupReadsByUmi -s adjacency -o /dev/stdout -i $bam1 | $fgbio CallMolecularConsensusReads  -M 1 --read-group-id=$patient_name  -o /dev/stdout -i /dev/stdin |  $fgbio FilterConsensusReads  -M 1  -N 30 -r $ref  -o $bam2 -i /dev/stdin};
	warn $cmd2;
	system($cmd2) ;#unless -e $bam2;
	#unlink $bam1;
	
	my $cmd3 = qq{samtools sort -n -@ 4 $bam2  -m 3G -T /data-beegfs/tmp/ | samtools fastq -1 $f1 -2 $f2 - && pigz -f -p 4 $f1 && pigz -f -p 4 $f2};
	#unlink $bam2;
	$f1.=".gz";
	$f2.=".gz";

	system($cmd3);# unless -e $f1;
	
	warn $f1;
	my $cmd4 = qq{$bwa  mem $index $f1 $f2 -t 4  | samtools view -S -b - > $bam3  2>/dev/null };
	warn $cmd4;
	system($cmd4) ;#unless -e $bam3 ;
	#unlink $f1;
	#unlink $f2;
	$pm->finish();
	
}



$pm->wait_all_children();

open (LIST,">$file");
foreach my $chr (@chrs){
	$chr ="star" if $chr eq "*";
	my $bam3 = $dir."/$patient_name.$chr.align.bam";
	die() unless -e $bam3;
	print LIST $bam3."\n";
}
close LIST;
exit(0);
