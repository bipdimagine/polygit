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
my $filein;
my $bamout;
GetOptions(
	'project=s' => \$project_name,
	'patient=s' => \$patient_name,
	'fork=s' => \$fork,
	'out=s'  => \$out,
	'filein=s' => \$filein,
	'bamout=s' => \$bamout,
);
	
#my $patient_name = "2006N080737_LEBMI";
my $project = $buffer->newProject(-name=>"$project_name");
my $patient = $project->getPatient($patient_name);
my $method = $patient->alignmentMethod();
my $dir = $project->getAlignmentPipelineDir($method."/".$patient->name);
my $final2 = $dir."/".$patient->name.".consensus.umi.bam";
$final2 = $bamout if $bamout;
if($out){
	print $final2;
	exit(0);
}
exit(0) if -e $final2;
my $samtools = $buffer->software("samtools");
my $sambamba = $buffer->software("sambamba");
die() unless $filein;

$| =1;
my $method = $patient->alignmentMethod();
my $dir = $project->getAlignmentPipelineDir($method."/".$patient->name);
die() unless -e $filein;


my $run = $patient->getRun();
my $machine = $run->machine;
my $run_name = $run->plateform_run_name();
my $run_date = $run->date;
my $constructor =  $run->machine_constructor();
my $constructormachine= $constructor."-".$machine ;
my $plateform =  $run->plateform();
my $bar_code= $patient->barcode();
my $project = $patient->getProject;
my $project_name= $project->name();
	
my $picard_path = $buffer->software('picard_path');
my $java = $buffer->software('java');
$java ="java" unless -e $java;
my $picard =  $java." -jar ".$buffer->software('picard_path');
my $filebai = $final2;
$filebai =~ s/bam/bai/;
my $cmd_picard  = $picard." AddOrReplaceReadGroups CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT  INPUT=/dev/stdin "." OUTPUT=".$final2." RGDS=".$project_name." RGLB=".$run_name." RGSM=".$patient_name." RGCN=".$plateform." RGID=".$run_name." RGPL=".$constructormachine." RGPU=1 " ;
#warn $picard." AddOrReplaceReadGroups CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT  INPUT=/dev/stdin "." OUTPUT=".$final2." RGDS=".$project_name." RGLB=".$run_name." RGSM=".$patient_name." RGCN=".$plateform." RGID=".$run_name." RGPL=".$constructormachine." RGPU=1 ";
#die();
#$cmd_picard = $cmd_picard."  && mv $bamout".".bai $final2.bai ";


my $final1 = $dir."/".$patient->name.".tmp.bam";


my $cmd1 =qq{$samtools cat -b $filein >$final1};
warn $cmd1;
system($cmd1) unless -e $final1;
my $file_touch = $dir."/ok.merge";
my $cmd2 = qq{cd /data-beegfs/tmp;$samtools sort $final1 -@ $fork | $cmd_picard && touch $file_touch};
#my $cmd2 = qq{$sambamba sort --tmpdir=/data-beegfs/tmp $final1 -t $fork -o $final2 && touch $file_touch};

warn $cmd2;
system($cmd2) unless -e $final2;

warn "$file_touch";
die() unless -e "$file_touch";
warn "ok";
#unlink "$dir/ok.merge";
exit(0);
