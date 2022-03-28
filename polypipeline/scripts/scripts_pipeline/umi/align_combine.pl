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
my $json_file;
GetOptions(
	'project=s' => \$project_name,
	'patient=s' => \$patient_name,
	'fork=s' => \$fork,
	'out=s'  => \$out,
	'filein=s'  => \$json_file,
);



#my $patient_name = "2006N080737_LEBMI";
my $project = $buffer->newProject(-name=>"$project_name");
my $patient = $project->getPatient($patient_name);

my $index = $project->getGenomeIndex("bwa");#
$index .= "/all.fa";#if $prog eq "bwa";
my $method = $patient->alignmentMethod();
my $dir = $project->getAlignmentPipelineDir($method."/".$patient->name);
my $final2 = "$dir/list_file.txt";
$fork=40 unless $fork;
my $pm = new Parallel::ForkManager(int($fork/4));
 my $reference = $project->genomeFasta();
if($out){
	print $final2;
	exit(0);
}

die("miss json file ".$dir."/files.json  vs $json_file") unless -e $json_file;
open (JSON,"$json_file");
my $files = decode_json(<JSON>);
	
	
	my $nb_files = scalar @$files;
	my @bamfiles;
	$pm->run_on_finish(
    	sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    		push(@bamfiles,$data->{file});
    }
    );		

	
my @bam_files;
my $hh;
my $nb =0;
foreach my $hfile (@$files){
	$nb ++;
	#next if exists $hh->{$hfile->{ubam}};
	$hh->{$hfile->{ubam}} ++;
	warn "-----************ ".$nb."/".$nb_files."**************";
	my $pid = $pm->start and next;
	
	my $bam = align_and_sort($hfile);
	warn "°°°°°°°°°°°°° ".$nb."/".$nb_files."°°°°°°°°°°°°°°";
	$pm->finish(0,{file=>$bam});
}
$pm->wait_all_children();
warn "end fork";
warn $final2;
warn "######################################";
open (LIST,">$final2");
print LIST join("\n",@bamfiles);
close(LIST);


exit(0);


sub align_and_sort {
	my ($hfile) = @_;
	my $fastq1 = $hfile->{fastq1};
	my $fastq2 =$hfile->{fastq2};
	my $ubam = $hfile->{ubam};
	my $bam = $ubam;
	$bam =~ s/u\.bam/align.merge\.bam/;
	return $bam if -e $bam;
	die() unless -e $fastq1;
	die() unless -e $fastq2;
	die() unless -e $ubam;
	my $picard = $buffer->software("picard");
	my $java = $buffer->software("java");
	my $bwa = $buffer->software("bwa");
	my $cmd = qq{$bwa  mem $index $fastq1 $fastq2 -t 4 2>/dev/null| $java -Djava.io.tmpdir=/data-beegfs/tmp -jar $picard SortSam -I /dev/stdin -O /dev/stdout -SORT_ORDER queryname  2>/dev/null |  $java -Djava.io.tmpdir=/data-beegfs/tmp -jar $picard MergeBamAlignment CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT ALIGNER_PROPER_PAIR_FLAGS=false R=$reference ALIGNED=/dev/stdin UNMAPPED=$ubam OUTPUT=$bam 2>/dev/null};
	warn $cmd;
	system($cmd);
	warn "end";
	return $bam;
}

exit(0);