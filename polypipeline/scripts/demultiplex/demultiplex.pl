#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use Getopt::Long;
use Data::Dumper;

my $dir ="";
my $run;
my $file2;
my $hiseq;
my $sample_sheet;
my $mask;
my $single;
my $mismatch;
my $dragen;
my $cellranger_type;

GetOptions(
	'dir=s' => \$dir,
	'run=s' => \$run,
	'hiseq=s' => \$hiseq,
	'sample_sheet=s' => \$sample_sheet,
	'mask=s' => \$mask,
	'single=s' => \$single,
	'mismatch=s' => \$mismatch,
	'dragen=s' => \$dragen,
	'cellranger_type=s' => \$cellranger_type,
);


my $basecall_dir =  "$dir/Data/Intensities/BaseCalls/";
die("no basecall dir $basecall_dir") unless -e  $basecall_dir;
exit(0) unless -e  $basecall_dir;
if ($sample_sheet){
	die("problem with sample sheet : $sample_sheet") unless -e $sample_sheet;
}
else {
	die "no default " unless -e $basecall_dir."/SampleSheet.csv";
	warn "use $basecall_dir/SampleSheet.csv";
}
my $cmd_cfg = qq{ /software/bin/bcl2fastq };
die($dir. "not exists") unless -e $dir;
die() unless  $run;
die("hiseq: HISEQ2500-1 or HISEQ2500-2 or NOVASEQ or ISEQ or NEXTSEQ500 or NOVASEQ-Abel or 10X or MISEQ") if ($hiseq ne "HISEQ2500-1" && $hiseq ne 'HISEQ2500-2' && $hiseq ne 'NOVASEQ' && $hiseq ne 'ISEQ' && $hiseq ne 'NEXTSEQ500' && $hiseq ne 'NOVASEQ-Abel' && $hiseq ne '10X' && $hiseq ne 'MISEQ');
my $root= "/data-isilon/data/sequences";
my $tmp = "/data-beegfs/tmp/$run";
my $reports_dir = "/data-isilon/sequencing/ngs/demultiplex/$run";

if (-d $tmp){
	die("$tmp already exists ");
	die();
}

my $final_dir = $root."/ILLUMINA/$hiseq/IMAGINE/$run";
$final_dir = $root."/ILLUMINA/$hiseq/LAVOISIER/$run" if $hiseq eq 'NEXTSEQ500';
$final_dir = $root."/ILLUMINA/$hiseq/HEMATO/$run" if $hiseq eq 'MISEQ';

if (-d $final_dir){
	die($final_dir. " already exists "); 
}

mkdir $tmp;
my $wtmp = "$tmp/Unaligned";
warn $dir;
my $complete = "$dir/RTAComplete.txt";
warn $complete;
my $checkComplete=1;
$checkComplete =0 if -f $complete;
while($checkComplete ==1){
	warn "sleep 3300s";
	sleep(3300);
	$checkComplete =0 if -f $complete;
}

my $basecall_dir =  "$dir/Data/Intensities/BaseCalls/";
die($basecall_dir ." doesnt exists") unless -e $basecall_dir;
my $cmd1 = "$cmd_cfg --input-dir $basecall_dir  --output-dir=$wtmp --ignore-missing-bcl --ignore-missing-positions --ignore-missing-filter --ignore-missing-control --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0  --reports-dir $reports_dir -R $dir";
if($hiseq eq "10X"){    
	$cmd1= "cd $tmp ; /software/distrib/$cellranger_type/latest mkfastq --id=$run --run=$dir --csv=$sample_sheet --output-dir=$wtmp";
	$cmd1 = $cmd1. " --barcode-mismatches=".$mismatch if $mismatch;
	$cmd1 = $cmd1. " --barcode-mismatches=0" unless $mismatch;
}
elsif($dragen == 1){
	my $dragen_output = "/staging/RUN/".$run;
	$cmd1="	dragen --bcl-input-directory $dir  --output-directory $dragen_output --bcl-conversion-only true";
	$cmd1 = $cmd1. " --sample-sheet=$sample_sheet" if $sample_sheet ;
warn $cmd1;
}

else{
	if ($mismatch){
		$cmd1 = $cmd1. " --barcode-mismatches=".$mismatch;
	}
	else{
		$cmd1 = $cmd1. " --barcode-mismatches=0";
	}
	$cmd1 = $cmd1. " --sample-sheet=$sample_sheet" if $sample_sheet ;
	$cmd1 = $cmd1. " --use-bases-mask=$mask" if $mask;S
	$cmd1 = $cmd1. " --ignore-dual-index" if $single;
	warn $cmd1;
}
if($dragen == 1){
	warn $cmd1;
	system ("ssh masson\@10.200.27.109 $cmd1");
}
else{
	warn $cmd1;
	system("cd $tmp &&"."/software/bin/run_cluster.pl -cpu=40 -cmd=\"$cmd1 && touch $tmp/done.txt\"") if (-e "/software/bin/run_cluster.pl");
	system("cd $tmp &&"."$cmd1 "." && touch done.txt") unless (-e "/software/bin/run_cluster.pl");
	if ($? == -1) {
 	 print "configure failed !!!!!!!! \n";
 	 rmdir $wtmp;
 	 rmdir $tmp;
  	 die();
	}
}


if (-e $tmp."/done.txt"){
mkdir $final_dir;
system("move_fastq.sh $wtmp  $final_dir");
system("rm $final_dir/*Undetermined");
my $fileReport= $run."_laneBarcode.html";
system("cp $tmp/Unaligned/Reports/html/*/all/all/all/laneBarcode.html /data-isilon/sequencing/ngs/demultiplex/$run/$fileReport");
#system("rm -rf $tmp");

warn "count_sequence.pl -dir=$final_dir";
system("count_sequence.pl -dir=$final_dir");
}
else {
	exit (1);
}
#echo find /data-ibrix/data/sequences/tmp/$2/Unaligned -name *.gz | xargs mv -t /data-ibrix/data/sequences/Hiseq_run/$2/
