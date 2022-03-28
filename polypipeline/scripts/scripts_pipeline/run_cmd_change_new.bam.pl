#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use colored;
#use Set::IntSpan;
use GBuffer;
use Data::Dumper;
use Getopt::Long;
use Carp;
use calling_target;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use String::ProgressBar;
use List::Util qw(sum);
use File::Temp;



my $project_name;
my $file2;
my $file1;
my $patient_name;
my $bc;
my $dir_out;
my $fork;
my $dir;
GetOptions(
	'file1=s'   => \$file1,
	'file2=s'   => \$file2,
	'patient=s' => \$patient_name,
	'bc=s' => \$bc,
	'fork=s' =>\$fork,
	'dir=s' => \$dir,
);
my @chromosomes = `samtools idxstats $file1 | cut -f 1,3`;
chomp(@chromosomes);
warn Dumper @chromosomes;
#die();
my $tdir  ="/tmp/";
my $bfile = File::Temp->new( TEMPLATE =>"Pipeline.XXXXXXX",
                        DIR => $dir,
                        SUFFIX => ".bam");

system("cp $file1 $bfile ");
system("cp $file1.bai $bfile.bai ");
#$bfile = $file1;
unless ($dir){
my @dirs = split("/",$file2);
pop(@dirs);
 $dir = join("/",@dirs);
}
my $pm = new Parallel::ForkManager($fork);
my $jobs;
$pm->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    	if (defined $data && exists $data->{job}){
    			
    		my $jobid = $data->{job};
    		delete $jobs->{$jobid};
    		warn "\t\ttfather " . scalar(keys %{$jobs});
    	}
    	
    }
  );

my @files;
my $id = time;
foreach my $l (@chromosomes) {
	my ($chr,$nb) = split(" ",$l);
	warn "skip :".$l if $nb == 0;
	next if $nb == 0;
	my $cmd;
	my $file;
	$id = $chr;
	if  ($chr eq "*"){
		next;
		 $file = $dir."/$patient_name.unaligned.MT.bam";
		 $cmd = qq{samtools view -h  -f 4 $bfile |  perl -lane '\$a =\$_;\$a=~s/chrMT/chrM/g; \$a=~s/:$bc/:$patient_name/;print \$a ' | samtools view -bS - > $file ; samtools index $file};
	}
	else {
		 $file = $dir."/$patient_name.$chr.MT.bam";
		 $cmd = qq{  samtools view -h $bfile $chr  | perl -lane '\$a =\$_;\$a=~s/chrMT/chrM/g; \$a=~s/:$bc/:$patient_name/;print \$a ' | samtools view -bS - > $file ; samtools index $file };
	}
	
	push(@files,$file);
	next if -e $file.".bai";
	$jobs->{$chr} = $cmd;
	
	my $pid = $pm->start and next;
	warn "start ".$chr;
	system($cmd);
	sleep(1);
	my $r ;
	$r->{job} = $chr if -e $file.".bai";
	$r->{chromosome} = $chr;# if -e $file.".bai";
	$pm->finish(0,$r);
	
}
$pm->wait_all_children();

if (scalar keys %{$jobs}){
	warn Dumper $jobs;
	die();
}
my $cmd_merge =  qq{sambamba merge -t $fork $file2 }.join(" ",@files);
warn $cmd_merge;
system("$cmd_merge");
die() unless -e $file2;

foreach my $file (@files) {
	#unlink $file;
	#unlink $file.".bai";
}
unlink $bfile;
unlink $bfile.".bai";