#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use POSIX;
use Term::Spinner;
use Data::Dumper;
use File::Temp qw/tempfile tempdir /;
use File::Slurp;
use String::ProgressBar;
use Term::StatusBar;
use Parallel::ForkManager;
 use Digest::MD5 qw(md5 md5_hex md5_base64);
use Term::ANSIColor;
 use warnings;
use IPC::Run3 'run3';
use Time::Piece;
use Time::Seconds qw/ ONE_DAY /;
use lib "$Bin/../../GenBo/lib/";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../../GenBo/lib/kyoto/";
use lib "$Bin/../../GenBo/lib/GenBoDB/writeDB";
use lib "$Bin/../packages";
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
use GBuffer;
use GenBoProject;
use colored; 
use Text::Table;
use file_util;
use File::Temp qw/ tempfile tempdir /;; 
use Term::Menus;
 use Proc::Simple;

my @running_jobs;

$SIG{INT} = \&tsktsk;
$SIG{KILL} = \&tsktsk;
sub tsktsk {
	warn "kill jobs";
  	foreach my $j (@running_jobs){
  		$j->kill() if $j->poll();
  	}
}

my $projectName;
my $patients_name;
my $limit;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	#'low_calling=s' => \$low_calling,
);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
my $dir_log = $buffer->config->{project_pipeline}->{bds}."/".$projectName.".dragen.".time;
system("mkdir $dir_log && chmod a+rwx $dir_log");

my $patients = $project->get_only_list_patients($patients_name);
my $script_perl = $Bin."/scripts/";

my $status = new Term::StatusBar (
                    label => 'jobs Done : ',
                   showTime=>1,
                   subTextAlign =>"center",
                    totalItems => scalar(@$patients),  ## Equiv to $status->setItems(10)
                   startRow =>1,
 );
 $status->{maxCol} = 200;
 $|=1;
 
system("clear");# unless $noprint;
print colored(['bright_blue on_black'],"2-  DRAGEN MOVE PROGRESS ")."\n";
$status->start();# unless $noprint
my $nb_jobs = scalar(@$patients);
my $pending = $nb_jobs;
my $ok = 0;
my $failed =0;
my $running = 0;
my $FH = $status->{fh};
print $FH "\033[5;0H", (' 'x(5))."total: ".colored(['bright_cyan on_black'],$nb_jobs)." Pending: ".colored(['bright_blue on_black'],$pending)." Running: ".colored(['bright_magenta on_black'],$running)." OK: ".colored(['bright_green on_black'],$ok)."  Failed: ".colored(['bright_red on_black'],$failed)."  ";
	
my $start_time = time;
my $jobs =[];



foreach my $patient (@$patients){
	my $dir_pipeline = $patient->getDragenDir("pipeline");
	my $prefix = $patient->name;
	my $bam_pipeline = $dir_pipeline."/".$prefix.".bam";
	
	next if -e $patient->getBamFileName("dragen-align");
	next if -e $bam_pipeline;
	$running = 1;
	$pending --;
	my $job;
	$job->{name} = $patient->name;
	$job->{cmd} = "perl $script_perl/dragen_align_calling.pl -project=$projectName -patient=".$patient->name;
	my $myproc = Proc::Simple->new(); 
	$myproc->redirect_output ($dir_log."/".$patient->name."align.log", $dir_log."/".$patient->name."align.err");
	push(@running_jobs,$myproc);
	my $t1 = time;
	$myproc->start($job->{cmd});
	
	my $FH = $status->{fh};
	print $FH "\033[5;0H", (' 'x(5))."total: ".colored(['bright_cyan on_black'],$nb_jobs)." Pending: ".colored(['bright_blue on_black'],$pending)." Running: ".colored(['bright_magenta on_black'],$running)." OK: ".colored(['bright_green on_black'],$ok)."  Failed: ".colored(['bright_red on_black'],$failed)."  ";
	
	#$myproc->start("sleep 5");
	while ($myproc->poll()){
		my $ttt = Time::Seconds->new( time -$start_time );
		$status->subText(colored(['bright_yellow on_black'],$ttt->pretty));
		sleep 5;
	}
	my $ttt = Time::Seconds->new( time - $t1 );
	$job->{elapse} = $ttt->pretty;
	$job->{time} = time - $t1;
	if ($myproc->exit_status() == 0){
		$ok ++;
		$job->{ok} ++;
	}
	else {
		$failed ++;
		$job->{failed} ++;
	}
	push(@$jobs,$job);
	 $FH = $status->{fh};
	print $FH "\033[4;0H", (' 'x(5))."total: ".colored(['bright_cyan on_black'],$nb_jobs)." Pending: ".colored(['bright_blue on_black'],$pending)." Running: ".colored(['bright_magenta on_black'],$running)." OK: ".colored(['bright_green on_black'],$ok)."  Failed: ".colored(['bright_red on_black'],$failed)."  ";
	$status->update();
		
}
system("clear");
my $nbo =0;
my $sum_time =0;
my (@error) = grep {exists $_->{failed}} @$jobs;
my (@ok) = grep {exists $_->{ok}} @$jobs;
my $sum = 0;
map {$sum+=$_->{time}} @ok;
my $mean =0;
$mean = int($sum/scalar(@ok)) if @ok; 
my $ttt = Time::Seconds->new($mean);
my $tt1 = Time::Seconds->new($sum);
my @text;
	
if (@error) {
	foreach my $job (@$jobs){
		print colored(['bright_red on_black'],$job->{name}." ".$job->{elapse})."\n";
		print "log : ".$dir_log."/".$job->{name}."align.err";
	}
	die();
}
if (@$jobs){
print colored(['bright_green on_black'],"1-  DRAGEN ALIGN MEAN TIME : ".$ttt->pretty." total ".$tt1->pretty)."\n";
}
else {
	print colored(['bright_cyan on_black'],"1-  DRAGEN ALIGN NONE ")."\n";
}
print colored(['bright_blue on_black'],"2-  DRAGEN MOVE PROGRESS ")."\n";

my $st_move = time;
my $cmds =[];
foreach my $patient (@$patients){
	#next if -e $patient->getBamFileName("dragen-align");
	my $dir_pipeline = $patient->getDragenDir("pipeline");
	push (@$cmds,"perl $script_perl/dragen_move.pl -project=$projectName -patient=".$patient->name);
}
if (@$cmds){
open (CLUSTER, " | run_cluster.pl -cpu=10  -limit=4 -row=4");
foreach my $cmd (@$cmds){
	print CLUSTER $cmd."\n";
	
}

close(CLUSTER);
}
my $ttt_move = Time::Seconds->new( time - $st_move);
system("clear");
print colored(['bright_green on_black'],"1-  DRAGEN ALIGN OK :".$tt1->pretty." MEAN :".$ttt->pretty)."\n";
print colored(['bright_green on_black'],"2-  DRAGEN MOVE OK  :".$ttt_move->pretty)."\n";
print colored(['bright_cyan on_black'],"3-  DRAGEN GENOTYPE RUNNING  : ");

my $spinner = Term::Spinner->new(output_handle => \*STDOUT,);

my $cmd_genotype = "perl $script_perl/dragen_genotype.pl -project=$projectName -patient=$patients_name ";
#$cmd_genotype = "sleep 5 ";
my $myproc1 = Proc::Simple->new(); 
my $st1 = time;

	push(@running_jobs,$myproc1);
	#$myproc1->signal_on_destroy("KILL");
	my $t1 = time;
	$myproc1->redirect_output ($dir_log."/genotype.log", $dir_log."/genotype.err");
	$myproc1->start($cmd_genotype);
	my $nb = 0;
	while ($myproc1->poll()){
		$spinner->draw();
		$spinner->advance();
		sleep 1;
	}
	undef $spinner;
 my $ttt_genotype = Time::Seconds->new( time - $st1);	
if ($myproc1->exit_status()!= 0){
	print "\033[3;0H\033[2K";
	print colored(['bright_red on_black'],"3-  DRAGEN GENOTYPE ERROR  ")."\n";
	die("genotype");
}

print "\033[3;0H\033[2K";
print colored(['bright_green on_black'],"3-  DRAGEN GENOTYPE   OK ".$ttt_genotype->pretty)."\n";
print colored(['bright_cyan on_black'],"4-  LMDB DEPTH  ")."\n";




my $cmd_pipeline = "$Bin/../bds_pipeline.pl -project=$projectName -config=$Bin/pipeline.cfg -pipeline=dragen -yes=1";
my $exit_pipeline = system($cmd_pipeline);
$t1 = time;
my $ttt_pipeline = Time::Seconds->new( time - $st1);
if ($exit_pipeline== 0){
	system("clear");
	print colored(['bright_green on_black'],"1-  DRAGEN ALIGN OK :".$tt1->pretty." MEAN :".$ttt->pretty)."\n";
	print colored(['bright_green on_black'],"2-  DRAGEN MOVE OK  :".$ttt->pretty)."\n";
	print colored(['bright_green on_black'],"3-  DRAGEN GENOTYPE   OK :".$ttt_genotype->pretty)."\n";
	print colored(['bright_green on_black'],"4-  USUAL PIPELINE  OK :".$ttt_pipeline->pretty)."\n";
	
}
else {
	print "======================================================================================\n";
	print colored(['bright_green on_black'],"1-  DRAGEN ALIGN OK :".$tt1->pretty." MEAN :".$ttt->pretty)."\n";
	print colored(['bright_green on_black'],"2-  DRAGEN MOVE OK  :".$ttt->pretty)."\n";
	print colored(['bright_green on_black'],"3-  DRAGEN GENOTYPE   OK ".$ttt_genotype->pretty)."\n";
	print colored(['bright_red on_black'],"4-  USUAL PIPELINE  ERROR".$ttt_pipeline->pretty)."\n";
	print "======================================================================================\n";
}
print colored(['bright_cyan on_black'],"5-  DRAGEN CNV  RUNNIG : ")."\n";;
my $status2 = new Term::StatusBar (
                    label => 'jobs Done : ',
                   showTime=>1,
                   subTextAlign =>"center",
                    totalItems => scalar(@$patients),  ## Equiv to $status->setItems(10)
                   startRow =>6,
 );
$status2->start();
my $exit;
my $tcnv = time;
foreach my $patient (@$patients){
	my $cmd = qq{perl $script_perl/dragen_cnv.pl -project=$projectName -patient=}.$patient->name;
	 $exit = system($cmd);
	last if $exit != 0;
	$status->update();
	
}
my $ttt_cnv = Time::Seconds->new( time - $tcnv);
system("clear");
if ($exit == 0) {
	print "======================================================================================\n";
	print colored(['bright_green on_black'],"1-  DRAGEN ALIGN OK :".$tt1->pretty." MEAN :".$ttt->pretty)."\n";
	print colored(['bright_green on_black'],"2-  DRAGEN MOVE OK  :".$ttt->pretty)."\n";
	print colored(['bright_green on_black'],"3-  DRAGEN GENOTYPE   OK ".$ttt_genotype->pretty)."\n";
	print colored(['bright_green on_black'],"4-  USUAL PIPELINE  OK ".$ttt_pipeline->pretty)."\n";
	print colored(['bright_green on_black'],"5-  DRAGEN CNV  OK  : ".$ttt_cnv->pretty)."\n";
	print "======================================================================================\n";
}
else {
	print "======================================================================================\n";
	print colored(['bright_green on_black'],"1-  DRAGEN ALIGN OK :".$tt1->pretty." MEAN :".$ttt->pretty)."\n";
	print colored(['bright_green on_black'],"2-  DRAGEN MOVE OK  :".$ttt->pretty)."\n";
	print colored(['bright_green on_black'],"3-  DRAGEN GENOTYPE   OK ".$ttt_genotype->pretty)."\n";
	print colored(['bright_green on_black'],"4-  USUAL PIPELINE  OK ".$ttt_pipeline->pretty)."\n";
	print colored(['bright_red on_black'],"5-  DRAGEN CNV  ERROR  : ".$ttt_cnv->pretty)."\n";
	print "======================================================================================\n";
}

#
#foreach my $job (@$jobs){
#	if (exists $job->{ok}){	
#		print colored(['bright_green on_black'],$job->{name}." ".$job->{elapse})."\n";
#		$sum_time += $job->{time};
#		$nbo ++;
#	}
#	else {
#		print colored(['bright_red on_black'],$job->{name}." ".$job->{elapse})."\n";	
#		#print colored(['bright_green on_black'].$job->name." ".$job->elapse)."\n";
#	}
#	
#}
#
#if ($nbo > 0) {
#	$sum_time = int($sum_time/$nbo);
#	my $ttt = Time::Seconds->new( $sum_time);
#	print "\n------------------------\n";
#	print colored(['bright_green on_black'],"mean job time : ".$ttt->pretty)."\n";
#	print "\n------------------------\n"
#}

sub run_cluster {
	
}
