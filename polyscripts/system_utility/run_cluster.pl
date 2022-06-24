#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use Getopt::Long;
#!/usr/bin/perl
use POSIX;
use Proc::Simple;
use Data::Dumper;
use File::Temp qw/tempfile tempdir /;
use File::Slurp;
use String::ProgressBar;
use Term::StatusBar;
#use Timer::Simple;
use Parallel::ForkManager;
 use Digest::MD5 qw(md5 md5_hex md5_base64);
use Term::ANSIColor;
 use warnings;
use IPC::Run3 'run3';
use Time::Piece;
use Time::Seconds qw/ ONE_DAY /;
#$qsub .= "-n 1 --ntasks-per-node=1 --cpus-per-task=$cpus " if( $cpus > 0 );


my $cpus;
my $cmd;
my $process = [];
#$| =1;
my $fork;
my $jobs;
$| =1;
my $limit;
my $noprint;
my $row =0;
GetOptions(
	'cpu=s' => \$cpus,
	'cmd=s' => \$cmd,
	'fork=s' => \$fork,
	'limit=s' => \$limit,
	'row=s' => \$row,
);
unless ($noprint){
for  (my $i=$row;$i<($row+5);$i++){
	print "\033[$i;0H\033[2K";
} 
system("clear") if $row ==0;
}

my $pm = new Parallel::ForkManager($fork);
my $running_jobs ={};
$SIG{'INT'} = sub {
	print "\n wait killing jobs !!!\n";
 	end();
exit(1);	
};
#die() unless $cmd;
die() unless $cpus;

my @cmds;
my $print_stdout;
unless ($cmd){
  @cmds = <STDIN>;
	chomp(@cmds);
	
}
else {
  push(@cmds,$cmd);	
  $print_stdout =1;
}

  $print_stdout =1 if scalar(@cmds) == 1;
my $qsub = "srun ";
my $bds = "bds-cluster";






 my %hcmd;
 my $nb_max = scalar(@cmds);

my $ok = 0;
my $nb = 0;  
my $limit_jobs = 400;
   $limit_jobs = $limit if $limit;  

my $commands;
 my $runnings_jobs;  
foreach my $cmd (@cmds){
	my $hcmd;
	if($cmd =~/:::/){
		my @c = split(/:::/,$cmd);
		$hcmd->{cmd} = $c[0];
		$hcmd->{cpus} = $c[1];
	}
	else {
	$hcmd->{cmd} = $cmd;
	$hcmd->{cpus} = $cpus;
	}
	
	push(@$commands,$hcmd);
#	$hcmd{$cmd} =undef;
}


my $file = "/tmp/tmp.".time."cmd";
	my $t =0;
my $status = new Term::StatusBar (
                    label => 'jobs Done : ',
                   showTime=>1,
                   subTextAlign =>"center",
                    totalItems => scalar(@$commands),  ## Equiv to $status->setItems(10)
                    startRow => $row,
                   
 );
 $status->{maxCol} = 200;
 $limit_jobs =10000 unless $limit_jobs;
 my $nb_jobs = scalar(@$commands);
while (@$commands) {
	
	last if ($t>=$limit_jobs); 
	my $cmd = shift(@$commands);
	$t++;
	run_cmd($cmd);
	
	
}
	$status->subText("Submitting Jobs ...." );

sleep(1);
#warn "end submitted";
$|=1;
#system("clear") unless $noprint;
$status->start() unless $noprint;
my $completed ;
my $failed ;	
my $cancel;
my $start_time = time;

while (keys %$runnings_jobs ){
	my $jobidline = join(",", keys %$runnings_jobs);	
	#print $jobidline."\n";
#	warn "/cm/shared/apps/slurm/16.05.8/bin/sacct -b -n -j $jobidline";
	my @t = `/cm/shared/apps/slurm/16.05.8/bin/sacct -b -n -j $jobidline`;

	chomp(@t);
	my $hstatus;
	my $running = 0;
	my $pending =0;
	foreach my $l (@t){
		my @case = split(" ",$l);
		my $jobid = $case[0];
		next unless exists $jobs->{$jobid};
		my $exit = pop(@case);
		my $stat = pop(@case);
		$jobs->{$jobid}->{status} = $stat;
		
		$hstatus->{$stat}->{$jobid} ++;
		if (lc($stat) eq "failed"){
			system("mv slurm-".$jobid."out"." slurm-".$jobid."failed")
		}
		if (lc($stat) eq "completed"){
			delete $runnings_jobs->{$jobid};
			unlink "slurm-".$jobid.".out";
			$completed->{$jobid} ++;
			unlink "slurm-".$jobid."out";
			if (@$commands){
				my $cmd = shift(@$commands);
				run_cmd($cmd);
			}
			$status->update() unless $noprint;
			
		}
		elsif (lc($stat) eq "pending"){
			$pending ++;
		}
		elsif (lc($stat) eq "running"){
			$running ++;
		}
		elsif (lc($stat) eq "cancel"){
			delete $runnings_jobs->{$jobid};
			$cancel->{$jobid} ++;
			unlink "slurm-".$jobid."out";
			$status->update()  ;
		}
		
		else {
			delete $runnings_jobs->{$jobid};
			if (@$commands){
				my $cmd = shift(@$commands);
				run_cmd($cmd);
			}
			$failed->{$jobid} ++;
			$status->update();
		}
		my $failed = scalar(keys %$failed) +scalar(keys %$cancel);
		my $ok = scalar(keys %$completed);
		my $ttt = Time::Seconds->new( time -$start_time );
			
			#warn
		#$status->subText("Job OK :".colored(['bright_green on_black'],$ok)."     Job Failed: ".colored(['bright_red on_black'],$failed) ) unless $noprint;
		$status->subText(colored(['bright_yellow on_black'],$ttt->pretty));
		
		my $FH = $status->{fh};
		my $r = 4;
		$r = ($row -1) + $r if $row > 0;
		print $FH "\033[$r;0H", (' 'x(5))."total: ".colored(['bright_cyan on_black'],$nb_jobs)." Pending: ".colored(['bright_blue on_black'],$pending)." Running: ".colored(['bright_magenta on_black'],$running)." OK: ".colored(['bright_green on_black'],$ok)."  Failed: ".colored(['bright_red on_black'],$failed)."  ";
		
		#print $FH "\033[5;0H", (' 'x(25)).colored(['bright_yellow on_black'],$ttt->pretty).color('reset');
		#
		#print  $FH " ".$ttt->pretty;
	}
	#sleep(1);
	sleep(5) if keys %$runnings_jobs;
	#system("clear");
}








	end();

DESTROY {
	print "\n wait killing jobs !!!\n";
	
	 end();
	
	
}



sub end {
	foreach my $c (keys %$runnings_jobs){
		$cancel->{$c} ++;
		unlink "slurm-".$c."out";
		system("scancel ".$c);
	}
	my $fail = scalar(keys %$failed) +scalar(keys %$cancel);
	if ($fail == 0){
	

	
	my $hstatus;
	my $running = 0;
	my $pending =0;
	print "\n\n";
	my $jobidline = join(",", keys %$completed);
	my @t = `/cm/shared/apps/slurm/16.05.8/bin/sacct --format=jobid,elapsed -j $jobidline`;
		chomp(@t);
		 print color('bold green');
		 my $nbs;
	foreach my $l (@t){
		my @case = split(" ",$l);
		my $jobid = $case[0];
		next unless exists $jobs->{$jobid};
		print  $jobid." ".$case[-1]."\n";
		my ($h,$m,$s) = split(":",$case[-1]);
		$nbs +=(3600*($h*1)) + (60*$m)+$s;
	}
	
	$nbs /= scalar(keys %$completed);
	$nbs = int($nbs);
	my $val = Time::Seconds->new($nbs);
	
	print " job's done \n ";
	print " mean job's time :  ".$val->pretty." \n";
	 print color('reset');
	exit(0);
}
else {
	my $t =time;
open(ERROR,">job.$t.error") ;

foreach  my $id ( keys %$failed){
	my $cmd = $jobs->{$id}->{cmd};
	
	  print colored(['red'],"$cmd => FAILED \n");
	    print color('reset');
	   print ERROR $cmd."\n";
}

#open(ERROR,">job.$t.error") ;

foreach  my $id (keys %$cancel){
	my $cmd = $jobs->{$id}->{cmd};
	  print colored(['red'],"$cmd => CANCEL \n");
	    print color('reset');
	   print ERROR $cmd."\n";
}
close (ERROR);

print  colored(['red'], scalar(keys %$cancel)." : CANCEL \n");
print  colored(['red '], scalar(keys %$failed)." : FAILED \n");
 print color('reset');
print "list of erro jobs are here job.$t.error\n\n";

exit(1);
}	
exit(0);
}	
	
sub run_cmd {
	my ($cmd) = @_;
	my $out;
	my $cpus = $cmd->{cpus};
	my $script = qq{#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=$cpus
};
	my $run = $script."\n  ".$cmd->{cmd}."\n";

	open(TOTO,">$file") or die();
	print TOTO $run;
	close TOTO;
	my $jobid = `/cm/shared/apps/slurm/16.05.8/bin/sbatch  --parsable $file`;
	unlink $file;
	chomp($jobid);
	#warn $file;
	#run3 ["sbatch  "],\$file,\$out ;
	$cmd->{jobid} = $jobid;
	$jobs->{$jobid}->{cmd} = $cmd->{cmd};
	
	$jobs->{$jobid}->{status} = "submitted";
	$runnings_jobs->{$jobid} ++;
	
}

