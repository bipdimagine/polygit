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
   
#$qsub .= "-n 1 --ntasks-per-node=1 --cpus-per-task=$cpus " if( $cpus > 0 );


my $cpus;
my $cmd;
my $process = [];
#$| =1;
my $fork;

$| =1;

GetOptions(
	'cpu=s' => \$cpus,
	'cmd=s' => \$cmd,
	'fork=s' => \$fork,
);


 $fork = int((40*16)/$cpus) unless $fork;

my $pm = new Parallel::ForkManager($fork);
my $running_jobs ={};
$SIG{'INT'} = sub {
	print "\n wait killing jobs !!!\n";
	 foreach my $h (keys %$running_jobs){
	 	system("scancel -n $h");
	 }
	 
	 my @pids = $pm->running_procs;
	 
	 foreach my $proc (@pids){
	 	system("kill $proc");
	}

	sleep(5);
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

$qsub .= " --cpus-per-task=$cpus -n 1 -N 1 " if( $cpus > 0 );



#$myproc->redirect_output ("toto.out", undef);
#my $syntax = "perl -e 'print \"#!/bin/sh \\n".$cmd."\"'| $qsub";
my $process;


 my %hcmd;
 my $nb_max = scalar(@cmds);
my $status = new Term::StatusBar (
                    label => 'jobs Done : ',
                   showTime=>1,
                   subTextAlign =>"left",
                    totalItems => $nb_max,  ## Equiv to $status->setItems(10)
                   
    );
my $ok = 0;
my $failed = 0;    
my $nb = 0;  

  system("clear");     
$pm->run_on_finish(
    sub { 
    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
  
    	unless (defined($h) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
	delete 		 $running_jobs->{$h->{digest}};
    
    	$nb ++;
    	warn "end jobs ".$h->{status};
    	if ($h->{status} ==0 ){
    		$ok ++;
    			delete $hcmd{$h->{cmd}};
    	}
    	else {
    		$failed ++;
    	}
    	
    	# update($status,1);
    	 my @pids = $pm->running_procs;
    	 my $nb_run = scalar(@pids);
    	 $status->update();
 		#$status->subText("\e[42mOK : ".$ok ."\e[0m  \e[41mFailed : ".$failed."\e[0m   running : $nb_run " );
 		$status->subText("OK : ".$ok ."  Failed : ".$failed."   total: $nb_max  running : $nb_run  " );
    	
    }
    
    );
foreach my $cmd (@cmds){
	$hcmd{$cmd} =undef;
}


my $hpid;
$status->start;
foreach my $cmd (@cmds){
	my $t = time;
	 my $digest = md5_hex(time."-".$cmd."-".rand(50000));
	 $running_jobs->{$digest}++;
	my $pid = $pm->start and next;

	my $syntax = "$qsub --job-name=$digest perl -e 'system(\"$cmd\") ==0  or die \"system  failed: $?\" ' ";#
	$syntax .= ">/dev/null 2>/dev/null" unless $print_stdout;
	my $return = system($syntax) ;
	$status->update();
	$pm->finish(0,{cmd=>$cmd,status=>$return,digest=>$digest});
	}
	$pm->wait_all_children();
	print "\n\n\n";
	end();

DESTROY {
	print "\n wait killing jobs !!!\n";
	
	 end();
	
	
}

sub update {
	my ($st,$nb) = @_;
	 for (my $i=0;$i<$nb;$i++){
 	 $st->update();
 	}
 
}

sub end {
my $failed = scalar keys(%hcmd); 	
 my @pids = $pm->running_procs;
 
 foreach my $pid (@pids){
 	system("kill $pid");
 }
 
if ($failed == 0){
	 print color('bold green');
	print " job's done \n ";
	 print color('reset');
	exit(0);
}
else {
	my $t =time;
open(ERROR,">job.$t.error") ;
foreach my $cmd (keys %hcmd){
	  print colored(['red '],"$cmd => FAILED \n");
	    print color('reset');
	   print ERROR $cmd."\n";
}
close (ERROR);
print "list of erro jobs are here job.$t.error";
exit(1);
}	
	
}

